"""Read-only ZIP audit tests using entirely fictional temporary archives."""
import json
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile
from research.oncology.improve_archive_audit import (
    BASE, REQUIRED_TABLES, STUDIES, audit_archive,
)


def archive(path, *, full_splits=False, change=None):
    payload = {
        BASE + name: "sample\tresponse_metric\nFAKE\t0.5\n"
        for name in REQUIRED_TABLES
    }
    payload[BASE + "x_data/cancer_gene_expression.tsv"] = (
        "sample\tgene_a\tgene_b\nSYN\t1\t2\n"
    )
    payload[BASE + "y_data/response.tsv"] = (
        "dataset\tcell_id\tdrug_id\tAUC\nSYN\t1\t2\t0.1\n"
    )
    for study in STUDIES:
        payload[BASE + f"splits/{study}_all.txt"] = "0\n"
        if full_splits:
            for i in range(10):
                for kind in ("train", "val", "test"):
                    payload[BASE + f"splits/{study}_split_{i}_{kind}.txt"] = "0\n"
    if change:
        change(payload)
    with ZipFile(path, "w") as output:
        for name, content in payload.items():
            output.writestr(name, content)


class ImproveArchiveAuditTests(unittest.TestCase):
    def test_metadata_and_split_names_only(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            archive(path, full_splits=True)
            result = audit_archive(path)
            self.assertEqual(result["missing_required_tables"], [])
            self.assertTrue(all(f["present"] == 31 for f in result["split_manifests"].values()))
            self.assertEqual(result["unrecognized_split_filenames"], [])
            self.assertIn("y_data/response.tsv", {f["table"] for f in result["table_headers"]})
            self.assertEqual(list(Path(d).iterdir()), [path])

    def test_incomplete_splits_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            archive(path)
            result = audit_archive(path)
            self.assertEqual(result["split_manifests"]["CCLE"]["present"], 1)
            self.assertEqual(len(result["split_manifests"]["CCLE"]["missing"]), 30)

    def test_no_sample_rows_printed_or_reported(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            archive(path)
            rendered = json.dumps(audit_archive(path))
            self.assertNotIn("FAKE", rendered)
            self.assertNotIn('"0.1"', rendered)

    def test_reject_unsafe_zip_member(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "bad.zip"
            with ZipFile(path, "w") as file:
                file.writestr("../escape.tsv", "x\n")
            with self.assertRaisesRegex(ValueError, "unsafe ZIP"):
                audit_archive(path)

    def test_duplicate_zip_member_rejected(self):
        import warnings
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "duplicate.zip"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                with ZipFile(path, "w") as file:
                    file.writestr(BASE + "x_data/drug_info.tsv", "a\n")
                    file.writestr(BASE + "x_data/drug_info.tsv", "b\n")
            with self.assertRaisesRegex(ValueError, "duplicate ZIP"):
                audit_archive(path)

    def test_duplicate_tsv_columns_flagged_without_shifting(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "columns.zip"
            archive(path, change=lambda data: data.update({
                BASE + "y_data/response.tsv": "drug\tdrug\nA\tB\n"
            }))
            result = audit_archive(path)
            entry = next(t for t in result["table_headers"]
                         if t["table"] == "y_data/response.tsv")
            self.assertEqual(entry["status"], "needs_review")
            self.assertEqual(entry["column_count"], 2)
            self.assertIn("duplicate_named_columns", entry["header_quality_flags"])
            self.assertEqual(entry["duplicate_named_column_count"], 1)
            self.assertEqual(entry["duplicate_named_columns_preview"], ["drug"])

    def test_untrusted_torch_file_not_deserialized(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            archive(path, change=lambda data: data.update({
                BASE + "processed/train_data.pt": b"NOT VALID PYTORCH CONTENT"
            }))
            result = audit_archive(path)
            self.assertIn("processed/train_data.pt",
                          result["serialized_python_artifacts_not_loaded"])

    def test_fictional_md5_does_not_pass_published_md5(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            archive(path)
            self.assertFalse(audit_archive(path, verify_md5=True)
                             ["integrity"]["matches_published_md5"])

    def test_wide_methylation_header_over_64kb(self):
        # Original defect: 64-KiB readline cut off a valid wide header.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "wide.zip"
            columns = [f"cg{i:08d}" for i in range(10000)]
            header = "sample\t" + "\t".join(columns) + "\n"
            self.assertGreater(len(header), 65536)
            archive(path, change=lambda data: data.update({
                BASE + "x_data/cancer_DNA_methylation.tsv": header + "SYN\t0\n"
            }))
            result = audit_archive(path)
            methylation = next(h for h in result["table_headers"]
                               if h["table"] == "x_data/cancer_DNA_methylation.tsv")
            self.assertEqual(methylation["status"], "ok")
            self.assertEqual(methylation["column_count"], 10001)
            self.assertEqual(methylation["columns_preview"][0], "sample")
            self.assertEqual(len(methylation["columns_preview"]), 16)

    def test_header_without_newline_at_end_of_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "header_only.zip"
            archive(path, change=lambda data: data.update({
                BASE + "x_data/drug_info.tsv": "drug_id\tsmiles"
            }))
            result = audit_archive(path)
            drug = next(h for h in result["table_headers"]
                        if h["table"] == "x_data/drug_info.tsv")
            self.assertEqual(drug["status"], "ok")
            self.assertEqual(drug["column_count"], 2)

    def test_oversized_header_is_reported_without_abort(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "oversize.zip"
            archive(path, change=lambda data: data.update({
                BASE + "x_data/cancer_DNA_methylation.tsv": "cg00000000\t" * 20 + "end\n"
            }))
            with patch("research.oncology.improve_archive_audit.MAX_HEADER_BYTES", 100):
                result = audit_archive(path)
            self.assertEqual(result["headers_exceeding_safety_cap"],
                             ["x_data/cancer_DNA_methylation.tsv"])
            oversized = next(h for h in result["table_headers"]
                             if h["table"] == "x_data/cancer_DNA_methylation.tsv")
            self.assertEqual(oversized["status"], "header_exceeds_safety_cap")
            self.assertIsNone(oversized["column_count"])
            self.assertTrue(result["response_table_present"])

    def test_leading_unnamed_methylation_index_is_preserved(self):
        # A wide molecular table may be exported with an empty index heading.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "leading_blank.zip"
            features = [f"cg{i:08d}" for i in range(10000)]
            header = "\t" + "\t".join(features) + "\n"
            self.assertGreater(len(header), 65536)
            archive(path, change=lambda data: data.update({
                BASE + "x_data/cancer_DNA_methylation.tsv": header + "SYN\t0\n"
            }))
            result = audit_archive(path)
            entry = next(t for t in result["table_headers"]
                         if t["table"] == "x_data/cancer_DNA_methylation.tsv")
            self.assertEqual(entry["status"], "needs_review")
            self.assertEqual(entry["column_count"], 10001)
            self.assertEqual(entry["unnamed_column_count"], 1)
            self.assertEqual(entry["unnamed_column_positions_preview"], [0])
            self.assertIsNone(entry["columns_preview"][0])
            self.assertEqual(entry["columns_preview"][1], features[0])
            self.assertIn("unnamed_columns", entry["header_quality_flags"])
            self.assertTrue(result["response_table_present"])

    def test_multiple_blank_positions_and_trailing_tab(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "multiple_blank.zip"
            archive(path, change=lambda data: data.update({
                BASE + "x_data/drug_info.tsv": "\tdrug_id\t\tSMILES\t\nSYN\t1\t2\t3\t4\n"
            }))
            result = audit_archive(path)
            entry = next(t for t in result["table_headers"]
                         if t["table"] == "x_data/drug_info.tsv")
            self.assertEqual(entry["column_count"], 5)
            self.assertEqual(entry["unnamed_column_positions_preview"], [0, 2, 4])
            self.assertEqual(entry["unnamed_column_count"], 3)
            self.assertEqual(entry["columns_preview"], [None, "drug_id", None, "SMILES", None])
            self.assertIn("x_data/drug_info.tsv", result["headers_needing_review"])

    def test_all_unnamed_columns_are_flagged_not_silently_assigned(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "all_blank.zip"
            archive(path, change=lambda data: data.update({
                BASE + "x_data/drug_info.tsv": "\t\t\nx\ty\tz\n"
            }))
            result = audit_archive(path)
            entry = next(t for t in result["table_headers"]
                         if t["table"] == "x_data/drug_info.tsv")
            self.assertEqual(entry["status"], "needs_review")
            self.assertEqual(entry["unnamed_column_count"], 3)
            self.assertEqual(entry["column_count"], 3)
            self.assertEqual(entry["columns_preview"], [None, None, None])

    def test_missing_zip_error(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "does not exist"):
                audit_archive(Path(d) / "absent.zip")


if __name__ == "__main__":
    unittest.main()
