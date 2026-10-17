"""Fictional repeat/metadata cases for IMPROVE Phase 1.3; no real data."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from research.oncology.improve_split_audit import BASE, STUDIES
from research.oncology.improve_repeat_provenance_audit import audit_repeat_provenance


def make_zip(path, *, change=None):
    rows = ["source\tstudy\timprove_sample_id\timprove_chem_id\tauc\n"]
    files = {}
    for study in STUDIES:
        start = len(rows) - 1
        triples = [("S1", "D1"), ("S1", "D1"), ("S2", "D2")]
        for sample, drug in triples:
            rows.append(f"{study}\t{study}\t{sample}\t{drug}\tSECRET_LABEL\n")
        files[BASE + f"splits/{study}_all.txt"] = "".join(
            f"{i}\n" for i in range(start, start + 3))
        for part, ix in zip(("train", "val", "test"), range(start, start + 3)):
            files[BASE + f"splits/{study}_split_0_{part}.txt"] = str(ix) + "\n"
    files[BASE + "y_data/response.tsv"] = "".join(rows)
    if change:
        change(files)
    with ZipFile(path, "w") as z:
        for key, value in files.items():
            z.writestr(key, value)


class RepeatProvenanceAuditTests(unittest.TestCase):
    def test_five_studies_and_repeats_are_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            make_zip(p)
            result = audit_repeat_provenance(p)
            self.assertEqual(result["response_rows"], 15)
            self.assertEqual(len(result["studies"]), 5)
            for v in result["studies"].values():
                self.assertEqual(v["repeated_pair_groups"], 1)
                self.assertEqual(v["repeated_pair_excess_rows"], 1)
                self.assertEqual(v["repeat_groups_spanning_partitions"], 1)
                self.assertFalse(v["confirmed_assay_occasion_identity"])

    def test_no_biological_ids_or_response_values_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            make_zip(p)
            output = str(audit_repeat_provenance(p))
            for secret in ("SECRET_LABEL", "'S1'", "'D1'", "'S2'", "'D2'"):
                self.assertNotIn(secret, output)

    def test_header_and_source_mismatch_exposed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            def change(files):
                key = BASE + "y_data/response.tsv"
                files[key] = files[key].replace("CCLE\tCCLE\tS1", "CCLE\tgCSI\tS1", 1)
            make_zip(p, change=change)
            r = audit_repeat_provenance(p)
            self.assertEqual(r["studies"]["CCLE"]["study_label_not_equal_to_source_name"], 1)
            self.assertEqual(r["studies"]["CCLE"]["source_label_mismatches"], 0)
            self.assertEqual(r["studies"]["CCLE"]["distinct_study_field_values"], 2)

    def test_study_identifier_values_not_in_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            def change(files):
                key = BASE + "y_data/response.tsv"
                files[key] = files[key].replace(
                    "CCLE\\tCCLE\\tS1", "CCLE\\tSTUDY_IDENTIFIER_SECRET_987\\tS1", 1)
            make_zip(p, change=change)
            report = audit_repeat_provenance(p)
            self.assertNotIn("STUDY_IDENTIFIER_SECRET_987", str(report))
            self.assertIn("study_field_value_multiplicity_histogram", report["studies"]["CCLE"])

    def test_missing_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            make_zip(p, change=lambda f: f.pop(BASE + "splits/CCLE_split_0_test.txt"))
            with self.assertRaises(KeyError):
                audit_repeat_provenance(p)

    def test_invalid_fold_rejected(self):
        with self.assertRaisesRegex(ValueError, "fold must"):
            audit_repeat_provenance(Path("none.zip"), 10)

    def test_pair_split_overlap_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "fictional.zip"
            def change(files):
                files[BASE + "splits/CCLE_split_0_test.txt"] = "0\n2\n"
            make_zip(p, change=change)
            with self.assertRaisesRegex(ValueError, "split row overlap"):
                audit_repeat_provenance(p)


if __name__ == "__main__":
    unittest.main()
