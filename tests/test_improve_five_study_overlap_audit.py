"""Synthetic five-source IMPROVE biological overlap tests. No patient records."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile
from research.oncology.improve_split_audit import BASE, STUDIES
from research.oncology.improve_five_study_overlap_audit import audit_five_study_overlap


def fictitious_archive(path: Path, alter=None):
    tsv = ["source\timprove_sample_id\timprove_chem_id\tstudy\tauc\n"]
    parts = {}
    for source in STUDIES:
        pos = len(tsv) - 1
        if source == "CCLE":
            data = [("S_A", "D_A"), ("CCLE_ONLY", "D_0"),
                    ("S_B", "D_B")]
        elif source == "gCSI":
            data = [("S_B", "D_B"), ("gCSI_ONLY", "D_0"),
                    ("S_A", "D_A")]
        else:
            data = [(source + "_A", "D_A"),
                    (source + "_V", "D_0"),
                    (source + "_T", "D_B")]
        for sample, drug in data:
            tsv.append(f"{source}\t{sample}\t{drug}\t{source}\tSECRET_OUTCOME\n")
        parts[BASE + f"splits/{source}_all.txt"] = (
            f"{pos}\n{pos+1}\n{pos+2}\n")
        parts[BASE + f"splits/{source}_split_0_train.txt"] = f"{pos}\n"
        parts[BASE + f"splits/{source}_split_0_val.txt"] = f"{pos+1}\n"
        parts[BASE + f"splits/{source}_split_0_test.txt"] = f"{pos+2}\n"
    parts[BASE + "y_data/response.tsv"] = "".join(tsv)
    if alter:
        alter(parts)
    with ZipFile(path, "w") as z:
        for key, value in parts.items():
            z.writestr(key, value)


class ImproveFiveStudyOverlapTests(unittest.TestCase):
    def test_all_ten_study_pairs_present(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path)
            r = audit_five_study_overlap(path)
            self.assertEqual(r["response_total_rows"], 15)
            self.assertEqual(len(r["studies"]), 5)
            self.assertEqual(r["study_pair_count"], 10)
            self.assertEqual(r["studies"]["CCLE"]["unique_model_drug_pairs"], 3)

    def test_cross_study_train_and_test_pair_repetition(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path)
            r = audit_five_study_overlap(path)
            overlap = r["study_pair_overlaps"]["CCLE_vs_gCSI"]
            self.assertEqual(overlap["all_vs_all"]["sample_drug_pairs"], 2)
            self.assertEqual(overlap["train_vs_foreign_test"]["sample_drug_pairs"], 1)
            self.assertEqual(overlap["foreign_train_vs_test"]["sample_drug_pairs"], 1)
            self.assertEqual(overlap["shared_raw_response_row_indices"], 0)
            self.assertEqual(
                r["pooled_partition_risks"]["cross_source_train_to_other_source_heldout_pairs"], 2)

    def test_row_disjoint_does_not_imply_cross_source_pair_disjoint(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path)
            r = audit_five_study_overlap(path)
            self.assertEqual(r["studies"]["CCLE"]["split_row_index_overlap"]["train_test"], 0)
            self.assertEqual(r["pooled_partition_risks"]["train_vs_test_shared_pairs"], 2)

    def test_within_source_repeated_pair_is_reported_not_removed(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            def alter(parts):
                part = BASE + "y_data/response.tsv"
                parts[part] = parts[part].replace(
                    "CTRPv2\tCTRPv2_T\tD_B", "CTRPv2\tCTRPv2_A\tD_A")
            fictitious_archive(path, alter)
            r = audit_five_study_overlap(path)
            self.assertEqual(r["studies"]["CTRPv2"]["repeated_assay_pair_excess_rows"], 1)
            self.assertEqual(r["studies"]["CTRPv2"]["pairs_with_multiple_rows"], 1)
            self.assertEqual(r["studies"]["CTRPv2"]["within_source_pair_intersections"]
                             ["train_test"]["sample_drug_pairs"], 1)

    def test_bad_index_is_rejected_before_overlap_calculation(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path, lambda parts: parts.update({
                BASE + "splits/CCLE_split_0_test.txt": "999\n"}))
            with self.assertRaisesRegex(ValueError, "index outside response"):
                audit_five_study_overlap(path)

    def test_duplicate_split_membership_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path, lambda parts: parts.update({
                BASE + "splits/CCLE_split_0_train.txt": "0\n0\n"}))
            with self.assertRaisesRegex(ValueError, "duplicate index"):
                audit_five_study_overlap(path)

    def test_no_names_or_labels_are_published(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path)
            text = str(audit_five_study_overlap(path))
            for word in ("SECRET_OUTCOME", "CCLE_ONLY", "gCSI_ONLY", "S_A", "S_B"):
                self.assertNotIn(word, text)

    def test_missing_manifest_and_invalid_fold_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fictitious_archive(path)
            with self.assertRaisesRegex(ValueError, "missing required source manifest"):
                audit_five_study_overlap(path, fold=1)
            with self.assertRaisesRegex(ValueError, "fold must"):
                audit_five_study_overlap(path, fold=10)


if __name__ == "__main__":
    unittest.main()
