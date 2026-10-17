"""Synthetic-only ZIP tests for zero-based IMPROVE fold and key overlap audits."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from research.oncology.improve_split_audit import (
    BASE, STUDIES, audit_split_indices,
)


def create_zip(path, *, modify=None):
    response = ["source\timprove_sample_id\timprove_chem_id\tstudy\tauc\n"]
    body = {}
    for study in STUDIES:
        beginning = len(response) - 1
        response += [
            f"{study}\tSHARED\tD1\t{study}\t0.1\n",
            f"{study}\tS2\tD2\t{study}\t0.2\n",
            f"{study}\tSHARED\tD1\t{study}\t0.3\n",
            f"{study}\tS4\tD4\t{study}\t0.4\n",
        ]
        body[BASE + f"splits/{study}_all.txt"] = "".join(
            f"{i}\n" for i in range(beginning, beginning + 4)
        )
        body[BASE + f"splits/{study}_split_0_train.txt"] = (
            f"{beginning}\n{beginning + 1}\n"
        )
        body[BASE + f"splits/{study}_split_0_val.txt"] = f"{beginning + 3}\n"
        body[BASE + f"splits/{study}_split_0_test.txt"] = f"{beginning + 2}\n"
    body[BASE + "y_data/response.tsv"] = "".join(response)
    if modify:
        modify(body)
    with ZipFile(path, "w") as archive:
        for name, content in body.items():
            archive.writestr(name, content)


class SplitAuditTests(unittest.TestCase):
    def test_five_complete_studies_and_entity_overlap(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path)
            report = audit_split_indices(path)
            self.assertEqual(report["response_row_count"], 20)
            self.assertEqual(report["rows_missing_sample_or_drug_key"], 0)
            self.assertEqual(len(report["split_audit"]), 5)
            item = report["split_audit"]["CCLE"]
            self.assertEqual(
                (item["all_rows"], item["train_rows"],
                 item["validation_rows"], item["test_rows"]), (4, 2, 1, 1)
            )
            self.assertEqual(item["row_overlap_train_test"], 0)
            self.assertEqual(
                item["train_test_entity_overlap"]["overlapping_sample_drug_pairs"], 1
            )
            self.assertEqual(
                item["train_test_entity_overlap"]["overlapping_sample_ids"], 1
            )
            self.assertEqual(item["study_all_indices_with_different_source_label"], 0)

    def test_fold_not_available_is_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path)
            with self.assertRaisesRegex(ValueError, "missing response/split file"):
                audit_split_indices(path, fold=1)

    def test_illegal_fold(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "fold must"):
                audit_split_indices(Path(folder) / "fictional.zip", fold=10)

    def test_duplicate_index(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path, modify=lambda payload: payload.update({
                BASE + "splits/CCLE_split_0_test.txt": "2\n2\n"
            }))
            with self.assertRaisesRegex(ValueError, "duplicate index"):
                audit_split_indices(path)

    def test_out_of_range(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path, modify=lambda payload: payload.update({
                BASE + "splits/CCLE_split_0_test.txt": "20\n"
            }))
            with self.assertRaisesRegex(ValueError, "index outside response"):
                audit_split_indices(path)

    def test_cross_partition_row_overlap_detected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path, modify=lambda payload: payload.update({
                BASE + "splits/CCLE_split_0_test.txt": "0\n2\n"
            }))
            report = audit_split_indices(path)
            self.assertEqual(report["split_audit"]["CCLE"]["row_overlap_train_test"], 1)

    def test_outside_original_study(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path, modify=lambda payload: payload.update({
                BASE + "splits/CCLE_split_0_test.txt": "2\n4\n"
            }))
            report = audit_split_indices(path)
            self.assertEqual(
                report["split_audit"]["CCLE"]["partition_indices_outside_study_all"], 1
            )

    def test_no_response_values_in_report(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fictional.zip"
            create_zip(path)
            report = audit_split_indices(path)
            self.assertNotIn("0.1", str(report))
            self.assertNotIn("SHARED", str(report))


if __name__ == "__main__":
    unittest.main()
