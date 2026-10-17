"""Fictional data only: feature-key audit checks exact joins without values."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from research.oncology.improve_feature_key_audit import audit_feature_keys
from research.oncology.improve_split_audit import BASE, STUDIES


def make_archive(path, change=None):
    data = {}
    response = ["source\timprove_sample_id\timprove_chem_id\tstudy\tauc\n"]
    for study in STUDIES:
        i = len(response) - 1
        response.extend([
            f"{study}\tS{i}\tD{i}\t{study}\t0.12345\n",
            f"{study}\tS{i+1}\tD{i+1}\t{study}\t0.87654\n",
            f"{study}\tS{i+2}\tD{i+2}\t{study}\t0.33333\n",
        ])
        data[BASE + f"splits/{study}_split_0_train.txt"] = f"{i}\n"
        data[BASE + f"splits/{study}_split_0_val.txt"] = f"{i+1}\n"
        data[BASE + f"splits/{study}_split_0_test.txt"] = f"{i+2}\n"
    data[BASE + "y_data/response.tsv"] = "".join(response)
    data[BASE + "x_data/cancer_gene_expression.tsv"] = (
        "\tENSG0000003\tENSG0000005\n"
        "S0\t9999999\t8888888\n"
        "S2\t7777777\t6666666\n"
    )
    data[BASE + "x_data/drug_SMILES.tsv"] = (
        "improve_chem_id\tcanSMILES\nD0\tC(C)O\nD1\tNCCO\n"
    )
    data[BASE + "x_data/drug_info.tsv"] = (
        "DrugID\tSMILES\tcanSMILES\tNAME\tPUBCHEM_ID\tCAS_ID\timprove_chem_id\n"
        "local01\tOC\tCO\tFictional\t100\t99\tD2\n"
    )
    data[BASE + "x_data/drug_ecfp4_nbits512.tsv"] = (
        "improve_chem_id\tecfp4.0\nD0\t1\n"
    )
    if change:
        change(data)
    with ZipFile(path, "w") as zipfile:
        for name, body in data.items():
            zipfile.writestr(name, body)


class ImproveFeatureKeyAuditTests(unittest.TestCase):
    def test_counts_and_partition_specific_coverage(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            make_archive(p)
            report = audit_feature_keys(p)
            gene = report["feature_tables"]["cancer_gene_expression.tsv"]
            self.assertEqual(gene["unique_matching_response_ids"], 2)
            self.assertEqual(gene["per_study"]["CCLE"]["train"]["exact_key_matched_rows"], 1)
            self.assertEqual(gene["per_study"]["CCLE"]["val"]["exact_key_matched_rows"], 0)
            self.assertEqual(gene["per_study"]["CCLE"]["test"]["exact_key_matched_rows"], 1)

    def test_drug_info_uses_named_seventh_key_not_drugid(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            make_archive(p)
            info = audit_feature_keys(p)["feature_tables"]["drug_info.tsv"]
            self.assertEqual(info["key_column_position"], 6)
            self.assertEqual(info["unique_matching_response_ids"], 1)
            self.assertEqual(info["per_study"]["CCLE"]["test"]["exact_key_matched_rows"], 1)

    def test_response_values_and_feature_values_not_in_report(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            make_archive(p)
            result = str(audit_feature_keys(p))
            for forbidden in ("9999999", "8888888", "0.12345", "S0", "D0", "C(C)O"):
                self.assertNotIn(forbidden, result)
            self.assertIn("cancer_mutation.parquet", result)

    def test_missing_features_and_unmatched_keys_are_reported(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            def change(data):
                data[BASE + "x_data/cancer_gene_expression.tsv"] = (
                    "\tENSG1\nUNKNOWN_SAMPLE\t1\n"
                )
            make_archive(p, change)
            report = audit_feature_keys(p)
            gene = report["feature_tables"]["cancer_gene_expression.tsv"]
            self.assertEqual(gene["unique_matching_response_ids"], 0)
            self.assertEqual(gene["join_status"], "no_exact_key_matches_needs_crosswalk")
            self.assertEqual(report["feature_tables"]["cancer_RPPA.tsv"]["status"],
                             "missing_file")

    def test_duplicate_and_empty_keys_are_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            def change(data):
                data[BASE + "x_data/cancer_gene_expression.tsv"] = (
                    "\tG1\nS0\t1\nS0\t2\n\t0\n"
                )
            make_archive(p, change)
            gene = audit_feature_keys(p)["feature_tables"]["cancer_gene_expression.tsv"]
            self.assertEqual(gene["duplicate_key_rows"], 1)
            self.assertEqual(gene["empty_id_rows"], 1)
            self.assertEqual(gene["unique_nonblank_ids"], 1)

    def test_header_mismatch_rejected_for_drug_fingerprints(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "f.zip"
            def change(data):
                data[BASE + "x_data/drug_ecfp4_nbits512.tsv"] = "WRONG\tFeature\nD0\t1\n"
            make_archive(p, change)
            with self.assertRaisesRegex(ValueError, "unexpected drug feature key header"):
                audit_feature_keys(p)

    def test_invalid_fold_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "fold must"):
                audit_feature_keys(Path(d) / "absent.zip", fold=10)

    def test_missing_archive_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "ZIP not found"):
                audit_feature_keys(Path(d) / "absent.zip")


if __name__ == "__main__":
    unittest.main()
