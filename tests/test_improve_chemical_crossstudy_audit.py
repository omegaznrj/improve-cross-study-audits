"""Fictional chemical aliases and cross-study biological-pair overlap tests."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from research.oncology.improve_chemical_crossstudy_audit import (
    audit_chemical_and_cross_study,
)
from research.oncology.improve_split_audit import BASE


def fixture(path, alter=None):
    data = {
        BASE + "y_data/response.tsv": (
            "source\timprove_sample_id\timprove_chem_id\tstudy\tauc\n"
            "CCLE\tS1\tD1\tCCLE\tHIDDEN\n"
            "CCLE\tS2\tD2\tCCLE\tHIDDEN\n"
            "CCLE\tS3\tD3\tCCLE\tHIDDEN\n"
            "gCSI\tS1\tD2\tgCSI\tHIDDEN\n"
            "gCSI\tS1\tD1\tgCSI\tHIDDEN\n"
            "gCSI\tS3\tD3\tgCSI\tHIDDEN\n"
        ),
        BASE + "splits/CCLE_split_0_train.txt": "0\n",
        BASE + "splits/CCLE_split_0_val.txt": "1\n",
        BASE + "splits/CCLE_split_0_test.txt": "2\n",
        BASE + "splits/gCSI_split_0_train.txt": "3\n",
        BASE + "splits/gCSI_split_0_val.txt": "4\n",
        BASE + "splits/gCSI_split_0_test.txt": "5\n",
        BASE + "x_data/drug_SMILES.tsv": (
            "improve_chem_id\tcanSMILES\n"
            "D1\tCCO\nD2\tN\nD3\tCC\n"
        ),
        BASE + "x_data/drug_info.tsv": (
            "DrugID\tSMILES\tcanSMILES\tNAME\tPUBCHEM_ID\tCAS_ID\timprove_chem_id\n"
            "a\tCCO\tCCO\tA\t21\t10\tD1\n"
            "b\tCCO\tC(C)O\tA alias\t21\t10\tD1\n"
            "c\tO\tO\tB\t22\t20\tD2\n"
            "d\tO\tO\tB alias\t23\t20\tD2\n"
        ),
        BASE + "x_data/drug_ecfp4_nbits512.tsv": (
            "improve_chem_id\t" +
            "\t".join(f"ecfp4.{i}" for i in range(512)) + "\n" +
            "D1\t" + "\t".join(["1"] + ["0"] * 511) + "\n" +
            "D2\t" + "\t".join(["1"] + ["0"] * 511) + "\n" +
            "D3\t" + "\t".join(["0"] + ["1"] + ["0"] * 510) + "\n"
        ),
    }
    if alter:
        alter(data)
    with ZipFile(path, "w") as archive:
        for name, content in data.items():
            archive.writestr(name, content)


class ImproveChemicalCrossStudyAuditTests(unittest.TestCase):
    def test_cross_source_shared_pairs_reported(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fixture(path)
            result = audit_chemical_and_cross_study(path)
            overlaps = result["cross_study_overlap"]
            self.assertEqual(overlaps["ccle_all_vs_gcsi_all"]["pairs"], 2)
            self.assertEqual(overlaps["ccle_all_vs_gcsi_all"]["drugs"], 3)
            self.assertEqual(overlaps["ccle_train_vs_gcsi_val"]["pairs"], 1)
            self.assertEqual(overlaps["ccle_train_vs_gcsi_test"]["pairs"], 0)
            self.assertEqual(overlaps["gcsi_train_vs_ccle_test"]["pairs"], 0)

    def test_canonical_reference_and_metadata_mismatches_counted(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fixture(path)
            report = audit_chemical_and_cross_study(path)
            chemical = report["chemical_text_and_fingerprint_diagnostics"]
            self.assertEqual(chemical["distinct_training_chemical_ids"], 2)
            self.assertEqual(chemical["drug_info_comparison"]["ids_with_repeated_metadata_rows"], 2)
            self.assertEqual(chemical["drug_info_comparison"]["ids_with_multiple_distinct_nonblank_alias_smiles"], 1)
            self.assertEqual(chemical["drug_info_comparison"]["ids_with_any_exact_alias_reference_match"], 1)
            self.assertEqual(chemical["drug_info_comparison"]["ids_with_alias_but_no_exact_reference_text_match"], 1)
            self.assertEqual(chemical["drug_info_comparison"]["ids_with_multiple_nonblank_pubchem_ids"], 1)

    def test_identical_ecfp_across_distinct_ids_not_silently_treated_unique(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fixture(path)
            chemical = audit_chemical_and_cross_study(path)[
                "chemical_text_and_fingerprint_diagnostics"]["fingerprints"]
            self.assertEqual(chemical["identical_512bit_fingerprint_groups_across_distinct_ids"], 1)
            self.assertEqual(chemical["distinct_ids_in_identical_fingerprint_groups"], 2)
            self.assertEqual(chemical["duplicate_key_rows"], 0)

    def test_duplicate_primary_smiles_rows_reported(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            def duplicate(data):
                data[BASE + "x_data/drug_SMILES.tsv"] += "D1\tCOC\n"
            fixture(path, duplicate)
            item = audit_chemical_and_cross_study(path)[
                "chemical_text_and_fingerprint_diagnostics"]["drug_smiles"]
            self.assertEqual(item["duplicate_reference_id_rows"], 1)
            self.assertEqual(item["ids_with_more_than_one_reference_smiles"], 1)

    def test_no_raw_drug_names_smiles_or_heldout_labels_emitted(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            fixture(path)
            report = str(audit_chemical_and_cross_study(path))
            for forbidden in ("HIDDEN", "C(C)O", "'D1'", "'S1'", "B alias"):
                self.assertNotIn(forbidden, report)
            self.assertTrue(audit_chemical_and_cross_study(path)["no_response_target_values_used"])

    def test_missing_reference_smiles_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fictional.zip"
            def remove(data):
                del data[BASE + "x_data/drug_SMILES.tsv"]
            fixture(path, remove)
            with self.assertRaisesRegex(ValueError, "missing chemical table"):
                audit_chemical_and_cross_study(path)

    def test_fingerprint_bit_count_must_be_exactly_512(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "short.zip"
            def make_short(data):
                data[BASE + "x_data/drug_ecfp4_nbits512.tsv"] = (
                    "improve_chem_id\tecfp4.0\nD1\t1\nD2\t1\n"
                )
            fixture(path, make_short)
            with self.assertRaisesRegex(ValueError, "512 ordered ECFP4 bits"):
                audit_chemical_and_cross_study(path)

    def test_bad_fold_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "fold must"):
                audit_chemical_and_cross_study(Path(d) / "none.zip", fold=10)


if __name__ == "__main__":
    unittest.main()
