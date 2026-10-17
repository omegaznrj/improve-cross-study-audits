"""Synthetic fixture for all-response-linked chemical text audit."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from research.oncology.improve_split_audit import BASE
from research.oncology.improve_full_chemical_identity_audit import audit_full_chemical_identity


def fixture(path):
    prefix = BASE + "x_data/"
    fp_header = "improve_chem_id\t" + "\t".join(f"ecfp4.{i}" for i in range(512)) + "\n"
    fp_one = "C1\t" + "\t".join("0" for _ in range(512)) + "\n"
    fp_two = "C2\t" + "\t".join("1" for _ in range(512)) + "\n"
    values = {
        BASE + "y_data/response.tsv":
            "source\timprove_sample_id\timprove_chem_id\tstudy\tauc\n"
            "CCLE\tS1\tC1\tB1\tSECRET_OUTCOME\n"
            "gCSI\tS2\tC2\tB2\tSECRET_OUTCOME\n",
        prefix + "drug_SMILES.tsv":
            "improve_chem_id\tcanSMILES\nC1\tC\nC2\tCC\n",
        prefix + "drug_info.tsv":
            "improve_chem_id\tcanSMILES\tPUBCHEM_ID\nC1\tC\t1\n"
            "C2\tCC\t2\nC2\tC-C\t2\n",
        prefix + "drug_ecfp4_nbits512.tsv": fp_header + fp_one + fp_two,
    }
    with ZipFile(path, "w") as z:
        for name, data in values.items():
            z.writestr(name, data)


class FullChemicalTests(unittest.TestCase):
    def test_all_drugs_and_flags_without_labels(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "mock.zip"
            fixture(p)
            r = audit_full_chemical_identity(p)
            self.assertEqual(r["response_linked_drug_ids"], 2)
            chem = r["chemical_diagnostics"]
            self.assertEqual(chem["drug_smiles"]["ids_with_rows"], 2)
            self.assertEqual(chem["drug_info_comparison"]["ids_with_repeated_metadata_rows"], 1)
            self.assertEqual(chem["drug_info_comparison"]["ids_with_multiple_distinct_nonblank_alias_smiles"], 1)
            self.assertEqual(chem["fingerprints"]["ids_with_rows"], 2)
            self.assertNotIn("SECRET_OUTCOME", str(r))
            self.assertFalse(r["structural_chemistry_adjudicated"])

    def test_missing_archive_and_duplicate_zip_members(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "mock.zip"
            with self.assertRaisesRegex(ValueError, "archive not found"):
                audit_full_chemical_identity(p)
            fixture(p)
            with ZipFile(p, "a") as z:
                z.writestr(BASE + "y_data/response.tsv", "duplicate")
            with self.assertRaisesRegex(ValueError, "duplicate archive"):
                audit_full_chemical_identity(p)


if __name__ == "__main__":
    unittest.main()
