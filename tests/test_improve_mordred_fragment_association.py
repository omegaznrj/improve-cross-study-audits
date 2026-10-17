"""Synthetic regression tests for fragment association; no source records."""
import io
import unittest
import zipfile

from research.oncology.improve_split_audit import BASE
from research.oncology.improve_structure_equivalence_audit import _rdkit
from research.oncology.improve_mordred_fragment_association import association, FEATURES


class FragmentAssociationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.Chem, _ = _rdkit()
        except RuntimeError as err:
            raise unittest.SkipTest(str(err))

    def fixture(self, descriptors):
        data=io.BytesIO()
        with zipfile.ZipFile(data,"w") as z:
            z.writestr(BASE+"x_data/drug_SMILES.tsv",
                       "improve_chem_id\tcanSMILES\nD1\tCCO\nD2\tCCO.[Na+]\n")
            z.writestr(BASE+"x_data/drug_mordred.tsv",
                       "improve_chem_id\t"+"\t".join(FEATURES)+"\n"+descriptors)
        return zipfile.ZipFile(io.BytesIO(data.getvalue()))

    def test_disconnected_association(self):
        with self.fixture("D1\t"+"\t".join(["1"]*6)+"\nD2\t"+"\t".join(["100000000"]*6)+"\n") as z:
            r=association(z,{"D1","D2"},self.Chem)
        self.assertEqual(r["linked_drug_counts_by_reference_fragment_status"],
                         {"multiple_fragments":1,"one_fragment":1})
        self.assertEqual(r["descriptor_counts_by_reference_fragment_status"]["mordred.Diameter"]
                         ["multiple_fragments_exact_1e8"],1)

    def test_duplicate_key_rejected(self):
        with self.fixture("D1\t"+"\t".join(["1"]*6)+"\nD1\t"+"\t".join(["1"]*6)+"\n") as z:
            with self.assertRaisesRegex(ValueError,"duplicate"):
                association(z,{"D1"},self.Chem)


if __name__=="__main__":
    unittest.main()
