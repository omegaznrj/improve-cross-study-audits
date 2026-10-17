"""Synthetic tests for the aggregate-only Mordred descriptor investigation."""
import io
import unittest
import zipfile
from research.oncology.improve_split_audit import BASE
from research.oncology.improve_mordred_descriptor_extremes import extreme_descriptors


class MordredDescriptorTests(unittest.TestCase):
    def make_zip(self, rows):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as z:
            z.writestr(BASE+"x_data/drug_mordred.tsv",rows)
        return zipfile.ZipFile(io.BytesIO(buf.getvalue()))

    def test_only_extreme_descriptors_and_aggregate_ranges(self):
        with self.make_zip("improve_chem_id\tA\tB\tC\nD1\t1000001\t-3\tNA\nD2\t2000000\t4\t3\nOTHER\t-1e12\t0\t9\n") as z:
            r=extreme_descriptors(z,{"D1","D2"})
        self.assertEqual(r["linked_drugs"],2)
        self.assertEqual([x["descriptor"] for x in r["extreme_descriptor_columns"]],["A"])
        self.assertEqual(r["extreme_descriptor_columns"][0]["minimum"],1000001)
        self.assertEqual(r["extreme_descriptor_columns"][0]["values_at_or_above_threshold"],2)
        self.assertNotIn("D1",str(r))

    def test_duplicate_linked_key_rejected(self):
        with self.make_zip("improve_chem_id\tA\nD1\t3\nD1\t4\n") as z:
            with self.assertRaisesRegex(ValueError,"duplicate"):
                extreme_descriptors(z,{"D1"})

    def test_invalid_threshold(self):
        with self.make_zip("improve_chem_id\tA\nD1\t3\n") as z:
            with self.assertRaises(ValueError):
                extreme_descriptors(z,{"D1"},0)


if __name__=="__main__":
    unittest.main()
