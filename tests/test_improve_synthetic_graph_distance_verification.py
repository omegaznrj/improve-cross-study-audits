"""Synthetic RDKit graph-distance tests. Mordred installation optional."""
import unittest
import json
from research.oncology.improve_synthetic_graph_distance_verification import run


class SyntheticGraphDistanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run()["cases"]

    def test_connected_graph_has_no_sentinel(self):
        d = self.result["connected_ethanol"]
        self.assertEqual(d["fragments"], 1)
        self.assertFalse(d["rdkit_distance_matrix_has_1e8"])

    def test_two_disconnected_graphs_have_large_distance(self):
        for key in ("disconnected_ethanol_sodium","disconnected_ethanol_water"):
            d = self.result[key]
            self.assertEqual(d["fragments"],2)
            self.assertTrue(d["rdkit_distance_matrix_has_1e8"])
            self.assertGreater(d["rdkit_off_diagonal_1e8_pairs"],0)

    def test_entire_report_is_json_serializable(self):
        report = run()
        self.assertIsInstance(json.dumps(report), str)

    def test_optional_mordred_results_are_explicit(self):
        for case in self.result.values():
            self.assertIn(case["mordred"]["status"],
                          {"not_installed","evaluated","evaluation_error"})


if __name__=="__main__":
    unittest.main()
