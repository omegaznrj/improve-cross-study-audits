"""Unit checks for isolated Mordred report integrity; no external data."""
import unittest
from research.oncology.improve_isolated_mordred_reproduction import TARGETS, verify_values


class IsolatedMordredTests(unittest.TestCase):
    def test_complete_numeric_descriptor_report(self):
        verify_values({"synthetic": {"descriptors": {
            name: {"status": "numeric", "value": 1.0} for name in TARGETS
        }}})

    def test_missing_descriptor_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Missing"):
            verify_values({"synthetic": {"descriptors": {"Radius": {"status": "numeric", "value": 1.0}}}})

    def test_nonfinite_numeric_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            verify_values({"synthetic": {"descriptors": {
                name: {"status": "numeric", "value": float("inf") if name == "Diameter" else 1.0}
                for name in TARGETS
            }}})


if __name__ == "__main__":
    unittest.main()
