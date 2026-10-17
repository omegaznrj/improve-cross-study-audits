"""Isolated Mordred replication on synthetic connected/disconnected molecules.

Run from a separate venv containing mordred==1.2.0 and RDKit. This script
fails rather than silently treating an unsuccessful import as 'not installed'.
It does not read IMPROVE, source records, outcomes, or private identifiers.
"""
from __future__ import annotations

import importlib.metadata
import json
import math

from research.oncology.improve_synthetic_graph_distance_verification import MOLECULES

TARGETS = ("Diameter", "Radius", "WPath", "ECIndex")


def verify_values(case_results):
    """Assert reporting integrity, not a presupposed scientific conclusion."""
    for label, case in case_results.items():
        if set(case["descriptors"]) != set(TARGETS):
            raise ValueError("Missing or duplicate Mordred targets: " + label)
        for name, value in case["descriptors"].items():
            if value["status"] == "numeric" and not math.isfinite(value["value"]):
                raise ValueError("Nonfinite numeric result")


def run():
    from rdkit import Chem, rdBase
    from mordred import Calculator, descriptors

    catalog = Calculator(descriptors, ignore_3D=True).descriptors
    selected = {str(d): d for d in catalog if str(d) in TARGETS}
    if set(selected) != set(TARGETS):
        raise RuntimeError("Mordred targets unavailable: " + str(sorted(set(TARGETS)-set(selected))))
    calc = Calculator([selected[name] for name in TARGETS], ignore_3D=True)
    cases = {}
    for label, smiles in MOLECULES.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RuntimeError("Cannot parse synthetic SMILES: " + label)
        result = {}
        for name, value in zip(TARGETS, calc(mol)):
            try:
                number = float(value)
                if math.isfinite(number):
                    result[name] = {"status": "numeric", "value": number}
                else:
                    result[name] = {"status": "nonfinite", "representation": str(value)}
            except (TypeError, ValueError, OverflowError):
                result[name] = {"status": "not_numeric", "representation": str(value)}
        cases[label] = {"smiles": smiles, "fragments": len(Chem.GetMolFrags(mol)),
                        "descriptors": result}
    verify_values(cases)
    return {
        "schema": "tumorlab-isolated-mordred-reproduction-v1",
        "rdkit_version": rdBase.rdkitVersion,
        "mordred_version": importlib.metadata.version("mordred"),
        "cases": cases,
        "reads_improve_dataset": False, "training_performed": False,
        "caveat": "Reproduces local package calculations only, not the undocumented IMPROVE generation pipeline.",
    }


def main():
    print(json.dumps(run(), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
