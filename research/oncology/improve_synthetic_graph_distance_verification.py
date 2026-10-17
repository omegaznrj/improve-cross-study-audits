"""Phase 1.5F: controlled disconnected-graph distance and Mordred verification.

No IMPROVE source data, chemical IDs, or clinical outcomes are read. Synthetic
molecules only. RDKit distance matrices are always tested; installed Mordred
is optional and its actual outputs are labelled independently.
"""
from __future__ import annotations
import json
import math
import importlib
from importlib import metadata

MOLECULES = {
    "connected_ethanol": "CCO",
    "disconnected_ethanol_sodium": "CCO.[Na+]",
    "disconnected_ethanol_water": "CCO.O",
    "connected_butane": "CCCC",
}


def run():
    from rdkit import Chem
    from rdkit.Chem import rdmolops
    from rdkit import rdBase
    try:
        from mordred import Calculator
        from mordred import TopologicalIndex, EccentricConnectivityIndex
        mordred_status = metadata.version("mordred")
    except (ImportError, metadata.PackageNotFoundError):
        Calculator = None
        mordred_status = "not_installed"
    result = {}
    for label, smi in MOLECULES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise RuntimeError("synthetic SMILES not parseable: " + label)
        distances = rdmolops.GetDistanceMatrix(mol)
        flat = [float(v) for row in distances for v in row]
        connected = len(Chem.GetMolFrags(mol)) == 1
        entry = {
            "synthetic_smiles": smi,
            "fragments": len(Chem.GetMolFrags(mol)),
            "rdkit_distance_matrix_max": float(max(flat)),
            "rdkit_distance_matrix_has_1e8": any(v == 1e8 for v in flat),
            "rdkit_off_diagonal_1e8_pairs": int(sum(
                distances[i,j] == 1e8 for i in range(mol.GetNumAtoms())
                for j in range(i+1,mol.GetNumAtoms())
            )),
            "mordred": {"status": mordred_status},
        }
        if Calculator is not None:
            try:
                # Keep names dynamic: pinned generators may offer different APIs.
                from mordred import descriptors
                wanted = {"Diameter","Radius","WPath","ECIndex"}
                chosen = [d for d in Calculator(descriptors, ignore_3D=True).descriptors
                          if str(d) in wanted]
                vals = Calculator(chosen, ignore_3D=True)(mol)
                entry["mordred"] = {
                    "status": "evaluated", "package_version": mordred_status,
                    "descriptor_values": {
                        str(d): (float(v) if isinstance(v,(int,float)) and math.isfinite(float(v))
                                 else str(v))
                        for d,v in zip(chosen,vals)
                    },
                }
            except Exception as exc:
                entry["mordred"] = {"status": "evaluation_error",
                                    "package_version": mordred_status,
                                    "error_type": type(exc).__name__}
        result[label] = entry
    return {
        "schema": "tumorlab-synthetic-graph-distance-verification-v1",
        "rdkit_version": rdBase.rdkitVersion,
        "cases": result,
        "reads_improve_dataset": False,
        "training_performed": False,
        "interpretation_limits": [
            "RDKit distance-matrix behavior does not prove IMPROVE's original descriptor generator and preprocessing settings.",
            "Optional local Mordred installation/version may differ from the one used to generate IMPROVE.",
            "The source values still require a documented generation provenance comparison.",
        ],
    }


def main():
    print(json.dumps(run(),indent=2))


if __name__ == "__main__":
    main()
