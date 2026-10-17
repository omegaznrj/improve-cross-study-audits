"""Phase 1.5E: test whether large Mordred topological values track disconnected SMILES.

Read-only, aggregate-only evidence. Assesses correspondence, not causation or
descriptor recomputation. No response outcomes or compound identifiers emitted.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
from pathlib import Path
import zipfile
from collections import Counter

from research.oncology.improve_split_audit import BASE, _read_response
from research.oncology.improve_structure_equivalence_audit import _rdkit
from research.oncology.improve_numeric_quality_audit import _value_kind

FEATURES = ("mordred.VR1_A", "mordred.VR2_A", "mordred.ECIndex",
            "mordred.Diameter", "mordred.Radius", "mordred.WPath")


def _records(z, filename, needed):
    with z.open(BASE + "x_data/" + filename) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t", strict=True)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(filename + ": invalid header")
        if not set(needed).issubset(reader.fieldnames):
            raise ValueError(filename + ": missing columns")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(filename + ": malformed row")
            yield row


def association(z, linked, Chem):
    structure = {}
    for r in _records(z, "drug_SMILES.tsv", ("improve_chem_id", "canSMILES")):
        k = r["improve_chem_id"].strip()
        if k not in linked:
            continue
        if k in structure:
            raise ValueError("duplicate linked reference SMILES")
        mol = Chem.MolFromSmiles(r["canSMILES"].strip())
        structure[k] = None if mol is None else len(Chem.GetMolFrags(mol)) > 1
    counts = {name: Counter() for name in FEATURES}
    totals = Counter()
    matched = set()
    for r in _records(z, "drug_mordred.tsv", ("improve_chem_id",) + FEATURES):
        k = r["improve_chem_id"].strip()
        if k not in linked:
            continue
        if k in matched:
            raise ValueError("duplicate linked Mordred key")
        matched.add(k)
        disconnected = structure.get(k)
        status = "unparseable_or_unmatched" if disconnected is None else (
            "multiple_fragments" if disconnected else "one_fragment")
        totals[status] += 1
        for name in FEATURES:
            kind, value = _value_kind(r[name])
            if kind != "valid":
                counts[name][status + "_invalid"] += 1
            else:
                counts[name][status + "_valid"] += 1
                if abs(value) >= 1e6:
                    counts[name][status + "_abs_ge_1e6"] += 1
                if abs(value) >= 1e9:
                    counts[name][status + "_abs_ge_1e9"] += 1
                if value == 100000000:
                    counts[name][status + "_exact_1e8"] += 1
    return {
        "linked_drug_counts_by_reference_fragment_status": dict(sorted(totals.items())),
        "descriptor_counts_by_reference_fragment_status": {
            k: dict(sorted(v.items())) for k, v in counts.items()
        },
        "notes": [
            "Associations are not proof of a specific descriptor-generation rule.",
            "Reference SMILES are not verified physical preparations.",
            "Descriptor generator/version and source preprocessing remain unverified.",
            "No compound records or response outcomes exported.",
        ],
    }


def audit(path):
    Chem, rdkit_version = _rdkit()
    with zipfile.ZipFile(path) as z:
        members = z.namelist()
        if len(members) != len(set(members)):
            raise ValueError("duplicate ZIP members")
        linked = set(_read_response(z)["drug_code_by_id"])
        return {"schema": "improve-mordred-fragment-association-v1",
                "rdkit_version": rdkit_version,
                **association(z, linked, Chem),
                "training_performed": False, "outcomes_read": False}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[2] / "data/raw/improve_cross_study/csa_data.zip")
    args=p.parse_args(argv)
    try:
        print(json.dumps(audit(args.archive),indent=2))
    except (OSError, ValueError, zipfile.BadZipFile, csv.Error, RuntimeError) as exc:
        p.exit(2, "IMPROVE fragment association: " + str(exc) + "\n")


if __name__ == "__main__":
    main()
