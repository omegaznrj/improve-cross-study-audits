"""Read-only, structure-aware IMPROVE SMILES concordance audit.

RDKit required. No merging, standardization in source, or response values.
Canonical isomeric SMILES equality is a computable representation-level result,
not proof that different formulations or assay preparations are identical.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import zipfile

from research.oncology.improve_archive_audit import DEFAULT_ARCHIVE
from research.oncology.improve_chemical_crossstudy_audit import _read_tsv
from research.oncology.improve_split_audit import _read_response


def _rdkit():
    try:
        from rdkit import Chem, rdBase
    except ImportError as exc:
        raise RuntimeError(
            "RDKit is required; install with .venv\\Scripts\\python.exe -m pip install rdkit"
        ) from exc
    return Chem, rdBase.rdkitVersion


def _signature(text, Chem):
    if not text:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    return (
        Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        len(Chem.GetMolFrags(mol)),
    )


def classify(reference, alias, Chem):
    """Exclusively describe structural-string relations; do not adjudicate formulations."""
    a, b = _signature(reference, Chem), _signature(alias, Chem)
    if a is None or b is None:
        return "invalid_or_blank_structure"
    if a[0] == b[0]:
        return "same_canonical_isomeric_structure"
    if a[1] == b[1]:
        return "same_connectivity_stereo_difference_or_unspecified"
    if a[2] != b[2]:
        return "different_fragment_counts_review_salt_or_mixture"
    return "different_connectivity_or_other_structure_difference"


def audit_structure_equivalence(path=DEFAULT_ARCHIVE):
    Chem, version = _rdkit()
    path = Path(path)
    if not path.is_file():
        raise ValueError("IMPROVE archive missing")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member names")
        rows = _read_response(z)
        drugs = set(rows["drug_code_by_id"]) - {""}
        refs = defaultdict(set)
        aliases = defaultdict(set)
        for row in _read_tsv(z, "drug_SMILES.tsv", ("improve_chem_id", "canSMILES")):
            key = row["improve_chem_id"].strip()
            if key in drugs:
                refs[key].add(row["canSMILES"].strip())
        for row in _read_tsv(z, "drug_info.tsv", ("improve_chem_id", "canSMILES", "PUBCHEM_ID")):
            key = row["improve_chem_id"].strip()
            if key in drugs and row["canSMILES"].strip():
                aliases[key].add(row["canSMILES"].strip())
        stats = Counter()
        affected_ids = defaultdict(set)
        for drug in drugs:
            if len(refs[drug]) != 1:
                stats["missing_or_nonunique_reference_ids"] += 1
                continue
            reference = next(iter(refs[drug]))
            if not aliases[drug]:
                stats["no_nonblank_alias_ids"] += 1
            for alias in aliases[drug]:
                category = classify(reference, alias, Chem)
                stats["alias_" + category] += 1
                affected_ids[category].add(drug)
            if len(aliases[drug]) > 1:
                stats["ids_with_multiple_alias_strings"] += 1
        report = {
            "schema": "tumorlab-improve-rdkit-structure-concordance-v1",
            "rdkit_version": version,
            "response_linked_drug_ids": len(drugs),
            "aggregate_alias_comparisons": dict(sorted(stats.items())),
            "distinct_drug_ids_by_category": {
                k: len(v) for k, v in sorted(affected_ids.items())
            },
            "outcome_values_read": False,
            "chemical_equivalence_fully_adjudicated": False,
            "training_performed": False,
            "limitations": [
                "Exact RDKit canonical isomeric structure match does not establish formulation or assay equivalence.",
                "Connectivity-only agreement may mask stereo differences or missing stereo assignments.",
                "Different SMILES structures may involve tautomer, protonation, or normalization conventions; not automatically distinct drugs.",
                "No PubChem external validation, environment chemistry, concentrations, or assay occasion resolved.",
                "All raw identifiers and source structure strings excluded from output.",
            ],
        }
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = ap.parse_args(argv)
    try:
        print(json.dumps(audit_structure_equivalence(args.archive), indent=2))
    except (RuntimeError, OSError, ValueError, zipfile.BadZipFile) as exc:
        ap.exit(2, "IMPROVE RDKit equivalence audit: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
