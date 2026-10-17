"""Read-only chemistry-aware reconciliation of IMPROVE response-linked SMILES.

RDKit is REQUIRED. Classifies each metadata representation against the sole
reference SMILES by exact graph, stereo-agnostic graph, largest-fragment-only
(salts/mixtures), or unresolved. No chemical identity is silently rewritten.
Only aggregate counts are emitted; no compound IDs or SMILES leave the ZIP.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import zipfile

from research.oncology.improve_archive_audit import DEFAULT_ARCHIVE
from research.oncology.improve_split_audit import _read_response
from research.oncology.improve_chemical_crossstudy_audit import _read_tsv


def _rdkit():
    try:
        from rdkit import Chem, rdBase
    except ImportError as e:
        raise RuntimeError("RDKit is required: install rdkit in the project .venv") from e
    return Chem, rdBase.rdkitVersion


def classify(reference: str, alias: str, Chem) -> str:
    if not reference or not alias:
        return "missing_string"
    try:
        a = Chem.MolFromSmiles(reference)
        b = Chem.MolFromSmiles(alias)
    except Exception:
        return "invalid_smiles"
    if a is None or b is None:
        return "invalid_smiles"
    if Chem.MolToSmiles(a, isomericSmiles=True) == Chem.MolToSmiles(b, isomericSmiles=True):
        return "same_isomeric_graph"
    if Chem.MolToSmiles(a, isomericSmiles=False) == Chem.MolToSmiles(b, isomericSmiles=False):
        return "same_connectivity_stereo_difference"
    def largest(mol):
        pieces = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        return sorted(pieces, key=lambda x: (x.GetNumHeavyAtoms(), Chem.MolToSmiles(x)))[-1]
    try:
        la, lb = largest(a), largest(b)
        if (len(Chem.GetMolFrags(a)) > 1 or len(Chem.GetMolFrags(b)) > 1):
            if Chem.MolToSmiles(la, isomericSmiles=True) == Chem.MolToSmiles(lb, isomericSmiles=True):
                return "same_largest_fragment_extra_components"
            if Chem.MolToSmiles(la, isomericSmiles=False) == Chem.MolToSmiles(lb, isomericSmiles=False):
                return "same_largest_fragment_stereo_difference"
    except Exception:
        return "unresolved_fragment_parsing"
    return "different_connectivity_or_unresolved_chemistry"


def audit_chemical_equivalence(path: Path = DEFAULT_ARCHIVE) -> dict:
    Chem, rdkit_version = _rdkit()
    if not path.is_file():
        raise ValueError("archive missing")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate ZIP member names")
        rows = _read_response(z)
        ids = set(rows["drug_code_by_id"]) - {""}
        reference = defaultdict(set)
        for item in _read_tsv(z, "drug_SMILES.tsv", ("improve_chem_id", "canSMILES")):
            key = (item["improve_chem_id"] or "").strip()
            if key in ids:
                reference[key].add((item["canSMILES"] or "").strip())
        aliases = defaultdict(set)
        for item in _read_tsv(z, "drug_info.tsv", ("improve_chem_id", "canSMILES", "PUBCHEM_ID")):
            key = (item["improve_chem_id"] or "").strip()
            if key in ids:
                aliases[key].add((item["canSMILES"] or "").strip())
        if any(len(reference[k]) != 1 or not next(iter(reference[k])) for k in ids):
            raise ValueError("reference SMILES missing or ambiguous for response-linked ID")
        tally = Counter()
        compound_flags = Counter()
        for key in ids:
            ref = next(iter(reference[key]))
            statuses = [classify(ref, other, Chem) for other in aliases[key]]
            for kind in statuses:
                tally[kind] += 1
            compound_flags["chemicals_with_nonidentical_alias_text"] += int(
                any(other != ref and bool(other) for other in aliases[key]))
            compound_flags["chemicals_with_unresolved_structure"] += int(
                any(k in {"invalid_smiles", "unresolved_fragment_parsing",
                          "different_connectivity_or_unresolved_chemistry"}
                    for k in statuses))
            compound_flags["chemicals_with_stereo_or_fragment_warning"] += int(
                any("stereo_difference" in k or "fragment" in k for k in statuses))
    return {
        "schema": "improve-rdkit-chem-equivalence-v1",
        "rdkit_version": rdkit_version,
        "response_linked_chemical_ids": len(ids),
        "distinct_alias_strings_compared": sum(tally.values()),
        "comparison_categories": dict(sorted(tally.items())),
        "chemical_flags": dict(sorted(compound_flags.items())),
        "chemistry_adjudicated_for_clinical_or_training_use": False,
        "no_raw_identifiers_or_molecules_exported": True,
        "source_files_modified": False,
        "notes": [
            "Same isomeric graph means RDKit normalization of input SMILES matches, not proof of formulation equivalence.",
            "Stereo-agnostic and largest-fragment matches are warning categories, not permission to merge chemical IDs.",
            "Tautomers, protonation, isotopes, normalization choices and authoritative compound registry links still require review.",
            "Raw metadata and compound descriptors remain untouched; no training or response values interpreted.",
        ],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = p.parse_args(argv)
    try:
        result = audit_chemical_equivalence(args.archive)
    except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as err:
        p.exit(2, "IMPROVE chemistry equivalence: " + str(err) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
