"""Phase 1.4 stereochemistry triage, not forced chemical equivalence.

Read-only original IMPROVE ZIP, RDKit canonical isomeric SMILES comparison.
Reports aggregate mutually specified stereo conflicts separately from
unspecified stereo and keeps ambiguous identities unresolved.
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
from research.oncology.improve_structure_equivalence_audit import _rdkit, _signature


def stereo_triage(reference: str, alias: str, Chem) -> str:
    a, b = _signature(reference, Chem), _signature(alias, Chem)
    if a is None or b is None:
        return "invalid_or_blank"
    if a[0] == b[0]:
        return "same_isomeric_structure"
    if a[1] != b[1]:
        return "other_structural_discrepancy"
    # Compare assigned tetrahedral stereo atom invariants, without assuming that
    # unspecified stereo agrees or that CIP atom indices are stable.
    ra, rb = Chem.MolFromSmiles(reference), Chem.MolFromSmiles(alias)
    def assigned(m):
        chiral = Chem.FindMolChiralCenters(m, includeUnassigned=True)
        return len(chiral), sum(cip != "?" for _, cip in chiral)
    n1, q1 = assigned(ra)
    n2, q2 = assigned(rb)
    if q1 == 0 or q2 == 0 or q1 != q2 or n1 != n2:
        return "stereo_annotation_incomplete_or_different"
    # When stereo assignments exist on both sides, canonical isomeric strings
    # differ, but we must not assert an enantiomer relationship blindly:
    # bond stereo and atom correspondence require separate mapping.
    return "specified_stereo_disagreement_needs_atom_mapping"


def audit_stereo(path=DEFAULT_ARCHIVE):
    Chem, version = _rdkit()
    path = Path(path)
    if not path.is_file():
        raise ValueError("archive not found")
    with zipfile.ZipFile(path) as z:
        members = z.namelist()
        if len(members) != len(set(members)):
            raise ValueError("duplicate archive members")
        drug_ids = set(_read_response(z)["drug_code_by_id"]) - {""}
        refs = defaultdict(set)
        aliases = defaultdict(set)
        for row in _read_tsv(z, "drug_SMILES.tsv", ("improve_chem_id", "canSMILES")):
            key = row["improve_chem_id"].strip()
            if key in drug_ids:
                refs[key].add(row["canSMILES"].strip())
        for row in _read_tsv(z, "drug_info.tsv", ("improve_chem_id", "canSMILES", "PUBCHEM_ID")):
            key = row["improve_chem_id"].strip()
            if key in drug_ids and row["canSMILES"].strip():
                aliases[key].add(row["canSMILES"].strip())
        classes = Counter()
        ids_by_class = defaultdict(set)
        for key in drug_ids:
            if len(refs[key]) != 1:
                classes["nonunique_or_missing_reference"] += 1
                continue
            reference = next(iter(refs[key]))
            if not aliases[key]:
                classes["missing_alias"] += 1
            for alias in aliases[key]:
                category = stereo_triage(reference, alias, Chem)
                classes[category] += 1
                ids_by_class[category].add(key)
    return {
        "schema": "improve-stereo-discrepancy-triage-v1",
        "rdkit_version": version,
        "response_linked_drug_ids": len(drug_ids),
        "comparison_counts": dict(sorted(classes.items())),
        "affected_drug_id_counts": {k: len(v) for k, v in sorted(ids_by_class.items())},
        "chemical_identical_assay_preparations_proven": False,
        "raw_ids_or_smiles_exported": False,
        "source_modified": False,
        "notes": [
            "Stereo ambiguity cannot be fixed by arbitrarily assigning or stripping stereochemistry.",
            "Specified stereo differences require atom-mapped review including double-bond stereo.",
            "A shared connectivity graph does not establish identical molecular intervention.",
            "Source paper, authoritative chemical identifiers and physical formulations remain future verification gates.",
        ],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = p.parse_args(argv)
    try:
        print(json.dumps(audit_stereo(args.archive), indent=2))
    except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as exc:
        p.exit(2, "IMPROVE stereo triage: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
