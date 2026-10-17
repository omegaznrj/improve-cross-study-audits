"""Read-only atom-correspondence stereochemistry audit of IMPROVE aliases.

Compares RDKit CIP atom assignments and E/Z bond stereo across *all* full
connectivity graph mappings (capped). Conservative on symmetry and missing
stereo; does not assert formulation equivalence. Aggregate reporting only.
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
from research.oncology.improve_structure_equivalence_audit import _rdkit, _signature

MAX_MATCHES = 256


def _stereo(m, Chem):
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in m.GetAtoms()
             if a.HasProp("_CIPCode")}
    bonds = {}
    for b in m.GetBonds():
        if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            bonds[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = str(b.GetStereo())
    return atoms, bonds


def mapped_stereo_relation(a: str, b: str, Chem) -> str:
    sig_a, sig_b = _signature(a, Chem), _signature(b, Chem)
    if sig_a is None or sig_b is None:
        return "invalid_structure"
    if sig_a[0] == sig_b[0]:
        return "same_canonical_isomeric_structure"
    if sig_a[1] != sig_b[1]:
        return "different_connectivity"
    ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
    # Full graph matches rather than merely shared substructure. Match chirality
    # disabled so stereo disagreements are observable.
    matches = mb.GetSubstructMatches(
        ma, useChirality=False, uniquify=False, maxMatches=MAX_MATCHES + 1
    )
    if not matches or len(matches) > MAX_MATCHES:
        return "mapping_unresolved_or_symmetry_cap"
    ca, ba = _stereo(ma, Chem)
    cb, bb = _stereo(mb, Chem)
    diagnoses = set()
    for mapping in matches:
        mapped_a_atoms = {mapping[i]: v for i, v in ca.items()}
        mapped_a_bonds = {
            frozenset(mapping[i] for i in key): v for key, v in ba.items()
        }
        atom_overlap = set(mapped_a_atoms) & set(cb)
        bond_overlap = set(mapped_a_bonds) & set(bb)
        if any(mapped_a_atoms[i] != cb[i] for i in atom_overlap) or any(
            mapped_a_bonds[i] != bb[i] for i in bond_overlap
        ):
            diagnoses.add("assigned_stereo_conflict")
        elif set(mapped_a_atoms) != set(cb) or set(mapped_a_bonds) != set(bb):
            diagnoses.add("incomplete_stereo_annotation")
        else:
            diagnoses.add("stereo_mapping_equivalent_but_serialization_differs")
    return next(iter(diagnoses)) if len(diagnoses) == 1 else "ambiguous_symmetric_mapping"


def audit_atom_mapped_stereo(path=DEFAULT_ARCHIVE):
    Chem, version = _rdkit()
    path = Path(path)
    if not path.is_file():
        raise ValueError("archive not found")
    with zipfile.ZipFile(path) as z:
        entries = z.namelist()
        if len(entries) != len(set(entries)):
            raise ValueError("duplicate archive members")
        drug_ids = set(_read_response(z)["drug_code_by_id"]) - {""}
        refs, aliases = defaultdict(set), defaultdict(set)
        for row in _read_tsv(z, "drug_SMILES.tsv", ("improve_chem_id", "canSMILES")):
            key = row["improve_chem_id"].strip()
            if key in drug_ids:
                refs[key].add(row["canSMILES"].strip())
        for row in _read_tsv(z, "drug_info.tsv", ("improve_chem_id", "canSMILES", "PUBCHEM_ID")):
            key = row["improve_chem_id"].strip()
            if key in drug_ids and row["canSMILES"].strip():
                aliases[key].add(row["canSMILES"].strip())
        counts, affected = Counter(), defaultdict(set)
        for key in drug_ids:
            if len(refs[key]) != 1:
                counts["missing_or_nonunique_reference"] += 1
                continue
            reference = next(iter(refs[key]))
            if not aliases[key]:
                counts["missing_alias"] += 1
            for alias in aliases[key]:
                category = mapped_stereo_relation(reference, alias, Chem)
                counts[category] += 1
                affected[category].add(key)
    return {
        "schema": "improve-atom-mapped-stereo-v1",
        "rdkit_version": version,
        "response_linked_drug_ids": len(drug_ids),
        "alias_comparison_counts": dict(sorted(counts.items())),
        "affected_ids_by_class": {k: len(v) for k, v in sorted(affected.items())},
        "chemical_formulations_verified": False,
        "source_records_modified": False,
        "raw_identifiers_exported": False,
        "limitations": [
            "CIP agreement is representation-level evidence, not validation of chemical provenance.",
            "Stereo bond labels and symmetry are conservatively treated; ambiguous matches stay unresolved.",
            "Unassigned stereochemistry requires verified source chemical records to resolve.",
            "No PubChem verification, assay concentrations, formulations, or human case evidence.",
        ],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = p.parse_args(argv)
    try:
        print(json.dumps(audit_atom_mapped_stereo(args.archive), indent=2))
    except (RuntimeError, OSError, ValueError, zipfile.BadZipFile) as exc:
        p.exit(2, "IMPROVE atom mapped stereo audit: " + str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
