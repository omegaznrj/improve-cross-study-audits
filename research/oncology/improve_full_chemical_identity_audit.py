"""Read-only IMPROVE all-response-chemical concordance check; not structural adjudication."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import zipfile

from research.oncology.improve_archive_audit import DEFAULT_ARCHIVE
from research.oncology.improve_split_audit import _read_response
from research.oncology.improve_chemical_crossstudy_audit import _chemical_records


def audit_full_chemical_identity(path: Path = DEFAULT_ARCHIVE) -> dict:
    if not path.is_file():
        raise ValueError("archive not found")
    with zipfile.ZipFile(path) as z:
        members = z.namelist()
        if len(members) != len(set(members)):
            raise ValueError("duplicate archive entries")
        rows = _read_response(z)
        selected = set(rows["drug_code_by_id"])
        selected.discard("")
        result = _chemical_records(z, selected)
    return {
        "schema": "tumorlab-improve-all-drugs-chemical-text-v1",
        "response_rows": rows["rows"],
        "response_linked_drug_ids": len(selected),
        "chemical_diagnostics": result,
        "outcome_values_read": False,
        "structural_chemistry_adjudicated": False,
        "training_performed": False,
        "limitations": [
            "Covers all response-linked IDs across five studies, not merely CCLE/gCSI training drugs.",
            "Repeated metadata keys are not joined to response rows.",
            "Distinct SMILES text can represent equivalent molecules; exact matching is not RDKit standardization.",
            "Different fingerprints do not prove distinct chemical identity, and identical fingerprints are not definitive equivalence.",
            "Salt, tautomer, stereochemistry, scaffold and PubChem validity remain unverified.",
        ],
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    args = p.parse_args(argv)
    try:
        print(json.dumps(audit_full_chemical_identity(args.archive), indent=2))
    except (OSError, ValueError, UnicodeDecodeError, zipfile.BadZipFile) as err:
        p.exit(2, "IMPROVE chemical identity audit: " + str(err) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
