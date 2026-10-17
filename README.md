# IMPROVE cross-study drug-response: reproducible data qualification

**Ziad Fawaz · Independent computational oncology research · October 2026**

This repository contains read-only Python audits and synthetic regression tests developed while investigating the IMPROVE Cross-Study Analysis (CSA) monotherapy drug-response benchmark. Its purpose is to help researchers distinguish what the dataset *records* from what a predictive or mechanistic model might *assume* about chemical identity, experimental independence, and molecular feature values.

The work is exploratory and has not undergone independent scientific review. No treatment-response predictor has been trained or evaluated here. These static **cell-line assays are not longitudinal patient outcomes**, and these scripts are not clinical decision-support software.

## Why these questions matter

A harmonized drug identifier is useful for joining data across screens, but it does not necessarily establish which exact stereoisomer or preparation generated a particular response. Similarly, repeated cell-line–drug identifiers may reflect legitimate repeats or source assays, rather than independent experimental occasions. Large computed descriptors may reflect graph conventions, not erroneous chemistry.

These are distinct issues with a common consequence: the scientific question determines which source distinctions must survive preprocessing. The audits therefore preserve raw provenance, separate uncertainty from error, and avoid silently creating labels or independent observations.

## Investigated source and findings

The original IMPROVE CSA v2 archive contains **587,709 response records**, **785 response-linked cell-line identifiers**, and **749 response-linked harmonized chemical identifiers** across CCLE, CTRPv2, GDSCv1, GDSCv2, and gCSI. These are *not* 587,709 independent patients.

- **Chemical representations:** 86 of 749 harmonized drug IDs have multiple distinct nonblank metadata SMILES. Iterative RDKit comparisons found five harmonized IDs with explicitly conflicting assigned stereochemical representations: `Drug_109`, `Drug_1094`, `Drug_192`, `Drug_538`, and `Drug_885`. Source-specific experimental preparations remain unverified. **The five-case adjudication also involved manual inspection of source metadata, not a single automated flag.**
- **Cross-study provenance:** repeated harmonized cell-line–compound pairs and overlaps between folds or studies can invalidate a naive *independent-pair* interpretation; they do not by themselves establish duplicated biological assays, mistakes in the source benchmark, or an inappropriate official task.
- **Numeric feature qualification:** linked methylation measurements include 625,792 missing cells (~4.07%); all 785 linked cell models contain at least one missing methylation value. Continuous and discretized copy number share the same 192 missing-value positions among all 785 linked models. Numerous mutation-count features are zero-heavy, which is not itself invalid.
- **Mordred/RDKit:** 23 of 749 linked reference SMILES are multi-fragment structures. All 23 have IMPROVE `Diameter` and `Radius` values of 100,000,000. Controlled synthetic molecules independently reproduce the same values with **Mordred 1.2.0 and RDKit 2026.09.1**; this supports a disconnected-graph explanation but does not confirm the generator versions used by IMPROVE.

The reports do not assert that any original dataset annotations are wrong, that experimental formulations were misidentified, or that published prediction performance is affected. The archived test and code checkpoint is noted in `RELEASE_MANIFEST.json`.

## Obtain the original data separately

Dataset: *Cross-Study Benchmark Dataset for Monotherapy Drug Response Prediction*, v2.

https://doi.org/10.5281/zenodo.15258883

Download the original `csa_data.zip` from the publisher and place it, **without extracting or editing its contents**, at:

```text
data/raw/improve_cross_study/csa_data.zip
```

The archive is not redistributed by this repository. Its source checksum and membership can be checked using `python -m research.oncology.improve_archive_audit`. Source licensing and permissions remain those of the original data provider.

## Reproduce with Python 3.11

Create a clean environment. On Windows, a **short path** can help with RDKit DLL loading:

```cmd
py -3.11 -m venv C:\TLVerify
C:\TLVerify\Scripts\activate.bat
python -m pip install -r requirements-audit.txt
```

On Linux/macOS, replace this with a standard Python 3.11 `venv` and activate it before installation. Exact package versions are recorded in `requirements-audit.txt`; platform-compatible builds are required.

From the repository root, test the **synthetic logic first**:

```cmd
python -m unittest discover -s tests -p "test_improve_*.py" -v
```

Then, after obtaining the original archive independently:

```cmd
python -m research.oncology.improve_full_chemical_identity_audit
python -m research.oncology.improve_structure_equivalence_audit
python -m research.oncology.improve_atom_mapped_stereo_audit
python -m research.oncology.improve_stereo_resolution_audit
python -m research.oncology.improve_five_study_overlap_audit
python -m research.oncology.improve_repeat_provenance_audit
python -m research.oncology.improve_numeric_quality_audit
python -m research.oncology.improve_missingness_geometry_audit
python -m research.oncology.improve_targeted_numeric_diagnostics
python -m research.oncology.improve_mordred_descriptor_extremes
python -m research.oncology.improve_mordred_fragment_association
```

To reproduce disconnected-graph behavior **without the IMPROVE archive**:

```cmd
python -m research.oncology.improve_synthetic_graph_distance_verification
python -m research.oncology.improve_isolated_mordred_reproduction
```

All audits are descriptive; no modeling or data correction is performed. Some reports output pooled diagnostic counts; **those counts are not preprocessing parameters to fit on held-out data**. The chemical-level five-case review additionally requires independent source-metadata interpretation. The complete 749-drug Mordred regeneration comparison was paused and is **not part of this release**.

## Boundaries, reproducibility and reuse

- No original study rows, raw SMILES export, clinical/patient records, downloaded ZIP, private review spreadsheets, credentials, model weights, or LLM conversation files are provided.
- Individual reported numbers are tied to the declared IMPROVE source release. Cross-study repeated keys do not prove the same laboratory occasion.
- This is an independent, **LLM-assisted** analytical and programming investigation. Results came from executing Python/RDKit checks, testing synthetic controls, and inspecting original source metadata; human oncology and chemical adjudication is still needed.
- Source references and methods are available for specialists to critique or extend in other scientific directions.
- The code is provided under the accompanying MIT license. The dataset and dependencies retain their own terms.

## References

- Partin et al., *Cross-Study Benchmark Dataset for Monotherapy Drug Response Prediction*, Zenodo v2 (2025), https://doi.org/10.5281/zenodo.15258883
- Partin et al., *Benchmarking community drug response prediction models: datasets, models, tools, and metrics for cross-dataset generalization analysis*, *Briefings in Bioinformatics* 27(1), bbaf667 (2026), https://doi.org/10.1093/bib/bbaf667
- Moriwaki et al., *Mordred: a molecular descriptor calculator*, *Journal of Cheminformatics* 10:4 (2018), https://doi.org/10.1186/s13321-018-0258-y
