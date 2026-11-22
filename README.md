# ECRA Rainbow Trout: Chemical-space and evidence audit

**Project status (22 November 2026): SCIENTIFIC WORKING RELEASE; QSAR DATA GATE: HOLD.**

This repository contains a reproducibility and evidence-review package for the *Molecular Diversity* working manuscript **“Evidence-gated chemical-space coverage and experimental comparability in historical rainbow trout acute-toxicity data: implications for similarity-based read-across.”**

## Dataset and scientific interpretation

- The frozen, retained EPA ECOTOX-derived *Oncorhynchus mykiss* 96-hour LC50 register comprises **169 observations, 107 distinct CAS identifiers and 59 reference IDs**.
- Fingerprint and chemical-space analyses are descriptive and provisional: **82 compounds structure-screenable; 25 excluded**.
- Earlier source-table checks documented **22 numerical matches across three original articles**, of which **20 were trout-compatible**. This does **not** constitute independent-study approval or model eligibility.
- **27 source-linked records** from references 12004 (10), 10775 (7) and 12665 (10) still require original numerical table authentication and scientific adjudication.
- **0 observations** have passed the controlled ECRA QSAR release gate. No valid ECRA empirical toxicity prediction or experimentally validated read-across performance is claimed.

## Repository layout

- `data/`: frozen analytical tables and machine-readable evidence snapshots.
- `results/`: descriptive results, comparisons, and visualizations.
- `provenance/`: dataset audit trail, source references, processing metadata and file checksums.
- `manuscript/`: unpublished scientific working documents (retain private until publication clearance).
- `code/`: executable processing scripts, when provenance and licenses are verified.
- `supplementary/`: journal supplementary assets.
- `docs/`: decision logs and analysis limitations.

## Important limitations

**Bibliographic verification is not numerical-table authentication. A matching LC50 does not prove an independent experiment.** Model release requires verified chemical identity, original source page/row/value/unit, exposure basis, study lineage, and two independent scientific reviewer decisions.

The current repository is **private**. Do not cite an internal GitHub commit as a public archived dataset. Public deposit, redistribution of ECOTOX-derived records, DOI minting, and final supplement publication require provenance/license review, authorship clearance, and release checks.

See [DATA_LIMITATIONS.md](DATA_LIMITATIONS.md) and [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md).

## Status of this deposit

This is an incremental deposit of files verifiably available during the current session; the contents manifest lists included and excluded material. A complete original ECOTOX extract, all scripts/environments, and all source-table PDFs were **not** assumed present or independently verified.
