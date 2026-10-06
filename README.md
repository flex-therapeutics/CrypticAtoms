# CrypticAtoms

Benchmarking protein cryptic pocket sampling at the atomic scale: 21 apo–holo reference pairs across 17 proteins.

This repository contains the benchmark annotations, exact prepared reference structures, archived sample scores, and code for the metrics and analyses presented in the manuscript. It includes nine sampling pipelines and two Rosetta controls. The manuscript's Figure 4 compares backbone and side-chain recovery for c-Abl, HIV RT (3V81), and GLTP.

## Install

Python 3.12 was used for release validation. From the repository root:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-tested.txt
```

The tested requirements pin the validation environment. The smaller `requirements.txt` supports structural metrics and table reproduction; `requirements-geometry.txt` adds RDKit for generated-ligand coordinate diagnostics. Historical ligand results used RDKit 2022.09.5. Their supplied tables reproduce in the tested environment; coordinate-level RDKit results should be checked when changing versions.

## Reproduce the manuscript analyses

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python reproduce.py
python -m unittest discover -s tests -v
```

The reproduction command copies inputs to `reproduced/`, runs the analyses there, and compares the resulting tables with the manuscript's archived values. It leaves source tables unchanged. It checks recovery tiers and pocket-map spread for 231 method–target combinations, joint apo/holo recovery, docking-improvement recovery, chronology, stable-residue retention, generated-ligand summaries, backbone/side-chain diagnostics, and paired statistical analyses. `reproduced/validation.json` records the comparisons. The coordinate tests check all 21 reference baselines, invariance to rigid-body transformations, and six worked examples from Figure 4.

The statistical scripts preserve the original 20,000 protein-cluster bootstrap resamples, exact cluster sign flips, fixed seeds, complete-case comparisons, and separate Holm families. Alternate holo references remain grouped by protein. Missing paired observations are not assigned zero scores. Sample-success percentages use the stated 30-sample budget. Docking-improvement counts treat pairs without eligible samples as unsuccessful and exclude reference pairs without positive holo improvement.

## Score a new conformation

```bash
python -m benchmark.score --target cAbl --sample examples/cAbl_learned.pdb --chain A
```

Use the relevant target name from `benchmark/targets.json`, a PDB or mmCIF file, and its protein chain. Residue numbers must match the supplied references. For sequentially numbered predictions, pass `--residue-map mapping.json`, containing a complete mapping from input residue numbers (string keys) to reference residue numbers (integer values). This command does not infer sequence correspondence. Resolve insertions, deletions, chain selection and construct differences before scoring; duplicate residue/atom keys are rejected.

The output includes distances to apo and holo, recovery tiers, matched alignment counts, missing pocket residues, stable-residue retention, pocket-residue clashes, and backbone/side-chain diagnostics. The archived metric uses available matched heavy atoms; inspect the reported coverage rather than interpreting a partially observed pocket as a complete reconstruction.

- `benchmark/metrics.py`: pocket-Cα alignment and element-constrained Hungarian atom correspondence. Local residue alignment selects correspondence; final errors are measured in the pocket frame.
- `benchmark/chemical_symmetry.py`: strict atom-name correspondence with permitted chemical symmetry swaps, used in the paired sensitivity audit.
- `benchmark/diagnostics.py`: complete-residue stability, backbone/side-chain errors and severe inter-residue contacts, excluding pairs linked by one or two reference bonds.
- `benchmark/ligand_metrics.py`: ligand radius, Morgan-fingerprint similarity and bond/clash geometry. A molecular graph and generated protein coordinates are required; missing or unparseable graphs are not passing outcomes.
- `analysis/compute_recovery.py`: holo recovery, joint apo/holo recovery and two-reference distance-map spread.
- `analysis/compute_docking_recovery.py`: recovery of 50%, 75% and 100% of the reference holo docking improvement.
- `analysis/compute_reference_chronology.py`: cumulative recovery and comparisons across documented date boundaries. Dates do not establish training membership.
- `analysis/compute_sensitivity.py`: the finalized 4,002-sample chemical-correspondence comparison.
- `analysis/compute_diagnostics.py`: stable-residue and ligand summaries and Figure 4 diagnostics from supplied sample-level measurements.
- `analysis/review_statistics.py` and `analysis/review_docking_statistics.py`: paired inference.

## Annotations and data

`benchmark/targets.json` is the authoritative annotation file: PDB records and chains, prepared-reference paths, ligand CCD identifiers, pocket residues, moving residues, and archived reference distances. Preserve those masks and baselines when comparing methods. Unarchived shell definitions remain null where appropriate.

`references/manifest.csv` records the 42 prepared reference files, public PDB source links and SHA-256 hashes. The supplied files preserve the original preparations and formats; a fresh raw PDB download may differ in assembly, chain, residue numbering or missing-atom treatment. References with alternate holo structures share proteins; the two HIV RT comparisons retain distinct apo preparations.

`analysis/` contains sample-level measurements and manuscript figure data. `controls/` contains archived Rosetta scores and optional generation/scoring protocols. `examples/` contains the six selected protein conformations shown in Figure 4 and a manifest with their reported metric values.

## Scope

The repository reproduces the reported analyses from supplied measurements and scores new, appropriately mapped protein conformations. It does not bundle all generated ensembles, model checkpoints, inference environments, generated ligand coordinates, docking poses, or PyMOL rendering assets. It therefore does not rerun all inference or redocking experiments. The archived geometry tables remain available, and the coordinate diagnostic functions can be applied when those inputs are supplied.

Optional Rosetta generation requires a separately installed, appropriately licensed PyRosetta environment. `python controls/run_rosetta.py --workers 8 --samples 30` prepares the bundled apo inputs and runs the two controls; `python controls/score_rosetta.py` scores those outputs. These commands are not run by the reproduction script. Their protocol sources are supplied; new Rosetta generation was not executed during release preparation.

## Release status and rights

The benchmark repository is hosted at https://github.com/flex-therapeutics/CrypticAtoms. Code is distributed under the MIT license supplied by Flex Therapeutics (see `LICENSE`). Third-party reference data retain their applicable terms and provenance. No Rosetta binaries, proprietary model weights, credentials, or full inference environments are included.
