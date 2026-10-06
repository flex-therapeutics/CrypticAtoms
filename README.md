# CrypticAtoms

## References

The benchmark contains 21 apo–holo reference pairs across 17 proteins. `benchmark/targets.json` records the PDB identifiers, chains, ligand CCD identifiers, pocket residues, moving residues, and reference distances.

`references/` contains the 42 prepared apo and holo structures. `references/manifest.csv` lists their paths, public PDB source links, and SHA-256 hashes. Use these prepared structures and annotations for consistent scoring; raw PDB downloads may differ in assembly, residue numbering, or missing-atom treatment. Alternate holo references share proteins, and the two HIV RT comparisons retain distinct apo preparations.

## Run

From the repository root with Python 3.12:

```bash
python -m pip install -r requirements-tested.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python reproduce.py
```
