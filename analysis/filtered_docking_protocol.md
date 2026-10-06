# Docking after a pocket-residue clash filter

Exclude a sampled conformation if more than 5% of its reference-pocket residues participate in the existing severe inter-residue heavy-atom clash criterion. Retain conformations at exactly 5%. This filters entire conformations; it does not delete residues, repair coordinates, graft pockets, or alter docking boxes.

Select up to ten eligible conformations per method/target in original sample-index order, without using docking scores for selection. Use every eligible conformation if fewer than ten remain. Reuse matched finite docking scores where available and dock the selected conformations that lack results using the established Vina 1.2.7, exhaustiveness 8, rigid-receptor protocol and reference ligand/box. Apo2Mol is the current 1,000-step run and uses direct generated-pocket docking. BBFlow uses the complete matching structure and docking set.

`clash_filtered_docking_per_sample.csv` records selection, clashes, source structure paths, reused scores, and jobs required. `new_docking_results.csv` records supplementary docking. New results and receptor/pose files are isolated below each method directory. Residue mapping is checked before docking, using existing sequence-verified mappings where available. Input coordinates are not repaired or optimized.

`clash_filtered_docking_per_target.csv` reports the best score among the selected eligible conformations minus the apo reference score, in kcal/mol; negative is favorable. Empty targets remain missing, not zero. `clash_filtered_docking_summary.csv` gives mean and median across available targets, target coverage, full-ten coverage, and the original scores on the same target subset. Aggregates over differing target sets and docking budgets should not be interpreted as a matched ranking. A full-ten-only summary is included as an additional diagnostic, not a replacement for reporting missing targets.

Figure 3c uses the best selected score per target. Figure A5 shows the individual scores of exactly the same selected conformations.

Prepared receptors and supplementary docking logs remain external to this analysis package.
