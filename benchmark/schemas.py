"""Pydantic schemas for the cryptic-pocket benchmark.

Two things live here:
- `TargetSpec`: the per-target benchmark input (one entry per protein in `targets.json`).
- `ResultRow`: the per-sample output schema every method's exporter writes into
  `results/<method>/<run_id>.parquet`.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import pandas as pd
from pydantic import BaseModel, Field, field_validator

# Repo root = parent of the benchmark/ package. Paths in targets.json may be
# stored repo-relative (e.g. "data/structures/apo/CB1_apo.pdb") so the repo
# stays portable; field validators below resolve them to absolutes at load time.
REPO_ROOT = Path(__file__).resolve().parent.parent


def _resolve_repo_path(p: str | None) -> str | None:
    if p is None or p == "":
        return p
    pp = Path(p)
    if pp.is_absolute():
        return str(pp)
    return str((REPO_ROOT / pp).resolve())


class ApoSpec(BaseModel):
    pdb_code: str
    chain: str
    contig: str | None = None  # e.g. "A96-292"
    path: str    # apo structure (PDB or CIF); repo-relative in targets.json, resolved on load

    @field_validator("path")
    @classmethod
    def _resolve_path(cls, v: str) -> str:
        return _resolve_repo_path(v)


class HoloSpec(BaseModel):
    pdb_code: str
    chain: str
    path: str    # holo structure (CIF); repo-relative in targets.json, resolved on load

    @field_validator("path")
    @classmethod
    def _resolve_path(cls, v: str) -> str:
        return _resolve_repo_path(v)


class LigandSpec(BaseModel):
    drug_ccd: str                         # primary ligand of interest
    transferred_ccds: list[str]           # ccds copied apo→holo frame (drug + structural cofactors)
    extra_cofactor_ccds: list[str] = []   # non-drug ligands present in holo (e.g. ZN, GDP)


class ShellSpec(BaseModel):
    # residues split by distance threshold from the drug heavy atoms.
    # `pocket` = within threshold; `far` = beyond. Used for shell-fix experiments.
    pocket: list[int]
    far: list[int] | None = None  # Unarchived shell complements remain unknown.


class PocketSpec(BaseModel):
    # canonical pocket: residues within 5 Å of drug heavy atoms in the holo frame.
    residues_5A: list[int]
    shells: dict[int, ShellSpec] | None = None  # keys: 10, 15 (Å thresholds)
    # subset of `residues_5A` whose heavy-atom set has any atom displaced more than
    # `Baselines.moving_threshold_A` between apo and holo, evaluated in the **pocket-Cα
    # alignment frame** (the canonical convention). Populated by `compute_baselines.py`.
    # This is the residue set used to compute `ResultRow.aa_rmsd_moving_vs_{apo,holo}`.
    moving_residues_2A: list[int] = []


class AlignmentBaselines(BaseModel):
    """apo↔holo distances measured in one specific superposition frame.

    The moving-residue set is taken from `PocketSpec.moving_residues_2A` (canonical,
    pocket-frame), so `aa_rmsd_moving` here is "how that set looks in *this* frame".
    """
    ca_rmsd_full: float | None = None
    ca_rmsd_pocket: float | None = None
    aa_rmsd_moving: float | None = None


class Baselines(BaseModel):
    # primary: superposition on pocket-Cα (residues_5A). Canonical for ResultRow metrics.
    align_pocket: AlignmentBaselines = Field(default_factory=AlignmentBaselines)
    # secondary: superposition on full-chain matched Cα. For reference / comparison.
    align_full: AlignmentBaselines = Field(default_factory=AlignmentBaselines)
    moving_threshold_A: float = 2.0
    computed_with: str | None = None  # short note on alignment/method


class TargetSpec(BaseModel):
    name: str
    apo: ApoSpec
    holo: HoloSpec
    ligands: LigandSpec
    pocket: PocketSpec
    baselines: Baselines = Field(default_factory=Baselines)
    apo_lig_path: str | None = None  # apo protein + transferred holo ligands (for rfd3); resolved on load

    @field_validator("apo_lig_path")
    @classmethod
    def _resolve_apo_lig_path(cls, v: str | None) -> str | None:
        return _resolve_repo_path(v)


class ResultRow(BaseModel):
    """One generated structure from one method/run/case/sample.

    Long-format: every method's exporter emits one of these per output structure.
    `hparams` is a free-form JSON dict — method-specific knobs (partial_t, rank, ddim_steps, ...).
    """
    target: str
    method: Literal["rfd3", "dynamicbind", "apo2mol", "bioemu", "confornet", "Flexgen"]
    run_id: str         # free-form identifier (e.g. "multi_target_2026-05-11")
    case: str           # method-defined scenario name (e.g. "fix10_fixlig", "default")
    sample_idx: int     # 0-based within (target, case, run_id)
    hparams: dict[str, Any] = Field(default_factory=dict)
    struct_path: str    # absolute path to generated structure

    # canonical metrics — null when not applicable
    rmsd_full_vs_apo: float | None = None
    rmsd_full_vs_holo: float | None = None
    rmsd_pocket_vs_apo: float | None = None
    rmsd_pocket_vs_holo: float | None = None
    # all-heavy-atom RMSD over the target's `pocket.moving_residues_2A` set
    aa_rmsd_moving_vs_apo: float | None = None
    aa_rmsd_moving_vs_holo: float | None = None

    # ligand metrics — null for methods that don't predict the ligand
    lig_rmsd: float | None = None
    lig_rmsd_pocket: float | None = None
    lig_centroid_dist: float | None = None

    status: Literal["ok", "failed"] = "ok"
    error: str | None = None

    @field_validator("hparams")
    @classmethod
    def _hparams_json_serializable(cls, v: dict[str, Any]) -> dict[str, Any]:
        json.dumps(v)  # raises if any value is not JSON-serializable
        return v


def rows_to_parquet(rows: list[ResultRow], path: str | Path) -> None:
    """Write a list of ResultRow to parquet. `hparams` is serialized to a JSON string column."""
    records = []
    for r in rows:
        d = r.model_dump()
        d["hparams"] = json.dumps(d["hparams"], sort_keys=True)
        records.append(d)
    df = pd.DataFrame.from_records(records)
    df.to_parquet(path, index=False)


def read_results(path: str | Path) -> pd.DataFrame:
    """Read a results parquet back, parsing `hparams` JSON into dicts."""
    df = pd.read_parquet(path)
    df["hparams"] = df["hparams"].apply(json.loads)
    return df
