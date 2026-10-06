"""Shared structural metrics used by CrypticAtoms.

The manuscript aligns structures on matched pocket C-alpha atoms. Heavy-atom
correspondence is element-constrained Hungarian assignment within each residue,
chosen after local backbone alignment; RMSD is measured in the pocket frame.
The atom-name and full-chain helpers remain available for diagnostic uses.
"""
from __future__ import annotations

from pathlib import Path

import biotite.structure as struc
import biotite.structure.io as bsio
import numpy as np


def parse_contig_residues(contig: str) -> dict[str, set[int]]:
    """Parse an rfd3-style contig string into chain → set of residue IDs.

    Example:
        'A123-200,15,A216-429' → {'A': {123,...,200, 216,...,429}}
    Numeric-only tokens are gap counts (residues missing) and are ignored.
    """
    out: dict[str, set[int]] = {}
    for tok in contig.split(","):
        tok = tok.strip()
        if not tok or tok.lstrip("-").isdigit():
            continue  # gap count, not a residue range
        chain = tok[0]
        rng = tok[1:]
        if "-" in rng:
            a, b = rng.split("-", 1)
            ids = range(int(a), int(b) + 1)
        else:
            ids = [int(rng)]
        out.setdefault(chain, set()).update(int(i) for i in ids)
    return out


def load_protein_chain(
    path: str | Path,
    chain: str,
    res_ids: set[int] | None = None,
) -> struc.AtomArray:
    """Heavy-atom protein chain, first altloc, first model.

    If `res_ids` is given, restrict to those residue IDs (use this to drop
    crystallization fusion partners — e.g. BRIL/T4L inserts with res_id >= 1000).
    """
    atoms = bsio.load_structure(str(path), model=1, altloc="first")
    atoms = atoms[atoms.chain_id == chain]
    atoms = atoms[struc.filter_amino_acids(atoms)]
    atoms = atoms[atoms.element != "H"]
    if res_ids is not None:
        atoms = atoms[np.isin(atoms.res_id, list(res_ids))]
    return atoms


def _key(arr: struc.AtomArray) -> list[tuple[int, str]]:
    return list(zip(arr.res_id.tolist(), arr.atom_name.tolist()))


def match_atoms(
    a: struc.AtomArray, b: struc.AtomArray
) -> tuple[struc.AtomArray, struc.AtomArray]:
    """Return subsets of a and b sharing (res_id, atom_name), in matching order."""
    a_idx_by_key = {k: i for i, k in enumerate(_key(a))}
    b_idx_by_key = {k: i for i, k in enumerate(_key(b))}
    shared = sorted(set(a_idx_by_key) & set(b_idx_by_key))
    if not shared:
        return a[np.array([], dtype=int)], b[np.array([], dtype=int)]
    a_idx = np.array([a_idx_by_key[k] for k in shared])
    b_idx = np.array([b_idx_by_key[k] for k in shared])
    return a[a_idx], b[b_idx]


def superimpose_on_ca(
    mobile: struc.AtomArray, target: struc.AtomArray
) -> tuple[struc.AtomArray, float, int]:
    """Kabsch-align `mobile` onto `target` using Cα of residues present in both chains.

    Returns the rigidly transformed `mobile` (all heavy atoms moved), the Cα RMSD over
    the matched set, and the number of matched Cα atoms.
    """
    m_ca = mobile[mobile.atom_name == "CA"]
    t_ca = target[target.atom_name == "CA"]
    m_match, t_match = match_atoms(m_ca, t_ca)
    if len(m_match) < 3:
        raise ValueError(f"need ≥3 matched Cα to align; got {len(m_match)}")
    _, transform = struc.superimpose(t_match, m_match)
    aligned_mobile = transform.apply(mobile)
    aligned_ca = transform.apply(m_match)
    ca_rmsd = float(struc.rmsd(t_match, aligned_ca))
    return aligned_mobile, ca_rmsd, len(m_match)


def superimpose_on_pocket_ca(
    mobile: struc.AtomArray, target: struc.AtomArray, pocket_res_ids: list[int]
) -> tuple[struc.AtomArray, float, int]:
    """Kabsch-align `mobile` onto `target` using only Cα of `pocket_res_ids` present in both.

    Returns (aligned_mobile, ca_rmsd over the matched pocket Cα, n matched Cα).
    Use this when you want pocket motion expressed *relative to the pocket frame*.
    """
    pocket_set = set(pocket_res_ids)
    m_ca = mobile[(mobile.atom_name == "CA") & np.isin(mobile.res_id, list(pocket_set))]
    t_ca = target[(target.atom_name == "CA") & np.isin(target.res_id, list(pocket_set))]
    m_match, t_match = match_atoms(m_ca, t_ca)
    if len(m_match) < 3:
        raise ValueError(
            f"need ≥3 matched pocket Cα to align; got {len(m_match)} "
            f"(pocket size {len(pocket_set)})"
        )
    _, transform = struc.superimpose(t_match, m_match)
    aligned_mobile = transform.apply(mobile)
    aligned_ca = transform.apply(m_match)
    ca_rmsd = float(struc.rmsd(t_match, aligned_ca))
    return aligned_mobile, ca_rmsd, len(m_match)


def _kabsch_align_residue(
    a_res: struc.AtomArray, b_res: struc.AtomArray
) -> struc.AtomArray:
    """Per-residue Kabsch alignment of a_res onto b_res using backbone (N, CA, C).
    Returns a copy of a_res with coords in b_res's per-residue frame. Falls back
    to no-op if fewer than 3 backbone atoms align."""
    bb = {"N", "CA", "C"}
    a_bb_mask = np.isin(a_res.atom_name, list(bb))
    b_bb_mask = np.isin(b_res.atom_name, list(bb))
    # Match by name on backbone
    name_to_b = {nm: i for i, nm in enumerate(b_res.atom_name) if nm in bb}
    a_pts, b_pts = [], []
    for i, nm in enumerate(a_res.atom_name):
        if nm in name_to_b:
            a_pts.append(a_res.coord[i])
            b_pts.append(b_res.coord[name_to_b[nm]])
    if len(a_pts) < 3:
        return a_res
    a_pts = np.array(a_pts); b_pts = np.array(b_pts)
    a_c = a_pts.mean(0); b_c = b_pts.mean(0)
    H = (a_pts - a_c).T @ (b_pts - b_c)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    out = a_res.copy()
    out.coord = (a_res.coord - a_c) @ R.T + b_c
    return out


def _hungarian_pair_within_residue(
    a_res: struc.AtomArray, b_res: struc.AtomArray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue Hungarian assignment constrained to same element.
    Atom-to-atom pairing is determined in a per-residue-aligned frame (Kabsch
    on backbone N/CA/C) so that backbone displacement in the global frame
    doesn't confuse the matching. Caller uses the returned indices on the
    ORIGINAL (globally-aligned) coords to compute RMSD."""
    from scipy.optimize import linear_sum_assignment
    a_local = _kabsch_align_residue(a_res, b_res)
    a_idx_all = []; b_idx_all = []
    a_el = np.array([str(e) for e in a_res.element])
    b_el = np.array([str(e) for e in b_res.element])
    for el in np.unique(np.concatenate([a_el, b_el])):
        ai = np.where(a_el == el)[0]
        bi = np.where(b_el == el)[0]
        if len(ai) == 0 or len(bi) == 0:
            continue
        # Hungarian on per-residue-aligned coordinates for clean shape matching
        D = np.linalg.norm(
            a_local.coord[ai][:, None, :] - b_res.coord[bi][None, :, :], axis=2
        )
        r, c = linear_sum_assignment(D)
        a_idx_all.extend(ai[r].tolist())
        b_idx_all.extend(bi[c].tolist())
    return np.array(a_idx_all, dtype=int), np.array(b_idx_all, dtype=int)


def per_residue_max_displacement(
    a_aligned: struc.AtomArray, b: struc.AtomArray, res_ids: list[int]
) -> dict[int, float]:
    """Max heavy-atom displacement per residue (symmetry-aware via per-residue
    Hungarian matching within element). Skips residues with zero matched atoms."""
    out: dict[int, float] = {}
    for rid in res_ids:
        a_res = a_aligned[a_aligned.res_id == rid]
        b_res = b[b.res_id == rid]
        if len(a_res) == 0 or len(b_res) == 0:
            continue
        ai, bi = _hungarian_pair_within_residue(a_res, b_res)
        if len(ai) == 0:
            continue
        disp = np.linalg.norm(a_res.coord[ai] - b_res.coord[bi], axis=1)
        out[rid] = float(disp.max())
    return out


def find_moving_residues(
    a_aligned: struc.AtomArray,
    b: struc.AtomArray,
    pocket_res_ids: list[int],
    threshold: float = 2.0,
) -> list[int]:
    """Subset of pocket_res_ids whose max heavy-atom displacement exceeds `threshold` Å."""
    disps = per_residue_max_displacement(a_aligned, b, pocket_res_ids)
    return sorted(rid for rid, d in disps.items() if d > threshold)


def all_atom_rmsd_over_residues(
    a_aligned: struc.AtomArray, b: struc.AtomArray, res_ids: list[int]
) -> float:
    """All-heavy-atom RMSD over the given residues (symmetry-aware via per-residue
    Hungarian matching constrained to same element). Avoids inflating RMSD from
    symmetric atom-name swaps."""
    if not res_ids:
        return float("nan")
    sq = 0.0; n = 0
    for rid in res_ids:
        a_res = a_aligned[a_aligned.res_id == rid]
        b_res = b[b.res_id == rid]
        if len(a_res) == 0 or len(b_res) == 0:
            continue
        ai, bi = _hungarian_pair_within_residue(a_res, b_res)
        if len(ai) == 0:
            continue
        d = np.linalg.norm(a_res.coord[ai] - b_res.coord[bi], axis=1)
        sq += float((d ** 2).sum()); n += len(d)
    if n == 0:
        return float("nan")
    return float(np.sqrt(sq / n))


def ca_rmsd_over_residues(
    a_aligned: struc.AtomArray, b: struc.AtomArray, res_ids: list[int]
) -> float:
    """Cα RMSD over the given residues, between the already-aligned a and b."""
    if not res_ids:
        return float("nan")
    a_ca = a_aligned[(a_aligned.atom_name == "CA") & np.isin(a_aligned.res_id, res_ids)]
    b_ca = b[(b.atom_name == "CA") & np.isin(b.res_id, res_ids)]
    a_match, b_match = match_atoms(a_ca, b_ca)
    if len(a_match) == 0:
        return float("nan")
    return float(struc.rmsd(b_match, a_match))
