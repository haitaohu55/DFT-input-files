#!/usr/bin/env python3
"""Calculate indirect BTBT A prefactor from extracted g_cv,nu(q)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = Path(__file__).resolve().parent
GC_TABLE = OUTDIR / "gcv_table.dat"
POSCAR = ROOT / "7_h5" / "POSCAR"
MEFF_REPORT = ROOT / "5_meff" / "OUTPUT" / "meff"

TEMPERATURE_K = 300.0
A_CONSTANT = 5.3464012e-4
EV_TO_J = 1.602176634e-19
HBAR_J_S = 1.054571817e-34
AMU_TO_KG = 1.66053906660e-27

# Atomic masses in atomic mass units. POSCAR for this run contains O and V.
ATOMIC_MASS_AMU = {
    "O": 15.9994,
    "V": 50.9415,
}


def parse_poscar_mass(poscar: Path) -> tuple[float, float, dict[str, int]]:
    lines = poscar.read_text().splitlines()
    scale = float(lines[1].split()[0])
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    lattice *= scale
    volume_a3 = abs(float(np.linalg.det(lattice)))

    symbols = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(symbols) != len(counts):
        raise ValueError(f"Cannot parse symbols/counts from {poscar}")

    composition = dict(zip(symbols, counts))
    missing = sorted(set(composition) - set(ATOMIC_MASS_AMU))
    if missing:
        raise ValueError(f"Missing atomic masses for {missing}")

    mass_amu = sum(ATOMIC_MASS_AMU[sym] * n for sym, n in composition.items())
    return mass_amu, volume_a3, composition


def read_f0_candidate(meff_report: Path) -> float | None:
    if not meff_report.exists():
        return None
    text = meff_report.read_text(errors="replace")
    match = re.search(r"m_r_parallel/m0\s*=\s*([0-9.+\-Ee]+)", text)
    if match is None:
        return None
    return float(match.group(1))


def main() -> None:
    data = np.loadtxt(GC_TABLE)
    if data.ndim != 2 or data.shape[1] < 10:
        raise ValueError(f"Unexpected table shape in {GC_TABLE}: {data.shape}")

    nu = data[:, 0].astype(int)
    phonon_mev = data[:, 1]
    phonon_ev = data[:, 2]
    nph = data[:, 3]
    re_g = data[:, 4]
    im_g = data[:, 5]
    abs_g_ev = data[:, 6]
    abs_g2_ev2 = data[:, 7]

    if len(nu) != 42:
        raise ValueError(f"Expected 42 phonon branches for 14 atoms, got {len(nu)}")
    if np.any(phonon_ev <= 0):
        raise ValueError("All phonon energies must be positive for A_ind")
    if np.any(~np.isfinite(data)):
        raise ValueError("Input gcv table contains non-finite values")

    mass_amu, volume_a3, composition = parse_poscar_mass(POSCAR)
    mass_kg = mass_amu * AMU_TO_KG

    # Convention used here:
    #   g = D * sqrt(hbar / (2 M_cell omega))
    # so D = |g| * sqrt(2 M_cell epsilon_ph) / hbar.
    # VASP g is in eV, epsilon_ph is in eV, and D is reported in eV/cm.
    d_j_per_m = abs_g_ev * EV_TO_J * np.sqrt(2.0 * mass_kg * phonon_ev * EV_TO_J) / HBAR_J_S
    d_ev_per_cm = (d_j_per_m / EV_TO_J) / 100.0

    bose_factor = 1.0 + 2.0 * nph
    d2_over_e = d_ev_per_cm**2 / phonon_ev
    a_coeff = A_CONSTANT * bose_factor * d2_over_e
    total_coeff = float(np.sum(a_coeff))

    max_a_idx = int(np.argmax(a_coeff))
    max_g_idx = int(np.argmax(abs_g2_ev2))

    out_table = OUTDIR / "A_ind_table.dat"
    header = "\n".join(
        [
            "A_ind mode table for VBM -> CBM intervalley transition",
            f"source_gcv_table = {GC_TABLE}",
            f"poscar = {POSCAR}",
            f"temperature_K = {TEMPERATURE_K:.1f}",
            f"composition = {composition}",
            f"cell_mass_amu = {mass_amu:.8f}",
            f"cell_mass_kg = {mass_kg:.12e}",
            f"cell_volume_A3 = {volume_a3:.12f}",
            "D convention: D = |g| * sqrt(2*M_cell*epsilon_ph) / hbar",
            "A_ind,nu = 5.3464012e-4 * (1+2Nph) * D^2/epsilon_ph * F0^(5/2)",
            "A_coeff_per_F0_5_2 means A_ind,nu / F0^(5/2)",
            "columns: nu phonon_meV phonon_eV Nph_300K Re_g_eV Im_g_eV abs_g_eV abs_g2_eV2 D_eV_cm bose_1p2N D2_over_eV A_coeff_per_F0_5_2 frac_of_total",
        ]
    )
    table = np.column_stack(
        [
            nu,
            phonon_mev,
            phonon_ev,
            nph,
            re_g,
            im_g,
            abs_g_ev,
            abs_g2_ev2,
            d_ev_per_cm,
            bose_factor,
            d2_over_e,
            a_coeff,
            a_coeff / total_coeff,
        ]
    )
    np.savetxt(
        out_table,
        table,
        fmt=[
            "%3d",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
            "%.12e",
        ],
        header=header,
    )

    f0_candidate = read_f0_candidate(MEFF_REPORT)
    f0_lines: list[str] = []
    if f0_candidate is not None:
        f0_factor = f0_candidate**2.5
        f0_lines = [
            "",
            "Optional numerical substitution:",
            f"F0_candidate_from_5_meff_m_r_parallel = {f0_candidate:.12g}",
            f"F0_candidate^(5/2) = {f0_factor:.12e}",
            f"A_total_if_F0_candidate = {total_coeff * f0_factor:.12e}",
            f"A_max_branch_if_F0_candidate = {a_coeff[max_a_idx] * f0_factor:.12e}",
        ]

    summary_lines = [
        "A_ind summary for VBM -> CBM intervalley transition",
        f"source_gcv_table = {GC_TABLE}",
        f"temperature_K = {TEMPERATURE_K:.1f}",
        f"phonon_branches = {len(nu)}",
        f"composition = {composition}",
        f"cell_mass_amu = {mass_amu:.8f}",
        f"cell_mass_kg = {mass_kg:.12e}",
        f"cell_volume_A3 = {volume_a3:.12f}",
        "",
        "Formula evaluated:",
        "A_ind,total = 5.3464012e-4 * F0^(5/2) * sum_nu[(1+2N_nu)*D_nu^2/epsilon_nu]",
        "D_nu(eV/cm) = |g_nu(eV)| * sqrt(2*M_cell*epsilon_nu(J)) / hbar / 100",
        "",
        f"sum_all_branches_A_coeff_per_F0_5_2 = {total_coeff:.12e}",
        (
            "max_A_contribution_branch = "
            f"nu={nu[max_a_idx]}, phonon_meV={phonon_mev[max_a_idx]:.12f}, "
            f"Nph={nph[max_a_idx]:.12e}, |g|={abs_g_ev[max_a_idx]:.12e} eV, "
            f"D={d_ev_per_cm[max_a_idx]:.12e} eV/cm, "
            f"A_coeff/F0^(5/2)={a_coeff[max_a_idx]:.12e}, "
            f"fraction={a_coeff[max_a_idx] / total_coeff:.12e}"
        ),
        (
            "max_abs_g_branch = "
            f"nu={nu[max_g_idx]}, phonon_meV={phonon_mev[max_g_idx]:.12f}, "
            f"Nph={nph[max_g_idx]:.12e}, |g|={abs_g_ev[max_g_idx]:.12e} eV, "
            f"D={d_ev_per_cm[max_g_idx]:.12e} eV/cm, "
            f"A_coeff/F0^(5/2)={a_coeff[max_g_idx]:.12e}, "
            f"fraction={a_coeff[max_g_idx] / total_coeff:.12e}"
        ),
        "",
        "Top branches by A contribution:",
    ]
    for idx in np.argsort(a_coeff)[::-1][:10]:
        summary_lines.append(
            f"nu={nu[idx]:2d}  phonon_meV={phonon_mev[idx]:12.6f}  "
            f"Nph={nph[idx]:.6e}  |g|={abs_g_ev[idx]:.6e} eV  "
            f"D={d_ev_per_cm[idx]:.6e} eV/cm  "
            f"A_coeff/F0^(5/2)={a_coeff[idx]:.6e}  frac={a_coeff[idx] / total_coeff:.6e}"
        )

    summary_lines.extend(f0_lines)
    summary_lines.extend(
        [
            "",
            f"table = {out_table}",
        ]
    )

    summary = OUTDIR / "A_ind_summary.txt"
    summary.write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
