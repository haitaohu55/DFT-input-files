#!/usr/bin/env python3
"""Extract g_cv,nu(q) for the VBM -> CBM transition from vaspelph.h5."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = Path(__file__).resolve().parent
H5_FILE = ROOT / "7_h5" / "OUTPUT" / "vaspelph.h5"
TABLE_FILE = BASE / "gcv_table.dat"
SUMMARY_FILE = BASE / "gcv_summary.txt"

# Transition requested in the notes:
# initial VBM: kv, final CBM: kc, q = kc - kv.
KV = np.array([0.4651162791, 0.5, 0.0])
KC = np.array([0.0, 0.0, 0.0])

VBAND = 56  # VASP 1-based VBM band index.
CBAND = 57  # VASP 1-based CBM band index.
SPIN = 0
TEMPERATURE = 300.0
KB_EV_PER_K = 8.617333262145e-5


def delta_mod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return a - b folded into [-0.5, 0.5)."""
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return d - np.round(d)


def fold_01(x: np.ndarray) -> np.ndarray:
    """Fold reduced coordinates into [0, 1)."""
    return np.asarray(x, dtype=float) - np.floor(np.asarray(x, dtype=float))


def nearest_kpoint(kpts: np.ndarray, target: np.ndarray) -> tuple[int, float]:
    distances = np.linalg.norm(delta_mod(kpts, target), axis=1)
    index = int(np.argmin(distances))
    return index, float(distances[index])


def local_band_index(vasp_band_1based: int, band_start: int, nbands: int) -> int:
    """Map a VASP 1-based band index into the local HDF5 band dimension."""
    index = vasp_band_1based - int(band_start)
    if not 0 <= index < nbands:
        raise ValueError(
            f"Cannot map VASP band {vasp_band_1based} with "
            f"band_start={band_start}, nbands={nbands}."
        )
    return index


def bose_occupation(phonon_e: np.ndarray, temperature: float) -> np.ndarray:
    occ = np.full_like(phonon_e, np.nan, dtype=float)
    mask = phonon_e > 1e-12
    occ[mask] = 1.0 / np.expm1(phonon_e[mask] / (KB_EV_PER_K * temperature))
    return occ


def main() -> None:
    with h5py.File(H5_FILE, "r") as h5:
        vk = np.array(h5["kpoints/vkpt_k"])
        vkp = np.array(h5["kpoints/vkpt_kp"])

        elph = h5["matrix_elements/elph"]
        phonon_e = np.array(h5["matrix_elements/phonon_eigenvalues"])

        _, nkp, nk, nmodes, nbands_kp, nbands_k, _ = elph.shape
        band_start_k = int(np.ravel(h5["matrix_elements/band_start_k"][()])[0])
        band_start_kp = int(np.ravel(h5["matrix_elements/band_start_kp"][()])[0])

        ik_initial, dk_initial = nearest_kpoint(vk, KV)
        ikp_final, dk_final = nearest_kpoint(vkp, KC)

        ib_v = local_band_index(VBAND, band_start_k, nbands_k)
        ib_c = local_band_index(CBAND, band_start_kp, nbands_kp)

        g_real = elph[SPIN, ikp_final, ik_initial, :, ib_c, ib_v, 0]
        g_imag = elph[SPIN, ikp_final, ik_initial, :, ib_c, ib_v, 1]
        g = np.array(g_real) + 1j * np.array(g_imag)

        # VASP stores these phonon eigenvalues here as meV-scale values for this run.
        # Convert to eV for Bose factors while keeping the raw meV column in the table.
        phonon_mev = phonon_e[ikp_final, ik_initial, :]
        eph = phonon_mev / 1000.0
        nph = bose_occupation(eph, TEMPERATURE)
        g_abs = np.abs(g)
        g2 = g_abs**2
        absorption_weight = g2 * nph
        emission_weight = g2 * (nph + 1.0)

        q = delta_mod(vkp[ikp_final], vk[ik_initial])
        q_folded = fold_01(q)

    header_lines = [
        "# g_cv,nu(q) extracted from vaspelph.h5",
        f"# H5_FILE = {H5_FILE}",
        f"# initial VBM kv = {vk[ik_initial].tolist()}, target = {KV.tolist()}, distance = {dk_initial:.6e}",
        f"# final   CBM kc = {vkp[ikp_final].tolist()}, target = {KC.tolist()}, distance = {dk_final:.6e}",
        f"# q = kc - kv folded [-0.5,0.5) = {q.tolist()}",
        f"# q folded [0,1) = {q_folded.tolist()}",
        f"# VBAND = {VBAND}, CBAND = {CBAND}, local ib_v = {ib_v}, local ib_c = {ib_c}",
        f"# elph dimensions used: nkp = {nkp}, nk = {nk}, nmodes = {nmodes}, nbands_kp = {nbands_kp}, nbands_k = {nbands_k}",
        f"# temperature_K = {TEMPERATURE}",
        "# columns: nu phonon_meV phonon_eV Nph_300K Re_g_eV Im_g_eV abs_g_eV abs_g2_eV2 abs_weight_eV2 em_weight_eV2",
    ]

    rows = np.column_stack(
        [
            np.arange(1, nmodes + 1),
            phonon_mev,
            eph,
            nph,
            g.real,
            g.imag,
            g_abs,
            g2,
            absorption_weight,
            emission_weight,
        ]
    )

    np.savetxt(
        TABLE_FILE,
        rows,
        header="\n".join(line[2:] for line in header_lines)
        + "\nnu phonon_meV phonon_eV Nph_300K Re_g_eV Im_g_eV abs_g_eV abs_g2_eV2 abs_weight_eV2 em_weight_eV2",
        comments="# ",
        fmt=["%3d"] + ["%.12e"] * 9,
    )

    max_g2 = int(np.nanargmax(g2))
    max_abs = int(np.nanargmax(absorption_weight))
    max_em = int(np.nanargmax(emission_weight))
    summary = [
        "g_cv,nu(q) summary",
        f"q [-0.5,0.5) = {q}",
        f"q [0,1) = {q_folded}",
        f"initial VBM k = {vk[ik_initial]}",
        f"final CBM k = {vkp[ikp_final]}",
        f"VBAND = {VBAND}, CBAND = {CBAND}",
        "",
        f"max |g|^2: nu={max_g2 + 1}, phonon_meV={phonon_mev[max_g2]:.6f}, |g|={g_abs[max_g2]:.12e} eV, |g|^2={g2[max_g2]:.12e} eV^2",
        f"max absorption |g|^2*Nph at {TEMPERATURE:.0f} K: nu={max_abs + 1}, weight={absorption_weight[max_abs]:.12e} eV^2",
        f"max emission |g|^2*(Nph+1) at {TEMPERATURE:.0f} K: nu={max_em + 1}, weight={emission_weight[max_em]:.12e} eV^2",
        "",
        f"table = {TABLE_FILE}",
    ]
    SUMMARY_FILE.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print("\n".join(header_lines))
    print(f"Wrote {TABLE_FILE}")
    print(f"Wrote {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
