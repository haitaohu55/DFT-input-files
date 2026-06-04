#!/usr/bin/env python3
"""
Generate VASP KPOINTS around CBM and VBM for effective-mass fitting.
Coordinates for CBM/VBM are VASP reciprocal fractional coordinates.
The k-space offsets are Cartesian, in Angstrom^{-1}; this is what the parabolic fit uses.
"""
import argparse, json, itertools
from pathlib import Path
import numpy as np


def read_poscar_lattice(path):
    lines = Path(path).read_text().splitlines()
    scale = float(lines[1].split()[0])
    A = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)], dtype=float) * scale
    # Rows are a1, a2, a3. Reciprocal rows b_i satisfy a_i dot b_j = 2*pi delta_ij.
    B = 2 * np.pi * np.linalg.inv(A).T
    return A, B


def frac_offset_from_cart_dk(B_rows, dk_cart):
    # dk_cart = sum_i df_i b_i = B_rows.T @ df_col
    return np.linalg.solve(B_rows.T, np.asarray(dk_cart, dtype=float))


def build_cloud(k0_frac, h, include_cross=True):
    vals = [-h, 0.0, h]
    cloud = []
    for dx, dy, dz in itertools.product(vals, vals, vals):
        cloud.append(np.array([dx, dy, dz], dtype=float))
    return cloud


def write_kpoints(out, points):
    with open(out, "w") as f:
        f.write("Effective-mass point cloud around CBM and VBM\n")
        f.write(f"{len(points)}\n")
        f.write("Reciprocal\n")
        for p in points:
            k = p["k_frac"]
            f.write(f"{k[0]: .12f} {k[1]: .12f} {k[2]: .12f} 1\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poscar", default="POSCAR")
    ap.add_argument("--h", type=float, default=0.01, help="Cartesian k offset in Angstrom^-1")
    ap.add_argument("--cbm", nargs=3, type=float, default=[0.0, 0.0, 0.0], help="CBM reciprocal fractional coordinate")
    ap.add_argument("--vbm", nargs=3, type=float, default=[0.42, 0.5, 0.0], help="VBM reciprocal fractional coordinate")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    A, B = read_poscar_lattice(args.poscar)
    cbm = np.array(args.cbm, dtype=float)
    vbm = np.array(args.vbm, dtype=float)
    h = args.h

    points = []
    for group, k0 in [("CBM", cbm), ("VBM", vbm)]:
        for dk in build_cloud(k0, h):
            df = frac_offset_from_cart_dk(B, dk)
            k = k0 + df
            label = f"dk=({dk[0]:+.5f},{dk[1]:+.5f},{dk[2]:+.5f})A^-1"
            points.append({
                "group": group,
                "k0_frac": k0.tolist(),
                "dk_cart_Ainv": dk.tolist(),
                "df_frac": df.tolist(),
                "k_frac": k.tolist(),
                "label": label,
            })

    stem = args.out or f"mass_h{h:.4f}".replace(".", "p")
    kp_file = f"{stem}"
    meta_file = f"{stem}.json"
    write_kpoints(kp_file, points)
    meta = {
        "poscar": str(args.poscar),
        "h_Ainv": h,
        "lattice_rows_A": A.tolist(),
        "reciprocal_rows_Ainv": B.tolist(),
        "cbm_frac": cbm.tolist(),
        "vbm_frac": vbm.tolist(),
        "points": points,
    }
    with open(meta_file, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote {kp_file} with {len(points)} explicit k-points")
    print(f"Wrote {meta_file}")
    print("Reciprocal vector lengths |b_i| [A^-1]:", np.linalg.norm(B, axis=1))
    print("Fractional offsets for +/-h along Cartesian x/y/z:")
    for name, dk in [("x", [h,0,0]), ("y", [0,h,0]), ("z", [0,0,h])]:
        print(name, frac_offset_from_cart_dk(B, dk))


if __name__ == "__main__":
    main()
