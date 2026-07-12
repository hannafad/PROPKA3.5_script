#!/usr/bin/env python3
"""Run PROPKA3 on every PDB in input_files/ and write reports to output_files/.

PROPKA predicts per-residue pKa values from heavy-atom geometry, so hydrogens
in the input PDB are ignored. Each run produces a <name>.pka report.

Two CSVs are also written, each with one row per input structure (so structures
are directly comparable):
  - output_files/free_energy_summary.csv -- the folding free energy at optimum
    stability, the pH it occurs at, and the full folding-free-energy-vs-pH curve
    (pH 0-14, neutral reference).
  - output_files/pka_summary.csv -- a pKa column per residue selected in
    config.yaml.
"""
import csv
import re
import subprocess
from pathlib import Path

import yaml

INPUT_DIR = Path("input_files")
OUTPUT_DIR = Path("output_files")
CONFIG_FILE = Path("config.yaml")
ENERGY_CSV = OUTPUT_DIR / "free_energy_summary.csv"
PKA_CSV = OUTPUT_DIR / "pka_summary.csv"


def load_config():
    """Read residue-selection settings; tolerate a missing/empty config."""
    if not CONFIG_FILE.exists():
        return {"positions": [], "chains": []}
    data = yaml.safe_load(CONFIG_FILE.read_text()) or {}
    residues = data.get("residues") or {}
    return {
        "positions": [int(p) for p in (residues.get("positions") or [])],
        "chains": [str(c) for c in (residues.get("chains") or [])],
    }


def parse_pka(pka_path):
    """Extract the summary pKa table and folding free-energy data from a .pka."""
    residues = {}      # (chain, resnum) -> {"label", "pKa"}
    dg_curve = {}      # pH (float) -> folding free energy (kcal/mol)
    optimal_pH = None
    optimal_dG = None

    section = None
    for line in pka_path.read_text().splitlines():
        stripped = line.strip()

        if stripped.startswith("SUMMARY OF THIS PREDICTION"):
            section = "summary"
            continue
        if "as a function of pH (using neutral reference)" in stripped:
            section = "dg"
            continue

        if section == "summary":
            if stripped.startswith("Group") or not stripped:
                continue
            if stripped.startswith("-"):
                section = None
                continue
            parts = stripped.split()
            # Group is "RESNAME RESNUM CHAIN", then pKa, model-pKa, [ligand...]
            try:
                resname, resnum, chain, pka = parts[0], int(parts[1]), parts[2], float(parts[3])
            except (IndexError, ValueError):
                continue
            residues[(chain, resnum)] = {"label": f"{resname}{resnum}{chain}", "pKa": pka}

        elif section == "dg":
            if not stripped:
                section = None
                continue
            parts = stripped.split()
            try:
                dg_curve[float(parts[0])] = float(parts[1])
            except (IndexError, ValueError):
                section = None

        if "pH of optimum stability" in stripped:
            m = re.search(
                r"optimum stability is\s+([-\d.]+).*?free energy is\s+([-\d.]+)",
                stripped,
            )
            if m:
                optimal_pH = float(m.group(1))
                optimal_dG = float(m.group(2))

    return residues, dg_curve, optimal_pH, optimal_dG


def select_residues(residues, cfg):
    """Filter a structure's residues by the configured positions/chains."""
    positions, chains = cfg["positions"], cfg["chains"]
    out = {}
    for (chain, resnum), info in residues.items():
        if chains and chain not in chains:
            continue
        if positions and resnum not in positions:
            continue
        out[(chain, resnum)] = info
    return out


def main():
    cfg = load_config()
    rows = []              # per-structure records
    residue_labels = {}    # (chain, resnum) -> column label, union across files

    for pdb in sorted(INPUT_DIR.glob("*.pdb")):
        print(f"==> {pdb.name}")
        # propka3 writes <name>.pka into the current dir, so run it from output_files/
        subprocess.run(["propka3", pdb.resolve()], cwd=OUTPUT_DIR)

        pka_path = OUTPUT_DIR / f"{pdb.stem}.pka"
        if not pka_path.exists():
            print(f"    (no .pka produced for {pdb.name}, skipping)")
            continue

        residues, dg_curve, optimal_pH, optimal_dG = parse_pka(pka_path)
        residues = select_residues(residues, cfg)
        for key, info in residues.items():
            residue_labels.setdefault(key, info["label"])

        rows.append({
            "file": pdb.stem,
            "optimal_pH": optimal_pH,
            "optimal_pH_dG_kcal_mol": optimal_dG,
            "dg_curve": dg_curve,
            "residues": {key: info["pKa"] for key, info in residues.items()},
        })

    if not rows:
        print("No structures processed; no CSV written.")
        return

    # Stable column order across structures.
    pH_points = sorted({pH for row in rows for pH in row["dg_curve"] if pH.is_integer()})
    residue_keys = sorted(residue_labels)

    # Free-energy CSV: optimum-stability pH/dG plus the full dG-vs-pH curve.
    energy_fields = (
        ["file", "optimal_pH", "optimal_pH_dG_kcal_mol"]
        + [f"dG_pH{int(pH)}" for pH in pH_points]
    )
    with ENERGY_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=energy_fields)
        writer.writeheader()
        for row in rows:
            record = {
                "file": row["file"],
                "optimal_pH": row["optimal_pH"],
                "optimal_pH_dG_kcal_mol": row["optimal_pH_dG_kcal_mol"],
            }
            for pH in pH_points:
                record[f"dG_pH{int(pH)}"] = row["dg_curve"].get(pH)
            writer.writerow(record)

    # pKa CSV: one pKa column per selected residue.
    pka_fields = ["file"] + [residue_labels[key] for key in residue_keys]
    with PKA_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=pka_fields)
        writer.writeheader()
        for row in rows:
            record = {"file": row["file"]}
            for key in residue_keys:
                record[residue_labels[key]] = row["residues"].get(key)
            writer.writerow(record)

    print(f"\nWrote {ENERGY_CSV} and {PKA_CSV} "
          f"({len(rows)} structures, {len(residue_keys)} residue columns).")


if __name__ == "__main__":
    main()
