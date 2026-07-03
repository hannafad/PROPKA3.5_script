#!/usr/bin/env python3
"""Run PROPKA3 on every PDB in input_files/ and write reports to output_files/.

PROPKA predicts per-residue pKa values from heavy-atom geometry, so hydrogens
in the input PDB are ignored. Each run produces a <name>.pka report.
"""
import subprocess
from pathlib import Path

INPUT_DIR = Path("input_files")
OUTPUT_DIR = Path("output_files")

for pdb in sorted(INPUT_DIR.glob("*.pdb")):
    print(f"==> {pdb.name}")
    # propka3 writes <name>.pka into the current dir, so run it from output_files/
    subprocess.run(["propka3", pdb.resolve()], cwd=OUTPUT_DIR)
