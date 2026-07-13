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

When detailed_analysis is enabled in config.yaml, one CSV per selected residue
is also written into output_files/residue_pka/, breaking each pKa down into the
determinants PROPKA used (desolvation, hydrogen bonds, Coulombic interactions).
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
RESIDUE_DIR = OUTPUT_DIR / "residue_pka"

# One determinant field: an energy value plus its partner group (NAME ID CHAIN),
# or a "XXX 0 X" placeholder when there is no interaction. ID is a residue number
# for amino acids but an atom name (e.g. "O1") for ligand groups, so it is not
# restricted to digits. Each determinant line ends with exactly three of these
# (sidechain H-bond, backbone H-bond, Coulombic), so findall yields those three.
DETERMINANT_RE = re.compile(r"([+-]?\d+\.\d+)\s+(XXX|[A-Za-z][\w+-]*)\s+(\S+)\s+(\S)")
INTERACTIONS = ("sidechain", "backbone", "coulombic")


def load_config():
    """Read config settings; tolerate a missing/empty config."""
    data = {}
    if CONFIG_FILE.exists():
        data = yaml.safe_load(CONFIG_FILE.read_text()) or {}
    residues = data.get("residues") or {}
    return {
        "positions": [int(p) for p in (residues.get("positions") or [])],
        "chains": [str(c) for c in (residues.get("chains") or [])],
        "detailed_analysis": bool(data.get("detailed_analysis", False)),
        "propka": data.get("propka") or {},
    }


def build_propka_args(popts):
    """Turn the config's `propka:` options into propka3 CLI arguments.

    Keeping these in the config (rather than hand-editing the command) means
    every run is reproducible and logged from config.yaml alone. Any option
    left null/blank falls back to PROPKA's own default.
    """
    args = []
    if popts.get("pH") is not None:
        args += ["--pH", str(popts["pH"])]
    if popts.get("reference"):
        args += ["--reference", str(popts["reference"])]
    if popts.get("chain"):
        args += ["--chain", str(popts["chain"])]
    if popts.get("titrate_only"):
        args += ["--titrate_only", str(popts["titrate_only"])]
    for opt in ("window", "grid"):
        if popts.get(opt):
            args += [f"--{opt}", *[str(v) for v in popts[opt]]]
    if popts.get("quiet"):
        args.append("--quiet")
    # Escape hatch for any flag not mapped above, recorded verbatim.
    args += [str(a) for a in (popts.get("extra_args") or [])]
    return args


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


def parse_determinants(pka_path):
    """Parse the detailed determinant table (how each pKa was built up).

    Returns {label: record}. Each record holds the desolvation terms and, for
    each interaction type, a {partner_label: energy} map. Partner labels differ
    between structures, which is why the CSV writer unions them into columns.
    """
    records = {}
    current = None
    in_table = False

    for line in pka_path.read_text().splitlines():
        if not in_table:
            # The determinant table opens with a "RESIDUE ... BURIED ..." header.
            if line.lstrip().startswith("RESIDUE") and "BURIED" in line:
                in_table = True
            continue

        stripped = line.strip()
        if not stripped:
            current = None            # blank line separates residue blocks
            continue
        if stripped.startswith(("Coupled", "SUMMARY")):
            break                     # end of the determinant table

        # Residue label sits in fixed columns regardless of continuation lines.
        # Dashed rule lines fall through here and are skipped by the guard below.
        resname, resnum, chain = line[0:3].strip(), line[3:8].strip(), line[8:9].strip()
        if not (resname and resnum.isdigit()):
            continue
        label = f"{resname}{resnum}{chain}"

        matches = DETERMINANT_RE.findall(line)
        if len(matches) != 3:
            continue                  # not a well-formed determinant line

        # A block's first line also carries pKa/buried/desolvation; continuation
        # lines start blank there and only add further interactions.
        head = line[9:DETERMINANT_RE.search(line).start()].replace("%", "").split()
        if head:
            current = {
                "resname": resname, "resnum": int(resnum), "chain": chain,
                "pKa": float(head[0].rstrip("*")),   # * marks a coupled residue
                "buried_pct": int(head[1]),
                "desolv_regular": float(head[2]), "desolv_regular_n": int(head[3]),
                "desolv_re": float(head[4]), "desolv_re_n": int(head[5]),
                "sidechain": {}, "backbone": {}, "coulombic": {},
            }
            records[label] = current
        if current is None:
            continue

        for slot, (val, pname, pnum, pchain) in zip(INTERACTIONS, matches):
            if pname == "XXX":
                continue              # placeholder: no interaction of this type
            partner = f"{pname}{pnum}{pchain}"
            current[slot][partner] = current[slot].get(partner, 0.0) + float(val)

    return records


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


def _partner_sort_key(partner):
    """Order interaction columns by the partner's residue number."""
    m = re.search(r"\d+", partner)
    return (int(m.group()) if m else 0, partner)


def write_residue_csvs(detailed, cfg):
    """Write one determinant-breakdown CSV per selected residue.

    `detailed` is an ordered list of (file_stem, {label: record}). Structures
    are rows. Interaction columns are the union of partners seen across
    structures: a structure that lacks a given interaction gets 0 (it neither
    raises nor lowers the pKa); a structure lacking the residue entirely gets
    blanks (nothing to say about a residue that is not there).
    """
    RESIDUE_DIR.mkdir(exist_ok=True)
    for stale in RESIDUE_DIR.glob("*.csv"):
        stale.unlink()

    # Residues to emit, applying the same position/chain selection as the pKa CSV.
    selected = {}                     # label -> resnum (for filename ordering)
    for _, records in detailed:
        for label, rec in records.items():
            if cfg["chains"] and rec["chain"] not in cfg["chains"]:
                continue
            if cfg["positions"] and rec["resnum"] not in cfg["positions"]:
                continue
            selected.setdefault(label, rec["resnum"])

    for label in selected:
        # Union of interaction partners across every structure that has this residue.
        partners = {slot: set() for slot in INTERACTIONS}
        for _, records in detailed:
            rec = records.get(label)
            if rec:
                for slot in INTERACTIONS:
                    partners[slot].update(rec[slot])
        for slot in INTERACTIONS:
            partners[slot] = sorted(partners[slot], key=_partner_sort_key)

        prefix = {"sidechain": "sidechain_hbond", "backbone": "backbone_hbond",
                  "coulombic": "coulombic"}
        fieldnames = (
            ["file", "pKa", "buried_pct", "desolv_regular", "desolv_regular_natoms",
             "desolv_re", "desolv_re_natoms"]
            + [f"{prefix[slot]}_{p}" for slot in INTERACTIONS for p in partners[slot]]
        )

        path = RESIDUE_DIR / (re.sub(r"[^A-Za-z0-9+_-]", "_", label) + ".csv")
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for stem, records in detailed:
                rec = records.get(label)
                if rec is None:
                    writer.writerow({"file": stem})   # residue absent: all blank
                    continue
                record = {
                    "file": stem,
                    "pKa": rec["pKa"],
                    "buried_pct": rec["buried_pct"],
                    "desolv_regular": rec["desolv_regular"],
                    "desolv_regular_natoms": rec["desolv_regular_n"],
                    "desolv_re": rec["desolv_re"],
                    "desolv_re_natoms": rec["desolv_re_n"],
                }
                for slot in INTERACTIONS:
                    for p in partners[slot]:
                        record[f"{prefix[slot]}_{p}"] = rec[slot].get(p, 0.0)
                writer.writerow(record)

    print(f"Wrote {len(selected)} per-residue determinant CSVs to {RESIDUE_DIR}/")


def main():
    cfg = load_config()
    rows = []              # per-structure records
    residue_labels = {}    # (chain, resnum) -> column label, union across files
    detailed = []          # (file_stem, {label: determinant record}), if enabled

    propka_args = build_propka_args(cfg["propka"])
    print(f"PROPKA options: {' '.join(propka_args) or '(defaults)'}\n")

    for pdb in sorted(INPUT_DIR.glob("*.pdb")):
        print(f"==> {pdb.name}")
        # propka3 writes <name>.pka into the current dir, so run it from output_files/
        subprocess.run(["propka3", *propka_args, pdb.resolve()], cwd=OUTPUT_DIR)

        pka_path = OUTPUT_DIR / f"{pdb.stem}.pka"
        if not pka_path.exists():
            print(f"    (no .pka produced for {pdb.name}, skipping)")
            continue

        residues, dg_curve, optimal_pH, optimal_dG = parse_pka(pka_path)
        residues = select_residues(residues, cfg)
        for key, info in residues.items():
            residue_labels.setdefault(key, info["label"])

        if cfg["detailed_analysis"]:
            detailed.append((pdb.stem, parse_determinants(pka_path)))

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

    if cfg["detailed_analysis"]:
        write_residue_csvs(detailed, cfg)


if __name__ == "__main__":
    main()
