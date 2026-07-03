# PROPKA3 CLI flags and how to add them to the script

The runner (`run_propka.py`) calls `propka3` on each PDB with no extra options:

```python
subprocess.run(["propka3", pdb.resolve()], cwd=OUTPUT_DIR)
```

To change PROPKA's behaviour, add flags to that list **before** `pdb.resolve()`.
Each flag (and its value, if any) is a separate string. Example — predict at
pH 7.4, using the low-pH reference state, and quietly:

```python
subprocess.run(
    ["propka3", "--pH", "7.4", "-r", "low-pH", "-q", pdb.resolve()],
    cwd=OUTPUT_DIR,
)
```

The flags apply to every structure in the run, since the same command is used
for each file in the loop.

## Useful flags

| Flag | Argument | Default | What it does |
|------|----------|---------|--------------|
| `-o`, `--pH` | float | `7.0` | pH used for the stability (folding free-energy) and summary sections. **Does not change the predicted pKa values** — those are pH-independent. |
| `-r`, `--reference` | `neutral` \| `low-pH` | `neutral` | Reference state for the stability calculation. |
| `-c`, `--chain` | chain id(s) | all | Build the protein from only the named chain(s). Use `" "` for a blank chain id. |
| `-i`, `--titrate_only` | `"A:10,A:11"` | all | Treat only these residues as titratable (comma-separated `chain:resnum`). |
| `-g`, `--grid` | `lo hi step` | `0.0 14.0 0.1` | pH grid over which stability/charge properties are computed. |
| `-w`, `--window` | `lo hi step` | `0.0 14.0 1.0` | pH window shown in the stability profile output. |
| `-d`, `--display-coupled-residues` | — | off | Also report alternative pKa values arising from coupling between titratable groups. |
| `-k`, `--keep-protons` | — | off | Keep the protons already in the input PDB instead of stripping/repositioning. Rarely needed (PROPKA ignores H by default). |
| `--protonate-all` | — | off | Protonate all atoms. Does **not** affect the pKa calculation; only the internal structure handling. |
| `-q`, `--quiet` | — | off | Suppress non-warning console messages. |
| `-p`, `--parameters` | path | bundled `propka.cfg` | Use a custom parameter file. |
| `--log-level` | `DEBUG`…`CRITICAL` | `INFO` | Console verbosity. |

Multi-value flags like `--grid` take three separate strings, e.g.
`["--grid", "0.0", "14.0", "0.5"]`.

## Full list

For everything (including advanced mutation/alignment options), run:

```sh
propka3 --help
```
