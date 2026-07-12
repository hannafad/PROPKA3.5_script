# PROPKA3 pKa prediction

Run PROPKA3 on protein PDB files to predict per-residue pKa values.

## Setup

```sh
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

(Python 3.12 is required — PROPKA 3.5.1 is broken on 3.14. See docs.)

## Use

1. Export your Maestro structures to **PDB** and place them in `input_files/`.
   (PROPKA can't read `.mae`; hydrogens are optional — it ignores them.)
2. Run:

   ```sh
   python run_propka.py
   ```

3. Read the reports in `output_files/<name>.pka` — see the
   **SUMMARY OF THIS PREDICTION** section for predicted pKa per residue.
4. Compare structures in `output_files/pka_summary.csv` — one row per input
   PDB, with the pH of optimum stability and its folding free energy, the full
   folding-free-energy-vs-pH curve (`dG_pH0`…`dG_pH14`, neutral reference), and
   a pKa column per selected residue.

### Choosing which residues appear in the CSV

`config.yaml` controls the per-residue pKa columns. It is not tracked in git —
copy the template once and edit your local copy:

```sh
cp config.yaml.template config.yaml
```

(If `config.yaml` is absent the script still runs, defaulting to all residues.)
Leave `positions: []` to include every titratable residue, or list residue
numbers to compare only those sites across structures:

```yaml
residues:
  positions: [24, 56, 207]   # residue numbers; [] = all residues
  chains: []                 # restrict to chain IDs, e.g. [A]; [] = all chains
```

The free-energy columns are always written regardless of this selection.

See [`docs/PROPKA.md`](docs/PROPKA.md) for what PROPKA does, why hydrogens don't
matter, and how to interpret the output, and
[`docs/CLI_FLAGS.md`](docs/CLI_FLAGS.md) for the available PROPKA3 flags (pH,
reference state, chain selection, …) and how to add them to the script.

## References / Citations

These scripts utilise PROPKA 3.5.1, which was developed by:
- Søndergaard, Chresten R., Mats H. M. Olsson, Michał Rostkowski, and Jan H.
  Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical
  Calculation and Rationalization of pKa Values." *Journal of Chemical Theory
  and Computation* 7, no. 7 (2011): 2284–2295.
- Olsson, Mats H. M., Chresten R. Søndergaard, Michał Rostkowski, and Jan H.
  Jensen. "PROPKA3: Consistent Treatment of Internal and Surface Residues in
  Empirical pKa Predictions." *Journal of Chemical Theory and Computation* 7,
  no. 2 (2011): 525–537.
