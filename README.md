# PROPKA3 pKa prediction

Run PROPKA3 on protein PDB files to predict per-residue pKa values.

This is designed for **comparing different preparations of the same protein** —
e.g. energy-minimised vs not, or different side-chain orientations/protonation
states of the same residues. By running each preparation through PROPKA and
lining the results up side by side, you can measure how those small structural
differences shift the predicted pKa values, and (with the per-residue analysis
below) deduce which factors — desolvation, hydrogen bonds, Coulombic
interactions — are driving each change.

## Setup

```sh
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

(Python 3.12 is required — PROPKA 3.5.1 is broken on 3.14. See docs.)

## Use

1. Add **PDB** structures to `input_files/`.
   (hydrogens are not saved in pdb files but prepping of proteins is still recommended as modelling the hydrogens in other software will reposition the amino acid residues.)

2. Run:

   ```sh
   python run_propka.py
   ```

3. Read the reports in `output_files/<name>.pka` — see the
   **SUMMARY OF THIS PREDICTION** section for predicted pKa per residue.
4. Compare structures in the two CSVs, each one row per input PDB:
   - `output_files/free_energy_summary.csv` — the pH of optimum stability, its
     folding free energy, and the full folding-free-energy-vs-pH curve
     (`dG_pH0`…`dG_pH14`, neutral reference).
   - `output_files/pka_summary.csv` — a pKa column per selected residue.

### Choosing which residues appear in the pKa CSV

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

This selection only affects `pka_summary.csv`; `free_energy_summary.csv` is
always written in full.

### Per-residue determinant breakdown

Set `detailed_analysis: true` in `config.yaml` to also write one CSV per selected
residue into `output_files/residue_pka/`. Each file has one row per structure and
columns for every factor PROPKA used to build that residue's pKa: buried
fraction, desolvation terms, and each sidechain/backbone hydrogen bond and
Coulombic interaction — so you can see exactly why the pKa differs between
preparations.

Because interaction partners differ between structures, each partner becomes its
own column (e.g. `coulombic_GLU62A`, `sidechain_hbond_UNKO1Z` for a ligand
contact). A structure that simply lacks a given interaction gets `0` in that
column (it neither raises nor lowers the pKa); a structure missing the residue
entirely gets blanks.

> **Warning:** with `positions: []`, turning this on writes a file for **every**
> polar residue in the protein. Choose your residues first, then enable it.

### PROPKA run options

PROPKA's own command-line options (pH, reference state, chain selection, …) are
set in the `propka:` block of `config.yaml`, not passed on the command line — so
every run is reproducible and logged from the config alone. The exact options
used are printed at the start of each run. Anything left `null`/blank uses
PROPKA's default:

```yaml
propka:
  pH: 7.0              # pH for stability & summary sections (pKa is pH-independent)
  reference: neutral   # neutral | low-pH
  chain: null          # restrict to a single chain id, e.g. A
  titrate_only: null   # only these residues titratable, e.g. "A:10,A:11"
  window: null         # stability-profile pH window [lo, hi, step]
  grid: null           # pH grid for stability properties [lo, hi, step]
  quiet: false         # suppress PROPKA's non-warning console output
  extra_args: []       # any other propka3 flag verbatim, e.g. ["--log-level", "WARNING"]
```

See [`docs/PROPKA.md`](docs/PROPKA.md) for what PROPKA does, why hydrogens don't
matter, and how to interpret the output.

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
