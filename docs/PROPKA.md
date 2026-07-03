# PROPKA3 — how it's used here

## What PROPKA does

PROPKA3 predicts the **pKa values of ionizable groups** (ASP, GLU, HIS, CYS,
TYR, LYS, ARG, plus N- and C-termini and, optionally, ligand groups) in a
protein structure. It's an *empirical, structure-based* predictor: given the 3D
coordinates it estimates how each group's pKa is perturbed from its intrinsic
"model" value by its local environment (desolvation, hydrogen bonds, charge–
charge interactions).

Reference: Olsson, Søndergaard, Rostkowski & Jensen, *J. Chem. Theory Comput.*
2011 (PROPKA3). Package: <https://github.com/jensengroup/propka>.

## The key facts for this project

- **PROPKA uses heavy atoms only.** It ignores hydrogens in the input and, by
  default, strips and re-derives any it needs. So a PDB **without** hydrogens
  gives the *same* prediction as one with them — the lack of H atoms in a plain
  PDB is not a problem here.
- **Input is PDB (or mmCIF), not Maestro.** The open-source `propka3` cannot
  read `.mae`. Export your Maestro structures to PDB first (Maestro:
  *File → Export Structures → .pdb*, or `$SCHRODINGER/utilities/structconvert
  in.mae out.pdb`) and drop the `.pdb` files into `input_files/`.
- **It's non-destructive analysis.** PROPKA3 does not re-protonate or otherwise
  modify your structure and does not read the protonation states Maestro
  assigned — it independently *predicts* pKa from the coordinates. Your
  Maestro-prepped `.mae` files are untouched. (Maestro's Protein Prep Wizard
  runs PROPKA internally, so this is effectively an independent cross-check of
  the same algorithm.)
- If you instead need **re-protonated structures at a given pH**, PROPKA3 alone
  won't do it — that's a job for PDB2PQR (`pdb2pqr --with-ph 7
  --titration-state-method propka ...`). Out of scope for this project.

## Running it

```sh
source .venv/bin/activate      # env created with:  uv venv --python 3.12 .venv
                               #                    uv pip install -r requirements.txt
python run_propka.py           # processes every *.pdb in input_files/
```

Each structure produces `output_files/<name>.pka`.

> Note: PROPKA 3.5.1 is incompatible with Python 3.14 (a PEP 649 annotation
> change breaks parameter-file parsing). The env is pinned to Python 3.12 via
> `.python-version`.

## Reading the `.pka` report

The most useful section is **SUMMARY OF THIS PREDICTION**:

```
       Group      pKa  model-pKa   ligand atom-type
   ASP  21 A     3.36       3.80
   HIS  68 A     6.00       6.50
   LYS  11 A    11.02      10.50
   ...
```

- **Group** — residue type, number, chain.
- **pKa** — the *predicted* pKa in this structure. This is the number you want.
- **model-pKa** — the intrinsic reference pKa of that group type; the shift from
  it tells you how much the environment perturbs the group.

Interpretation at a chosen pH: a group is predominantly **protonated** when
pH < pKa and **deprotonated** when pH > pKa. So at pH 7, an ASP/GLU (pKa ~4) is
deprotonated (–1), a LYS/ARG (pKa ~10–12) is protonated (+1), and a HIS near
pKa 6–7 is ambiguous and worth inspecting.

Above the summary, PROPKA lists the **determinants** for each group (the
specific interactions driving each pKa shift). Below it are pH-dependent
**folding free-energy** and **net-charge / isoelectric-point** profiles,
controlled by the `pH`, `reference`, `window` and `grid` options.
