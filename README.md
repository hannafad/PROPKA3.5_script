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

See [`docs/PROPKA.md`](docs/PROPKA.md) for what PROPKA does, why hydrogens don't
matter, and how to interpret the output.
