# Bioinformatics - Lab 5 (Python)

Two independent exercises, based on the lab tasks, using the **original FASTA filenames** you supplied. No C++ compiler needed.

## Run on Windows (CMD or VS Code terminal)

Open this folder in the terminal and execute:

```cmd
py -m pip install -r requirements.txt
py ex1.py
py ex2.py
```

Results are written to `outputs/`:

- `gel_electrophoresis.png` - 10 DNA sample lanes and a reference ladder
- `gel_results.txt` - fragment lengths, estimated molecular weights and migration distances
- `motif_frequencies.txt` - overlapping 3-, 4- and 5-base patterns occurring 6-10 times

The results are **already generated and included** in this ZIP.

## Input files (included)

- `data/480540219_genomic.fna` - Influenza A sequence; used by `ex1.py`
- `data/GCF_000919155.1_ViralProj244477_genomic.fna.gz` - Duck picornavirus sequence; used by `ex2.py`

The files are loaded relative to the Python scripts, so running the programs from any current directory works. To swap inputs, use `py ex1.py --fasta path\to\file.fna` or `py ex2.py --fasta path\to\file.fna.gz`.

## Exercise 1

A seeded random generator selects 10 fragments from the first 2,000 DNA bases. Fragment molecular weight is estimated as `617.96 * bp + 36.04` Da. Migration distance is modeled as a logarithmic function of fragment size across a 10 cm gel. This is an **educational simulation**, not a physical gel measurement.

## Exercise 2

The program counts overlapping 3-, 4-, and 5-nucleotide substrings across the whole input DNA and selects those occurring between 6 and 10 times inclusive. Ambiguous `N`-containing substrings are excluded. Since the input sequence differs from your colleague's, the text results will also differ.

## Optional checks

```cmd
py -m unittest -v test_lab5
```

## GitHub

Upload the entire folder (including `data` and `outputs`). If your repository is already set up:

```cmd
git add .
git commit -m "Complete Bioinformatics Lab 5 in Python"
git push
```
