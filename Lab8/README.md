# Needleman–Wunsch DNA Alignment Lab

Python 3 Tkinter desktop application. Implements the Needleman–Wunsch **global alignment algorithm from scratch**, with no BioPython and no external packages. The VB6 reference project is not copied or translated.

## Run

1. Install Python 3 with Tkinter (on Windows, standard Python from python.org includes Tkinter).
2. Open a terminal in this folder.
3. Run `python app.py` (or `py app.py` on Windows).

## Example

- S1: `ACCGTGAAGCCAATAC`
- S2: `AGCGTGCAGCCAATAC`
- Match = `+1`, mismatch = `-1`, gap = `-2`

Press **Align sequences**. The scoring matrix has 17 × 17 cells. Gold cells mark the optimal traceback path. Green is the current traceback step. Use **Next/Previous traceback step** to follow it, or select any step in the traceback tab.

## FASTA files

Use **Load FASTA…** next to either sequence. Supported: `.fna`, `.fasta`, `.fa`, `.fas`, and gzipped FASTA `.gz`. The program reads the first record, displays its header and length, and inserts a selected region into the editable sequence box. `Start` is 1-based, `Length` defaults to 120 bases. To load another region, change Start/Length and press **Load FASTA…** again. To compare the first 60 bases, enter Length = 60 before loading.

A maximum of 120 bases per input is enforced for a readable matrix and to avoid accidentally allocating a huge full-genome DP table. Full Needleman–Wunsch requires **O(n × m) time and memory** when displaying the complete scoring matrix. FASTA sequences may contain IUPAC ambiguity letters; only identical letters count as matches. `U` is normalized to `T`.

## Algorithm

1. Create score matrix `F` of `(len(S1)+1) × (len(S2)+1)`.
2. Initialize `F[0,0]=0`, `F[i,0]=i×gap`, `F[0,j]=j×gap`.
3. For each `(i,j)`, calculate:
   - diagonal: `F[i-1,j-1] + match` if characters equal, else `+ mismatch`
   - up: `F[i-1,j] + gap`
   - left: `F[i,j-1] + gap`
   - store the maximum and its direction
4. Trace back from bottom-right to `(0,0)`. Diagonal aligns both characters; up adds a gap to S2; left adds a gap to S1. Ties prefer diagonal, then up, then left.
5. Display the aligned strings and score.

## Files

- `app.py` — Tkinter interface, matrix visualization, traceback stepping, FASTA loading and export.
- `needleman_wunsch.py` — standalone alignment and FASTA parsing functions.
- `test_alignment.py` — standard-library unit tests.
- `480540219_genomic.fna` — supplied Influenza A FASTA.
- `GCF_000919155.1_ViralProj244477_genomic.fna.gz` — supplied compressed Duck picornavirus FASTA.
