
# pspA Detection & Classification Pipeline

A reproducible pipeline for detecting **pspA** from genome assemblies, rescuing missing cases from short reads, classifying BoxB by Hollingshead et al 2000 pspA family defenition, and (optionally) building trees.

This repository contains five scripts that work together:

- `orchestrate_pspA_workflow.py` — end-to-end orchestration of de novo detection, rescue, BoxB classification, and optional tree.
- `denovo_pspA_detect.py` — BLASTX-based de novo detection with 200‑bp upstream validation. 
- `core_boxb_classifier.py` — per-sequence BoxB classifier (BLASTX vs BoxB protein references). 
- `rescue_refguided.py` — reference-guided rescue using minimap2 + samtools; trims the upstream 200 nt for rescued candidates.
- `pspA_utils.py` — shared helpers (BLAST DB checks, ID normalization, upstream trimming, breadth/depth parsing, row writers).

---

## Contents

```
.
├── orchestrate_pspA_workflow.py
├── denovo_pspA_detect.py
├── core_boxb_classifier.py
├── rescue_refguided.py
└── pspA_utils.py
```

## Dependencies

### Python (>=3.8)
- **Biopython** (`Bio`) for FASTA/sequence utilities. 
- **pandas** for tabular parsing/writing. 
- **ETE3** for tree rendering (optional; only if `--build-tree` or de novo NT tree used).

### External tools
- **NCBI BLAST+** — `blastx`, `blastn`, `blastdbcmd` (DB/subject modes used by de novo & BoxB steps). 
- **MAFFT** — multiple sequence alignment (AA/NT alignments). 
- **FastTree** — fast phylogeny (AA `-wag` / NT `-gtr`). 
- **minimap2** — short-read mapping in the rescue step. 
- **samtools** — BAM conversion/sorting/indexing and depth calculation. 

> **Note:** All binaries must be on your `PATH`. See the official installation docs for each tool. If you use Conda, the Bioconda channel provides BLAST+, MAFFT, FastTree, minimap2, and samtools.

### Quick environment setup (optional)
Using Conda (recommended for reproducibility):

```bash
conda create -n pspA python=3.10 biopython pandas ete3 -c conda-forge
conda activate pspA
conda install -c bioconda blast mafft fasttree minimap2 samtools
```

---

## Input layout & reference files

- **Assemblies directory**: contains `<assembly_id>.fasta` files (IDs are basenames). Used by de novo detection.
- **Reads directory**: per-strain subfolders, each with `R1_001.fastq.gz` and `R2_001.fastq.gz`. Used by the rescue step. 
- **pspA protein references**: Hollingshead protein DB/FASTA for BLASTX (DB preferred; subject FASTA fallback).
- **Rx1 upstream (~200 bp) references**: nucleotide DB/FASTA used to validate upstream for de novo hits. 
- **BoxB protein references**: DB/FASTA for BoxB classification.

> The rescue step expects a combined FASTA with pspA **plus ~200 bp upstream**; rescued sequences are trimmed to remove the first 200 nt prior to classification.

---

## Typical usage

### Orchestrated workflow

```bash
python3 orchestrate_pspA_workflow.py \
  --assemblies-dir /strains \
  --reads-dir /reads \
  --rescue-ref-fasta 01_fasta/05_Hollingshead_pspAandRx1Upstream/Rx1upstream_and_pspAref.fasta \
  --boxb-db 01_fasta/03_Hollingshead_onlyboxB/CladeDefineRegion_prot_DB \
  --boxb-subject-fasta 01_fasta/03_Hollingshead_onlyboxB/boxBs_hollingsheadAndRx1.fasta \
  --protein-db 01_fasta/02_Hollingshead_pspA/hollingshead_proteinDB \
  --upstream-subject-fasta 01_fasta/05_Hollingshead_pspAandRx1Upstream/Rx1upstream_and_pspAref.fasta \
  --work-dir PathToOutput \
  --threads X \
  --build-tree
```

Key flags (selected):
- `--protein-db` **or** `--protein-subject-fasta` (fallback) control BLASTX target for de novo detection. 
- `--upstream-db` **or** `--upstream-subject-fasta` (fallback) control upstream validation. 
- `--boxb-db` / `--boxb-subject-fasta` control BoxB classification.
- `--rescue-ref-fasta` is the combined upstream+pspA reference set used by rescue.
- `--build-tree` renders a BoxB AA tree with references colored by Hollingshead family. 

### De novo detection (standalone)
```bash
python3 denovo_pspA_detect.py \
  --assemblies-dir <assemblies_dir> \
  --protein-db <pspA_protein_db_basename> \
  --upstream-subject-fasta <rx1_upstream_fasta> \
  --work-dir <work_dir> \
  --threads 12
```
Outputs: `01_denovo/filtered_pspA_hits.fasta`, `01_denovo/summary_fixed_200bp_upstream.csv`, `01_denovo/denovo_summary.csv`; optional NT tree if `--denovo-nt-tree` is used.

### Rescue (standalone)
```bash
python3 rescue_refguided.py \
  --strains <id1> <id2> ... \
  --reads-dir <reads_dir> \
  --ref-fasta-with-200bp <combined_ref.fa> \
  --work-dir <work_dir> \
  --threads 12
```
Outputs are under `04_best/` (minimap2 index, sorted/indexed BAMs, per-reference depth, trimmed FASTA, and `rescue_summary.csv`).

### BoxB classifier (standalone)
```bash
python3 core_boxb_classifier.py \
  --input <coding_nt.fa> \
  --boxb-db <boxb_db_basename> \
  --boxb-subject-fasta <boxb_refs.fasta> \
  --out-dir <out_dir>
```
Writes per-record JSON and CSV with the best AA BoxB hit and family call.

---

## Output overview (orchestrated run)
- `01_denovo/` — filtered coding hits + upstream validation and summary.
- `04_best/` — rescue indexes, BAMs, per-ref depth, trimmed FASTAs, rescue summary.
- `05_boxb_per_sample/` — per-sample BoxB classification outputs (CSV/JSON).
- `06_boxb_reports/` — consolidated `boxB_summary.csv` and `pspA_boxB_compact.csv`.
- `07_tree/` — optional BoxB AA tree (`boxB_tree.nwk`, PNG, legend CSV).

---

## Notes & thresholds
- De novo BLASTX filters: `pident ≥ 75`, `AA length ≥ 200` by default. Upstream BLASTN gate keeps hits with `pident ≥ 75` and `alignment length ≥ 100`; `HighConfidence` if `pident ≥ 90` and `length ≥ 150`. 
- Non-overlapping hit selection per assembly, with distance and hard caps to avoid overlapping HSPs while allowing separate loci.
- Rescue qualification: breadth at depth ≥ *min_depth* (default 10) must be ≥ 0.80 to be deemed `Rescue_Qualified`; otherwise `Rescue_LowBreadth`. 
- BoxB status labels are harmonized to `Matched` or `No_BoxB_match`.

---


## Citation & licensing
Susan K. Hollingshead, Robert Becker, David E. Briles, "Diversity of PspA: Mosaic Genes and Evidence for Past Recombination in Streptococcus pneumoniae" 2000 Infection and Imunity

---

## Acknowledgements
- Hollingshead pspA reference set

