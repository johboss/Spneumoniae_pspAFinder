
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

### Five-sample test run

Run this from the repository directory after activating the `ete3` conda environment. The test inputs contain assemblies and GFF files for five samples; no short-read rescue inputs are required unless a sample has no de-novo candidate.

```bash
conda activate ete3

python3 orchestrate_pspA_workflow.py \
  --assemblies-dir /Users/john.boss/30_pspA/11_pspAFinder/06_test \
  --reads-dir /Users/john.boss/30_pspA/11_pspAFinder/06_test \
  --rescue-ref-fasta /Users/john.boss/30_pspA/11_pspAFinder/01_fasta/05_Hollingshead_pspAandRx1Upstream/Rx1upstream_and_pspAref.fasta \
  --boxb-db /Users/john.boss/30_pspA/11_pspAFinder/01_fasta/03_Hollingshead_onlyboxB/CladeDefineRegion_prot_DB \
  --boxb-subject-fasta /Users/john.boss/30_pspA/11_pspAFinder/01_fasta/03_Hollingshead_onlyboxB/CladeDefineRegion_ManualTrim20240513.fasta \
  --protein-db /Users/john.boss/30_pspA/11_pspAFinder/01_fasta/02_Hollingshead_pspA/hollingshead_proteinDB \
  --upstream-subject-fasta /Users/john.boss/30_pspA/11_pspAFinder/01_fasta/05_Hollingshead_pspAandRx1Upstream/Rx1upstream_and_pspAref.fasta \
  --work-dir /Users/john.boss/30_pspA/11_pspAFinder/07_testOut \
  --gff-root /Users/john.boss/30_pspA/11_pspAFinder/06_test \
  --metadata-csv /Users/john.boss/39_Metadata/Metadata_ALL.csv \
  --max-upstream-length 500 \
  --orf-window 2000 \
  --threads 2 \
  --jobs 1 \
  --build-pspa-alignment \
  --build-upstream-tree
```

The run writes BoxB reports under `07_testOut/06_boxb_reports/`, including `confirmed_pspA_sequences.fasta`. Confirmed coordinate-bearing hits are refined with Prodigal and blastp in a default ±2000 bp window before upstream extraction. Only predicted ORFs containing the BLASTX start or end coordinate are eligible, preventing a separate nearby gene from being selected. Short windows use Prodigal `-p meta` because `-p single` requires at least 20 kb of sequence. Use `--orf-window` to change that window. If Prodigal cannot refine a locus, the pipeline records a Low-confidence fallback using the original BLASTX coordinates instead of stopping the run. Add `--build-pspa-alignment` to also write `confirmed_pspA_alignment.fasta` using MAFFT. Upstream sequences, alignment, tree, metadata, and iTOL datasets are written under `07_testOut/08_upstream_phylogeny/`. Upstream regions are capped at 500 bp by default; if a GFF implies a longer region, the 500 bp nearest pspA is retained. Use `--max-upstream-length` to change this limit. To start from a clean test output, remove or rename `07_testOut` before running the command.

### De novo detection (standalone)
```bash
python3 denovo_pspA_detect.py \
  --assemblies-dir <assemblies_dir> \
  --protein-db <pspA_protein_db_basename> \
  --upstream-subject-fasta <rx1_upstream_fasta> \
  --work-dir <work_dir> \
  --threads 12
```
Outputs: `01_denovo/filtered_pspA_hits.fasta`, `01_denovo/filtered_pspA_hit_coordinates.csv`, `01_denovo/summary_fixed_200bp_upstream.csv`, `01_denovo/denovo_summary.csv`; optional NT tree if `--denovo-nt-tree` is used.

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
- `06_boxb_reports/` — consolidated `boxB_summary.csv`, `pspA_boxB_compact.csv`, authoritative BoxB-confirmed loci in `confirmed_pspA_hits.csv` (including BLASTX and refined ORF coordinates/warnings), and all confirmed pspA nucleotide sequences in `confirmed_pspA_sequences.fasta` (plus its optional MAFFT alignment).
- `07_tree/` — optional BoxB AA tree (`boxB_tree.nwk`, PNG, legend CSV).
- `08_upstream_phylogeny/` — BoxB-confirmed upstream sequences, metadata, `phandango_metadata.csv` with exact tree-tip IDs, MAFFT alignment, FastTree tree, categorical/gradient iTOL datasets, converter-generated per-strain `DATASET_RANGE` alternatives, and `node_popup.txt` metadata for iTOL node selection when `--build-upstream-tree` is enabled.

The upstream analysis is driven exclusively by `06_boxb_reports/confirmed_pspA_hits.csv`; candidate loci in the de-novo upstream summary are not used. Configure annotation and metadata inputs with `--gff-root` and `--metadata-csv` (both have defaults for the standard T7 layout).

### Standalone promoter architecture analysis

After manually reviewing and, if needed, editing `08_upstream_phylogeny/upstream_alignment.fasta`, run the downstream promoter analysis separately:

```bash
conda activate ete3
python promoter_architecture.py \
  --alignment 07_testOut/08_upstream_phylogeny/upstream_alignment.fasta \
  --metadata 07_testOut/08_upstream_phylogeny/phandango_metadata.csv \
  --output-dir 07_testOut/11_promoter_analysis
```

The module searches each reviewed upstream sequence for sigma70-like `-35` (`TTGACA`) and `-10` (`TATAAT`) motifs, allowing up to two mismatches per motif and a 14–20 bp spacer. It writes promoter calls, merged Phandango metadata, and association results without changing the main pspAFinder workflow.

The same run also writes `promoter_annotation.gff`. This is a positional Phandango annotation track, not a metadata column: coordinates use the selected reference record's alignment columns, including gaps. By default the first FASTA record is the reference; select another with `--reference-id`. The track includes `minus35`, `spacer`, `minus10`, and `pspA_start`. If the alignment coordinate system includes the pspA coding region, add its alignment-coordinate length to create the complete ORF feature:

```bash
python promoter_architecture.py \
  --alignment 08_upstream_phylogeny/upstream_alignment_trimmed.fasta \
  --metadata 09_phandango/phandango_metadata.csv \
  --output-dir 11_promoter_analysis \
  --reference-id 18-1__hit3 \
  --pspa-length 2300
```

Keep SNP/indel data in a separate Phandango variant file such as `upstream_variants.tsv`; the GFF contains only positional promoter and pspA annotations.

The same run writes `upstream_variants.tsv`, `variation_annotation.gff`, and `consensus.fasta`. These compare every alignment column with the majority consensus while retaining gap columns: `SNP` records contain base substitutions, `Deletion` records contain consensus-base-to-gap events, and `Insertion` records contain bases in columns where the consensus is a gap. Load `variation_annotation.gff` as the positional variation track; do not merge these events into categorical Phandango metadata.

### Grouped upstream plus pspA variation comparison

To compare selected closely related strains without comparing between groups, run:

```bash
python compare_promotor_variation.py \
  --upstream-fasta /Volumes/T7/pspA/02_JPIAMR_B/08_upstream_phylogeny/upstream_sequences.fasta \
  --gene-fasta /Volumes/T7/pspA/02_JPIAMR_B/06_boxb_reports/confirmed_pspA_sequences.fasta \
  --output-dir /Volumes/T7/pspA/02_JPIAMR_B/12_promotor_variation_analysis \
  --threads 4
```

The script matches records by their `j###` prefix, concatenates upstream and confirmed pspA sequences, and runs MAFFT separately for each requested group. It writes `group1_combined.fasta` through `group4_combined.fasta`, corresponding aligned FASTA files, `pairwise_difference_summary.csv`, and `comparison_parameters.txt`. `Total_differences` includes both nucleotide substitutions and alignment gap/indel columns; the summary also reports these as `SNP_differences` and `Indel_differences`.

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
