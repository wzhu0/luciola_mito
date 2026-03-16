# Luciola Mitogenome Phylogeny Pipeline

# Overview

Phylogenetic inference of *Luciola* fireflies from whole mitogenomes across 12 populations, using 13 protein-coding genes (PCGs) + 2 rRNA genes.

| Step | Directory | Description |
|------|-----------|-------------|
| 00 | `00_raw_reads/` | Copy raw reads from hpc_exchange |
| 01 | `01_assembly/` | Mitogenome assembly (GetOrganelle) |
| 02 | `02_annotation/` | Gene annotation (MITOS2) |
| 03 | `03_align/` | Gene extraction, alignment (MACSE/MAFFT), trimming (trimAl) |
| 04 | `04_iqtree/` | Concatenated and partitioned ML tree (IQ-TREE) |
| 05 | `05_revbayes/` | Time-calibrated Bayesian inference (RevBayes) |

## Sample information
We have 122 DNAseq samples and 12 RNAseq samples. This pipeline covers the DNAseq samples only. 
A different pipeline might be needed for RNAseq samples.

DNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/shared/data_luciola/{2020,2022,2025}/`
RNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/wzhu/data_luciola/RNASeq/`

**RNAseq samples**
Each RNAseq sample has two tissue libraries (head and body) that need to be combined before assembly.
Mitochondrial transcript levels vary greatly across genes which might make GetOrganelle struggle.

Worth checking [MITGARD](https://github.com/pedronachtigall/MITGARD) - a reference-guided assembly tool to recover the mitochondrial genome from RNA-seq data of any Eukaryote species.

## Step 00: Copy raw reads
Run step 00 and step 01 in batches if home directory space is limited.

```bash
# Edit the utils/sample_list.txt file to copy reads files you need
# Run on login node inside a screen session
screen -S cp_reads
bash utils/cp_reads.sh
# detach with Ctrl+A D, reattach with: screen -r cp_reads
```

---

## Step 01: Mitogenome assembly

### GetOrganelle Installation & Initialization
```bash
conda create -n getorganelle -c bioconda -c conda-forge getorganelle bowtie2
conda activate getorganelle
# download and initialize the database
get_organelle_config.py --add animal_mt
```

### Run assembly
```bash
mkdir -p logs
bash 01_assembly/scripts/01_run_getOrganelle.sh
```
One SLURM job is submitted per sample. Jobs that are already complete or queued are skipped automatically

### Check results
```bash
bash 01_assembly/scripts/02_check_assembly.sh
```

All samples should show `circular genome` or `1 scaffold(s)` with length ~16,000–17,000 bp.
A circular genome is ideal. A scaffold of ~16,000+ bp is typically missing only the control region and is usually fine for downstream analysis.
Samples marked `INCOMPLETE` need manual inspection.

### Collect assemblies

Once all samples are complete, copy final assemblies to `01_assembly/assemblies/`

```bash
bash 01_assembly/scripts/03_collect_assemblies.sh
```

---

## Step 02: Gene annotation using mitos2

### Installation
```bash
conda create -n mitos2 -c bioconda -c conda-forge mitos=2.1.10
conda activate mitos2
# mitos2 uses its own curated reference database
mkdir -p ~/luciola/mito/utils/mitos2_db
cd ~/luciola/mito/utils/mitos2_db
wget https://zenodo.org/records/4284483/files/refseq63m.tar.bz2
tar -xjf refseq63m.tar.bz2
rm refseq63m.tar.bz2
cd ~/luciola/mito
```

MITOS2 doesn't use our own reference .gb file. It works purely with its built-in RefSeq database (refseq63m) to do de novo annotation using BLAST and HMM models trained on thousands of metazoan mitogenomes. 

### Run annotation

```bash
bash 02_annotation/scripts/01_run_mitos2.sh
```
One SLURM job is submitted per sample. Output GFFs are written to `02_annotation/output/<sample>/result.gff`.

### QC

Check gene counts and assembly statistics for all samples:

```bash
bash 02_annotation/scripts/02_mito_curation_qc.sh
```

check annotation quality `02_annotation/curation_qc.tsv` 

**Interpreting the output:**
assembly length, contig count, PCG count (expect 13), rRNA count (expect 2), tRNA count (expect 22), duplicated genes.
- `TRNA=21` is acceptable since tRNAs are not used in the phylogeny inference. 
- All other deviations (missing PCGs, missing rRNAs, duplicate genes) require manual inspection of the GFF before proceeding.

---

## Step 03: Gene extraction, alignment, and trimming

### 03a: Extract gene sequences

```bash
sbatch 03_align/scripts/01_extract.slurm
```

Extracts 13 PCGs + rrnS + rrnL from each assembly. 

- Genes that wrap around the origin of circular genomes are handled automatically.
- rrnL and rrnS boundaries are recomputed from flanking tRNA positions rather than taken directly from the GFF, because MITOS2 frequently truncates the 5' end of rRNA genes. The logic is:
  - **+ strand:** rrnL 3' end = trnL1_start + 19; rrnS 3' end = trnV_start − 1
  - **− strand:** rrnL 3' end = trnL1_end − 19; rrnS 3' end = trnV_end + 1
  - Expected lengths: rrnL = 1251 bp, rrnS = 747 bp

### Check extracted sequences

```bash
bash 03_align/scripts/01a_check_extracted.sh
```

**Interpreting the output:**


### 03b: Align PCGs (MACSE)

```bash
sbatch 03_align/scripts/02_align_pcg.slurm
```

Uses MACSE v2 with invertebrate mitochondrial code (`-gc_def 5`). Outputs codon-aware nucleotide and amino acid alignments. 
Frameshifts are represented as `!` in MACSE output and are replaced with `-` after alignment.

### 03c: Align rRNAs (MAFFT)

```bash
sbatch 03_align/scripts/02_align_rrna.slurm
```

Uses MAFFT `--auto` for rrnS and rrnL.

### Check alignments

```bash
bash 03_align/scripts/02a_check_alignments.sh
```

Recommended to do eye inspection with AliView as well.

### 03d: Trim alignments (trimAl)

```bash
sbatch 03_align/scripts/03_trim.slurm
```

Uses `trimAl -automated1` for all 15 genes. This removes gap-rich columns while retaining informative variation. 

Recommended to do another round of eye inspection.

### 03e: Convert to NEXUS (for RevBayes)

```bash
sbatch 03_align/scripts/04_fasta2nexus.slurm
```

Converts trimmed FASTA alignments to NEXUS format and writes them to `05_revbayes/data/`. 

---

## Step 04: ML phylogeny (IQ-TREE)

```bash
sbatch 04_iqtree/scripts/run_iqtree.slurm
```

### What the script does

1. **Concatenation:** all 15 trimmed gene alignments are concatenated into a single FASTA (`04_iqtree/data/concat.fasta`) 

2. **Partition file:** written to `04_iqtree/data/partitions.txt`. PCGs are split into three codon position partitions (pos1, pos2, pos3). rRNAs are a single partition each. Total: 13×3 + 2 = 41 partitions.

3. **IQ-TREE run:** ModelFinder (`-m MFP`) selects the best substitution model per partition. Ultrafast bootstrap (`-B 1000`) provides branch support values.


### Output

All IQ-TREE output files are written to `04_iqtree/output/` with prefix `luciola_mito`.
The key file is `04_iqtree/output/luciola_mito.treefile`.

---

## Step 05: Time-calibrated phylogeny (RevBayes)

```bash
sbatch 05_revbayes/scripts/run_revbayes_timetree.slurm
```