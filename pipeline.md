# Luciola Mitogenome Phylogeny Pipeline

# Overview

Phylogenetic inference of *Luciola* fireflies from whole mitogenomes across 11 populations, using 13 protein-coding genes + 2 rRNA genes.

| Step | Directory | Description |
|------|-----------|-------------|
| 00 | `00_raw_reads/` | Copy raw reads from hpc_exchange |
| 01 | `01_assembly/` | Mitogenome assembly (GetOrganelle) |
| 02 | `02_annotation/` | Gene annotation (MITOS2) |
| 03 | `03_align/` | Gene extraction, alignment (MAFFT), trimming (trimAl) |
| 04 | `04_revbayes/` | Time-calibrated Bayesian inference (RevBayes) |
| 05 | `05_iqtree/` | Concatenated and partitioned ML tree (IQ-TREE) |


## Step 00: Copy raw reads
Run step 00 and step 01 in batches if home directory space is limited.

```bash
# Change the utils/sample_list.txt file to copy reads files you need
# Run on login node inside a screen session
screen -S cp_reads
bash utils/cp_reads.sh
# detach with Ctrl+A D, reattach with: screen -r cp_reads
```

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
mkdir logs
bash 01_assembly/scripts/01_run_getOrganelle.sh
```

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

### QC

Check gene counts and assembly statistics for all samples:

```bash
bash 02_annotation/scripts/02_mito_curation_qc.sh
```

check annotation quality `02_annotation/curation_qc.tsv` 