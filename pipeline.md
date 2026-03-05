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
```bash
# Change the utils/sample_list.txt file to copy files you need
# Run on login node inside a screen session
screen -S cp_reads
bash utils/cp_reads.sh
# detach with Ctrl+A D, reattach with: screen -r cp_reads
```
### GetOrganelle Installation & Initialization
```bash
conda create -n getorganelle -c bioconda -c conda-forge getorganelle bowtie2
conda activate getorganelle
# download and initialize the database
get_organelle_config.py --add animal_mt
```

### Run GetOrganelle assembly
```bash
bash scripts/run_getOrganelle.sh
```
