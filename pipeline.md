# Luciola Mitogenome Phylogeny Pipeline

# Overview

Phylogenetic inference of *Luciola* fireflies from whole mitogenomes across 12 populations, using 13 protein-coding genes (PCGs) + 2 rRNA genes.

| Step | Directory | Description |
|------|-----------|-------------|
| 00 | `00_raw_reads/` | Copy raw reads from hpc_exchange |
| 01 | `01_assembly/` | Mitogenome assembly (GetOrganelle / MITGARD) |
| 02 | `02_annotation/` | Gene annotation (MITOS2) |
| 03 | `03_align/` | Gene extraction, alignment (MACSE/MAFFT), trimming (trimAl) |
| 04 | `04_iqtree/` | Concatenated and partitioned ML tree (IQ-TREE) |
| 05 | `05_revbayes/` | Time-calibrated Bayesian inference (RevBayes) |

## Sample information

We have 122 DNAseq samples and 12 RNAseq samples (ItCol population), plus 25 downloaded GenBank mitogenomes and 1 Luciola singapura with only 15 genes from GenBank.

DNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/shared/data_luciola/{2020,2022,2025}/`
RNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/wzhu/data_luciola/RNASeq/`

**RNAseq samples (ItCol)**
Each RNAseq sample has two tissue libraries (head and body) that need to be combined before assembly.
Assembled using [MITGARD](https://github.com/pedronachtigall/MITGARD) (reference-guided) rather than GetOrganelle, since mitochondrial transcript levels vary greatly across genes. Intergenic regions may be unreliable for RNA-seq assemblies (low/no transcription coverage) and should be treated with caution in downstream analyses.

---

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

### DNAseq samples (GetOrganelle)

#### Installation & Initialization
```bash
conda create -n getorganelle -c bioconda -c conda-forge getorganelle bowtie2
conda activate getorganelle
get_organelle_config.py --add animal_mt
```

#### Run assembly
```bash
mkdir -p logs
bash 01_assembly/scripts/01_run_getOrganelle.sh
```
One SLURM job is submitted per sample. Jobs that are already complete or queued are skipped automatically

#### Check results
```bash
bash 01_assembly/scripts/02_check_assembly.sh
```

All samples should show `circular genome` or `1 scaffold(s)` with length ~16,000–17,000 bp.
A circular genome is ideal. A scaffold of ~16,000+ bp is typically missing only the control region and is usually fine for downstream analysis.
Samples marked `INCOMPLETE` need manual inspection.

### RNAseq samples (MITGARD)

```bash
# merge head + body libraries
bash 01_assembly/scripts/00_merge_rnaseq.sh
# reference-guided assembly
sbatch 01_assembly/scripts/01_run_mitgard.slurm  
bash 01_assembly/scripts/02_check_rnaseq_assembly.sh
bash 01_assembly/scripts/03_collect_rnaseq_assemblies.sh
```

### GenBank assemblies

Public mitogenome sequences from GenBank are added directly to `01_assembly/assemblies/` as `.fasta` files. We can run a quick sanity check on these:

```bash
bash 01_assembly/scripts/02_check_genebank_assembly.sh
```

All samples should have `NSEQS=1` and a length of ~16,000–17,000 bp.

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
# DNAseq samples
sbatch 02_annotation/scripts/01_run_mitos2.slurm
# RNAseq samples
bash 02_annotation/scripts/01_run_mitos2_rnaseq.sh
```

Output GFFs are written to `02_annotation/output/<sample>/result.gff`.

### QC

Check gene counts and assembly statistics for all samples:

```bash
bash 02_annotation/scripts/02_mito_curation_qc.sh
```

check annotation quality `02_annotation/curation_qc.tsv` 
**Interpreting the output:**
 Expected: 13 PCGs, 2 rRNAs, 22 tRNAs.

- `TRNA=21` is acceptable (tRNAs are not used in phylogeny inference)
- Missing PCGs, missing rRNAs, or duplicate genes require manual inspection before proceeding
- **SwLa-LuSpp-22-f-2020_pooled_101:** Assembly is a scaffold (non-circular) missing rrnS, trnV, and the control region (OL) despite extensive parameter tuning (varied `-R`, `-w`, `-k`). The 15,622 bp scaffold covers all 13 PCGs and rrnL. rrnS is added as all `?` in the alignment to maintain consistent taxon sampling across partitions.

---

## Step 03: Gene extraction, alignment, and trimming

### 01: Extract gene sequences

```bash
sbatch 03_align/scripts/01_extract.slurm
```

Extracts 13 PCGs + rrnS + rrnL from each assembly. 

- Genes that wrap around the origin of circular genomes are handled automatically.
- rrnL and rrnS boundaries are recomputed from flanking tRNA positions rather than taken directly from the GFF, because MITOS2 frequently truncates the 5' end of rRNA genes. The logic is:
  - **+ strand:** rrnL 3' end = trnL1_start + 19; rrnS 3' end = trnV_start − 1
  - **− strand:** rrnL 3' end = trnL1_end − 19; rrnS 3' end = trnV_end + 1
  - Expected lengths: rrnL ≈ 1251 bp, rrnS ≈ 747 bp

### 02: Add GenBank samples to extracted fastas

```bash
bash 03_align/scripts/02_add_gb_samples.sh
```

Append them to the per-gene fasta files in `03_align/01_extracted/`.

### Check extracted sequences

```bash
bash 03_align/scripts/03_check_extracted.sh
```


### 04a: Align PCGs (MACSE)

```bash
sbatch 03_align/scripts/04_align_pcg.slurm
```

Uses MACSE v2 with invertebrate mitochondrial code (`-gc_def 5`). Outputs codon-aware nucleotide and amino acid alignments. 
Frameshifts are represented as `!` in MACSE output and are replaced with `-` after alignment.

### 04b: Align rRNAs (MAFFT)

```bash
sbatch 03_align/scripts/04_align_rrna.slurm
```

Uses MAFFT `--auto` for rrnS and rrnL.

### 05: Trim alignments (trimAl)

```bash
sbatch 03_align/scripts/05_trim.slurm
```

**PCGs:** gap-rich columns (gap fraction > 0.7, i.e. `-gt 0.3`) are masked as `?` rather than removed, preserving alignment length and codon reading frame for partitioned analysis. This avoids frame-shifting that would result from column removal.

**rRNAs:** trimmed with `trimAl -gt 0.3`, removing columns where more than 70% of sequences have a gap.

### 06: Convert to NEXUS (for RevBayes)

```bash
sbatch 03_align/scripts/06_fasta2nexus.slurm
```

Converts trimmed FASTA alignments to NEXUS format and writes them to `05_revbayes/data/`. 

---

## Step 04: ML phylogeny (IQ-TREE)

```bash
sbatch 04_iqtree/scripts/run_iqtree.slurm
```

### What the script does

1. **Concatenation:** all 15 trimmed gene alignments are concatenated into a single FASTA (`04_iqtree/data/concat.fasta`)

2. **Partition file:** written to `04_iqtree/data/partitions.txt`. PCG codon positions are pooled across all 13 genes into 3 partitions (pos1, pos2, pos3); rrnS and rrnL are combined into a single rRNA partition, matching the RevBayes model.

3. **IQ-TREE run:** GTR+G4+I model applied independently per partition. Ultrafast bootstrap (`-B 1000`) with 10 independent tree search runs (`--runs 10`) to reduce dependence on starting tree.

### Output

All IQ-TREE output files are written to `04_iqtree/output/` with prefix `luciola_mito`.
The key file is `04_iqtree/output/luciola_mito.treefile`.

---

## Step 05: Time-calibrated phylogeny (RevBayes)

```bash
sbatch 05_revbayes/scripts/run_revbayes_timetree.slurm
```

Partitioned GTR+G4+I model with 4 independent site categories: PCG codon positions 1, 2, 3, and rRNA. Strict clock.

---

## Future work

Intergenic regions could be added as a fifth partition to increase phylogenetic signal, particularly for population-level relationships where PCGs may be too conservative. However, this could be tricky as 1. gene order might differ among species; 2. RNA-seq samples don't have eliable intergenic assembly.
