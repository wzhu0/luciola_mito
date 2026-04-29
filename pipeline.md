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
| 05 | `05_revbayes/` | Bayesian tree inference (RevBayes) |
| 06 | `06_mptp/` | Species delimitation (mPTP) |

## Sample information

We have 122 DNAseq samples and 12 RNAseq samples (ItCol population), plus 23 downloaded GenBank mitogenomes (used to be 25 but 2 were excluded due to ambiguous identification) and 1 *Luciola singapura* with per-gene accessions from GenBank. 
Total: 158 samples after exclusions.

DNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/shared/data_luciola/{2020,2022,2025}/`
RNAseq reads: `/mnt/netvolumes/srva229/bayes/hpc_exchange/wzhu/data_luciola/RNASeq/`

**RNAseq samples (ItCol)**
Each RNAseq sample has two tissue libraries (head and body) that need to be combined before assembly.
Assembled using [MITGARD](https://github.com/pedronachtigall/MITGARD) (reference-guided) rather than GetOrganelle, since mitochondrial transcript levels vary greatly across genes. Intergenic regions may be unreliable for RNA-seq assemblies (low/no transcription coverage) and should be treated with caution in downstream analyses.

**GeneBank mitogenome accession numbers**

25 mitogenomes and 1 Per-gene sample *Luciola singapura* downloaded from GenBank, comprising outgroups and ingroups. Note that *Nipponoluciola cruciata* corresponds to *Luciola cruciata* in Gene Bank (the genus name change Ballantyne et al., 2022). Similarly, *Luciola substriata* was transferred to *Sclerotia substriata* by Ballantyne et al., 2016.

**Outgroups:** Aquatica (3 sequences), *Nipponoluciola cruciata* (6 sequences), and *Sclerotia substriata* (2 sequences).
| Sample | Accession |
|--------|-----------|
| Aquatica_lateralis_1 | LC306678.1 |
| Aquatica_lateralis_2 | NC_035755.1 |
| Aquatica_leii | NC_025276.1 |
| Luciola_unmunsana_2 | NC_050947.1 |
| Nipponoluciola_cruciata_1 | NC_022472.1 | 
| Nipponoluciola_cruciata_2 | LC677170.1 | 
| Nipponoluciola_cruciata_3 | LC306677.1 | 
| Nipponoluciola_cruciata_4 | OM718718.1 | 
| Nipponoluciola_cruciata_5 | OM718717.1 | 
| Nipponoluciola_cruciata_6 | AB849456.1 | 
| Luciola_curtithorax_1 | NC_038225.1 | 
| Luciola_curtithorax_2 | MG770613.1 | 
| Luciola_filiformis | PX289843.1 | 
| Luciola_kagiana_1 | OQ184181.2 | 
| Luciola_kagiana_2 | NC_072664.1 | 
| Luciola_parvula_1 | LC677171.1 | 
| Luciola_parvula_2 | NC_067969.1 | 
| Luciola_parvula_3 | OL944082.1 | 
| Luciola_sp_1 | OP747315.1 | 
| Luciola_sp_2 | OP747314.1 | 
| Sclerotia_substriata_1 | NC_027176.1 | 
| Sclerotia_substriata_2 | KP313820.1 | 
| Luciola_unmunsana_1 | MT134039.1 | 

| Luciola_singapura | Per-gene accessions (MW620428–MW620442): 12S=MW620428, 16S=MW620429, ATP6=MW620430, ATP8=MW620431, COX1=MW620432, COX2=MW620433, COX3=MW620434, CYTB=MW620435, ND1=MW620436, ND2=MW620437, ND3=MW620438, ND4=MW620439, ND4L=MW620440, ND5=MW620441, ND6=MW620442 | 


**Excluded samples:**

| Sample | Accession | Reason |
|--------|-----------|--------|
| Luciola_italica | CM142147.1 | annotation failure |
| Luciola_kagiana_3 | MW260619.1 | UNVERIFIED status in GenBank |
| Luciola_filiformis_2 | MW260625.1 | UNVERIFIED status in GenBank |

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

### 05b: Add missing rrnS for SwLa-22

```bash
bash 03_align/scripts/05b_add_missing_rrnS.sh
```

SwLa-LuSpp-22-f-2020_pooled_101 has no rrnS due to incomplete assembly. This script appends the sample to the trimmed rrnS alignment as an all `?` sequence of the correct alignment length. 

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

2. **Partition file:** written to `04_iqtree/data/partitions.txt`. PCG codon positions are pooled across all 13 genes into 3 partitions (pos1, pos2, pos3); rrnS and rrnL are combined into a single rRNA partition. Total: 4 partitions, matching the RevBayes model structure.

3. **IQ-TREE run:** GTR+G4+I model applied independently per partition. Ultrafast bootstrap (`-B 1000`) with 10 independent tree search runs (`--runs 10`) to reduce dependence on starting tree.

### Output

All IQ-TREE output files are written to `04_iqtree/output/` with prefix `luciola_mito`.
The key file is `04_iqtree/output/luciola_mito.treefile`.

---

## Step 05: Bayesian phylogenetic inference (RevBayes)

We run three separate RevBayes analyses, all using the same partitioned GTR+Gamma+I model (4 partitions: PCG codon positions 1, 2, 3, and rRNA). Each partition has its own GTR exchangeability rates, stationary frequencies, gamma shape parameter, and proportion of invariable sites.

| Analysis | Script | Output |
|----------|--------|--------|
| Time-calibrated tree | `run_revbayes_timetree.slurm` | Branch lengths in absolute time units (strict clock, rate from Hoehna et al. 2025) |
| Unrooted non-clock tree | `run_revbayes_unrooted_non-clock.slurm` | Branch lengths in substitutions per site, unrooted |
| Rooted non-clock tree | `run_revbayes_rooted_non-clock.slurm` | Branch lengths in substitutions per site, rooted on outgroup |

The unrooted and rooted non-clock trees are needed for mPTP species delimitation, which requires branch lengths in substitutions per site. The rooted version uses an outgroup constraint (one representative each from *Aquatica*, *Nipponoluciola*, and *Sclerotia*) 
```bash
sbatch 05_revbayes/scripts/run_revbayes_timetree.slurm
sbatch 05_revbayes/scripts/run_revbayes_unrooted_non-clock.slurm
sbatch 05_revbayes/scripts/run_revbayes_rooted_non-clock.slurm
```

### Output and convergence
 
MCMC trace files and the MCC trees are written to `05_revbayes/output/`. Convergence should be checked before proceeding.

**Future work for timetree analysis**

1. Mitocondrial intergenic regions could be added as a fifth partition to increase phylogenetic signal, particularly for population-level relationships where PCGs may be too conservative. However, this could be tricky as 1. gene order might differ among species; 2. RNA-seq samples don't have eliable intergenic assembly.
2. Relaxed clock for RevBayes time tree inference

### Post-process and visualization
 
After the MCMC runs are complete and convergence have been achieved, post process the outputs, download the MCC tree from the cluster and plot it locally.


```
 
### Tip name substitution
 
Raw RevBayes output uses internal sample codes. Run the rename script to replace them with neat, human-readable names before plotting. The mapping is defined in `utils/rename_map.txt`.
 
```bash
bash utils/apply_rename.sh 05_revbayes/output/Luciola_mito_timetree_MCC.tree
bash utils/apply_rename.sh 05_revbayes/output/Luciola_mito_unrooted_non-clock_MCC.tree
bash utils/apply_rename.sh 05_revbayes/output/Luciola_mito_rooted_non-clock_MCC.tree
```

### Download MCC trees
 
```bash
scp palmuc1:/home/wzhu/luciola/mito/05_revbayes/output/Luciola_mito_timetree_MCC.tree  05_revbayes/output/
scp palmuc1:/home/wzhu/luciola/mito/05_revbayes/output/Luciola_mito_unrooted_non-clock_MCC.tree  05_revbayes/output/
scp palmuc1:/home/wzhu/luciola/mito/05_revbayes/output/Luciola_mito_rooted_non-clock_MCC.tree    05_revbayes/output/

### Visualisation
**Timetree**
Produces an A4 PDF with population colour coding, 95% HPD age bars, alternating time bar background, and a geological timescale. Output goes to `05_revbayes/plots/`.
 
```bash
Rscript 05_revbayes/scripts/plot_rb_tree.R Luciola_mito_timetree_MCC.tree
```

**Unrooted and rooted non-clock trees**
First visualizing them in FigTree.

---

## Step 06: Species delimitation (mPTP)
 
mPTP (multi-rate Poisson Tree Processes) infers species boundaries directly from a rooted ultrametric tree by fitting separate Poisson processes to within-species and between-species branching events. We use the RevBayes MCC tree as input.
 
### Installation
 
```bash
cd ~/software
git clone https://github.com/Pas-Kapli/mptp.git
cd mptp
./autogen.sh
./configure --prefix=$HOME/software/mptp
make
make install
```

Add to PATH in `~/.bashrc`:

```bash
echo 'export PATH="$HOME/software/mptp/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
mptp --version
```
 
### Run mPTP
mPTP requires a rooted tree with branch lengths in substitutions per site.
---