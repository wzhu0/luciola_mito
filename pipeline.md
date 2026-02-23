### GetOrganelle Installation & Initialization
```bash
conda create -n getorganelle -c bioconda -c conda-forge getorganelle bowtie2
conda activate getorganelle
# download and initialize the database
get_organelle_config.py --add animal_mt
```

### Copy datasets
```bash
bash scripts/cp.reads.sh
```

### Run GetOrganelle assembly
```bash
sbatch scripts/run_getOrganelle.slurm
```