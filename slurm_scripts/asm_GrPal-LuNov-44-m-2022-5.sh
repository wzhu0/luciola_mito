#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/asm_GrPal-LuNov-44-m-2022-5.out
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=wenjie.zhu@lmu.de
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J asm_GrPal-LuNov-

module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate getorganelle

get_organelle_from_reads.py \
    -1 /home/wzhu/luciola/mito/data/GrPal-LuNov-44-m-2022-5_1_val_1.fq.gz \
    -2 /home/wzhu/luciola/mito/data/GrPal-LuNov-44-m-2022-5_2_val_2.fq.gz \
    -R 10 \
    -k 21,45,65,85,105 \
    -F animal_mt \
    -s /home/wzhu/luciola/mito/utils/Luciola_lusitanica.fasta \
    --genes /home/wzhu/luciola/mito/utils/Luciola_lusitanica.fasta \
    -o /home/wzhu/luciola/mito/getorganelle_out/GrPal-LuNov-44-m-2022-5 \
    -t 2
