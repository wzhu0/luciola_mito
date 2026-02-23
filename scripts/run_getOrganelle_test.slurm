#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/assemble_GrPa-LuSpp-3-m-2020_pooled_105.out
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=wenjie.zhu@lmu.de
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J asse_GrPa_3

module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate getorganelle

get_organelle_from_reads.py \
    -1 data/GrPa-LuSpp-3-m-2020_pooled_105_1_val_1.fixed.fq.gz \
    -2 data/GrPa-LuSpp-3-m-2020_pooled_105_2_val_2.fixed.fq.gz \
    -R 10 \
    -k 21,45,65,85,105 \
    -F animal_mt \
    -s utils/Luciola_lusitanica.fasta \
    --genes utils/Luciola_lusitanica.gb \
    -o getorganelle_out/GrPa-LuSpp-3-m-2020_pooled_105 \
    -t 2