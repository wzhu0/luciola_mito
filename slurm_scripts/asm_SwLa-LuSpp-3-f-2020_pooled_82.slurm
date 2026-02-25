#!/bin/bash
#SBATCH --partition=krypton
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/wzhu/luciola/mito/logs/asm_SwLa-LuSpp-3-f-2020_pooled_82.out
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=wenjie.zhu@lmu.de
#SBATCH --qos=normal_prio
#SBATCH -D /home/wzhu/luciola/mito
#SBATCH -J asm_SwLa-LuSpp-3

module load gnu/12
module load prebin/kry

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate getorganelle

get_organelle_from_reads.py \
    -1 /home/wzhu/luciola/mito/data/SwLa-LuSpp-3-f-2020_pooled_82_1_val_1.fixed.fq.gz \
    -2 /home/wzhu/luciola/mito/data/SwLa-LuSpp-3-f-2020_pooled_82_2_val_2.fixed.fq.gz \
    -R 10 \
    -k 21,45,65,85,105 \
    -F animal_mt \
    -s /home/wzhu/luciola/mito/utils/Luciola_lusitanica.fasta \
    --genes /home/wzhu/luciola/mito/utils/Luciola_lusitanica.fasta \
    -o /home/wzhu/luciola/mito/getorganelle_out/SwLa-LuSpp-3-f-2020_pooled_82 \
    -t 2
