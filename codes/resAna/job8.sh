#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p genacc_q
#SBATCH -t 7:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cwu3@fsu.edu
#SBATCH --mem 32G

module load R

cd /gpfs/research/chongwu/Chong/MWAS/resAna


R CMD BATCH  /gpfs/research/chongwu/Chong/MWAS/resAna/best_cis_IGAP.R /gpfs/research/chongwu/Chong/MWAS/logfile/IGAP.txt
