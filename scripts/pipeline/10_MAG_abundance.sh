#!/bin/bash
#SBATCH --job-name=coverm_abund
#SBATCH --partition=cu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/coverm_abund_%j.log
source ~/miniconda3/etc/profile.d/conda.sh
conda activate coverm_env
cd /data/habib/metagenome/biosurfactant
export TMPDIR=/data/habib/metagenome/biosurfactant/16_mag_abundance/tmp; mkdir -p "$TMPDIR"
R=/data/habib/metagenome/ncbi_upload; OUT=16_mag_abundance
coverm genome --genome-fasta-directory $OUT/all_MAGs_renamed -x fa \
  --coupled $R/C3_1.fq.gz $R/C3_2.fq.gz $R/C6.R1.fq.gz $R/C6.R2.fq.gz \
            $R/H2.R1.fq.gz $R/H2.R2.fq.gz $R/H6.R1.fq.gz $R/H6.R2.fq.gz \
  -m relative_abundance mean covered_fraction --min-covered-fraction 10 \
  -t ${SLURM_CPUS_PER_TASK} -o $OUT/MAG_relabund.tsv
echo "EXIT=$?"; wc -l $OUT/MAG_relabund.tsv
