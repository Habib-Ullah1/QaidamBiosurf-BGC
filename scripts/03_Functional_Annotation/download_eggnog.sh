#!/bin/bash
#SBATCH --job-name=eggnog_download
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/eggnog_download_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_env
cd /data/habib/metagenome/biosurfactant/databases/eggnog/

wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz \
    && gunzip eggnog.db.gz

wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz \
    && tar -zxf eggnog.taxa.tar.gz \
    && rm eggnog.taxa.tar.gz

wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz \
    && gunzip eggnog_proteins.dmnd.gz

echo "ALL DONE"
EOF
