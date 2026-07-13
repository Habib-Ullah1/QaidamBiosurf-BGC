#!/bin/bash
#SBATCH --job-name=phase5_eggnog
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/habib/metagenome/biosurfactant/logs/phase5_eggnog_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_env

PROJECT=/data/habib/metagenome
EGGNOG_DB=$PROJECT/biosurfactant/databases/eggnog
SAMPLES="C6 H2 H6"

for s in $SAMPLES; do
  echo "=== eggNOG-mapper for sample: $s ==="
  echo "Started: $(date)"
  mkdir -p $PROJECT/biosurfactant/$s/05_eggnog
  PROTEINS=$PROJECT/biosurfactant/$s/01_prodigal/${s}_proteins.faa
  OUTDIR=$PROJECT/biosurfactant/$s/05_eggnog
  emapper.py \
    -i $PROTEINS \
    --output ${s}_eggnog \
    --output_dir $OUTDIR \
    --data_dir $EGGNOG_DB \
    --itype proteins \
    --cpu 16 \
    --dmnd_db $EGGNOG_DB/eggnog_proteins.dmnd \
    --override \
    2> $PROJECT/biosurfactant/logs/eggnog_${s}.log
  echo "  Annotations: $(grep -v '^#' $OUTDIR/${s}_eggnog.emapper.annotations | wc -l)"
  echo "Finished: $(date)"
done

echo "Phase 5 COMPLETE"
