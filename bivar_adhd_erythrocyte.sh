#!/bin/bash
#SBATCH --job-name=bivar_adhd_erythrocyte
#SBATCH --error=/home/thales/MiXer/error/bivar_adhd_erythrocyte_%A_%a.err
#SBATCH --array=1-22%10
#SBATCH --time=192:00:00
#SBATCH --cpus-per-task=7
#SBATCH --nodes=1

set -o errexit

# Base & Reference Paths
BASE_DIR="/home/thales/paper_carlos/MiXer"
REF_PLINK="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink"
REF_SNPS="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_snps"
REF_LD="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_ld"

# Reference Files
BIM_FILE="$REF_PLINK/1000G.EUR.QC.@.bim"
LD_FILE="$REF_LD/1000G.EUR.QC.@.run4.ld"
SNPS_FILE="$REF_SNPS/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps"

# Traits & Params Files
TRAIT1_FILE="$BASE_DIR/jobs/ADHD_mixer.csv.gz"
TRAIT2_FILE="/home/thales/MiXer/traitfolder/erythrocyte_count_GSA.sumstats.gz"

FIT2_FILE="/home/thales/MiXer/results/fit2/erythrocyte_adhd/erythrocyte_adhd_modelfit2_${SLURM_ARRAY_TASK_ID}.json"

# Output File
OUT_FILE="/home/thales/MiXer/results/bivar/erythrocyte_adhd/erythrocyte_adhd_bivar_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$(dirname "$OUT_FILE")"

echo "🚀 Running MiXer test2 (bivar) for ADHD vs Erythrocyte Count (Chromosome ${SLURM_ARRAY_TASK_ID})..."

singularity exec --cleanenv /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif \
    python /tools/mixer/precimed/mixer.py test2 \
    --out "$OUT_FILE" \
    --lib /home/nico/Mixer2/mixer/src/build/lib/libbgmg.so \
    --bim-file "$BIM_FILE" \
    --ld-file "$LD_FILE" \
    --trait1-file "$TRAIT1_FILE" \
    --trait2-file "$TRAIT2_FILE" \
    --extract "$SNPS_FILE" \
    --load-params-file "$FIT2_FILE"
