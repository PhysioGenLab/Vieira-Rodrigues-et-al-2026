#!/bin/bash
#SBATCH --job-name=fit2_hb_adhd
#SBATCH --error=/home/thales/MiXer/error/fit2_hb_adhd_%A_%a.err
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
TRAIT2_FILE="/home/thales/MiXer/traitfolder/Hb_conc_GSA.sumstats.gz"

TRAIT1_PARAMS="$BASE_DIR/results/fit1/adhd/adhd_modelfit1_${SLURM_ARRAY_TASK_ID}.json"
TRAIT2_PARAMS="/home/thales/MiXer/results/fit1/hb/hb_modelfit1_${SLURM_ARRAY_TASK_ID}.json"

# Output File
OUT_FILE="/home/thales/MiXer/results/fit2/hb_adhd/hb_adhd_modelfit2_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$(dirname "$OUT_FILE")"

echo "🚀 Running MiXer fit2 for Hb vs ADHD (Chromosome ${SLURM_ARRAY_TASK_ID})..."

singularity exec --cleanenv /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif \
    python /tools/mixer/precimed/mixer.py fit2 \
    --out "$OUT_FILE" \
    --lib /home/nico/Mixer2/mixer/src/build/lib/libbgmg.so \
    --bim-file "$BIM_FILE" \
    --ld-file "$LD_FILE" \
    --trait1-file "$TRAIT1_FILE" \
    --trait2-file "$TRAIT2_FILE" \
    --extract "$SNPS_FILE" \
    --trait1-params-file "$TRAIT1_PARAMS" \
    --trait2-params-file "$TRAIT2_PARAMS"
