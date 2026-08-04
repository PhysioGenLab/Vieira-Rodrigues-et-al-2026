#!/bin/bash
#SBATCH --job-name=fit2_adhd_ferritin_moksnes_2022
#SBATCH --error=/home/thales/MiXer/error/fit2_adhd_ferritin_moksnes_2022_%A_%a.err
#SBATCH --array=1-22%10
#SBATCH --time=192:00:00
#SBATCH --cpus-per-task=7
#SBATCH --nodes=1

set -o errexit

# Base Directories
BASE_DIR="/home/thales/paper_carlos/MiXer"
REF_PLINK="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink"
REF_SNPS="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_snps"
REF_LD="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_ld"

# Traits
TRAIT1_FILE="$BASE_DIR/jobs/ADHD_mixer.csv.gz"
TRAIT2_FILE="/home/thales/MiXer/traitfolder/Ferritina_Moksnes_2022_mixer.csv"

# Fit1 Parameters (per-chromosome)
TRAIT1_PARAMS="$BASE_DIR/results/fit1/adhd/adhd_modelfit1_${SLURM_ARRAY_TASK_ID}.json"
TRAIT2_PARAMS="/home/thales/MiXer/results/fit1/Ferritina_Moksnes_2022_fit1_chr${SLURM_ARRAY_TASK_ID}.json"

# Reference Files
BIM_FILE="$REF_PLINK/1000G.EUR.QC.@.bim"
LD_FILE="$REF_LD/1000G.EUR.QC.@.run4.ld"
SNPS_FILE="$REF_SNPS/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps"

# Output File
OUT_FILE="/home/thales/MiXer/results/fit2/fit2_adhd_Ferritina_Moksnes_2022_${SLURM_ARRAY_TASK_ID}"

# Ensure output directory exists
mkdir -p "$(dirname "$OUT_FILE")"

# Quick check for critical files
for f in "$TRAIT1_PARAMS" "$TRAIT2_PARAMS" "$SNPS_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ File not found: $f"
        exit 1
    fi
done

echo "🚀 Running MiXer FIT2: ADHD vs Ferritin (Chromosome ${SLURM_ARRAY_TASK_ID})..."

singularity exec --cleanenv /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif \
    python /tools/mixer/precimed/mixer.py fit2 \
    --trait1-file "$TRAIT1_FILE" \
    --trait2-file "$TRAIT2_FILE" \
    --trait1-params-file "$TRAIT1_PARAMS" \
    --trait2-params-file "$TRAIT2_PARAMS" \
    --out "$OUT_FILE" \
    --extract "$SNPS_FILE" \
    --bim-file "$BIM_FILE" \
    --ld-file "$LD_FILE" \
    --lib "/home/nico/Mixer2/mixer/src/build/lib/libbgmg.so"

sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,MaxVMSize,Elapsed,TotalCPU
