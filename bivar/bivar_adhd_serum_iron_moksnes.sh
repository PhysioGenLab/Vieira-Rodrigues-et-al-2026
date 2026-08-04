#!/bin/bash
#SBATCH --job-name=bivar_adhd_serum_iron_moksnes
#SBATCH --error=/home/thales/MiXer/error/bivar_adhd_serum_iron_moksnes_2022_%A_%a.err
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

# Traits (z-score formatted sumstats)
TRAIT1_FILE="/home/thales/MiXer/traitfolder/FerroSerico_Moksnes_2022_mixer.csv"
TRAIT2_FILE="$BASE_DIR/jobs/ADHD_mixer.csv.gz"

# Fit2 Parameters (per-chromosome)
TRAIT1_PARAMS="/home/thales/MiXer/results/fit2/fit2_adhd_FerroSerico_Moksnes_2022_${SLURM_ARRAY_TASK_ID}.json"

# Reference Files
BIM_FILE="$REF_PLINK/1000G.EUR.QC.@.bim"
LD_FILE="$REF_LD/1000G.EUR.QC.@.run4.ld"
SNPS_FILE="$REF_SNPS/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps"

# Output File
OUT_FILE="/home/thales/MiXer/results/bivar/bivar_adhd_MOKSNES/bivar_adhd_FerroSerico_Moksnes_2022_${SLURM_ARRAY_TASK_ID}"

# Ensure output directory exists
mkdir -p "$(dirname "$OUT_FILE")"

# Quick check for critical files
for f in "$TRAIT1_PARAMS" "$TRAIT1_FILE" "$TRAIT2_FILE" "$SNPS_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ File not found: $f"
        exit 1
    fi
done

echo "🚀 Running MiXer test2 (bivar) for Serum Iron vs ADHD (Chromosome ${SLURM_ARRAY_TASK_ID})..."

singularity exec --cleanenv /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif \
    python /tools/mixer/precimed/mixer.py test2 \
    --trait1-file "$TRAIT1_FILE" \
    --trait2-file "$TRAIT2_FILE" \
    --load-params-file "$TRAIT1_PARAMS" \
    --out "$OUT_FILE" \
    --extract "$SNPS_FILE" \
    --bim-file "$BIM_FILE" \
    --ld-file "$LD_FILE" \
    --lib "/home/nico/Mixer2/mixer/src/build/lib/libbgmg.so"

sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,MaxVMSize,Elapsed,TotalCPU
