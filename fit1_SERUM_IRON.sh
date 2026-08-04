#!/bin/bash
#SBATCH --job-name=mixer_fit1_FerroSerico_Moksnes_2022
#SBATCH --error=/home/thales/MiXer/error/mixer_fit1_FerroSerico_Moksnes_2022_%A_%a.err
#SBATCH --array=1-22%10
#SBATCH --time=192:00:00
#SBATCH --cpus-per-task=7
#SBATCH --nodes=1

BASE_DIR="/home/thales/paper_carlos/MiXer"
REF_PLINK="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink"
REF_SNPS="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_snps"
REF_LD="$BASE_DIR/reference_files/1000G_EUR_Phase3_plink_ld"
TRAIT_DIR="/home/thales/MiXer/traitfolder"

BIM_FILE="$REF_PLINK/1000G.EUR.QC.@.bim"
LD_FILE="$REF_LD/1000G.EUR.QC.@.run4.ld"
SNPS_FILE="$REF_SNPS/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps"

if [[ ! -f "$SNPS_FILE" ]]; then
    echo "❌ File not found: $SNPS_FILE"
    exit 1
fi

TRAIT_FILE="$TRAIT_DIR/FerroSerico_Moksnes_2022_mixer.csv"
TRAIT_NAME="FerroSerico_Moksnes_2022"

OUT_FILE="/home/thales/MiXer/results/fit1/${TRAIT_NAME}_fit1_chr${SLURM_ARRAY_TASK_ID}"
mkdir -p "$(dirname "$OUT_FILE")"

echo "🚀 Running MiXer fit1 for trait ${TRAIT_NAME}, chromosome ${SLURM_ARRAY_TASK_ID}..."

singularity exec --cleanenv /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif \
    python /tools/mixer/precimed/mixer.py fit1 \
    --trait1-file "$TRAIT_FILE" \
    --out "$OUT_FILE" \
    --extract "$SNPS_FILE" \
    --bim-file "$BIM_FILE" \
    --ld-file "$LD_FILE" \
    --lib /home/nico/Mixer2/mixer/src/build/lib/libbgmg.so
