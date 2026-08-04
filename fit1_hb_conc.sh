#!/bin/bash
#SBATCH --job-name=fit1_hb
#SBATCH --error=/home/thales/MiXer/error/fit1_error/hb_%A_%a.err
#SBATCH --array=1-22%10
#SBATCH --time=192:00:00
#SBATCH --cpus-per-task=7
#SBATCH --nodes=1

set -o errexit

singularity exec --cleanenv --bind /home/nico,/home/thales /home/nico/Mixer2/mixer/mixer-1.3.0/singularity/mixer_newer2.sif python /tools/mixer/precimed/mixer.py fit1 \
      --trait1-file /home/thales/MiXer/traitfolder/Hb_conc_GSA.sumstats.gz \
      --out /home/thales/MiXer/results/fit1/hb/hb_modelfit1_${SLURM_ARRAY_TASK_ID} \
      --extract /home/nico/Mixer2/mixer/reference_files/1000G_EUR_Phase3_plink_snps/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file /home/nico/Mixer2/mixer/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/nico/Mixer2/mixer/reference_files/1000G_EUR_Phase3_plink_ld/1000G.EUR.QC.@.run4.ld \
      --lib /home/nico/Mixer2/mixer/src/build/lib/libbgmg.so

sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,MaxVMSize,Elapsed,TotalCPU
