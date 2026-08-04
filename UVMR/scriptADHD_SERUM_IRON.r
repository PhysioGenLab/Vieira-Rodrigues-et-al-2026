library(TwoSampleMR)
library(dplyr)
library(mr.raps)               # For MR-RAPS
library(MendelianRandomization) # For ConMix (mr_input function)
library(ieugwasr)              # For local ld_clump

# 1. ADHD is now the EXPOSURE
message("--- Loading Exposure: ADHD ---")
exposure_dat <- read_exposure_data(
  filename = "ADHD_DEMONTIS2022.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",      # ADHD usually has EAF
  pval_col = "P",
  samplesize_col = "N"
)

# Significance filter for ADHD (ADHD instruments)
exposure_dat <- exposure_dat %>% filter(pval.exposure < 5e-08)

# Local Clumping with ADHD SNPs
tmp_for_clump <- data.frame(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure)
clumped_snps <- ld_clump(tmp_for_clump,
                         plink_bin = "/home/meatavares/.conda/envs/plink_env/bin/plink",
                         bfile = "/scratch/meatavares/MVMR/EUR")
exposure_dat <- subset(exposure_dat, SNP %in% clumped_snps$rsid)

# 2. IRON is now the OUTCOME
message("--- Loading Outcome: Iron ---")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "FerroSerico_Moksnes_2022_qc.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "MAF",      # Iron file uses MAF
  pval_col = "P",
  samplesize_col = "N"
)

# 3. Harmonization
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

# 4. N adjustment for Steiger filtering
dat$samplesize.outcome[is.na(dat$samplesize.outcome)] <- 236612

# 5. MR and Diagnostics (Heterogeneity and Pleiotropy)
message("--- Running MR and Diagnostics ---")
res <- mr(dat)
write.csv(res, "res_univariada_ADHD_to_serum_iron.csv", row.names=F)

het <- mr_heterogeneity(dat)
write.csv(het, "heterogeneidade_ADHD_to_serum_iron.csv", row.names=F)

pleio <- mr_pleiotropy_test(dat)
write.csv(pleio, "pleiotropia_ADHD_to_serum_iron.csv", row.names=F)

# Original F-statistic
dat$rsq.exposure <- (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure))
dat$F_stat <- ((dat$samplesize.exposure - 2) * dat$rsq.exposure) / (1 - dat$rsq.exposure)
write.csv(dat, "dat_final_univariada_ADHD_to_serum_iron.csv", row.names=F)

# 5. MR-RAPS and ConMix
raps_res <- mr.raps.overdispersed.robust(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
cat(capture.output(print(raps_res)), file="mr_raps_ADHD_to_serum_iron_full.txt")

mr_input_con <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome)
conmix_res <- mr_conmix(mr_input_con)
cat(capture.output(print(conmix_res)), file="mr_conmix_ADHD_to_serum_iron.txt")

# 6. Steiger Filtering (Error handling included)
message("--- Running Steiger Filtering ---")
try({
    # Force calculation of r prior to running function to avoid numeric(0)
    dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
    dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)

    dat_steiger <- steiger_filtering(dat)
    res_steiger <- mr(dat_steiger)
    write.csv(res_steiger, "res_ADHD_to_serum_iron_post_steiger.csv", row.names=F)
}, silent = FALSE)

message("--- Process Completed ---")
