library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(MendelianRandomization)
library(mr.raps)

# 1. Exposure and Clumping
message("--- Loading Exposure ---")
exposure_dat <- read_exposure_data(
  filename = "TIBC_Moksnes_2022_qc.csv",
  sep = ",", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  eaf_col = "MAF", pval_col = "P", samplesize_col = "N"
)

exposure_dat <- exposure_dat %>% filter(pval.exposure < 5e-08)
tmp_for_clump <- data.frame(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure)
clumped_snps <- ld_clump(tmp_for_clump,
                         plink_bin = "/home/meatavares/.conda/envs/plink_env/bin/plink",
                         bfile = "/scratch/meatavares/MVMR/EUR")
exposure_dat <- subset(exposure_dat, SNP %in% clumped_snps$rsid)

# 2. Outcome
message("--- Loading Outcome ---")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "ADHD_DEMONTIS2022.csv",
  sep = ",", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  eaf_col = "EAF", pval_col = "P", samplesize_col = "N"
)

# 3. Harmonization
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

# --- STEIGER FIX ---
# Steiger requires the outcome sample size N
dat$samplesize.outcome[is.na(dat$samplesize.outcome)] <- 225534
dat$samplesize.exposure[is.na(dat$samplesize.exposure)] <- mean(dat$samplesize.exposure, na.rm=T)

# 4. MR and Diagnostics (Heterogeneity and Pleiotropy)
message("--- Running MR and Diagnostics ---")
res <- mr(dat)
write.csv(res, "res_univariada_TIBC.csv", row.names=F)

het <- mr_heterogeneity(dat)
write.csv(het, "heterogeneidade_TIBC.csv", row.names=F)

pleio <- mr_pleiotropy_test(dat)
write.csv(pleio, "pleiotropia_TIBC.csv", row.names=F)

# Original F-statistic
dat$rsq.exposure <- (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure))
dat$F_stat <- ((dat$samplesize.exposure - 2) * dat$rsq.exposure) / (1 - dat$rsq.exposure)
write.csv(dat, "dat_final_univariada_TIBC.csv", row.names=F)

# 5. MR-RAPS and ConMix
raps_res <- mr.raps.overdispersed.robust(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
cat(capture.output(print(raps_res)), file="mr_raps_TIBC_full.txt")

mr_input_con <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome)
conmix_res <- mr_conmix(mr_input_con)
cat(capture.output(print(conmix_res)), file="mr_conmix_TIBC.txt")

# 6. Steiger Filtering (With error protection)
message("--- Running Steiger Filtering ---")
try({
    # Force calculation of r before calling function to avoid numeric(0) error
    dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
    dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)

    dat_steiger <- steiger_filtering(dat)
    res_steiger <- mr(dat_steiger)
    write.csv(res_steiger, "res_TIBC_post_steiger.csv", row.names=F)
}, silent = FALSE)

message("--- Process Completed ---")
