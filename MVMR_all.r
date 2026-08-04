# 1. Environment and Library Setup
.libPaths(c("~/Rlibs", .libPaths()))

if (!require("MVMR", quietly = TRUE)) {
    message("Installing MVMR...")
    if (!require("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
    remotes::install_github("WSpiller/MVMR", upgrade = "never")
}

library(MVMR)
library(TwoSampleMR)
library(tibble)
library(dplyr)
library(MendelianRandomization)

# 2. Extract Exposures (Local)
message("--- Extracting Exposures ---")
mv_exposure_dat_rev <- mv_extract_exposures_local(
  filenames_exposure = c(
    "hb_conc_qc_.csv",
    "erythrocyte_count_qc.csv",
    "reticulocyte_count_qc.csv",
    "Ferritina_Moksnes_2022_qc.csv",
    "FerroSerico_Moksnes_2022_qc.csv",
    "TIBC_Moksnes_2022_qc.csv",
    "TSP_Moksnes_2022_qc.csv"
  ),
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  pval_threshold = 5e-08,
  plink_bin = "/home/meatavares/.conda/envs/plink_env/bin/plink",
  bfile = "/scratch/meatavares/MVMR/EUR",
  clump_r2 = 0.001,
  clump_kb = 10000,
  pop = "EUR",
  harmonise_strictness = 2
)

# Remove any potential duplicates to avoid matrix conflicts
mv_exposure_dat_rev <- mv_exposure_dat_rev %>%
  distinct(SNP, id.exposure, .keep_all = TRUE)

# 3. Read Outcome (ADHD)
message("--- Reading Outcome Data ---")
mv_outcome_dat_rev <- read_outcome_data(
  snps = mv_exposure_dat_rev$SNP,
  filename = "ADHD_DEMONTIS2022.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  pval_col = "P",
  samplesize_col = "N"
)

# 4. Harmonization and Primary Analysis
message("--- Harmonizing and Running MVMR ---")
mvdat_rev <- mv_harmonise_data(mv_exposure_dat_rev, mv_outcome_dat_rev)

# Analysis via TwoSampleMR
res_rev <- mv_multiple(mvdat_rev)
saveRDS(res_rev, file = "res_MVMR_Iron_ADHD.rds")

# 5. Generate Diagnostics (MVMR Package)
message("--- Generating Diagnostics (F-stat and Q-stat) ---")

# Safe extraction of harmonized matrices
beta_exp  <- as.matrix(mvdat_rev$exposure_beta)
se_exp    <- as.matrix(mvdat_rev$exposure_se)
beta_out  <- as.numeric(mvdat_rev$outcome_beta)
se_out    <- as.numeric(mvdat_rev$outcome_se)
snps_list <- rownames(mvdat_rev$exposure_beta)

no_exp <- ncol(beta_exp)
exposures_names <- colnames(mvdat_rev$exposure_beta)

# Save name mappings (ID vs Real Phenotype)
mapping <- data.frame(ID = paste0("betaX", 1:no_exp), Phenotype = exposures_names)
write.csv(mapping, "Mapping_Exposures_Names.csv", row.names = FALSE)

# Format data for the MVMR package
XGs_df <- as.data.frame(beta_exp)
colnames(XGs_df) <- paste0("betaX", 1:no_exp)
XGs_df$SNP <- snps_list

SEs_df <- as.data.frame(se_exp)
colnames(SEs_df) <- paste0("seX", 1:no_exp)
SEs_df$SNP <- snps_list

# Join to ensure uniform SNP ordering
mvmr_input_df <- inner_join(XGs_df, SEs_df, by = "SNP")

# Final formatting for MVMR package
mvmr_out_rev <- format_mvmr(
  BXGs   = mvmr_input_df %>% select(starts_with("betaX")),
  BYG    = beta_out,
  seBXGs = mvmr_input_df %>% select(starts_with("seX")),
  seBYG  = se_out,
  RSID   = mvmr_input_df$SNP
)

# 6. Model Validity Diagnostics
# F-statistics (Instrument Strength)
sres <- strength_mvmr(r_input = mvmr_out_rev, gencov = 0)

# Horizontal Pleiotropy (Q-statistic)
pres <- pleiotropy_mvmr(r_input = mvmr_out_rev, gencov = 0)

# Save diagnostics
write.csv(sres, "MVMR_Iron_Strength_Fstat.csv")
write.csv(pres, "MVMR_Iron_Pleiotropy_Qstat.csv")

message("--- Process completed! Check generated CSV files. ---")
