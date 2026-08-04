# 1. Functions and Libraries
source("RunLCV.R")
source("MomentFunctions.R")
library(data.table)
library(dplyr)

# 2. LD Scores
print("Loading LD Scores...")
ldscore_files <- list.files("eur_w_ld_chr/", pattern = "*.l2.ldscore.gz", full.names = TRUE)
ldscore_list <- lapply(ldscore_files, fread)
ldscore_all <- rbindlist(ldscore_list)
ldscore_all <- ldscore_all[, .(SNP, L2)]

# 3. Sumstats
print("Reading Sumstats...")
d1 <- fread("ADHD.sumstats.sumstats.gz") # Trait 1
d2 <- fread("Ferro.sumstats.sumstats.gz") # Trait 2

# 4. Merge and Allele Alignment
print("Merging and Aligning Alleles...")
data <- ldscore_all %>%
  inner_join(d1, by = "SNP") %>%
  inner_join(d2, by = "SNP", suffix = c(".adhd", ".iron"))

# Initial NA cleaning
data <- data %>% filter(!is.na(Z.adhd) & !is.na(Z.iron) & !is.na(L2))

# --- ALLELE ALIGNMENT LOGIC ---
data <- data %>%
  mutate(Z.iron.adj = case_when(
    A1.adhd == A1.iron ~ Z.iron,
    A1.adhd == A2.iron ~ Z.iron * -1,
    TRUE ~ as.numeric(NA)
  )) %>%
  filter(!is.na(Z.iron.adj))

# Order by position (Required for Jackknife)
data <- data %>% arrange(SNP)

# 5. Run LCV
print(paste("Running LCV with", nrow(data), "aligned SNPs..."))
results <- RunLCV(ell = data$L2,
                  z.1 = data$Z.adhd,
                  z.2 = data$Z.iron.adj,
                  sig.threshold = 30)

# 6. Display Results
cat("\n--- LCV RESULTS (ALIGNED) ---\n")
if(!is.null(results)){
  cat(sprintf("Estimated posterior gcp = %.2f (SE: %.2f)\n", results$gcp.pm, results$gcp.pse))
  cat(sprintf("P-value for GCP=0: %.3e\n", results$pval.gcpzero.2tailed))
  cat(sprintf("Genetic Correlation (rho): %.3f (SE: %.2f)\n", results$rho.est, results$rho.err))
}
