# 1. Load functions and libraries
source("RunLCV.R")
source("MomentFunctions.R")
library(data.table)
library(dplyr)

# 2. Load LD Scores
print("Loading LD Scores...")
ldscore_files <- list.files("eur_w_ld_chr/", pattern = "*.l2.ldscore.gz", full.names = TRUE)
ldscore_list <- lapply(ldscore_files, fread)
ldscore_all <- rbindlist(ldscore_list)
ldscore_all <- ldscore_all[, .(SNP, L2)] # Keep only essential columns to save RAM

# 3. Read Summary Statistics
print("Reading Sumstats...")
d1 <- fread("ADHD.sumstats.sumstats.gz")
d2 <- fread("reticulocitos.sumstats.gz")

# 4. Merge and Data Cleaning
print("Merging and cleaning data...")
data <- ldscore_all %>%
  inner_join(d1, by = "SNP", suffix = c("", ".adhd")) %>%
  inner_join(d2, by = "SNP", suffix = c(".adhd", ".retic"))

# REMOVE INCOMPLETE ROWS
# The is.na() function will drop empty fields
data <- data %>%
  filter(!is.na(Z.adhd) & !is.na(Z.retic) & !is.na(L2)) %>%
  filter(is.finite(Z.adhd) & is.finite(Z.retic) & is.finite(L2))

# Order by genomic position (Required for LCV Jackknife)
# If sumstats lack CHR and BP, sorting by SNP is the minimum fallback
data <- data %>% arrange(SNP) # If CHR/BP are present in data, use: arrange(CHR, BP)

# 5. Run LCV
print(paste("Running LCV with", nrow(data), "cleaned SNPs..."))
# Set sig.threshold to 30 to prevent errors from SNPs with extremely large effects
results <- RunLCV(ell = data$L2,
                  z.1 = data$Z.adhd,
                  z.2 = data$Z.retic,
                  sig.threshold = 30)

# 6. Display results
cat("\n--- LCV RESULTS ---\n")
if(!is.null(results)){
  cat(sprintf("Estimated posterior gcp = %.2f (SE: %.2f)\n", results$gcp.pm, results$gcp.pse))
  cat(sprintf("P-value for GCP=0: %.3e\n", results$pval.gcpzero.2tailed))
  cat(sprintf("Genetic Correlation (rho): %.2f (SE: %.2f)\n", results$rho.est, results$rho.err))
} else {
  print("Error calculating LCV.")
}
