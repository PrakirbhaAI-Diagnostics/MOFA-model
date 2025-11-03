# Step 1: Update xfun
install.packages("xfun")

# Step 2: Restart R session
# Step 3: Update all packages
update.packages(ask = FALSE)

# Step 4: Reinstall problematic packages
install.packages("miRBaseConverter")
install.packages("dplyr")

# Step 5: Verify xfun version
packageVersion("xfun")  # Should return >= 0.51


# Load all processed datasets
load("BRCA-preprocessed/rna_processed.rda")       # RNA
load("BRCA-preprocessed/mirna_processed_final.rda") # miRNA
load("BRCA-preprocessed/methyl_processed_final.rda") # Methylation
load("BRCA-preprocessed/CNV_processed_final.rda")   # CNV
load("BRCA-preprocessed/protein_processed.rda")    # Protein

# Create standardized 15-character barcodes for each modality
standardize_barcodes <- function(barcodes) {
  substr(barcodes, 1, 15)
}

# Apply to all datasets
colnames(cnv_filtered) <- standardize_barcodes(colnames(cnv_filtered))
colnames(m_values) <- standardize_barcodes(colnames(m_values))
colnames(normalized_counts) <- standardize_barcodes(colnames(normalized_counts))
colnames(prot_final) <- standardize_barcodes(colnames(prot_final))
colnames(rna_vst) <- standardize_barcodes(colnames(rna_vst))

# Get sample lists from all modalities
sample_lists <- list(
  cnv = colnames(cnv_filtered),
  methyl = colnames(m_values),
  mirna = colnames(normalized_counts),
  prot = colnames(prot_final),
  rna = colnames(rna_vst)
)

# Find intersection
common_samples <- Reduce(intersect, sample_lists)
cat("Final common samples:", length(common_samples), "\n")

# Subset each modality to common samples
cnv_common <- cnv_filtered[, common_samples]
methyl_common <- m_values[, common_samples]
mirna_common <- normalized_counts[, common_samples]
prot_common <- prot_final[, common_samples]
rna_common <- rna_vst[, common_samples]

# Check dimensions
cat("CNV dimensions:", dim(cnv_common), "\n")
cat("Methylation dimensions:", dim(methyl_common), "\n")
cat("miRNA dimensions:", dim(mirna_common), "\n")
cat("Protein dimensions:", dim(prot_common), "\n")
cat("RNA dimensions:", dim(rna_common), "\n")

# Verify sample order (should be TRUE)
all(
  identical(colnames(cnv_common), common_samples),
  identical(colnames(methyl_common), common_samples),
  identical(colnames(mirna_common), common_samples),
  identical(colnames(prot_common), common_samples),
  identical(colnames(rna_common), common_samples)
)

# For datasets that might have duplicate columns after truncation
remove_duplicates <- function(mat) {
  mat[, !duplicated(colnames(mat))]
}
cnv_filtered <- remove_duplicates(cnv_filtered)

# View first 5 samples across modalities
data.frame(
  CNV = colnames(cnv_common)[1:5],
  Methylation = colnames(methyl_common)[1:5],
  miRNA = colnames(mirna_common)[1:5],
  Protein = colnames(prot_common)[1:5],
  RNA = colnames(rna_common)[1:5]
)
all(colnames(rna_common) == colnames(mirna_common)) &
  all(colnames(rna_common) == colnames(methyl_common)) &
  all(colnames(rna_common) == colnames(cnv_common)) &
  all(colnames(rna_common) == colnames(prot_common))

rna_common <- as.matrix(rna_common)
mirna_common <- as.matrix(mirna_common)
methyl_common <- as.matrix(methyl_common)
cnv_common <- as.matrix(cnv_common)
prot_common <- as.matrix(prot_common)

#---------- clinical_samples---------
load("C:/Users/BITS/Documents/BRCA-preprocessed/RNA.rda")

library(TCGAbiolinks)
library(SummarizedExperiment)
# Load the clinical data from the RNA file
load("C:/Users/BITS/Documents/BRCA-preprocessed/RNA.rda")

# Store clinical data
clin <- data
clin.sample.info <- as.data.frame(colData(clin))

# Standardize the barcodes in the clinical data to 15 characters
clin.sample.info$standardized_barcode <- substr(clin.sample.info$barcode, 1, 15)

# Filter the clinical data to keep only the patients in the common samples
clin_common <- clin.sample.info[clin.sample.info$standardized_barcode %in% common_samples, ]

# View summary of the clinical data (e.g., vital status and pathologic stage)
cat("Vital Status:\n")
print(table(clin_common$vital_status))

cat("Pathologic Stage:\n")
print(table(clin_common$ajcc_pathologic_stage))

# View the filtered clinical data
head(clin_common)

# Define early and late stages
early_stages <- c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
late_stages  <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")

# Create a new column classifying patients as 'Early' or 'Late'
clin_common$stage_classification <- ifelse(
  clin_common$ajcc_pathologic_stage %in% early_stages, "Early",
  ifelse(clin_common$ajcc_pathologic_stage %in% late_stages, "Late", NA)
)

# View the classification
table(clin_common$stage_classification)

# View the first few rows of the clinical data with stage classification
head(clin_common[, c("ajcc_pathologic_stage", "stage_classification")])

# Remove patients with undefined stage classification
clin_common_filtered <- clin_common[!is.na(clin_common$stage_classification), ]

# View the distribution of early and late stages after filtering
table(clin_common_filtered$stage_classification)

# Save the clinical data with early/late classification
write.csv(clin_common, "BRCA_clinical_stage_classification.csv", row.names = FALSE)


# Load required libraries
install.packages("devtools")
devtools::install_github("r-lib/conflicted")
BiocManager::install("MOFA2")
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)


rna_common <- rna_common[order(apply(rna_common, 1, var), decreasing = TRUE)[1:10000], ]
mirna_common <- mirna_common[order(apply(mirna_common, 1, var), decreasing = TRUE)[1:500], ]
methyl_common <- methyl_common[order(apply(methyl_common, 1, var), decreasing = TRUE)[1:10000], ]
cnv_common <- cnv_common[order(apply(cnv_common, 1, var), decreasing = TRUE)[1:10000], ]

# Combine data into a list
data_list <- list(
  mRNA = rna_common,
  miRNA = mirna_common,
  Methylation = methyl_common,
  CNV = cnv_common,
  Protein = prot_common
)


# Create MOFA object
MOFAobject <- create_mofa(data_list)

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE  # Scale views to unit variance

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10  # Number of latent factors

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"  # Faster convergence
train_opts$seed <- 42  # For reproducibility

# Prepare MOFA object
MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

# Calculate variance explained
r2 <- calculate_variance_explained(MOFAobject)

# Plot variance explained
plot_variance_explained(MOFAobject)

# Plot factors (e.g., Factor 1 vs. Factor 2)
plot_factors(MOFAobject, factors = c(1, 6), color_by = "group")

calculate_variance_explained_per_sample(
  MOFAobject,
  views = "all",
  groups = "all",
  factors = "all")

calculate_variance_explained_per_sample(
  MOFAobject,
  views = "all",
  groups = "all",
  factors = "all")

factors_names(MOFAobject)
features_metadata(MOFAobject)


# UMAP visualization
MOFAobject <- run_umap(MOFAobject)
plot_dimred(MOFAobject, method = "UMAP", color_by = "group")

# Plot top weights for mRNA view, Factor 6
plot_top_weights(MOFAobject, view = "mRNA", factor = 6, nfeatures = 10)
plot_variance_explained(MOFAobject, factor = 6)

total_features <- methyl_common
factors <- get_factors(MOFAobject, factors = 6)[[1]]
factors <- as.data.frame(factors)
plot(factors$Factor6, total_features, xlab=' F6 val', ylab='features',main = "Factor 6 vs Total Features")


# Plot correlation between factors
plot_factor_cor(MOFAobject)

# Save
save_model(MOFAobject, "mofa_model.hdf5")

# Load
MOFAobject <- load_model("mofa_model.hdf5")