library(here)
library(SummarizedExperiment)
library(dplyr)

# Directory where the RDS files are stored
rds_dir <- here("data", "TCGA_data")

# Get a list of all RDS files in the directory
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Lists for storing expression and clinical data
expr_list <- list()  # Gene expression data
clinical_list <- list()  # Clinical data

# Loop through each RDS file
for (file in rds_files) {
    message("Loading: ", file)
    
    # Load the RDS file
    data <- readRDS(file)
    
    # Extract gene expression data (rows: genes, columns: samples)
    expr_matrix <- as.data.frame(assay(data))
    
    # Extract clinical information (colData)
    clinical_info <- as.data.frame(colData(data))
    
    # Add cancer type information (extract from the file name)
    cancer_type <- gsub(".rds", "", basename(file))  # Extract cancer type from the file name
    clinical_info$cancer_type <- cancer_type  # Add cancer type to clinical data
    
    # Add the data to the corresponding lists
    expr_list[[cancer_type]] <- expr_matrix
    clinical_list[[cancer_type]] <- clinical_info
}

# Merge all expression data by columns
merged_expr <- do.call(cbind, expr_list)

# Save the merged expression data as an RDS file
saveRDS(merged_expr, file = here("data", "TCGA_data", "merged_expression.rds"))

## Merge clinical data across all cancer types

# Convert all columns in each clinical data frame to character type
# This prevents type mismatch errors when combining
clinical_list <- lapply(clinical_list, function(df) {
    df %>% mutate_all(as.character)
})

# Combine all clinical data frames by rows, even if their columns differ
# The `.id` argument adds a new column indicating the cancer type (from the list name)
merged_clinical <- bind_rows(clinical_list, .id = "cancer_type")

# Save the merged clinical data to an RDS file
saveRDS(merged_clinical, file = here("data", "TCGA_data", "merged_clinical.rds"))


##merge experimental and clinical data
# Transpose expression matrix to make samples rows and genes columns
expr_t <- as.data.frame(t(merged_expr))

# Extract SampleID by removing cancer type prefix (e.g., "TCGA-ACC.")
expr_t$SampleID <- sub("^[^.]+\\.", "", rownames(expr_t))

# Add SampleID column to clinical data as well
merged_clinical$SampleID <- rownames(merged_clinical)

# Merge expression and clinical data by SampleID
merged_data <- merge(expr_t, merged_clinical, by = "SampleID", all = TRUE)

# Optionally reassign rownames to keep track by sample
rownames(merged_data) <- merged_data$SampleID

# Save merged data
saveRDS(merged_data, file = here("data", "TCGA_data", "merged_expr_clinical.rds"))
