library(here)
library(DESeq2)
library(ggplot2)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

# Load merged raw count expression matrix (genes x samples)
merged_expr <- readRDS(file = here("data", "TCGA_data", "merged_expression.rds"))

# Get sample IDs from column names
sample_ids <- colnames(merged_expr)

# Create a dummy sample information table (all samples assigned to "control")
# This is used when there is no group/condition information for DESeq2
sample_info <- data.frame(
  row.names = sample_ids,
  condition = rep("control", length(sample_ids))
)

# Create a DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = merged_expr,
  colData = sample_info,
  design = ~ 1
)

# Estimate size factors for normalization
dds <- estimateSizeFactors(dds)

# Perform variance stabilizing transformation (VST)
# Setting blind = TRUE assumes no design/batch effect
vst_data <- vst(dds, blind = TRUE)

# Extract the transformed expression matrix
vst_matrix <- assay(vst_data)

# Save the VST-transformed matrix as RDS
saveRDS(vst_matrix, file = here("data", "TCGA_data", "merged_VST_data.rds"))

two_genes <- subset(vst_matrix, rownames(vst_matrix) %in% c("ENSG00000177311.11","ENSG00000101966.13"))
two_genes_t <- t(two_genes)
two_genes_t <- data.frame(two_genes_t)
correlation <- cor.test(two_genes_t$ENSG00000177311.11, two_genes_t$ENSG00000101966.13, method = "pearson")
r_value <- correlation$estimate  # Pearson
p_value <- correlation$p.value
formatted_r <- formatC(r_value, digits = 3, format = "f")
formatted_p <- format.pval(p_value, digits = 3, eps = 2.2e-16)

g <- ggplot(two_genes_t, aes(x = ENSG00000177311.11, y = ENSG00000101966.13)) +
    geom_point(alpha = 0.3, color = "#444455") +
    geom_smooth(method = "lm", color = "#CC4444", se = FALSE) +
    annotate("text",
            x = 6,
            y = 13.7,
           label = paste0("r = ", formatted_r, "\np = ", formatted_p),
           size = 10, hjust = 0) +
    labs(x = "ZBTB38 Expression",
        y = "XIAP Expression"
    ) +
    theme_bw() + ylim(c(9,14)) +
    theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24)
    )

pdf(file = here("figures", "figureA.pdf"))
print(g)
dev.off()
