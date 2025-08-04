library(here)
library(ggplot2)
library(scales)
library(dplyr)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

vst_matrix <- readRDS(file = here("data", "TCGA_data", "merged_VST_data.rds"))
col_ids <- colnames(vst_matrix)
#Extract cancer type
cancer_types <- sapply(strsplit(col_ids, "\\."), `[`, 1)
two_genes <- subset(vst_matrix, rownames(vst_matrix) %in% c("ENSG00000177311.11", "ENSG00000101966.13"))

two_genes_types <- data.frame(
  SampleID = col_ids,
  cancer_type = cancer_types,
  t(two_genes)
)

gene1 <- "ENSG00000177311.11"  # ZBTB38
gene2 <- "ENSG00000101966.13"  # XIAP

cor_results <- data.frame(
  cancer_type = character(),
  r_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (cancer in unique(two_genes_types$cancer_type)) {
  subset_data <- two_genes_types[two_genes_types$cancer_type == cancer, ]
  
  if (nrow(subset_data) >= 3) {  # minimum sample numbers
    x <- as.numeric(subset_data[[gene1]])
    y <- as.numeric(subset_data[[gene2]])
    
    cor_test <- cor.test(x, y, method = "pearson")
    
    cor_results <- rbind(cor_results, data.frame(
      cancer_type = cancer,
      r_value = round(cor_test$estimate, 3),
      p_value = signif(cor_test$p.value, 3)
    ))
  }
}

# -log10(p-value)
cor_results$minus_log10_p <- -log10(cor_results$p_value)
cor_results <- cor_results[order(cor_results$minus_log10_p), ]

cor_results$cancer_type <- factor(cor_results$cancer_type, levels = cor_results$cancer_type)

g <- ggplot(cor_results, aes(x = cancer_type)) +
  geom_bar(aes(y = minus_log10_p), stat = "identity", fill = "steelblue") +
  geom_line(aes(y = r_value * max(minus_log10_p), group = 1), color = "#BB5555", size = 0.8) +
  geom_point(aes(y = r_value * max(minus_log10_p)), color = "#BB5555", size = 1) +
  scale_y_continuous(
    name = expression(-log[10](p~value)),
    sec.axis = sec_axis(~ . / max(cor_results$minus_log10_p), name = "Pearson correlation (r)")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 13),
    axis.text.y = element_text(hjust = 1, size = 12)) +
  labs(
    x = "Cancer Type"
  ) + coord_flip()

pdf(file = here("figures", "figureB.pdf"), width = 7, height = 6)
print(g)
dev.off()
