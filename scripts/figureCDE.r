library(here)
library(irlba)
library(Rtsne)
library(ggplot2)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

vst_matrix <- readRDS(file = here("data", "TCGA_data", "merged_VST_data.rds"))
col_ids <- colnames(vst_matrix)


tsne_input <- as.matrix(t(vst_matrix))

# PCA
set.seed(111)
pca_result <- irlba(tsne_input, nv = 50)
pca_scores <- pca_result$u %*% diag(pca_result$d)

# t-SNE
set.seed(111)
tsne_result <- Rtsne(pca_scores, dims = 2, perplexity = 30, verbose = TRUE)
sample_ids <- rownames(tsne_input)

tsne_df <- as.data.frame(tsne_result$Y)
rownames(tsne_df) <- sample_ids
#Extract cancer type
cancer_types <- sapply(strsplit(rownames(tsne_df), "\\."), `[`, 1)
tsne_df <- cbind(tsne_df, cancer_types)

colnames(tsne_df) <- c("tSNE1", "tSNE2", "tumor_type")
saveRDS(tsne_df, file = here("data", "TCGA_data", "tSNE_data.rds"))

#FigureC
g <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = tumor_type)) +
  geom_point(alpha = 0.7, size = 0.3) +
  labs(x = "t-SNE 1", y = "t-SNE 2") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "right",
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12)
) +
  scale_color_manual(
    values = rainbow(length(unique(tsne_df$tumor_type))),
    guide = guide_legend(override.aes = list(size = 3))
  )

pdf(file = here("figures","figureC.pdf"),height = 5, width = 7.5)
print(g)
dev.off()

#FigureD

two_genes <- subset(vst_matrix, rownames(vst_matrix) %in% c("ENSG00000177311.11", "ENSG00000101966.13"))
two_t <- t(two_genes)

two_genes_tsne <- data.frame(tsne_df,
    t(two_genes)
    )

colnames(two_genes_tsne) <- c("tSNE1", "tSNE2", "tumor_type", "XIAP", "ZBTB38")

#ZBTB38
g <- ggplot(two_genes_tsne, aes(x = tSNE1, y = tSNE2, color = ZBTB38)) +
  geom_point(alpha = 0.5, size = 0.3) +
  labs(x = "tSNE1", y = "tSNE2") +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = "right",
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12)) +
  scale_color_viridis_c(option = "magma",  limits = c(4, 15))

pdf(file = here("figures","figureD.pdf"),height = 5, width = 5.5)
print(g)
dev.off()

#FigureE (XIAP)
g <- ggplot(two_genes_tsne, aes(x = tSNE1, y = tSNE2, color = XIAP)) +
  geom_point(alpha = 0.5, size = 0.3) +
  labs(x = "tSNE1", y = "tSNE2") +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = "right",
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12)) +
  scale_color_viridis_c(option = "magma",  limits = c(7, 14))

pdf(file = here("figures","figureE.pdf"), height = 5, width = 5.5)
print(g)
dev.off()
