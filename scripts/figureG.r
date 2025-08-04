library(here)
library(ggplot2)
library(ggpubr)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

vst_matrix <- readRDS(file = here("data", "TCGA_data", "merged_VST_data.rds"))
merged_clinical <- readRDS(file = here("data", "TCGA_data", "merged_clinical.rds"))

two_genes <- subset(vst_matrix, rownames(vst_matrix) %in% c("ENSG00000177311.11","ENSG00000101966.13"))
two_genes_t <- t(two_genes)
two_genes_t <- data.frame(two_genes_t)
two_genes2 <- data.frame(barcode = sub(".*\\.", "", rownames(two_genes_t)), two_genes_t)


tumor_grade <- merged_clinical$tumor_grade
grade_major <- factor(tumor_grade, levels = c("G1", "G2", "G3", "G4"))
grade_major <- data.frame(barcode = rownames(merged_clinical), grade_major)
merged_grade <- merge(grade_major, two_genes2, by = 1)
colnames(merged_grade) <- c("barcode", "Tumor_grade", "XIAP", "ZBTB38")
saveRDS(merged_grade, file = here("data", "TCGA_data", "grade_two_genes.rds"))

#removing NA
merged_grade_sub <- subset(merged_grade, Tumor_grade %in% c("G1", "G2", "G3", "G4"))
#for statistic comparison
my_comparisons <- list(
  c("G1", "G2"),
  c("G3", "G4"),
  c("G2", "G4"),
  c("G1", "G4")
)


g <- ggplot(merged_grade_sub, aes(x = Tumor_grade, y = ZBTB38)) +
  geom_violin(fill = "#4488BB", alpha = 0.5) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5, color = "#005588") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5, color = "black") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "black") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = my_comparisons,
                     label = "p.format",
                     tip.length = 0.05,
                     label.y = c(15.5, 16, 16.5, 17),
                     size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12)
  )

pdf(file = here("figures", "figureG.pdf"), height = 4, width = 5)
print(g)
dev.off()
