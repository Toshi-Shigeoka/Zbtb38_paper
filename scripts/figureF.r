library(here)
library(ggplot2)
library(ggrepel)
library(dplyr)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

vst_matrix <- readRDS(file = here("data", "TCGA_data", "merged_VST_data.rds"))
target_gene <- "ENSG00000177311.11"
target_expr <- vst_matrix[target_gene, ]
cor_values <- apply(vst_matrix, 1, function(x) cor(x, target_expr, method = "pearson", use = "pairwise.complete.obs"))

#ranking
cor_df <- tibble(Gene = names(cor_values), Correlation = cor_values) %>%
  arrange(desc(Correlation))

#save
write.table(cor_df, file = here("data", "TCGA_data", "all_cancer_correlation_withZBTB38.txt"), sep = "\t", quote = F)

#"all_cancer_correlation_withZBTB38.txt" was modified manually to mark apoptotic genes

data <- read.table(file = here("data", "TCGA_data", "all_cancer_APOPTOTIC_correlation_withZBTB38_modified.txt"),sep = "\t", header = T)
data <- data.frame(rownames(data), data)
colnames(data) <- c("rank", "symbol", "score", "gene", "dot","class")
class(data$rank) <- "numeric"


g <- ggplot(data, aes(x=rank, y=score, label = gene)) + geom_abline(slope=0, intercept=0, colour='#BBBBBB') + geom_vline(xintercept=1, linetype="dashed") +  geom_vline(xintercept = 60660, linetype="dashed")+ 
  geom_line(color = "#4466CC", size = 1.5) +
  ylim(c(-0.6, 0.6))+ theme_bw() +
  geom_point(aes(x = rank, y = dot), color = "#999999", size = 0.3) +
  scale_x_continuous(limits = c(-10000, 70000),breaks = c(1,20000,40000,60660)) +
  geom_text_repel(aes(label = gene, color = class),force= 30,force_pull = 1,fontface = "bold",max.overlaps = Inf, point.padding = 0.1,box.padding = 0.5,min.segment.length = 0.1, size = 5.2) +
  theme(legend.position = "none",
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 16)) +
scale_color_manual(values = c("Apoptotic" = "#BB3333", "Anti_apoptotic" = "#11AA33"))
pdf(file = here("figures", "figureF.pdf"), width = 7, height = 5)
print(g)
dev.off()

