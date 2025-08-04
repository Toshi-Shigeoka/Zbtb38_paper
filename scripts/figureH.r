library(here)
library(readr)
library(ggplot2)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

#TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes should be downloaded from
#UCSC Xena Browser: https://xenabrowser.net/datapages/

#cnv data
cnv <- read_tsv(file = here("data","TCGA_data","TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"), col_names = TRUE)
cnv <- as.matrix(cnv)
rownames(cnv) <- cnv[,1]
cnv <- cnv[, -1]
class(cnv) <- "numeric"
cnv_t <- t(cnv)

#grade data
grade <- readRDS(file = here("data","TCGA_data","grade_two_genes.rds"))
grade <- grade[,1:2]
grade$id <- substr(grade$barcode, 1, 15)
merged <- merge(grade, cnv_t, by.x = 3, by.y = 0)

G34 <- subset(merged, Tumor_grade %in% c("G3", "G4"))
G34_cnv <- G34[,6:ncol(G34)]

#fix order
score_levels <- c(-2, -1, 0, 1, 2)


##CNV score of ZBTB38
zbt_vec <- as.numeric(G34_cnv[, "ZBTB38"])
zbt_table <- table(factor(zbt_vec, levels = score_levels))
zbt_frac <- zbt_table / sum(zbt_table)
df_zbt <- data.frame(
  Group = "ZBTB38",
  CNV = names(zbt_frac),
  Fraction = as.numeric(zbt_frac)
)

##CNV score of the other genes
all_vec <- as.numeric(as.matrix(G34_cnv))
all_table <- table(factor(all_vec, levels = score_levels))
other_table <- all_table - zbt_table
other_frac <- other_table / sum(other_table)
df_other <- data.frame(
  Group = "Others",
  CNV = names(other_frac),
  Fraction = as.numeric(other_frac)
)

## merge dfs
df_plot <- rbind(df_zbt, df_other)
df_plot$CNV <- factor(df_plot$CNV, levels = score_levels)
df_plot$Group <- factor(df_plot$Group, levels = c("Others", "ZBTB38"))


# Chi-squared test
contingency <- rbind(
  ZBTB38 = as.numeric(zbt_table),
  Others = as.numeric(other_table)
)
colnames(contingency) <- score_levels

chisq.test(contingency)


#    Pearson's Chi-squared test

#data:  contingency
#X-squared = 336.82, df = 4, p-value < 2.2e-16



g <- ggplot(df_plot, aes(x = Group, y = Fraction, fill = CNV)) +
  geom_bar(stat = "identity", width = 0.5) +           # ← 棒を細く
  scale_fill_manual(values = c("-2" = "#2222DD", "-1" = "#7777FF",
                               "0" = "#CCCCCC", "1" = "#FF7777", "2" = "#DD2222")) +
  labs(title = "      p < 2.2e-16",
       x = "", y = "Fraction") +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
  geom_hline(yintercept = 1, linetype = "solid", color = "gray30") +
  theme_minimal() +
  theme(
  axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1, size = 12),
  axis.text.y = element_text(size = 11)
  )

pdf(file = here("figures", "figureH.pdf"), width = 3, height = 4)
print(g)
dev.off()
