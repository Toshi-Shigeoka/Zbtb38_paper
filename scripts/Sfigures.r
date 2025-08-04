library(here)
library(readr)
library(ggplot2)
library(ggrepel)
library(Rtsne)
library(dplyr)
library(viridis)

if(!dir.exists(here("figures"))) dir.create(here("figures"))

data <- read_tsv(here("data", "MERAV_data", "all_genes_batch_adjusted_RAW.txt"))
data <- as.matrix(data)
tumor_data <- data[, grepl("Primary\\.Tumor", colnames(data))]
rownames(tumor_data) <- data[,1]


expr_matrix <- apply(tumor_data[,-1], 2, as.numeric)
rownames(expr_matrix) <- rownames(tumor_data)

expr_matrix_log <-  log((expr_matrix+1),2)
expr_matrix_t <- t(expr_matrix_log)
expr_matrix_t <- expr_matrix_t[, colSums(is.na(expr_matrix_t)) == 0]



two_genes_df <- data.frame(
  Sample = rownames(expr_matrix_t),
  ZBTB38 = expr_matrix_t[, "ZBTB38"],
  XIAP = expr_matrix_t[, "XIAP"]
)

two_genes_df$Type <- sapply(strsplit(two_genes_df$Sample, "_"), `[`, 2)
two_genes_df$Tissue <- sapply(strsplit(two_genes_df$Sample, "_"), `[`, 1)
colon <- subset(two_genes_df, Tissue == "Colon")
breast <- subset(two_genes_df, Tissue == "Breast")
haematopoietic <- subset(two_genes_df, Tissue == "Haematopoietic.And.Lymphoid.Tissue")

correlation <- cor.test(two_genes_df$ZBTB38, two_genes_df$XIAP, method = "pearson")


r_value <- correlation$estimate  # Pearson
p_value <- correlation$p.value
formatted_r <- formatC(r_value, digits = 3, format = "f")
formatted_p <- format.pval(p_value, digits = 3, eps = 2.2e-16)

g <- ggplot(two_genes_df, aes(x = ZBTB38, y = XIAP)) +
    geom_point(alpha = 0.3, color = "#444455") +
    geom_smooth(method = "lm", color = "#CC4444", se = FALSE) +
    annotate("text",
            x = 5,
            y = 9.9,
           label = paste0("r = ", formatted_r, "\np = ", formatted_p),
           size = 10, hjust = 0) +
    labs(title = "All",
        x = "ZBTB38 Expression",
        y = "XIAP Expression"
    ) +
    theme_bw()  + xlim(5, 10) + ylim(6.5, 10) +
    theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 25)
    )

pdf(file = here("figures", "SfigureA.pdf"))
print(g)
dev.off()

#Colon

correlation <- cor.test(colon$ZBTB38, colon$XIAP, method = "pearson")

r_value <- correlation$estimate  # Pearson
p_value <- correlation$p.value
formatted_r <- formatC(r_value, digits = 3, format = "f")
formatted_p <- format.pval(p_value, digits = 3, eps = 2.2e-16)

g <- ggplot(colon, aes(x = ZBTB38, y = XIAP)) +
    geom_point(alpha = 0.3, color = "#444455") +
    geom_smooth(method = "lm", color = "#CC4444", se = FALSE) +
    annotate("text",
            x = 6.5,
            y = 9.9,
           label = paste0("r = ", formatted_r, "\np = ", formatted_p),
           size = 10, hjust = 0) +
    labs(title = "Colon",
        x = "ZBTB38 Expression",
        y = "XIAP Expression"
    ) +
    theme_bw()  + ylim(7, 10) +
    theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 25)
    )

pdf(file = here("figures", "SfigureB.pdf"))
print(g)
dev.off()

#haematopoietic

correlation <- cor.test(haematopoietic$ZBTB38, haematopoietic$XIAP, method = "pearson")

r_value <- correlation$estimate  # Pearson
p_value <- correlation$p.value
formatted_r <- formatC(r_value, digits = 3, format = "f")
formatted_p <- format.pval(p_value, digits = 3, eps = 2.2e-16)

g <- ggplot(haematopoietic, aes(x = ZBTB38, y = XIAP)) +
    geom_point(alpha = 0.3, color = "#444455") +
    geom_smooth(method = "lm", color = "#CC4444", se = FALSE) +
    annotate("text",
            x = 5,
            y = 9.3,
           label = paste0("r = ", formatted_r, "\np = ", formatted_p),
           size = 10, hjust = 0) +
    labs(title = "Haematopoietic and Lymphoid",
        x = "ZBTB38 Expression",
        y = "XIAP Expression"
    ) +
    theme_bw()  +  ylim(6.5, 9.5) +
    theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 25)
    )

pdf(file = here("figures", "SfigureC.pdf"))
print(g)
dev.off()

#breast

correlation <- cor.test(breast$ZBTB38, breast$XIAP, method = "pearson")

r_value <- correlation$estimate  # Pearson
p_value <- correlation$p.value
formatted_r <- formatC(r_value, digits = 3, format = "f")
formatted_p <- format.pval(p_value, digits = 3, eps = 2.2e-16)

g <- ggplot(breast, aes(x = ZBTB38, y = XIAP)) +
    geom_point(alpha = 0.3, color = "#444455") +
    geom_smooth(method = "lm", color = "#CC4444", se = FALSE) +
    annotate("text",
            x = 6,
            y = 9.3,
           label = paste0("r = ", formatted_r, "\np = ", formatted_p),
           size = 10, hjust = 0) +
    labs(title = "Breast",
        x = "ZBTB38 Expression",
        y = "XIAP Expression"
    ) +
    theme_bw()  +  ylim(7, 9.5) +
    theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 25)
    )

pdf(file = here("figures", "SfigureD.pdf"))
print(g)
dev.off()

