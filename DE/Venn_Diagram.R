library(tximport)
library(readr)
library(tximportData)
library(rhdf5)
library(BiocManager)
library(DESeq2)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplot2)


# загрузка векторов в файл
load("sig_genes_AtK3_AtF3.RData")
load("sig_genes_AtK5_AtF5.RData")
load("sig_genes_AtF3_AtF5.RData")
load("up_sig_genes_AtK3_AtF3.RData")
load("up_sig_genes_AtK5_AtF5.RData")
load("up_sig_genes_AtF3_AtF5.RData")
load("down_sig_genes_AtK3_AtF3.RData")
load("down_sig_genes_AtK5_AtF5.RData")
load("down_sig_genes_AtF3_AtF5.RData")


### Построение диаграмм Венна
# список для диаграммы Вена
gene_lists <- list(
  "AtK3_vs_AtF3" = sig_genes_lfc,
  "AtK5_vs_AtF5" = sig_genes_lfc_5,
  "AtF3_vs_AtF5" = sig_genes_lfc3_5
)

png("Venn_diagram_all.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()




### для up генов
# список для диаграммы Вена
gene_lists_up <- list(
  "AtK3_vs_AtF3" = up_genes_lfc,
  "AtK5_vs_AtF5" = up_genes_lfc_5,
  "AtF3_vs_AtF5" = up_genes_lfc3_5
)

png("Venn_diagram_up.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_up) +
  scale_fill_gradient(low = "grey90", high = "red") +
  ggtitle("Up-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()


### для down генов
# список для диаграммы Вена
gene_lists_down <- list(
  "AtK3_vs_AtF3" = down_genes_lfc,
  "AtK5_vs_AtF5" = down_genes_lfc_5,
  "AtF3_vs_AtF5" = down_genes_lfc3_5
)

png("Venn_diagram_down.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_down) +
  scale_fill_gradient(low = "grey90", high = "blue") +
  ggtitle("Down-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()


