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

# Задаём директорию, где находятся файлы
dir <- "/media/eternus1/projects/milina/kallisto_output_de"

### КОД ДЛЯ ВСЕХ СЕМПЛОВ ###
# Создаем вектор с именами образцов (если у вас есть файл samples, его можно использовать)
samples <- c("AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtK3_1", "AtK3_2", "AtK3_3", "AtK5_1", "AtK5_2", "AtK5_3")

files <- file.path(dir, paste0(samples, ".k31"), "abundance.h5")

# Добавляем имена к каждому файлу для идентификации
names(files) <- samples

# Загружаем данные с помощью tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
# Проверяем загруженные данные
head(txi.kallisto$counts)

# Создание таблицы с информацией об условиях эксперимента
sampleTable <- data.frame(
  sampleName = names(files),
  #  group = factor(c("infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control", "control", "infected", "infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control")),
  day = factor(c("3", "3", "3", "5", "5", "3", "3", "3", "5", "5", "5")),
  variety = factor(c("AtF", "AtF", "AtF", "AtF", "AtF", "AtK", "AtK", "AtK", "AtK", "AtK", "AtK"))
  #  variety = factor(c("At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM"))
)

###ПОСТРОЕНИЕ ДИАГРАММ ВЕННА ДЛЯ НЕНОРМАЛИЗОВАННЫХ ЗНАЧЕНИЙ 
# Создание объекта DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = sampleTable, design = ~ 1 + day + variety + day:variety)
# Создание нового фактора для комбинации variety и day
dds$group <- factor(paste0(dds$variety, dds$day))
# Обновление дизайна эксперимента
design(dds) <- ~ group
# Проведение анализа
dds <- DESeq(dds)

# Сравнение групп AtK3 и AtF5
res_AtK3_AtF5 <- results(dds, contrast = c("group", "AtK3", "AtF5"))
significant_genes_AtK3_AtF5 <- subset(res_AtK3_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtK3_AtF5 <- rownames(significant_genes_AtK3_AtF5[rowSums(counts(dds)[rownames(significant_genes_AtK3_AtF5), ] != 0) > 0, ])

# Сравнение групп AtK5 и AtF5
res_AtK5_AtF5 <- results(dds, contrast = c("group", "AtK5", "AtF5"))
significant_genes_AtK5_AtF5 <- subset(res_AtK5_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtK5_AtF5 <- rownames(significant_genes_AtK5_AtF5[rowSums(counts(dds)[rownames(significant_genes_AtK5_AtF5), ] != 0) > 0, ])

# Сравнение групп AtF3 и AtF5
res_AtF3_AtF5 <- results(dds, contrast = c("group", "AtF3", "AtF5"))
significant_genes_AtF3_AtF5 <- subset(res_AtF3_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtF3_AtF5 <- rownames(significant_genes_AtF3_AtF5[rowSums(counts(dds)[rownames(significant_genes_AtF3_AtF5), ] != 0) > 0, ])

# Сравнение групп AtK3 и AtK5
res_AtK3_AtK5 <- results(dds, contrast = c("group", "AtK3", "AtK5"))
significant_genes_AtK3_AtK5 <- subset(res_AtK3_AtK5, padj <= 0.05)
expressed_genes_AtK3_AtK5 <- rownames(significant_genes_AtK3_AtK5[rowSums(counts(dds)[rownames(significant_genes_AtK3_AtK5), ] != 0) > 0, ])


# список для диаграммы Вена
gene_lists <- list(
  "AtK3_vs_AtF5" = expressed_genes_AtK3_AtF5,
  "AtK5_vs_AtF5" = expressed_genes_AtK5_AtF5,
  "AtF3_vs_AtF5" = expressed_genes_AtF3_AtF5,
  "AtK3_vs_AtK5" = expressed_genes_AtK3_AtK5
)

png("Venn_diagram_all.png", width = 12, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "grey90", high = "red")  # Настройка цветовой шкалы
dev.off()

