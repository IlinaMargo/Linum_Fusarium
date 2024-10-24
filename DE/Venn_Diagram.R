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



###ИЗВЛЕЧЕНИЕ ГЕНОВ ИЗ КОНТРОЛЬНЫХ ОБРАЗЦОВ, ЧТОБЫ ВЫЧЕСТЬ ИЗ ОСТАЛЬНЫХ ХИТМАПОВ
samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]

# Извлечение нормализованных данных для AtK3 и AtK5
samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]
counts_AtK3_AtK5 <- counts(dds)[, samples_AtK3_AtK5]

# Фильтрация генов НЕНОРМАЛИЗОВАННЫХ ДАННЫХ, у которых есть экспрессия хотя бы в одном образце (ненулевое значение)
expressed_genes_AtK3_AtK5 <- rownames(counts_AtK3_AtK5[rowSums(counts_AtK3_AtK5 != 0) > 0, ])

# Например, исключаем эти гены из объекта DESeq2
dds <- dds[!rownames(dds) %in% expressed_genes_AtK3_AtK5, ]


### Построение диаграмм Венна
# Сравнение групп AtK3 и AtF5
res_AtK3_AtF3 <- results(dds, contrast = c("group", "AtK3", "AtF3"))
significant_genes_AtK3_AtF3 <- subset(res_AtK3_AtF3, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtK3_AtF3 <- rownames(significant_genes_AtK3_AtF3[rowSums(counts(dds)[rownames(significant_genes_AtK3_AtF3), ] != 0) > 0, ])

# Сравнение групп AtK5 и AtF5
res_AtK5_AtF5 <- results(dds, contrast = c("group", "AtK5", "AtF5"))
significant_genes_AtK5_AtF5 <- subset(res_AtK5_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtK5_AtF5 <- rownames(significant_genes_AtK5_AtF5[rowSums(counts(dds)[rownames(significant_genes_AtK5_AtF5), ] != 0) > 0, ])

# Сравнение групп AtF3 и AtF5
res_AtF3_AtF5 <- results(dds, contrast = c("group", "AtF3", "AtF5"))
significant_genes_AtF3_AtF5 <- subset(res_AtF3_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5)
expressed_genes_AtF3_AtF5 <- rownames(significant_genes_AtF3_AtF5[rowSums(counts(dds)[rownames(significant_genes_AtF3_AtF5), ] != 0) > 0, ])


# список для диаграммы Вена
gene_lists <- list(
  "AtK3_vs_AtF3" = expressed_genes_AtK3_AtF3,
  "AtK5_vs_AtF5" = expressed_genes_AtK5_AtF5,
  "AtF3_vs_AtF5" = expressed_genes_AtF3_AtF5
)

png("Venn_diagram_all.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()




### для up генов
## AtF3 vs AtK3
up_significant_genes_AtK3_AtF3 <- subset(res_AtK3_AtF3, padj <= 0.05 & log2FoldChange >= 1.5)
up_expressed_genes_AtK3_AtF3 <- rownames(up_significant_genes_AtK3_AtF3[rowSums(counts(dds)[rownames(up_significant_genes_AtK3_AtF3), ] != 0) > 0, ])

## AtF5 vs AtK5
up_significant_genes_AtK5_AtF5 <- subset(res_AtK5_AtF5, padj <= 0.05 & log2FoldChange >= 1.5)
up_expressed_genes_AtK5_AtF5 <- rownames(up_significant_genes_AtK5_AtF5[rowSums(counts(dds)[rownames(up_significant_genes_AtK5_AtF5), ] != 0) > 0, ])


## AtF3 vs AtF5
up_significant_genes_AtF3_AtF5 <- subset(res_AtF3_AtF5, padj <= 0.05 & log2FoldChange >= 1.5)
up_expressed_genes_AtF3_AtF5 <- rownames(up_significant_genes_AtF3_AtF5[rowSums(counts(dds)[rownames(up_significant_genes_AtF3_AtF5), ] != 0) > 0, ])

# список для диаграммы Вена
gene_lists_up <- list(
  "AtK3_vs_AtF3" = up_expressed_genes_AtK3_AtF3,
  "AtK5_vs_AtF5" = up_expressed_genes_AtK5_AtF5,
  "AtF3_vs_AtF5" = up_expressed_genes_AtF3_AtF5
)

png("Venn_diagram_up.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_up) +
  scale_fill_gradient(low = "grey90", high = "red") +
  ggtitle("Up-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()


### для down генов
## AtF3 vs AtK3
down_significant_genes_AtK3_AtF3 <- subset(res_AtK3_AtF3, padj <= 0.05 & log2FoldChange <= - 1.5)
down_expressed_genes_AtK3_AtF3 <- rownames(down_significant_genes_AtK3_AtF3[rowSums(counts(dds)[rownames(down_significant_genes_AtK3_AtF3), ] != 0) > 0, ])

## AtF5 vs AtK5
down_significant_genes_AtK5_AtF5 <- subset(res_AtK5_AtF5, padj <= 0.05 & log2FoldChange <= - 1.5)
down_expressed_genes_AtK5_AtF5 <- rownames(down_significant_genes_AtK5_AtF5[rowSums(counts(dds)[rownames(down_significant_genes_AtK5_AtF5), ] != 0) > 0, ])

## AtF3 vs AtF5
down_significant_genes_AtF3_AtF5 <- subset(res_AtF3_AtF5, padj <= 0.05 & log2FoldChange <= - 1.5)
down_expressed_genes_AtF3_AtF5 <- rownames(down_significant_genes_AtF3_AtF5[rowSums(counts(dds)[rownames(down_significant_genes_AtF3_AtF5), ] != 0) > 0, ])

# список для диаграммы Вена
gene_lists_down <- list(
  "AtK3_vs_AtF3" = down_expressed_genes_AtK3_AtF3,
  "AtK5_vs_AtF5" = down_expressed_genes_AtK5_AtF5,
  "AtF3_vs_AtF5" = down_expressed_genes_AtF3_AtF5
)

png("Venn_diagram_down.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_down) +
  scale_fill_gradient(low = "grey90", high = "blue") +
  ggtitle("Down-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()



# Проверка почему up+down не равно all:
status_all_genes <- data.frame(
  Gene = intersection_all,
  AtK3_vs_AtF3 = ifelse(res_AtK3_AtF3[intersection_all, "log2FoldChange"] >= 1.5, "up", ifelse(res_AtK3_AtF3[intersection_all, "log2FoldChange"] <= -1.5, "down", "none")),
  AtK5_vs_AtF5 = ifelse(res_AtK5_AtF5[intersection_all, "log2FoldChange"] >= 1.5, "up", ifelse(res_AtK5_AtF5[intersection_all, "log2FoldChange"] <= -1.5, "down", "none")),
  AtF3_vs_AtF5 = ifelse(res_AtF3_AtF5[intersection_all, "log2FoldChange"] >= 1.5, "up", ifelse(res_AtF3_AtF5[intersection_all, "log2FoldChange"] <= -1.5, "down", "none"))
)

print(status_all_genes)
