# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = '3.18')
# BiocManager::install("tximport")
# BiocManager::install(c("rhdf5", "DESeq2", "tximportData"))
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("ashr")
# BiocManager::install("vsn")
# BiocManager::install("ComplexHeatmap")
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
library(gplots)

# Задаём директорию, где находятся файлы
dir <- "/media/eternus1/projects/milina/kallisto_LU_DE/"

### КОД ДЛЯ ВСЕХ СЕМПЛОВ ###
# Создаем вектор с именами образцов (если у вас есть файл samples, его можно использовать)
samples <- c("AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtF5_3", "AtK3_1", "AtK3_2", "AtK3_3", "AtK5_1", "AtK5_2", "AtK5_3")

files <- file.path(dir, paste0(samples, ".k31"), "abundance.h5")

# Добавляем имена к каждому файлу для идентификации
names(files) <- samples

# Загружаем данные с помощью tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance = "lengthScaledTPM", varReduce = FALSE) 
# Проверяем загруженные данные
txi.kallisto$counts

# Создание таблицы с информацией об условиях эксперимента
sampleTable <- data.frame(
  sampleName = names(files),
  #  group = factor(c("infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control", "control", "infected", "infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control")),
  day = factor(c("3", "3", "3", "5", "5", "3", "3", "3", "5", "5", "5")),
  variety = factor(c("AtF", "AtF", "AtF", "AtF", "AtF", "AtK", "AtK", "AtK", "AtK", "AtK", "AtK"))
  #  variety = factor(c("At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM"))
)

# Создание объекта DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = sampleTable, design = ~ 1 + day + variety + day:variety)
# Создание нового фактора для комбинации variety и day
dds$group <- factor(paste0(dds$variety, dds$day))
# Обновление дизайна эксперимента
design(dds) <- ~ group
# Проведение анализа
dds <- DESeq(dds)

# Нормализация данных с использованием Variance Stabilizing Transformation 
vsd <- vst(dds, blind=FALSE)
# Нормализация данных с использованием Regularized Log Transformation
rld <- rlog(dds, blind = FALSE)

# Извлечение нормализованных данных для всех генов
all_counts <- assay(rld)

# ## посмотрим на распределение lfc во всех образцах
# hmcol <- colorRampPalette(brewer.pal(9, "RdPu"))(100)
# rlogMat <- assay(rld)
# distsRL <- dist(t(assay(rld)))
# mat <- as.matrix(distsRL)
# hc <- hclust(distsRL)
# 
# png("heatmap_exp_through_all_samp.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# heatmap.2(
#   mat,
#   Rowv = as.dendrogram(hc),
#   symm = TRUE,
#   trace = "none",
#   col = rev(hmcol),
#   margin = c(13, 13)
# )
# dev.off()


# ###ИЗВЛЕЧЕНИЕ ГЕНОВ ИЗ КОНТРОЛЬНЫХ ОБРАЗЦОВ, ЧТОБЫ ВЫЧЕСТЬ ИЗ ОСТАЛЬНЫХ ХИТМАПОВ
# samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]
# 
# # Извлечение нормализованных данных для AtK3 и AtK5
# samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]
# counts_AtK3_AtK5 <- counts(dds)[, samples_AtK3_AtK5]
# 
# # Фильтрация генов НЕНОРМАЛИЗОВАННЫХ ДАННЫХ, у которых есть экспрессия хотя бы в одном образце (ненулевое значение)
# expressed_genes_AtK3_AtK5 <- rownames(counts_AtK3_AtK5[rowSums(counts_AtK3_AtK5) >= 200, ])
# nrow(dds)
# # Например, исключаем эти гены из объекта DESeq2
# dds <- dds[!rownames(dds) %in% expressed_genes_AtK3_AtK5, ]
# rld <- rld[!rownames(rld) %in% expressed_genes_AtK3_AtK5, ]
# nrow(dds)
# # Извлечение нормализованных данных для всех генов
# all_counts <- assay(rld)

### ХИТМАПЫ ДЛЯ НЕФИЛЬТРОВАННОГО
# ## Выбор образцов для AtF3 и AtK3
# samples_AtF3_AtK3 <- colnames(dds)[dds$group %in% c("AtF3", "AtK3")]

# # Центрирование значений по строкам (генам)
# centered_all_counts <- t(scale(t(all_counts[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))
# 
#  # Плавный переход от синего к красному через белый
# color_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
# 
# df <- as.data.frame(colData(dds)[,c("variety", "day")])
# 
# # Построение тепловой карты для всех генов
# png("heatmap_NF_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_all_counts, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtK3, ],
#          color = color_palette,
#          main = "Non-filtered genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
# )
# dev.off()




# ## Выбор образцов для AtF5 и AtK5
# samples_AtF5_AtK5 <- colnames(dds)[dds$group %in% c("AtF5", "AtK5")]
# 
# # Центрирование значений по строкам (генам)
# centered_all_counts_5 <- t(scale(t(all_counts[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для всех генов
# png("heatmap_NF_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_all_counts_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF5_AtK5, ],
#          color = color_palette,
#          main = "Non-filtered genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
# )
# dev.off()
# 
# ## Выбор образцов для AtF3 и AtF5
# samples_AtF3_AtF5 <- colnames(dds)[dds$group %in% c("AtF3", "AtF5")]
# 
# # Центрирование значений по строкам (генам)
# centered_all_counts_5_3 <- t(scale(t(all_counts[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для всех генов
# png("heatmap_NF_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_all_counts_5_3, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtF5, ],
#          color = color_palette,
#          main = "Non-filtered genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
# )
# dev.off()


# plotMA(res, ylim=c(-2,2), alpha = 0.5, main = "MA plot, p-val = 0.5")
# 
# #resultsNames(dds)
# 
# # сохраняем идентификаторы генов в индексы
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]
# 
# count_plot <- plotCounts(dds, gene=which.min(res$padj), intgroup="variety", returnData=TRUE)
# # Создание графика
# p <- ggplot(count_plot, aes(x=variety, y=count)) + 
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   scale_y_log10(breaks=c(25,100,400))
# 
# mcols(res)$description
# 
# # Построение графика
# meanSdPlot(assay(ntd))
# 
# ## рисуем pca plot главные компоненты для матрицы
# pcaData <- plotPCA(rld, intgroup=c("variety", "day"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# ggplot(pcaData, aes(PC1, PC2, color=variety, shape=day)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()



### ХИТМАПЫ ФИЛЬТРОВАННЫХ ЗНАЧЕНИЙ 
## Сравнение групп AtF3 и AtK3
res_AtF3_AtK3 <- results(dds, contrast = c("group", "AtF3", "AtK3"))
saveRDS(res_AtF3_AtK3, file = "res_AtF3_AtK3.rds")

# Преобразование объекта результатов в датафрейм
df_AtF3_AtK3 <- as.data.frame(res_AtF3_AtK3)

png("hist_val_AtF3_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_AtF3_AtK3 , aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) AtF3 vs AtK3",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_AtF3_AtK3 $log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_AtF3_AtK3 $log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 2
sig_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
# Сохранение вектора в файл
save(sig_genes_lfc, file = "sig_genes_AtK3_AtF3.RData")

# # Извлечение нормализованных данных для значимых генов
# filtered_counts_lfc <- assay(rld)[sig_genes_lfc, ]
# 
# # Выбор образцов для AtF3 и AtK3
# samples_AtF3_AtK3 <- colnames(dds)[dds$group %in% c("AtF3", "AtK3")]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_lfc <- t(scale(t(filtered_counts_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для значимых генов с |LFC| > 1.5
# png("heatmap_F_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_lfc, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtK3, ],
#          color = color_palette,
#          main = "Differential expressed genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
# )
# dev.off()





## Сравнение групп AtF5 и AtK5
res_AtF5_AtK5 <- results(dds, contrast = c("group", "AtF5", "AtK5"))
res_AtF3_AtK3 <- readRDS("res_AtF3_AtK3.rds")
# Преобразование объекта результатов в датафрейм
df_AtF5_AtK5 <- as.data.frame(res_AtF5_AtK5)

png("hist_val_AtF5_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_AtF5_AtK5, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) AtF5 vs AtK5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_AtF5_AtK5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_AtF5_AtK5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj <= 0.05)) # & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc_5, file = "sig_genes_AtK5_AtF5.RData")

# # Извлечение нормализованных данных для значимых генов
# filtered_counts_lfc_5 <- assay(rld)[sig_genes_lfc_5, ]
# 
# # Выбор образцов для AtF5 и AtK5
# samples_AtF5_AtK5 <- colnames(dds)[dds$group %in% c("AtF5", "AtK5")]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_lfc_5 <- t(scale(t(filtered_counts_lfc[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для значимых генов с |LFC| > 1.5
# png("heatmap_F_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_lfc_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF5_AtK5, ],
#          color = color_palette,
#          main = "Differential expressed genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
# )
# dev.off()

## Сравнение групп AtF5 и AtK5
res_AtF5_AtK5 <- results(dds, contrast = c("group", "AtF5", "AtK5"))
# Преобразование объекта результатов в датафрейм
saveRDS(res_AtF5_AtK5, file = "res_AtF5_AtK5.rds")

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc_5, file = "sig_genes_AtK5_AtF5.RData")

# # Извлечение нормализованных данных для значимых генов
# filtered_counts_lfc_5 <- assay(rld)[sig_genes_lfc_5, ]
# 
# # Выбор образцов для AtF5 и AtK5
# samples_AtF5_AtK5 <- colnames(dds)[dds$group %in% c("AtF5", "AtK5")]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_lfc_5 <- t(scale(t(filtered_counts_lfc[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для значимых генов с |LFC| > 1.5
# png("heatmap_F_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_lfc_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF5_AtK5, ],
#          color = color_palette,
#          main = "Differential expressed genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
# )
# dev.off()


## Сравнение групп AtF3 и AtF5
res_AtF3_AtF5 <- results(dds, contrast = c("group", "AtF3", "AtF5"))
saveRDS(res_AtF3_AtF5, file = "res_AtF3_AtF5.rds")

# Преобразование объекта результатов в датафрейм
df_AtF3_AtF5 <- as.data.frame(res_AtF3_AtF5)
head(df_AtF3_AtF5)
png("hist_val_AtF3_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_AtF3_AtF5, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) AtF3 vs AtF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_AtF3_AtF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_AtF3_AtF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj <= 0.05  & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc3_5, file = "sig_genes_AtF3_AtF5.RData")

# # Извлечение нормализованных данных для значимых генов
# filtered_counts_lfc3_5 <- assay(rld)[sig_genes_lfc3_5, ]
# 
# # Выбор образцов для AtF3 и AtK3
# samples_AtF3_AtF5 <- colnames(dds)[dds$group %in% c("AtF3", "AtF5")]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_lfc3_5 <- t(scale(t(filtered_counts_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для значимых генов с |LFC| > 1.5
# png("heatmap_F_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_lfc3_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtF5, ],
#          color = color_palette,
#          main = "Differential expressed genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
# )
# dev.off()





###ХИТМАПЫ ДЛЯ UP ГЕНОВ
## AtF3_vs_AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc, file = "up_sig_genes_AtK3_AtF3.RData")

# # Извлечение нормализованных данных для up-regulated генов
# filtered_counts_up_lfc <- assay(rld)[up_genes_lfc, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_up_lfc <- t(scale(t(filtered_counts_up_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для up-regulated генов с LFC > 1.5
# png("heatmap_UP_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_up_lfc, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtK3, ],
#          color = color_palette,
#          main = "Up-regulated Genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
# )
# dev.off()




## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc_5, file = "up_sig_genes_AtK5_AtF5.RData")

# # Извлечение нормализованных данных для up-regulated генов
# filtered_counts_up_lfc_5 <- assay(rld)[up_genes_lfc_5, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_up_lfc_5 <- t(scale(t(filtered_counts_up_lfc_5[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для up-regulated генов с LFC > 1.5
# png("heatmap_UP_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_up_lfc_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF5_AtK5, ],
#          color = color_palette,
#          main = "Up-regulated Genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
# )
# dev.off()


## Сравнение групп AtF3 и AtF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc3_5, file = "up_sig_genes_AtF3_AtF5.RData")

# # Извлечение нормализованных данных для up-regulated генов
# filtered_counts_up_lfc3_5 <- assay(rld)[up_genes_lfc3_5, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_up_lfc3_5 <- t(scale(t(filtered_counts_up_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для up-regulated генов с LFC > 1.5
# png("heatmap_UP_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_up_lfc3_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtF5, ],
#          color = color_palette,
#          main = "Up-regulated Genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
# )
# dev.off()





###ХИТМАПЫ ДЛЯ DOWN ГЕНОВ

### ХИТМАПЫ ДЛЯ DOWN ГЕНОВ
## AtF3_vs_AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj < 0.05 & log2FoldChange <= 1.5))
save(down_genes_lfc, file = "down_sig_genes_AtK3_AtF3.RData")

# # Извлечение нормализованных данных для down-regulated генов
# filtered_counts_down_lfc <- assay(rld)[down_genes_lfc, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_down_lfc <- t(scale(t(filtered_counts_down_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для down-regulated генов с LFC < -1.5
# png("heatmap_DOWN_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_down_lfc, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtK3, ],
#          color = color_palette,
#          main = "Down-regulated Genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
# )
# dev.off()




## Сравнение групп AtF5 и AtK5 !!! НЕТ НИ ОДНОГО ЗНАЧЕНИЯ !!! (из за фильтрации контрольных генов)
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc_5, file = "down_sig_genes_AtK5_AtF5.RData")

# # Извлечение нормализованных данных для down-regulated генов
# filtered_counts_down_lfc_5 <- assay(rld)[down_genes_lfc_5, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_down_lfc_5 <- t(scale(t(filtered_counts_down_lfc_5[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для down-regulated генов с LFC < -1.5
# png("heatmap_DOWN_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_down_lfc_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF5_AtK5, ],
#          color = color_palette,
#          main = "Down-regulated Genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
# )
# dev.off()




## Сравнение групп AtF3 и AtF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc3_5, file = "down_sig_genes_AtF3_AtF5.RData")

# # Извлечение нормализованных данных для down-regulated генов
# filtered_counts_down_lfc3_5 <- assay(rld)[down_genes_lfc3_5, ]
# 
# # Центрирование значений по строкам (генам)
# centered_filtered_counts_down_lfc3_5 <- t(scale(t(filtered_counts_down_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))
# 
# # Построение тепловой карты для down-regulated генов с LFC < -1.5
# png("heatmap_DOWN_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
# pheatmap(centered_filtered_counts_down_lfc3_5, 
#          cluster_rows = TRUE, 
#          show_rownames = FALSE, 
#          cluster_cols = FALSE, 
#          annotation_col = df[samples_AtF3_AtF5, ],
#          color = color_palette,
#          main = "Down-regulated Genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
# )
# dev.off()






### Список генов background для GO

# Извлечение таблицы каунтов
all_colnames <- colnames(dds)[dds$group]
counts_all <- counts(dds)[, all_colnames]

# гены - бэкграунд для GO (получается 9к элементов) - итог: вектор с названиями
exp_g_go_background <- rownames(counts_all[rowSums(counts_all != 0) > 1, ])

# Сохранение вектора в файл
save(exp_g_go_background, file = "go_background.RData")


### ОБЪЕДИНЕННЫЙ ХИТМАП
# объединить гены
# посчитать с повторностями и объединить их в единый table какой-то 
## значимые гены: sig_genes_lfc (F3/K3), sig_genes_lfc3_5 (F3/F5), sig_genes_lfc_5 (F5/K5)
# Объединяем списки значимых генов из трёх контрастов
all_significant_genes <- unique(c(sig_genes_lfc, sig_genes_lfc3_5, sig_genes_lfc_5))
# кол-во уникальных генов
length(all_significant_genes)

# Извлечение нормализованных данных для этих генов из rlog матрицы
subset_rlog <- assay(vsd)[all_significant_genes, ]
dim(subset_rlog)
# Центрирование значений по строкам (генам)
subset_rlog_centered <- t(scale(t(subset_rlog), center = TRUE, scale = FALSE))

# Задаем нужный порядок образцов
sample_order <- c("AtF3_1", "AtF3_2", "AtF3_3", 
                  "AtK3_1", "AtK3_2", "AtK3_3", 
                  "AtF5_1", "AtF5_2", 
                  "AtK5_1", "AtK5_2", "AtK5_3")

# Переупорядочиваем столбцы в нормализованной матрице по этому порядку
subset_rlog_centered <- subset_rlog_centered[, sample_order]


# Определяем минимальное и максимальное значения
min_value <- -max(abs(subset_rlog_centered))
max_value <- max(abs(subset_rlog_centered))

# Создаем breaks от min_value до max_value
breaks <- seq(min_value, max_value, length.out = 100)

# Определение цветовой палитры с белым цветом для значения 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Построение тепловой карты с заданными breaks
png("heatmap_signif_all.png", width = 14, height = 8, units = "in", res = 150, type = "cairo")
pheatmap(subset_rlog_centered, 
         cluster_rows = TRUE,   # Кластеризация по строкам (гены)
         cluster_cols = FALSE,  # Кластеризация по столбцам (образцы)
         main = "Heatmap of Differentially Expressed Genes in all Samples",
         show_rownames = FALSE, # Можно убрать названия генов, если их слишком много
         show_colnames = TRUE,  # Показываем названия образцов
         color = color_palette,  # Используем цветовую палитру
         breaks = breaks,        # Задаем заданные breaks
         distance = "manhattan",        # Используем Манхэттенское расстояние
         clustering_method = "centroid", # Метод кластеризации Ward.D2
         treeheight_row = 0             # Убираем дендрограмму строк
)
dev.off()



