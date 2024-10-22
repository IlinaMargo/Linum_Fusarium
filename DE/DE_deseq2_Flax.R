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


###ИЗВЛЕЧЕНИЕ ГЕНОВ ИЗ КОНТРОЛЬНЫХ ОБРАЗЦОВ, ЧТОБЫ ВЫЧЕСТЬ ИЗ ОСТАЛЬНЫХ ХИТМАПОВ
samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]

# Извлечение нормализованных данных для AtK3 и AtK5
samples_AtK3_AtK5 <- colnames(dds)[dds$group %in% c("AtK3", "AtK5")]
counts_AtK3_AtK5 <- counts(dds)[, samples_AtK3_AtK5]

# Фильтрация генов НЕНОРМАЛИЗОВАННЫХ ДАННЫХ, у которых есть экспрессия хотя бы в одном образце (ненулевое значение)
expressed_genes_AtK3_AtK5 <- rownames(counts_AtK3_AtK5[rowSums(counts_AtK3_AtK5 != 0) > 0, ])

# Например, исключаем эти гены из объекта DESeq2
dds <- dds[!rownames(dds) %in% expressed_genes_AtK3_AtK5, ]
rld <- rld[!rownames(rld) %in% expressed_genes_AtK3_AtK5, ]



### ХИТМАПЫ ДЛЯ НЕФИЛЬТРОВАННОГО
## Выбор образцов для AtF3 и AtK3
samples_AtF3_AtK3 <- colnames(dds)[dds$group %in% c("AtF3", "AtK3")]

# Центрирование значений по строкам (генам)
centered_all_counts <- t(scale(t(all_counts[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))

 # Плавный переход от синего к красному через белый
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))

df <- as.data.frame(colData(dds)[,c("variety", "day")])

# Построение тепловой карты для всех генов
png("heatmap_NF_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtK3, ],
         color = color_palette,
         main = "Non-filtered genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
)
dev.off()



## Выбор образцов для AtF5 и AtK5
samples_AtF5_AtK5 <- colnames(dds)[dds$group %in% c("AtF5", "AtK5")]

# Центрирование значений по строкам (генам)
centered_all_counts_5 <- t(scale(t(all_counts[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для всех генов
png("heatmap_NF_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF5_AtK5, ],
         color = color_palette,
         main = "Non-filtered genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
)
dev.off()

## Выбор образцов для AtF3 и AtF5
samples_AtF3_AtK5 <- colnames(dds)[dds$group %in% c("AtF3", "AtK5")]

# Центрирование значений по строкам (генам)
centered_all_counts_5_3 <- t(scale(t(all_counts[, samples_AtF3_AtK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для всех генов
png("heatmap_NF_AtF3_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts_5_3, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtK5, ],
         color = color_palette,
         main = "Non-filtered genes (AtF3 vs AtK5)"  # Заголовок тепловой карты
)
dev.off()


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
res_AtF3_AtK3 <- results(dds, contrast = c("group", "AtF3", "AtK3"), alpha = 0.5)

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj < 0.05 & abs(log2FoldChange) > 1.5))

# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc <- assay(rld)[sig_genes_lfc, ]

# Выбор образцов для AtF3 и AtK3
samples_AtF3_AtK3 <- colnames(dds)[dds$group %in% c("AtF3", "AtK3")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc <- t(scale(t(filtered_counts_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))

# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("heatmap_F_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtK3, ],
         color = color_palette,
         main = "Differential expressed genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF5 и AtK5
res_AtF5_AtK5 <- results(dds, contrast = c("group", "AtF5", "AtK5"), alpha = 0.5)

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj < 0.05 & abs(log2FoldChange) > 1.5))

# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc_5 <- assay(rld)[sig_genes_lfc_5, ]

# Выбор образцов для AtF5 и AtK5
samples_AtF5_AtK5 <- colnames(dds)[dds$group %in% c("AtF5", "AtK5")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc_5 <- t(scale(t(filtered_counts_lfc[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("heatmap_F_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF5_AtK5, ],
         color = color_palette,
         main = "Differential expressed genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF3 и AtF5
res_AtF3_AtF5 <- results(dds, contrast = c("group", "AtF3", "AtF5"), alpha = 0.5)

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj < 0.05 & abs(log2FoldChange) > 1.5))

# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc3_5 <- assay(rld)[sig_genes_lfc3_5, ]

# Выбор образцов для AtF3 и AtK3
samples_AtF3_AtF5 <- colnames(dds)[dds$group %in% c("AtF3", "AtF5")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc3_5 <- t(scale(t(filtered_counts_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("heatmap_F_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc3_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtF5, ],
         color = color_palette,
         main = "Differential expressed genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
)
dev.off()





###ХИТМАПЫ ДЛЯ UP ГЕНОВ
## AtF3_vs_AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj < 0.05 & log2FoldChange > 1.5))

# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc <- assay(rld)[up_genes_lfc, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc <- t(scale(t(filtered_counts_up_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtK3, ],
         color = color_palette,
         main = "Up-regulated Genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj < 0.05 & log2FoldChange > 1.5))

# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc_5 <- assay(rld)[up_genes_lfc_5, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc_5 <- t(scale(t(filtered_counts_up_lfc_5[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF5_AtK5, ],
         color = color_palette,
         main = "Up-regulated Genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
)
dev.off()


## Сравнение групп AtF3 и AtF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj < 0.05 & log2FoldChange > 1.5))

# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc3_5 <- assay(rld)[up_genes_lfc3_5, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc3_5 <- t(scale(t(filtered_counts_up_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc3_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtF5, ],
         color = color_palette,
         main = "Up-regulated Genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
)
dev.off()





###ХИТМАПЫ ДЛЯ DOWN ГЕНОВ

### ХИТМАПЫ ДЛЯ DOWN ГЕНОВ
## AtF3_vs_AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc <- rownames(subset(res_AtF3_AtK3, padj < 0.05 & log2FoldChange < -1.5))

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc <- assay(rld)[down_genes_lfc, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc <- t(scale(t(filtered_counts_down_lfc[, samples_AtF3_AtK3]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_AtF3_vs_AtK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtK3, ],
         color = color_palette,
         main = "Down-regulated Genes (AtF3 vs AtK3)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF5 и AtK5 !!! НЕТ НИ ОДНОГО ЗНАЧЕНИЯ !!! (из за фильтрации контрольных генов)
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc_5 <- rownames(subset(res_AtF5_AtK5, padj < 0.05 & log2FoldChange < -1.5))

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc_5 <- assay(rld)[down_genes_lfc_5, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc_5 <- t(scale(t(filtered_counts_down_lfc_5[, samples_AtF5_AtK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_AtF5_vs_AtK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF5_AtK5, ],
         color = color_palette,
         main = "Down-regulated Genes (AtF5 vs AtK5)"  # Заголовок тепловой карты
)
dev.off()
length(down_genes_lfc_5)



## Сравнение групп AtF3 и AtF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc3_5 <- rownames(subset(res_AtF3_AtF5, padj < 0.05 & log2FoldChange < -1.5))

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc3_5 <- assay(rld)[down_genes_lfc3_5, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc3_5 <- t(scale(t(filtered_counts_down_lfc3_5[, samples_AtF3_AtF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_AtF3_vs_AtF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc3_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_AtF3_AtF5, ],
         color = color_palette,
         main = "Down-regulated Genes (AtF3 vs AtF5)"  # Заголовок тепловой карты
)
dev.off()









# # Для гистограммы
# 
# hist(norm_counts, breaks = 50, main = "Histogram of normalized counts")
# 
# # Для плотности
# plot(density(norm_counts), main = "Density plot of normalized counts")
# 
# # Построение Q-Q графика для нормализованных данных
# qqnorm(norm_counts, main = "Q-Q Plot of normalized counts")
# qqline(norm_counts, col = "red")
# 
# ##проверка на биномиальное распределение 
# 
# library(MASS)
# # Вычисление среднего и дисперсии для каждого гена
# mean_values <- rowMeans(norm_counts)
# variance_values <- apply(norm_counts, 1, var)
# 
# # Построение графика: среднее против дисперсии
# plot(mean_values, variance_values, log="xy", main="Mean vs Variance for all genes",
#      xlab="Mean (log scale)", ylab="Variance (log scale)", pch=16, col="blue")
# abline(0, 1, col="red")  # Для Пуассоновского распределения среднее = дисперсия (тут скорее всего отрицательное биномиальное распределение)
# 
# ## сохранение результатов
# install.packages("openxlsx")
# library(openxlsx)
# #нефильтрованная дифэкспрессия
# # Преобразование результатов в data.frame
# res_df <- as.data.frame(res)
# 
# # Сохранение в формате .xlsx
# write.xlsx(res_df, "DE_Flax_unfiltered.xlsx", rowNames = TRUE)
# 
# #фильтрованная дифэкспрессия
# # Преобразование фильтрованных данных в data.frame
# significant_genes_df <- as.data.frame(significant_genes)
# 
# # Сохранение в формате .xlsx
# write.xlsx(significant_genes_df, "differential_expression_filtered.xlsx", rowNames = TRUE)
# 
# 
# 
# 
