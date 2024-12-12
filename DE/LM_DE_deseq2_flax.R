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
dir <- "/media/eternus1/projects/milina/kallisto_LU_DE/"

### КОД ДЛЯ ВСЕХ СЕМПЛОВ ###
## убрали из рассмотрения LMF5_3 и LMF3_1 из-за картинки экспрессии и схожести образцов
# Создаем вектор с именами образцов (если у вас есть файл samples, его можно использовать)
samples_LM <- c("LMF3_2", "LMF3_3", "LMF5_1", "LMF5_2", "LMK3_1", "LMK3_2", "LMK3_3", "LMK5_1", "LMK5_2", "LMK5_3")

files_LM <- file.path(dir, paste0(samples_LM, ".k31"), "abundance.h5")

# Добавляем имена к каждому файлу для идентификации
names(files_LM) <- samples_LM

# Загружаем данные с помощью tximport
txi.kallisto_LM <- tximport(files_LM, type = "kallisto", txOut = TRUE, countsFromAbundance = "lengthScaledTPM", varReduce = FALSE) 
# Проверяем загруженные данные
txi.kallisto_LM$counts

# Создание таблицы с информацией об условиях эксперимента
sampleTable_LM <- data.frame(
  sampleName = names(files_LM),
  #  group = factor(c("infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control", "control", "infected", "infected", "infected", "infected", "infected", "infected", "control", "control", "control", "control", "control")),
  day = factor(c("3", "3", "5", "5", "3", "3", "3", "5", "5", "5")),
  variety = factor(c("LMF", "LMF", "LMF", "LMF", "LMK", "LMK", "LMK", "LMK", "LMK", "LMK"))
  #  variety = factor(c("At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "At", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM", "LM"))
)

# Создание объекта DESeq2
dds_LM <- DESeqDataSetFromTximport(txi.kallisto_LM, colData = sampleTable_LM, design = ~ 1 + day + variety + day:variety)
# Создание нового фактора для комбинации variety и day
dds_LM$group <- factor(paste0(dds_LM$variety, dds_LM$day))
# Обновление дизайна эксперимента
design(dds_LM) <- ~ group
# Проведение анализа
dds_LM <- DESeq(dds_LM)

# Нормализация данных с использованием Variance Stabilizing Transformation 
vsd_LM <- vst(dds_LM, blind=FALSE)
# Нормализация данных с использованием Regularized Log Transformation
rld_LM <- rlog(dds_LM, blind = FALSE)

# Извлечение нормализованных данных для всех генов
all_counts_LM <- assay(rld_LM)

## посмотрим на все образцы
hmcol <- colorRampPalette(brewer.pal(9, "RdPu"))(100)
rlogMat_LM <- assay(rld_LM)
distsRL_LM <- dist(t(assay(rld_LM)))
mat_LM <- as.matrix(distsRL_LM)
hc_LM <- hclust(distsRL_LM)

png("heatmap_exp_through_all_samp_LM.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
heatmap.2(
  mat_LM,
  Rowv = as.dendrogram(hc_LM),
  symm = TRUE,
  trace = "none",
  col = rev(hmcol),
  margin = c(13, 13)
)
dev.off()



### ХИТМАПЫ ДЛЯ НЕФИЛЬТРОВАННОГО
## Выбор образцов для LMF3 и LMK3
samples_LMF3_LMK3 <- colnames(dds_LM)[dds_LM$group %in% c("LMF3", "LMK3")]

# Центрирование значений по строкам (генам)
centered_all_counts_LM <- t(scale(t(all_counts_LM[, samples_LMF3_LMK3]), center=TRUE, scale=FALSE))

# Плавный переход от синего к красному через белый
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))

df_LM <- as.data.frame(colData(dds_LM)[,c("variety", "day")])

# Построение тепловой карты для всех генов
png("LM_heatmap_NF_LMF3_vs_LMK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMK3, ],
         color = color_palette,
         main = "Non-filtered genes LM (LMF3 vs LMK3)"  # Заголовок тепловой карты
)
dev.off()




## Выбор образцов для LMF5 и LMK5
samples_LMF5_LMK5 <- colnames(dds_LM)[dds_LM$group %in% c("LMF5", "LMK5")]

# Центрирование значений по строкам (генам)
centered_all_counts_5_LM <- t(scale(t(all_counts_LM[, samples_LMF5_LMK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для всех генов
png("LM_heatmap_NF_LMF5_vs_LMK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF5_LMK5, ],
         color = color_palette,
         main = "Non-filtered genes (LMF5 vs LMK5)"  # Заголовок тепловой карты
)
dev.off()

## Выбор образцов для LMF3 и LMF5
samples_LMF3_LMF5 <- colnames(dds_LM)[dds_LM$group %in% c("LMF3", "LMF5")]

# Центрирование значений по строкам (генам)
centered_all_counts_5_3_LM <- t(scale(t(all_counts_LM[, samples_LMF3_LMF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для всех генов
png("LM_heatmap_NF_LMF3_vs_LMF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_all_counts_5_3_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMF5, ],
         color = color_palette,
         main = "Non-filtered genes LM (LMF3 vs LMF5)"  # Заголовок тепловой карты
)
dev.off()




### ХИТМАПЫ ФИЛЬТРОВАННЫХ ЗНАЧЕНИЙ 
## Сравнение групп LMF3 и LMK3
res_LMF3_LMK3 <- results(dds_LM, contrast = c("group", "LMF3", "LMK3"))

# Преобразование объекта результатов в датафрейм
df_LMF3_LMK3 <- as.data.frame(res_LMF3_LMK3)


png("hist_val_LMF3_LMK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_LMF3_LMK3, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) LMF3 vs LMK3",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_LMF3_LMK3$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_LMF3_LMK3$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()


# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 2
sig_genes_lfc_LM <- rownames(subset(res_LMF3_LMK3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc_LM, file = "sig_genes_LMF3_LMK3.RData")
# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc_LM <- assay(rld_LM)[sig_genes_lfc_LM, ]

# Выбор образцов для AtF3 и AtK3
samples_LMF3_LMK3 <- colnames(dds_LM)[dds_LM$group %in% c("LMF3", "LMK3")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc_LM <- t(scale(t(filtered_counts_lfc_LM[, samples_LMF3_LMK3]), center=TRUE, scale=FALSE))


# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("LM_heatmap_F_LMF3_vs_LMK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMK3, ],
         color = color_palette,
         main = "Differential expressed genes LM (LMF3 vs LMK3)"  # Заголовок тепловой карты
)
dev.off()





## Сравнение групп LMF5 и LMK5
res_LMF5_LMK5 <- results(dds_LM, contrast = c("group", "LMF5", "LMK5"))

# Преобразование объекта результатов в датафрейм
df_LMF5_LMK5 <- as.data.frame(res_LMF5_LMK5)

png("hist_val_LMF5_LMK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_LMF5_LMK5, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) LMF5 vs LMK5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_LMF5_LMK5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_LMF5_LMK5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()


# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc_5_LM <- rownames(subset(res_LMF5_LMK5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc_5_LM, file = "sig_genes_LMF5_LMK5.RData")
# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc_5_LM <- assay(rld_LM)[sig_genes_lfc_5_LM, ]

# Выбор образцов для AtF5 и AtK5
samples_LMF5_LMK5_LM <- colnames(dds_LM)[dds_LM$group %in% c("LMF5", "LMK5")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc_5_LM <- t(scale(t(filtered_counts_lfc_LM[, samples_LMF5_LMK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("heatmap_F_LMF5_vs_LMK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF5_LMK5, ],
         color = color_palette,
         main = "Differential expressed genes (LMF5 vs LMK5)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп LMF3 и LMF5
res_LMF3_LMF5 <- results(dds_LM, contrast = c("group", "LMF3", "LMF5"))

# Преобразование объекта результатов в датафрейм
df_LMF3_LMF5 <- as.data.frame(res_LMF3_LMF5)

png("hist_val_LMF3_LMF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_LMF3_LMF5, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) LMF3 vs LMF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_LMF3_LMF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_LMF3_LMF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()


# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc3_5_LM <- rownames(subset(res_LMF3_LMF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc3_5_LM, file = "sig_genes_LMF3_LMF5.RData")
# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc3_5_LM <- assay(rld_LM)[sig_genes_lfc3_5_LM, ]

# Выбор образцов для AtF3 и AtK3
samples_LMF3_LMF5 <- colnames(dds_LM)[dds_LM$group %in% c("LMF3", "LMF5")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc3_5_LM <- t(scale(t(filtered_counts_lfc3_5_LM[, samples_LMF3_LMF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для значимых генов с |LFC| > 1.5
png("heatmap_F_LMF3_vs_LMF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_lfc3_5, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df[samples_LMF3_LMF5, ],
         color = color_palette,
         main = "Differential expressed genes (LMF3 vs LMF5)"  # Заголовок тепловой карты
)
dev.off()


## Сравнение групп LMF3 и LMF5
res_LMK3_LMK5 <- results(dds_LM, contrast = c("group", "LMK3", "LMK5"))

# Преобразование объекта результатов в датафрейм
df_LMK3_LMK5 <- as.data.frame(res_LMK3_LMK5)

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_lfc_cont_LM <- rownames(subset(res_LMK3_LMK5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
save(sig_genes_lfc_cont_LM, file = "sig_genes_LMK3_LMK5.RData")
# Извлечение нормализованных данных для значимых генов
filtered_counts_lfc_cont_LM <- assay(rld_LM)[sig_genes_lfc_cont_LM, ]

# Выбор образцов для AtF3 и AtK3
samples_LMK3_LMK5 <- colnames(dds_LM)[dds_LM$group %in% c("LMK3", "LMK5")]

# Центрирование значений по строкам (генам)
centered_filtered_counts_lfc_cont_LM <- t(scale(t(filtered_counts_lfc_cont_LM[, samples_LMK3_LMK5]), center=TRUE, scale=FALSE))





###ХИТМАПЫ ДЛЯ UP ГЕНОВ
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc_LM <- rownames(subset(res_LMF3_LMK3, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc_LM, file = "up_genes_LMF3_LMK3.RData")
# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc_LM <- assay(rld_LM)[up_genes_lfc_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc_LM <- t(scale(t(filtered_counts_up_lfc_LM[, samples_LMF3_LMK3]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_LMF3_vs_LMK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMK3, ],
         color = color_palette,
         main = "Up-regulated Genes (LMF3 vs LMK3)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc_5_LM <- rownames(subset(res_LMF5_LMK5, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc_5_LM, file = "up_genes_LMF5_LMK5.RData")
# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc_5_LM <- assay(rld_LM)[up_genes_lfc_5_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc_5_LM <- t(scale(t(filtered_counts_up_lfc_5_LM[, samples_LMF5_LMK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_LMF5_vs_LMK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF5_LMK5, ],
         color = color_palette,
         main = "Up-regulated Genes LM (LMF5 vs LMK5)"  # Заголовок тепловой карты
)
dev.off()


## Сравнение групп LMF3 и LMF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc3_5_LM <- rownames(subset(res_LMF3_LMF5, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc3_5_LM, file = "up_genes_LMF3_LMF5.RData")
# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc3_5_LM <- assay(rld_LM)[up_genes_lfc3_5_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc3_5_LM <- t(scale(t(filtered_counts_up_lfc3_5_LM[, samples_LMF3_LMF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для up-regulated генов с LFC > 1.5
png("heatmap_UP_LMF3_vs_LMF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_up_lfc3_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMF5, ],
         color = color_palette,
         main = "Up-regulated Genes LM (LMF3 vs LMF5)"  # Заголовок тепловой карты
)
dev.off()

## Сравнение групп LMK3 и LMK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_lfc_cont_LM <- rownames(subset(res_LMK3_LMK5, padj <= 0.05 & log2FoldChange >= 1.5))
save(up_genes_lfc_cont_LM, file = "up_genes_LMK3_LMK5.RData")
# Извлечение нормализованных данных для up-regulated генов
filtered_counts_up_lfc_cont_LM <- assay(rld_LM)[up_genes_lfc_cont_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_up_lfc_cont_LM <- t(scale(t(filtered_counts_up_lfc_cont_LM[, samples_LMK3_LMK5]), center=TRUE, scale=FALSE))




###ХИТМАПЫ ДЛЯ DOWN ГЕНОВ
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc_LM <- rownames(subset(res_LMF3_LMK3, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc_LM, file = "down_sig_genes_LMF3_LMK3.RData")
# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc_LM <- assay(rld_LM)[down_genes_lfc_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc_LM <- t(scale(t(filtered_counts_down_lfc_LM[, samples_LMF3_LMK3]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_LMF3_vs_LMK3.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMK3, ],
         color = color_palette,
         main = "Down-regulated Genes (LMF3 vs LMK3)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп AtF5 и AtK5 !!! НЕТ НИ ОДНОГО ЗНАЧЕНИЯ !!! (из за фильтрации контрольных генов)
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc_5_LM <- rownames(subset(res_LMF5_LMK5, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc_5_LM, file = "down_sig_genes_LMK5_LMF5.RData")

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc_5_LM <- assay(rld_LM)[down_genes_lfc_5_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc_5_LM <- t(scale(t(filtered_counts_down_lfc_5_LM[, samples_LMF5_LMK5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_LMF5_vs_LMK5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF5_LMK5, ],
         color = color_palette,
         main = "Down-regulated Genes (LMF5 vs LMK5)"  # Заголовок тепловой карты
)
dev.off()




## Сравнение групп LMF3 и LMF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc3_5_LM <- rownames(subset(res_LMF3_LMF5, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc3_5_LM, file = "down_sig_genes_LMF3_LMF5.RData")

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc3_5_LM <- assay(rld_LM)[down_genes_lfc3_5_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc3_5_LM <- t(scale(t(filtered_counts_down_lfc3_5_LM[, samples_LMF3_LMF5]), center=TRUE, scale=FALSE))

# Построение тепловой карты для down-regulated генов с LFC < -1.5
png("heatmap_DOWN_LMF3_vs_LMF5.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
pheatmap(centered_filtered_counts_down_lfc3_5_LM, 
         cluster_rows = TRUE, 
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         annotation_col = df_LM[samples_LMF3_LMF5, ],
         color = color_palette,
         main = "Down-regulated Genes LM (LMF3 vs LMF5)"  # Заголовок тепловой карты
)
dev.off()


# Фильтрация значимых генов с порогом padj < 0.05 и LFC < -1.5
down_genes_lfc_cont_LM <- rownames(subset(res_LMK3_LMK5, padj <= 0.05 & log2FoldChange <= - 1.5))
save(down_genes_lfc_cont_LM, file = "down_sig_genes_LMK3_LMF5.RData")

# Извлечение нормализованных данных для down-regulated генов
filtered_counts_down_lfc_cont_LM <- assay(rld_LM)[down_genes_lfc_cont_LM, ]

# Центрирование значений по строкам (генам)
centered_filtered_counts_down_lfc_cont_LM <- t(scale(t(filtered_counts_down_lfc_cont_LM[, samples_LMK3_LMK5]), center=TRUE, scale=FALSE))




### Список генов background для GO

# Извлечение таблицы каунтов
all_colnames_LM <- colnames(dds_LM)[dds_LM$group]
counts_all_LM <- counts(dds_LM)[, all_colnames_LM]

# гены - бэкграунд для GO (получается 9к элементов) - итог: вектор с названиями
exp_g_go_background_LM <- rownames(counts_all_LM[rowSums(counts_all_LM != 0) > 1, ])

# Сохранение вектора в файл
save(exp_g_go_background_LM, file = "go_background_LM.RData")


### ОБЪЕДИНЕННЫЙ ХИТМАП
# объединить гены
# посчитать с повторностями и объединить их в единый table какой-то 
## значимые гены: sig_genes_lfc (F3/K3), sig_genes_lfc3_5 (F3/F5), sig_genes_lfc_5 (F5/K5)
# Объединяем списки значимых генов из трёх контрастов
all_significant_genes_LM <- unique(c(sig_genes_lfc_LM, sig_genes_lfc3_5_LM, sig_genes_lfc_5_LM, sig_genes_lfc_cont_LM))
# кол-во уникальных генов
length(all_significant_genes_LM)

# Извлечение нормализованных данных для этих генов из rlog матрицы
subset_rlog_LM <- assay(vsd_LM)[all_significant_genes_LM, ]
dim(subset_rlog_LM)
# Центрирование значений по строкам (генам)
subset_rlog_centered_LM <- t(scale(t(subset_rlog_LM), center = TRUE, scale = FALSE))

# Задаем нужный порядок образцов
sample_order_LM <- c("LMF3_2", "LMF3_3", 
                  "LMK3_1", "LMK3_2", "LMK3_3", 
                  "LMF5_1", "LMF5_2",
                  "LMK5_1", "LMK5_2", "LMK5_3")

# Переупорядочиваем столбцы в нормализованной матрице по этому порядку
subset_rlog_centered_LM <- subset_rlog_centered_LM[, sample_order_LM]


# Определяем минимальное и максимальное значения
min_value_LM <- -max(abs(subset_rlog_centered_LM))
max_value_LM <- max(abs(subset_rlog_centered_LM))

# Создаем breaks от min_value до max_value
breaks_LM <- seq(min_value_LM, max_value_LM, length.out = 100)

# Определение цветовой палитры с белым цветом для значения 0
color_palette_LM <- colorRampPalette(c("blue", "white", "red"))(length(breaks_LM) - 1)

# Построение тепловой карты с заданными breaks
png("LM_heatmap_signif_all.png", width = 14, height = 8, units = "in", res = 150, type = "cairo")
pheatmap(subset_rlog_centered_LM, 
         cluster_rows = TRUE,   # Кластеризация по строкам (гены)
         cluster_cols = FALSE,  # Кластеризация по столбцам (образцы)
         main = "Heatmap of Differentially Expressed Genes in all Samples LM",
         show_rownames = FALSE, # Можно убрать названия генов, если их слишком много
         show_colnames = TRUE,  # Показываем названия образцов
         color = color_palette_LM,  # Используем цветовую палитру
         breaks = breaks_LM,        # Задаем заданные breaks
         distance = "manhattan",        # Используем Манхэттенское расстояние
         clustering_method = "ward.D2", # Метод кластеризации Ward.D2
         treeheight_row = 0             # Убираем дендрограмму строк
)
dev.off()




### Построение диаграмм Венна
# список для диаграммы Вена
gene_lists_LM <- list(
  "LMK3_vs_LMF3" = sig_genes_lfc_LM,
  "LMK5_vs_LMF5" = sig_genes_lfc_5_LM,
  "LMF3_vs_LMF5" = sig_genes_lfc3_5_LM,
  "LMK3_vs_LMK5" = sig_genes_lfc_cont_LM
)

png("LM_Venn_diagram_all.png", width = 10, height = 6, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_LM) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.margin = margin(20, 50, 20, 50)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.1, 0.9))
dev.off()




### для up генов
# список для диаграммы Вена
gene_lists_up_LM <- list(
  "LMK3_vs_LMF3" = up_genes_lfc_LM,
  "LMK5_vs_LMF5" = up_genes_lfc_5_LM,
  "LMF3_vs_LMF5" = up_genes_lfc3_5_LM,
  "LMK3_vs_LMK5" = up_genes_lfc_cont_LM
)

png("LM_Venn_diagram_up.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_up_LM) +
  scale_fill_gradient(low = "grey90", high = "red") +
  ggtitle("Up-regulated genes LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.margin = margin(20, 50, 20, 50)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.1, 0.9))
dev.off()


### для down генов
# список для диаграммы Вена
gene_lists_down_LM <- list(
  "LMK3_vs_LMF3" = down_genes_lfc_LM,
  "LMK5_vs_LMF5" = down_genes_lfc_5_LM,
  "LMF3_vs_LMF5" = down_genes_lfc3_5_LM,
  "LMK3_vs_LMK5" = down_genes_cont_LM
)

png("LM_Venn_diagram_down.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_down_LM) +
  scale_fill_gradient(low = "grey90", high = "blue") +
  ggtitle("Down-regulated genes LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.margin = margin(20, 50, 20, 50)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.1, 0.9))
dev.off()




### Анализ GO

# Загрузка базы данных GO
go_annotations <- read.delim("/media/eternus1/projects/milina/db/Lu_2_go", 
                             header = FALSE, 
                             sep = "\t",
                             col.names = c("GeneID", "PFAM", "GO_term_description", "GO_ID"))

# Разделяем и извлекаем значения
go_annotations <- go_annotations %>%
  separate(GO_term_description, into = c("GO_term_description", "GO_ID"), sep = " ; ", extra = "merge") %>%
  mutate(GO_ID = trimws(GO_ID))  # Удаляем пробелы из GO_ID

head(go_annotations)
#теперь значения на своем месте в столбцах

# Фильтруем go_annotations, оставляя только строки, где GeneID есть в exp_g_go_background
filt_go_annotations_LM <- go_annotations %>%
  filter(GeneID %in% exp_g_go_background_LM)


# Выполнение анализа обогащения GO с использованием базы данных для льна
# Проведение анализа обогащения GO
# enricher чтобы установить свой собственный сет background genes


# Подготовка TERM2GENE и TERM2NAME
TERM2GENE <- go_annotations %>%
  select(GO_ID, GeneID) %>%
  distinct()

TERM2NAME <- go_annotations %>%
  select(GO_ID, GO_term_description) %>%
  distinct()

# Создаем список генов для каждого контраста
contrast_genes_LM <- list(
  'LMK3_LMF3' = sig_genes_lfc_LM,
  'LMK5_LMF5' = sig_genes_lfc_5_LM,
  'LMF3_LMF5' = sig_genes_lfc3_5_LM,
  'LMK3_LMK5' = sig_genes_lfc_cont_LM,
  'up_LMK3_LMF3' = up_genes_lfc_LM,
  'up_LMK5_LMF5' = up_genes_lfc_5_LM,
  'up_LMF3_LMF5' = up_genes_lfc3_5_LM,
  'up_LMK3_LMK5' = up_genes_lfc_cont_LM,
  'down_LMK3_LMF3' = down_genes_lfc_LM,
  'down_LMK5_LMF5' = down_genes_lfc_5_LM,
  'down_LMF3_LMF5' = down_genes_lfc3_5_LM,
  'down_LMK3_LMK5' = down_genes_cont_LM
)

# Объединение результатов анализа обогащения для всех контрастов
compare_results_LM <- compareCluster(
  geneClusters = contrast_genes_LM,
  fun = "enricher",
  pAdjustMethod = "BH",
  universe = filt_go_annotations_LM$GeneID,
  minGSSize = 10,
  maxGSSize = 1000,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)


# Извлекаем результаты в виде data.frame для каждого контраста
compare_results_df_LM <- as.data.frame(compare_results_LM)
dim(compare_results_df_LM)
# Фильтрация строк, где qvalue >= 0.01
compare_results_df_LM <- compare_results_df_LM %>%
  dplyr::filter(qvalue <= 0.01)

head(compare_results_df_LM)
# Построение dot plot с помощью ggplot2
library(ggplot2)

# Подготовка данных: выделение чисел из GeneRatio и BgRatio для точного вычисления
compare_results_df_LM <- compare_results_df_LM %>%
  mutate(
    BgRatio_numeric = as.numeric(sapply(strsplit(as.character(BgRatio), "/"), `[`, 1)) /
      as.numeric(sapply(strsplit(as.character(BgRatio), "/"), `[`, 2))
  )


png("LM_gg_GO_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(compare_results_df_LM, aes(x = Cluster, y = Description, size = BgRatio_numeric, color = qvalue)) +
  geom_point() +
  scale_size_continuous(range = c(3, 10)) +  # Устанавливаем диапазон размера точек
  scale_color_gradient(low = "red", high = "blue") +  # Градиент цвета для qvalue
  ggtitle("GO Enrichment Analysis Across Contrasts LM") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms", size = "BgRatio", color = "qvalue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Поворот подписей на оси X
dev.off()


png("LM_GO_dotplot.png", width = 16, height = 12, units = "in", res = 150, type = "cairo")
# Построение dot plot для объединенного результата
dotplot(compare_results_LM, showCategory = 20, color = 'qvalue') +
  ggtitle("GO Enrichment Analysis Across Contrasts LM") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## если гнать со всеми генами, а не только дифэкспрессирующимися

# Объединение результатов анализа обогащения для всех контрастов
compare_results_all_LM <- compareCluster(
  geneClusters = contrast_genes_LM,
  fun = "enricher",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 1000,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

compare_results_df_all_LM <- as.data.frame(compare_results_all_LM)
dim(compare_results_df_all_LM)
# Фильтрация строк, где qvalue >= 0.01
compare_results_df_all_LM <- compare_results_df_all_LM %>%
  dplyr::filter(qvalue <= 0.01)

head(compare_results_df_all_LM)

# Подготовка данных: выделение чисел из GeneRatio и BgRatio для точного вычисления
compare_results_df_all_LM <- compare_results_df_all_LM %>%
  mutate(
    BgRatio_numeric_LM = as.numeric(sapply(strsplit(as.character(BgRatio), "/"), `[`, 1)) /
      as.numeric(sapply(strsplit(as.character(BgRatio), "/"), `[`, 2))
  )


png("LM_GO_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(compare_results_df_all_LM, aes(x = Cluster, y = Description, size = BgRatio_numeric_LM, color = qvalue)) +
  geom_point() +
  scale_size_continuous(range = c(3, 10)) +  # Устанавливаем диапазон размера точек
  scale_color_gradient(low = "red", high = "blue") +  # Градиент цвета для qvalue
  ggtitle("GO Enrichment Analysis Across Contrasts LM") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms", size = "BgRatio", color = "qvalue") +
  theme(axis.text.y = element_text(size = 8))
dev.off()

###только для всех контрастов без up and down
contrast_genes_onlyall_LM <- list(
  "LMK3_vs_LMF3" = sig_genes_lfc_LM,
  "LMK5_vs_LMF5" = sig_genes_lfc_5_LM,
  "LMF3_vs_LMF5" = sig_genes_lfc3_5_LM
)


# Объединение результатов анализа обогащения для всех контрастов
compare_results_onlyall_LM <- compareCluster(
  geneClusters = contrast_genes_onlyall_LM,
  fun = "enricher",
  pAdjustMethod = "bonferroni",
  universe = filt_go_annotations_LM$GeneID,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

png("LM_GO_dotplot_only_all.png", width = 12, height = 8, units = "in", res = 150, type = "cairo")
# Построение dot plot для объединенного результата
dotplot(compare_results_onlyall_LM, showCategory = 30, color = 'qvalue') +
  ggtitle("GO Enrichment Analysis Across Contrasts LM") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




### Реактом анализ
library(XGR)
library(ggplot2)

## база данных reactome
file_path <- "/media/eternus1/projects/milina/R/RScripts/reactome_joint.xlsx"

# Загрузка данных из Excel-файла
reactome_data <- read_excel(file_path)
head(reactome_data)


# Список наборов генов
gene_sets_LM <- list(
  'LMK3_LMF3' = sig_genes_lfc_LM,
  'LMK5_LMF5' = sig_genes_lfc_5_LM,
  'LMF3_LMF5' = sig_genes_lfc3_5_LM,
  'LMK3_LMK5' = sig_genes_lfc_cont_LM,
  'up_LMK3_LMF3' = up_genes_lfc_LM,
  'up_LMK5_LMF5' = up_genes_lfc_5_LM,
  'up_LMF3_LMF5' = up_genes_lfc3_5_LM,
  'up_LMK3_LMK5' = up_genes_lfc_cont_LM,
  'down_LMK3_LMF3' = down_genes_lfc_LM,
  'down_LMK5_LMF5' = down_genes_lfc_5_LM,
  'down_LMF3_LMF5' = down_genes_lfc3_5_LM,
  'down_LMK3_LMK5' = down_genes_cont_LM
)

# Проведение анализа обогащения для каждого набора генов
enrichment_results_LM <- lapply(gene_sets_LM, function(genes) {
  xEnricherYours(
    data.file = genes,
    annotation.file = reactome_data,
    background.file = exp_g_go_background_LM,
    size.range = c(10, 2000),
    min.overlap = 5,
    test = "fisher",
    background.annotatable.only = TRUE,
    p.tail = "one-tail",
    p.adjust.method = "BH",
    verbose = TRUE,
    silent = FALSE
  )
})


# Вектор с названиями контрастов
contrast_names_LM <- c(
  "LMK3_LMF3",
  "LMK5_LMF5",
  "LMF3_LMF5",
  'LMK3_LMK5',
  "up_LMK3_LMF3",
  "up_LMK5_LMF5",
  "up_LMF3_LMF5",
  'up_LMK3_LMK5',
  "down_LMK3_LMF3",
  "down_LMK5_LMF5",
  "down_LMF3_LMF5",
  'down_LMK3_LMK5'
)

# Применение функции ко всем объектам в enrichment_results и сбор результатов в единый датафрейм
all_results_LM <- do.call(rbind, lapply(seq_along(enrichment_results_LM), function(i) {
  eTerm_obj_LM <- enrichment_results_LM[[i]]
  bg_ratios_LM <- calculate_BgRatio(eTerm_obj_LM)
  
  # Проверка, что bg_ratios содержит данные
  if (!is.null(bg_ratios_LM) && length(bg_ratios_LM) > 0) {
    # Создание датафрейма для текущего контраста
    dotplot_data_LM <- data.frame(
      TermID = eTerm_obj_LM$term_info[, "id"],
      TermName = eTerm_obj_LM$term_info[, "name"],
      adjp = eTerm_obj_LM$adjp,
      BgRatio = bg_ratios_LM,
      GeneSet = contrast_names_LM[i],  # Используем название контраста
      stringsAsFactors = FALSE
    )
    return(dotplot_data_LM)
  } else {
    # Если данных нет, возвращаем NULL, чтобы пропустить этот контраст
    return(NULL)
  }
}))

# Фильтрация результирующего датафрейма по padj < 0.01
filtered_results_LM <- all_results_LM %>%
  filter(adjp < 0.01)


# Определяем порядок для оси X
contrast_order <- c(
  'LMK3_LMF3',
  'up_LMK3_LMF3',
  'down_LMK3_LMF3',
  'LMF3_LMF5',
  'up_LMF3_LMF5',
  'down_LMF3_LMF5',
  'LMK5_LMF5',
  'up_LMK5_LMF5',
  'down_LMK5_LMF5',
  'LMK3_LMK5',
  'down_LMK3_LMK5'
)

# Преобразуем переменную Contrast в фактор с заданным порядком уровней
filtered_results_LM$GeneSet <- factor(filtered_results_LM$GeneSet, levels = contrast_order)

## Рисуем dotplot
png("LM_reactome_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(filtered_results_LM, aes(x = GeneSet, y = TermName, size = BgRatio, color = adjp)) +
  geom_point() +
  scale_size_continuous(range = c(3, 10)) +  # Настройка диапазона размера точек
  scale_color_gradient(low = "red", high = "blue") +  # Градиент цвета для adjp
  ggtitle("Reactome Enrichment Analysis Across Contrasts") +
  theme_minimal() +
  labs(
    x = "Contrast",
    y = "Reactome Pathways",
    size = "BgRatio",
    color = "Adjusted p-value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Поворот подписей на оси X для удобства чтения
  )
dev.off()




### рисуем связь путей, чтобы показать какие то значимые взаимосвязи
library(piano)
data("gsa_input")
# Форматируем gsc
gsc <- as.matrix(go_annotations[, c("GeneID", "GO_term_description")])
# Преобразуем в объект gene set collection
geneSets <- loadGSC(gsc)

## AtF3_AtK3
# # Извлекаем p-values
# pvals <- res_LMF3_LMK3$pvalue
# # Присваиваем имена для векторов
# names(pvals) <- rownames(res_LMF3_LMK3)
# 
# ## # Определяем направление изменений
# directions <- ifelse(res_LMF3_LMK3$log2FoldChange > 0, 1, -1)
# # Присваиваем имена для векторов
# names(directions) <- rownames(res_LMF3_LMK3)

# Извлекаем p-values
pvals <- res_LMF5_LMK5$pvalue
# Присваиваем имена для векторов
names(pvals) <- rownames(res_LMF5_LMK5)

## # Определяем направление изменений
directions <- ifelse(res_LMF5_LMK5$log2FoldChange > 0, 1, -1)
# Присваиваем имена для векторов
names(directions) <- rownames(res_LMF5_LMK5)

# Собираем gsa_input
gsa_input <- list(
  pvals = pvals,              # p-values
  directions = directions,    # Направление изменений
  gsc = gsc                   # Коллекция наборов генов
)

# Выявляем индексы, где p-values не NA
valid_indices <- !is.na(gsa_input$pvals)
# Отбираем только корректные данные
gsa_input$pvals <- gsa_input$pvals[valid_indices]
gsa_input$directions <- gsa_input$directions[valid_indices]


gsares <- runGSA(
  geneLevelStats = gsa_input$pvals,       # Передаем p-values
  directions = gsa_input$directions,     # Направления изменений
  gsc = geneSets,                   # Коллекция наборов генов
  nPerm = 10000                          # Количество перестановок
)

networkPlot2(gsares, class = "non", significance = 0.1)

