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

### Анализ для гриба во льне устойчивого сорта LM98 
# Задаём директорию, где находятся файлы
dir_FO <- "/media/eternus1/projects/milina/R/kallisto_FO_3/"

### КОД ДЛЯ ВСЕХ СЕМПЛОВ ###
## убрали из рассмотрения LMF5_3 и LMF3_1 из-за картинки экспрессии и схожести образцов
# Создаем вектор с именами образцов (если у вас есть файл samples, его можно использовать)
samples_FO <- c("LMF3_2", "LMF3_3", "LMF5_1", "LMF5_2", "LMF5_3", "AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtF5_3", "Fo_1", "Fo_2", "Fo_3")

files_FO <- file.path(dir_FO, paste0(samples_FO, ".k31"), "abundance.h5")

# Добавляем имена к каждому файлу для идентификации
names(files_FO) <- samples_FO

# Загружаем данные с помощью tximport
txi.kallisto_FO <- tximport(files_FO, type = "kallisto", txOut = TRUE, countsFromAbundance = "lengthScaledTPM", varReduce = FALSE) 
# Проверяем загруженные данные
txi.kallisto_FO$counts

# Создание таблицы с информацией об условиях эксперимента
sampleTable <- data.frame(
  sampleName = names(files_FO),
  variety = factor(c("LMF3", "LMF3", "LMF5", "LMF5", "LMF5", "AtF3", "AtF3", "AtF3", "AtF5", "AtF5", "AtF5", "Fo", "Fo", "Fo" ))
  )

# Создание объекта DESeq2
dds_FO <- DESeqDataSetFromTximport(txi.kallisto_FO, colData = sampleTable, design = ~ 1 + variety)
# Проведение анализа
dds_FO <- DESeq(dds_FO)

# Нормализация данных с использованием Regularized Log Transformation
rld_FO <- rlog(dds_FO, blind = FALSE)

# Извлечение нормализованных данных для всех генов
all_counts_FO <- assay(rld_FO)

# посмотрим на все образцы
hmcol <- colorRampPalette(brewer.pal(9, "RdPu"))(100)
rlogMat_FO <- assay(rld_FO)
distsRL_FO <- dist(t(assay(rld_FO)))
mat_FO <- as.matrix(distsRL_FO)
hc_FO <- hclust(distsRL_FO)

png("heatmap_exp_through_all_FO.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
heatmap.2(
  mat_FO,
  Rowv = as.dendrogram(hc_FO),
  symm = TRUE,
  trace = "none",
  col = rev(hmcol),
  margin = c(13, 13)
)
dev.off()



### Создание контрастов Atalante ###
## Сравнение групп Fo AtF3
res_FO_AtF3 <- results(dds_FO, contrast = c("variety", "Fo", "AtF3"))

# Преобразование объекта результатов в датафрейм
df_FO_AtF3 <- as.data.frame(res_FO_AtF3)

png("df_FO_AtF3_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_AtF3, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) Fo vs AtF3",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_AtF3$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_AtF3$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 2
sig_genes_lfc_FO_At3 <- rownames(subset(res_FO_AtF3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
#save(sig_genes_lfc_FO_At3, file = "sig_genes_lfc_FO_At3.RData")




## Сравнение групп LMF5 и LMK5
res_FO_AtF5 <- results(dds_FO, contrast = c("variety", "Fo", "AtF5"))

# Преобразование объекта результатов в датафрейм
df_FO_AtF5 <- as.data.frame(res_FO_AtF5)
##смотрим гистаграмму распределения по lFC
png("df_FO_AtF5_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_AtF5 , aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) Fo vs AtF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_AtF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_AtF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_FO_At5 <- rownames(subset(res_FO_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
#save(sig_genes_FO_At5, file = "sig_genes_FO_LMF5_LMK5.RData")



## Сравнение групп LMF3 и LMF5
res_FO_AtF3_AtF5 <- results(dds_FO, contrast = c("variety", "AtF3", "AtF5"))
# Преобразование объекта результатов в датафрейм
df_FO_AtF3_AtF5 <- as.data.frame(res_FO_AtF3_AtF5)
##смотрим гистаграмму распределения по lFC
png("df_FO_AtF3_AtF5_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_AtF3_AtF5 , aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) AtF3 vs AtF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_AtF3_AtF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_AtF3_AtF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_At3_5_FO <- rownames(subset(res_FO_AtF3_AtF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
### здесь получается только 1 дифэкспрессированный ген даже только по p-value. так что выкидываем из рассмотрения




### UP-regulated гены
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_FO_AtF3 <- rownames(subset(res_FO_AtF3, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_FO_LM, file = "up_genes_LMF3_LMK3.RData")

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_FO_AtF5 <- rownames(subset(res_FO_AtF5, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")

# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_At3_5_FO <- rownames(subset(res_FO_AtF3_AtF5, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")


### DOWN-regulated гены
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_FO_AtF3 <- rownames(subset(res_FO_AtF3, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_FO_LM, file = "up_genes_LMF3_LMK3.RData")

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_FO_AtF5 <- rownames(subset(res_FO_AtF5, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")

# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_At3_5_FO <- rownames(subset(res_FO_AtF3_AtF5, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")




### КОНТРАСТЫ LM98
res_FO_LMF3 <- results(dds_FO, contrast = c("variety", "Fo", "LMF3"))

# Преобразование объекта результатов в датафрейм
df_FO_LMF3 <- as.data.frame(res_FO_LMF3)

png("df_FO_LMF3_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_LMF3, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) Fo vs LMF3",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_LMF3$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_LMF3$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 2
sig_genes_lfc_FO_LM3 <- rownames(subset(res_FO_LMF3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
#save(sig_genes_lfc_FO_At3, file = "sig_genes_lfc_FO_At3.RData")


## Сравнение групп LMF5 и LMK5
res_FO_LMF5 <- results(dds_FO, contrast = c("variety", "Fo", "LMF5"))

# Преобразование объекта результатов в датафрейм
df_FO_LMF5 <- as.data.frame(res_FO_LMF5)
##смотрим гистаграмму распределения по lFC
png("df_FO_LMF5_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_LMF5 , aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) Fo vs LMF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_LMF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_LMF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_FO_LM5 <- rownames(subset(res_FO_LMF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
#save(sig_genes_FO_At5, file = "sig_genes_FO_LMF5_LMK5.RData")



## Сравнение групп LMF3 и LMF5
res_FO_LMF3_LMF5 <- results(dds_FO, contrast = c("variety", "LMF3", "LMF5"))
# Преобразование объекта результатов в датафрейм
df_FO_LMF3_LMF5 <- as.data.frame(res_FO_LMF3_LMF5)
##смотрим гистаграмму распределения по lFC
png("df_FO_LMF3_LMF5_hist_lfc.png", width = 8, height = 6, units = "in", res = 150, type="cairo")
ggplot(df_FO_LMF3_LMF5 , aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Log2 Fold Changes (LFC) LMF3 vs LMF5",
       x = "Log2 Fold Change (LFC)",
       y = "Frequency") +
  scale_x_continuous(breaks = seq(floor(min(df_FO_LMF3_LMF5$log2FoldChange, na.rm = TRUE)),
                                  ceiling(max(df_FO_LMF3_LMF5$log2FoldChange, na.rm = TRUE)), by = 1)) +
  theme_minimal()
dev.off()

# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_LM3_5_FO <- rownames(subset(res_FO_LMF3_LMF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))
### здесь получается только 1 дифэкспрессированный ген даже только по p-value. так что выкидываем из рассмотрения




### UP-regulated гены
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_FO_LMF3 <- rownames(subset(res_FO_LMF3, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_FO_LM, file = "up_genes_LMF3_LMK3.RData")

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_FO_LMF5 <- rownames(subset(res_FO_LMF5, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")

# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_LM3_5_FO <- rownames(subset(res_FO_LMF3_LMF5, padj <= 0.05 & log2FoldChange >= 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")


### DOWN-regulated гены
## LMF3_vs_LMK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_FO_LMF3 <- rownames(subset(res_FO_LMF3, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_FO_LM, file = "up_genes_LMF3_LMK3.RData")

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_FO_LMF5 <- rownames(subset(res_FO_LMF5, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")

# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_LM3_5_FO <- rownames(subset(res_FO_LMF3_LMF5, padj <= 0.05 & log2FoldChange <= - 1.5))
#save(up_genes_lfc_5_FO_LM, file = "up_genes_LMF5_LMK5.RData")




### Список генов background для GO At
# Извлечение таблицы каунтов
all_colnames_FO_At <- c("AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtF5_3", "Fo_1", "Fo_2", "Fo_3")
counts_all_FO_At <- counts(dds_FO)[, all_colnames_FO_At]

# гены - бэкграунд для GO (получается 9к элементов) - итог: вектор с названиями
exp_g_go_background_FO_At <- rownames(counts_all_FO_At[rowSums(counts_all_FO_At != 0) > 1, ])

# Сохранение вектора в файл
#save(exp_g_go_background_FO_LM, file = "go_background_FO_LM.RData")

### Список генов background для GO LM
# Извлечение таблицы каунтов
all_colnames_FO_LM <- c("LMF3_2", "LMF3_3", "LMF5_1", "LMF5_2", "LMF5_3", "Fo_1", "Fo_2", "Fo_3")
counts_all_FO_LM <- counts(dds_FO)[, all_colnames_FO_LM]

# гены - бэкграунд для GO (получается 9к элементов) - итог: вектор с названиями
exp_g_go_background_FO_LM <- rownames(counts_all_FO_LM[rowSums(counts_all_FO_LM != 0) > 1, ])


### ОБЪЕДИНЕННЫЙ ХИТМАП
# объединить гены
# посчитать с повторностями и объединить их в единый table какой-то 
## значимые гены: sig_genes_lfc (F3/K3), sig_genes_lfc3_5 (F3/F5), sig_genes_lfc_5 (F5/K5)
# Объединяем списки значимых генов из трёх контрастов
all_significant_genes_FO <- unique(c(sig_genes_lfc_FO_At3,
                                     sig_genes_FO_At5,
                                     sig_genes_At3_5_FO,
                                     sig_genes_lfc_FO_LM3,
                                     sig_genes_FO_LM5,
                                     sig_genes_LM3_5_FO))
# кол-во уникальных генов
length(all_significant_genes_FO)

# Извлечение нормализованных данных для этих генов из rlog матрицы
subset_rlog_FO <- assay(rld_FO)[all_significant_genes_FO, ]
dim(subset_rlog_FO)
# Центрирование значений по строкам (генам)
subset_rlog_centered_FO <- t(scale(t(subset_rlog_FO), center = TRUE, scale = FALSE))

# Задаем нужный порядок образцов
sample_order_FO <- c( "LMF3_2", "LMF3_3",
                        "LMF5_1", "LMF5_2", "LMF5_3", 
                        "Fo_1", "Fo_2", "Fo_3",
                        "AtF3_1","AtF3_2", "AtF3_3", 
                        "AtF5_1", "AtF5_2", "AtF5_3" 
                        )


# Переупорядочиваем столбцы в нормализованной матрице по этому порядку
subset_rlog_centered_FO <- subset_rlog_centered_FO[, sample_order_FO]


# Определяем минимальное и максимальное значения
min_value_FO <- -max(abs(subset_rlog_centered_FO))
max_value_FO <- max(abs(subset_rlog_centered_FO))

# Создаем breaks от min_value до max_value
breaks_FO <- seq(min_value_FO, max_value_FO, length.out = 100)

# Определение цветовой палитры с белым цветом для значения 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks_FO) - 1)

# Построение тепловой карты с заданными breaks
png("FO_heatmap_signif_all.png", width = 14, height = 8, units = "in", res = 150, type = "cairo")
pheatmap(subset_rlog_centered_FO, 
         cluster_rows = TRUE,   # Кластеризация по строкам (гены)
         cluster_cols = FALSE,  # Кластеризация по столбцам (образцы)
         main = "Heatmap of Differentially Expressed Genes in all Samples FO",
         show_rownames = FALSE, # Можно убрать названия генов, если их слишком много
         show_colnames = TRUE,  # Показываем названия образцов
         color = color_palette,  # Используем цветовую палитру
         breaks = breaks_FO,        # Задаем заданные breaks
         distance = "manhattan",        # Используем Манхэттенское расстояние
         clustering_method = "ward.D2", # Метод кластеризации Ward.D2
         treeheight_row = 0             # Убираем дендрограмму строк
)
dev.off()


### Построение диаграмм Венна At
# список для диаграммы Вена
gene_lists_FO<- list(
  "FO_vs_AtF3" = sig_genes_lfc_FO_At3,
  "FO_vs_AtF5" = sig_genes_FO_At5,
  "AtF3_vs_AtF5" = sig_genes_At3_5_FO 
)

png("FO_Venn_all.png", width = 10, height = 6, units = "in", res = 150, type="cairo")
ggVennDiagram(gene_lists_FO) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes At") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  coord_cartesian(xlim = c(-8, 12), ylim = c(-8, 5))
dev.off()


### Построение диаграмм Венна LM
# список для диаграммы Вена
gene_lists_FO<- list(
  "FO_vs_LMF3" = sig_genes_lfc_FO_LM3,
  "FO_vs_LMF5" = sig_genes_FO_LM5,
  "LMF3_vs_LMF5" = sig_genes_LM3_5_FO 
)

png("FO_Venn_all.png", width = 10, height = 6, units = "in", res = 150, type="cairo")
ggVennDiagram(gene_lists_FO) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  coord_cartesian(xlim = c(-8, 12), ylim = c(-8, 5))
dev.off()





### Анализ GO
###topGO
# Загрузка базы данных GO
go_annotations_FO <- read.delim("/media/eternus1/projects/milina/R/F39.go", 
                                header = FALSE, 
                                sep = "\t",
                                col.names = c("GeneID", "GO_term_description", "GO_ID"))

# Разделяем и извлекаем значения
go_annotations_FO <- go_annotations_FO %>%
  separate(GO_term_description, into = c("GO_term_description", "GO_ID"), sep = "_;_", extra = "merge") %>%
  mutate(GO_ID = trimws(GO_ID))  # Удаляем пробелы из GO_ID

go_map_FO <- split(go_annotations_FO$GO_ID, go_annotations_FO$GeneID)
#теперь значения на своем месте в столбцах

# Фильтруем go_annotations, оставляя только строки, где GeneID есть в exp_g_go_background
filt_go_annotations_At <- go_annotations_FO %>%
  filter(GeneID %in% exp_g_go_background_FO_At)


# Выполнение анализа обогащения GO с использованием базы данных для льна
# Проведение анализа обогащения GO

# Создаем список генов для каждого контраста
contrast_genes_At <- list(
  'Fo_AtF3' = sig_genes_lfc_FO_At3,
  'up_Fo_AtF3' = up_genes_FO_AtF3,
  'down_Fo_AtF3' = down_genes_FO_AtF3,
  'AtF3_AtF5' = sig_genes_At3_5_FO,
  'up_AtF3_AtF5' = up_genes_At3_5_FO,
  'down_AtF3_AtF5' = down_genes_At3_5_FO,
  'Fo_AtF5' = sig_genes_FO_At5,
  'up_Fo_AtF5'= up_genes_FO_AtF5,
  'down_Fo_AtF5' = down_genes_FO_AtF5
)

#### Для LM
results_list <- list()

# Цикл по всем контрастам для LM
for (contrast_name in names(contrast_genes_At)) {
  message("Processing contrast: ", contrast_name)
  
  # Текущий список генов
  sig_genes <- contrast_genes_At[[contrast_name]]
  
  # Создаем geneList для данного контраста
  geneList <- ifelse(exp_g_go_background_FO_At %in% sig_genes, 1, 0)
  names(geneList) <- exp_g_go_background_FO_At
  
  # Создаем объект topGOdata
  GOdata <- new("topGOdata",
                   description = paste("GO Analysis for", contrast_name),
                   ontology = "BP",  # Можно заменить на "MF" или "CC"
                   allGenes = geneList,
                   geneSel = function(x) x == 1,
                   annot = annFUN.gene2GO,
                   gene2GO = go_map_FO,
                   nodeSize = 10)
  
  # Выполняем тесты
  result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  result_ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  # Сбор результатов
  all_results <- GenTable(GOdata,
                             classicFisher = result_fisher,
                             elimKS = result_ks,
                             orderBy = "classicFisher",
                             topNodes = 20)
  # # Добавляем поправку на множественные сравнения
  all_results$adj_pval <- p.adjust(all_results$classicFisher, method = "BH")

  # Добавляем название контраста в результаты
  all_results$Contrast <- contrast_name
  
  # Сохраняем результаты в список
  results_list[[contrast_name]] <- all_results
}


# Объединяем все результаты в один датафрейм
final_results <- do.call(rbind, results_list)

# Сохраняем результаты в CSV
# write.csv(final_results, "GO_enrichment_results_all_contrasts.csv", row.names = FALSE)

# Визуализация для всех контрастов

#final_results_LM$Significance <- -log10(as.numeric(final_results_LM$adj))
final_results$GO.ID <- factor(final_results$GO.ID, levels = rev(unique(final_results$GO.ID)))

# Преобразование данных для общего dotplot
final_results_long <- final_results
#final_results_long_LM$Significance <- -log10(as.numeric(final_results_long_LM$classicFisher))

# Сортировка GO.ID для удобства визуализации
final_results_long$GO.ID <- factor(final_results_long$GO.ID, levels = rev(unique(final_results_long$GO.ID)))

# Добавляем вычисление BgRatio
total_genes_in_background_FO <- length(exp_g_go_background_FO_At)

final_results_long <- final_results_long %>%
  mutate(BgRatio = Annotated / total_genes_in_background_FO)

filtered_results <- final_results_long %>%
  filter(as.numeric(classicFisher) <= 0.05)
# filtered_results <- final_results_long %>%
#   filter(as.numeric(adj_pval) <= 0.05)
# Преобразуем classicFisher в числовой
filtered_results$classicFisher <- as.numeric(filtered_results$classicFisher)

# Построение dotplot с заданным порядком
png("topGO_dotplot_At_FO.png", width = 18, height = 12, units = "in", res = 150, type = "cairo")
ggplot(filtered_results, aes(x = Contrast, y = Term)) +
  geom_point(aes(size = BgRatio, color = adj_pval)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Contrast", 
    y = "GO Terms", 
    size = "BgRatio", 
    color = "adj p-value", 
    title = "GO Enrichment Analysis Across Contrasts Atalante FO"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14))
dev.off()




### LM98 ###
# Фильтруем go_annotations, оставляя только строки, где GeneID есть в exp_g_go_background
filt_go_annotations_LM <- go_annotations_FO %>%
  filter(GeneID %in% exp_g_go_background_FO_LM)


# Выполнение анализа обогащения GO с использованием базы данных для льна
# Проведение анализа обогащения GO

# Создаем список генов для каждого контраста
contrast_genes_LM <- list(
  'Fo_LMF3' = sig_genes_lfc_FO_LM3,
  'up_Fo_LMF3' = up_genes_FO_LMF3,
  'down_Fo_LMF3' = down_genes_FO_LMF3,
  'LMF3_LMF5' = sig_genes_LM3_5_FO,
  'up_LMF3_LMF5' = up_genes_LM3_5_FO,
  'down_LMF3_LMF5' = down_genes_LM3_5_FO,
  'Fo_LMF5' = sig_genes_FO_LM5,
  'up_Fo_LMF5'= up_genes_FO_LMF5,
  'down_Fo_LMF5' = down_genes_FO_LMF5
)

#### Для LM
results_list <- list()

# Цикл по всем контрастам для LM
for (contrast_name in names(contrast_genes_LM)) {
  message("Processing contrast: ", contrast_name)
  
  # Текущий список генов
  sig_genes <- contrast_genes_LM[[contrast_name]]
  
  # Создаем geneList для данного контраста
  geneList <- ifelse(exp_g_go_background_FO_LM %in% sig_genes, 1, 0)
  names(geneList) <- exp_g_go_background_FO_LM
  
  # Создаем объект topGOdata
  GOdata <- new("topGOdata",
                description = paste("GO Analysis for", contrast_name),
                ontology = "BP",  # Можно заменить на "MF" или "CC"
                allGenes = geneList,
                geneSel = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = go_map_FO,
                nodeSize = 10)
  
  # Выполняем тесты
  result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  result_ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  # Сбор результатов
  all_results <- GenTable(GOdata,
                          classicFisher = result_fisher,
                          elimKS = result_ks,
                          orderBy = "classicFisher",
                          topNodes = 20)
  # # Добавляем поправку на множественные сравнения
  all_results$adj_pval <- p.adjust(all_results$classicFisher, method = "BH")
  
  # Добавляем название контраста в результаты
  all_results$Contrast <- contrast_name
  
  # Сохраняем результаты в список
  results_list[[contrast_name]] <- all_results
}


# Объединяем все результаты в один датафрейм
final_results <- do.call(rbind, results_list)

# Сохраняем результаты в CSV
# write.csv(final_results, "GO_enrichment_results_all_contrasts.csv", row.names = FALSE)

# Визуализация для всех контрастов

#final_results_LM$Significance <- -log10(as.numeric(final_results_LM$adj))
final_results$GO.ID <- factor(final_results$GO.ID, levels = rev(unique(final_results$GO.ID)))

# Преобразование данных для общего dotplot
final_results_long <- final_results
#final_results_long_LM$Significance <- -log10(as.numeric(final_results_long_LM$classicFisher))

# Сортировка GO.ID для удобства визуализации
final_results_long$GO.ID <- factor(final_results_long$GO.ID, levels = rev(unique(final_results_long$GO.ID)))

# Добавляем вычисление BgRatio
total_genes_in_background_FO <- length(exp_g_go_background_FO_LM)

final_results_long <- final_results_long %>%
  mutate(BgRatio = Annotated / total_genes_in_background_FO)

filtered_results <- final_results_long %>%
  filter(as.numeric(classicFisher) <= 0.05)
# Преобразуем classicFisher в числовой
filtered_results$classicFisher <- as.numeric(filtered_results$classicFisher)

# Построение dotplot с заданным порядком
png("topGO_dotplot_LM_FO.png", width = 18, height = 12, units = "in", res = 150, type = "cairo")
ggplot(filtered_results, aes(x = Contrast, y = Term)) +
  geom_point(aes(size = BgRatio, color = adj_pval)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Contrast", 
    y = "GO Terms", 
    size = "BgRatio", 
    color = "adj p-value", 
    title = "GO Enrichment Analysis Across Contrasts LM98 FO"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14))
dev.off()


### рисуем связь путей, чтобы показать какие то значимые взаимосвязи
library(piano)
# Форматируем gsc
gsc <- as.matrix(go_annotations_FO[, c("GeneID", "GO_term_description")])
# Преобразуем в объект gene set collection
geneSets <- loadGSC(gsc)

# # LMF3_Fo
# # Извлекаем p-values
# pvals <- res_FO_LMF3$pvalue
# # Присваиваем имена для векторов
# names(pvals) <- rownames(res_FO_LMF3)
# 
# ## # Определяем направление изменений
# directions <- ifelse(res_FO_LMF3$log2FoldChange > 0, 1, -1)
# # Присваиваем имена для векторов
# names(directions) <- rownames(res_FO_LMF3)

# # LMF5_Fo
# # Извлекаем p-values
# pvals <- res_FO_LMF5$pvalue
# # Присваиваем имена для векторов
# names(pvals) <- rownames(res_FO_LMF5)
# 
# ## # Определяем направление изменений
# directions <- ifelse(res_FO_LMF5$log2FoldChange > 0, 1, -1)
# # Присваиваем имена для векторов
# names(directions) <- rownames(res_FO_LMF5)

# AtF3_Fo
# Извлекаем p-values
pvals <- res_FO_AtF3$pvalue
# Присваиваем имена для векторов
names(pvals) <- rownames(res_FO_AtF3)

## # Определяем направление изменений
directions <- ifelse(res_FO_AtF3$log2FoldChange > 0, 1, -1)
# Присваиваем имена для векторов
names(directions) <- rownames(res_FO_AtF3)

# # AtF5_Fo
# # Извлекаем p-values
# pvals <- res_FO_AtF5$pvalue
# # Присваиваем имена для векторов
# names(pvals) <- rownames(res_FO_AtF5)
# 
# ## # Определяем направление изменений
# directions <- ifelse(res_FO_AtF5$log2FoldChange > 0, 1, -1)
# # Присваиваем имена для векторов
# names(directions) <- rownames(res_FO_AtF5)



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
  gsc = geneSets,                        # Коллекция наборов генов
  nPerm = 10000                          # Количество перестановок
)

# pdf("networkplot_LMF5_LMK5_FO.pdf", width = 18, height = 12)
# networkPlot2(gsares, class = "non", significance = 0.95, adjusted = FALSE)
# dev.off()

networkPlot2(gsares, class = "non", significance = 0.05, adjusted = FALSE)

