### сравнения меджду сортами 
# это директория для анализа льна
# dir <- "/media/eternus1/projects/milina/kallisto_LU_DE/"
# директория для анализа гриба
dir <- "/media/eternus1/projects/milina/R/kallisto_FO_3/"
#samples <- c("AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtK3_1", "AtK3_2", "AtK3_3", "AtK5_1", "AtK5_2", "AtK5_3", "LMF3_2", "LMF3_3", "LMF5_1", "LMF5_2", "LMK3_1", "LMK3_2", "LMK3_3", "LMK5_1", "LMK5_2", "LMK5_3")
samples <- c("LMF3_1", "LMF3_2", "LMF5_1", "LMF5_2", "LMF5_3", "LMK3_1", "LMK3_2", "LMK3_3", "LMK5_1", "LMK5_2", "LMK5_3", "AtF3_1", "AtF3_2", "AtF3_3", "AtF5_1", "AtF5_2", "AtF5_3", "AtK3_1", "AtK3_2", "AtK3_3", "AtK5_1", "AtK5_2")

files <- file.path(dir, paste0(samples, ".k31"), "abundance.h5")

# Добавляем имена к каждому файлу для идентификации
names(files) <- samples

# Загружаем данные с помощью tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance = "lengthScaledTPM", varReduce = FALSE) 
# Проверяем загруженные данные
txi.kallisto$counts

# Создание таблицы с информацией об условиях эксперимента лен
# sampleTable <- data.frame(
#   sampleName = names(files),
#   day = factor(c("3", "3", "3", "5", "5", "3", "3", "3", "5", "5", "5", "3", "3", "5", "5", "3", "3", "3", "5", "5", "5")),
#   variety = factor(c("AtF", "AtF", "AtF", "AtF", "AtF", "AtK", "AtK", "AtK", "AtK", "AtK", "AtK", "LMF", "LMF", "LMF", "LMF", "LMK", "LMK", "LMK", "LMK", "LMK", "LMK"))
# )
# гриб
sampleTable <- data.frame(
  sampleName = names(files),
  day = factor(c("3", "3", "5", "5", "5", "3", "3", "3", "5", "5", "5", "3", "3", "3", "5", "5", "5", "3", "3", "3", "5", "5")),
  variety = factor(c("LMF", "LMF", "LMF", "LMF", "LMF", "LMK", "LMK", "LMK", "LMK", "LMK", "LMK", "AtF", "AtF", "AtF", "AtF", "AtF", "AtF", "AtK", "AtK", "AtK", "AtK", "AtK"))
)


# Создание объекта DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = sampleTable, design = ~ 1 + day + variety + day:variety)
# Создание нового фактора для комбинации variety и day
dds$group <- factor(paste0(dds$variety, dds$day))
# Обновление дизайна эксперимента
design(dds) <- ~ group
# Проведение анализа
dds <- DESeq(dds)
# Нормализация данных с использованием Regularized Log Transformation
rld <- rlog(dds, blind = FALSE)

# Извлечение нормализованных данных для всех генов
all_counts <- assay(rld)

### ХИТМАПЫ ФИЛЬТРОВАННЫХ ЗНАЧЕНИЙ 
## Сравнение групп AtF3 и LMK3
res_AtF3_LMF3 <- results(dds, contrast = c("group", "AtF3", "LMF3"))
# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 2
sig_genes_LM_At_3 <- rownames(subset(res_AtF3_LMF3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))


## Сравнение групп AtF5 и LMK5
res_AtF5_LMF5 <- results(dds, contrast = c("group", "AtF5", "LMF5"))
# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_LM_At_5 <- rownames(subset(res_AtF5_LMF5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))


## Сравнение групп LMK3 и AtK3
res_LMK3_AtK3 <- results(dds, contrast = c("group", "LMK3", "AtK3"))
# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_LM_At_3K <- rownames(subset(res_LMK3_AtK3, padj <= 0.05 & abs(log2FoldChange) >= 1.5))


## Сравнение групп AtK5 и LMK5
res_AtK5_LMK5 <- results(dds, contrast = c("group", "AtK5", "LMK5"))
# Фильтрация значимых генов с порогом padj < 0.05 и |log2FoldChange| > 1.5
sig_genes_LM_At_5K <- rownames(subset(res_AtK5_LMK5, padj <= 0.05 & abs(log2FoldChange) >= 1.5))


### значимые гены для каждого контраста 
## AtF3_LMF3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_LM_At_3 <- rownames(subset(res_AtF3_LMF3, padj <= 0.05 & log2FoldChange >= 1.5))

## Сравнение групп AtF5 и LMF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_LM_At_5 <- rownames(subset(res_AtF5_LMF5, padj <= 0.05 & log2FoldChange >= 1.5))

## Сравнение групп LMK3 и AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_LM_At_3K <- rownames(subset(res_LMK3_AtK3, padj <= 0.05 & log2FoldChange >= 1.5))

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
up_genes_LM_At_5K <- rownames(subset(res_AtK5_LMK5, padj <= 0.05 & log2FoldChange >= 1.5))



##down-regulated genes
## AtF3_LMF3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_LM_At_3 <- rownames(subset(res_AtF3_LMF3, padj <= 0.05 & log2FoldChange <= - 1.5))

## Сравнение групп AtF5 и LMF5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_LM_At_5 <- rownames(subset(res_AtF5_LMF5, padj <= 0.05 & log2FoldChange <= - 1.5))

## Сравнение групп LMK3 и AtK3
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_LM_At_3K <- rownames(subset(res_LMK3_AtK3, padj <= 0.05 & log2FoldChange <= - 1.5))

## Сравнение групп AtF5 и AtK5
# Фильтрация значимых генов с порогом padj < 0.05 и LFC > 1.5
down_genes_LM_At_5K <- rownames(subset(res_AtK5_LMK5, padj <= 0.05 & log2FoldChange <= - 1.5))


### Список генов background для GO

# Извлечение таблицы каунтов
all_colnames <- colnames(dds)[dds$group]
counts_all <- counts(dds)[, all_colnames]

# гены - бэкграунд для GO (получается 9к элементов) - итог: вектор с названиями
exp_g_go_background <- rownames(counts_all[rowSums(counts_all != 0) > 1, ])


### ОБЪЕДИНЕННЫЙ ХИТМАП
# объединить гены
# посчитать с повторностями и объединить их в единый table какой-то 
# Объединяем списки значимых генов из трёх контрастов
all_significant_genes <- unique(c(sig_genes_LM_At_3, sig_genes_LM_At_5, sig_genes_LM_At_3K, sig_genes_LM_At_5K))
# кол-во уникальных генов
length(all_significant_genes)

# Извлечение нормализованных данных для этих генов из rlog матрицы
subset_rlog <- assay(rld)[all_significant_genes, ]
dim(subset_rlog)
# Центрирование значений по строкам (генам)
subset_rlog_centered <- t(scale(t(subset_rlog), center = TRUE, scale = FALSE))

# Задаем нужный порядок образцов
sample_order <- c("AtF3_1", "AtF3_2", "AtF3_3", 
                  "LMF3_2", "LMF3_3",
                  "AtK3_1", "AtK3_2", "AtK3_3",
                  "LMK3_1", "LMK3_2", "LMK3_3",
                  "AtF5_1", "AtF5_2",
                  "LMF5_1", "LMF5_2",
                  "AtK5_1", "AtK5_2", "AtK5_3",
                  "LMK5_1", "LMK5_2", "LMK5_3"
                  )
# Переупорядочиваем столбцы в нормализованной матрице по этому порядку
subset_rlog_centered <- subset_rlog_centered[, sample_order]


# Определяем минимальное и максимальное значения
min_value <- -max(abs(subset_rlog_centered))
max_value <- max(abs(subset_rlog_centered))

# Создаем breaks от min_value до max_value
breaks<- seq(min_value, max_value, length.out = 100)

# Определение цветовой палитры с белым цветом для значения 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Построение тепловой карты с заданными breaks
png("At_LM_heatmap_signif_all.png", width = 14, height = 8, units = "in", res = 150, type = "cairo")
pheatmap(subset_rlog_centered, 
         cluster_rows = TRUE,   # Кластеризация по строкам (гены)
         cluster_cols = FALSE,  # Кластеризация по столбцам (образцы)
         main = "Heatmap of Differentially Expressed Genes in all Samples FO",
         show_rownames = FALSE, # Можно убрать названия генов, если их слишком много
         show_colnames = TRUE,  # Показываем названия образцов
         color = color_palette,  # Используем цветовую палитру
         breaks = breaks,        # Задаем заданные breaks
         distance = "manhattan",        # Используем Манхэттенское расстояние
         clustering_method = "ward.D2", # Метод кластеризации Ward.D2
         treeheight_row = 0             # Убираем дендрограмму строк
)
dev.off()




### Построение диаграмм Венна
# список для диаграммы Вена
gene_lists_LM_At <- list(
  "LMF3_vs_AtF3" = sig_genes_LM_At_3,
  "LMF5_vs_AtF5" = sig_genes_LM_At_5,
  "LMK3_vs_AtK3" = sig_genes_LM_At_3K,
  "LMF5_vs_AtK5" = sig_genes_LM_At_5K
)

png("At_LM_Venn_all.png", width = 10, height = 6, units = "in", res = 150, type="cairo")
ggVennDiagram(gene_lists_LM_At) +
  scale_fill_gradient(low = "grey90", high = "purple4") +
  ggtitle("All different expressed genes LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  coord_cartesian(xlim = c(-0.2, 1.2), ylim = c(-0.1, 0.9))
dev.off()




### для up генов
# список для диаграммы Вена
gene_lists_LM_At_up <- list(
  "up_LMF3_vs_AtF3" = up_genes_LM_At_3,
  "up_LMF5_vs_AtF5" = up_genes_LM_At_5,
  "up_LMK3_vs_AtK3" = up_genes_LM_At_3K,
  "up_LMF5_vs_AtK5" = up_genes_LM_At_5K
)

png("At_LM_Venn_diagram_up.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
# Построение диаграммы Венна
ggVennDiagram(gene_lists_LM_At_up) +
  scale_fill_gradient(low = "grey90", high = "red") +
  ggtitle("Up-regulated genes between At and LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 0.9))
dev.off()


### для down генов
# список для диаграммы Вена
gene_lists_LM_At_down <- list(
  "down_LMF3_vs_AtF3" = down_genes_LM_At_3,
  "down_LMF5_vs_AtF5" = down_genes_LM_At_5,
  "down_LMK3_vs_AtK3" = down_genes_LM_At_3K,
  "down_LMF5_vs_AtK5" = down_genes_LM_At_5K
)

png("At_LM_Venn_down.png", width = 10, height = 8, units = "in", res = 150, type="cairo")
ggVennDiagram(gene_lists_LM_At_down) +
  scale_fill_gradient(low = "grey90", high = "blue") +
  ggtitle("Down-regulated genes At LM") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 0.9))
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
go_map <- split(go_annotations$GO_ID, go_annotations$GeneID)

# Фильтруем go_annotations, оставляя только строки, где GeneID есть в exp_g_go_background
filt_go_annotations <- go_annotations %>%
  filter(GeneID %in% exp_g_go_background)


# Выполнение анализа обогащения GO с использованием базы данных для льна

# Создаем список генов для каждого контраста
contrast_genes <- list(
  "LMF3_vs_AtF3" = sig_genes_LM_At_3,
  "LMF5_vs_AtF5" = sig_genes_LM_At_5,
  "LMK3_vs_AtK3" = sig_genes_LM_At_3K,
  "LMF5_vs_AtK5" = sig_genes_LM_At_5K,
  "up_LMF3_vs_AtF3" = up_genes_LM_At_3,
  "up_LMF5_vs_AtF5" = up_genes_LM_At_5,
  "up_LMK3_vs_AtK3" = up_genes_LM_At_3K,
  "up_LMF5_vs_AtK5" = up_genes_LM_At_5K,
  "down_LMF3_vs_AtF3" = down_genes_LM_At_3,
  "down_LMF5_vs_AtF5" = down_genes_LM_At_5,
  "down_LMK3_vs_AtK3" = down_genes_LM_At_3K,
  "down_LMF5_vs_AtK5" = down_genes_LM_At_5K
)

results_list <- list()

# Цикл по всем контрастам для LM
for (contrast_name in names(contrast_genes)) {
  message("Processing contrast: ", contrast_name)
  
  # Текущий список генов
  sig_genes <- contrast_genes[[contrast_name]]
  
  # Создаем geneList для данного контраста
  geneList <- ifelse(exp_g_go_background %in% sig_genes, 1, 0)
  names(geneList) <- exp_g_go_background
  
  # Создаем объект topGOdata
  GOdata <- new("topGOdata",
                   description = paste("GO Analysis for", contrast_name),
                   ontology = "BP",  # Можно заменить на "MF" или "CC"
                   allGenes = geneList,
                   geneSel = function(x) x == 1,
                   annot = annFUN.gene2GO,
                   gene2GO = go_map,
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
  
  # Добавляем название контраста в результаты
  all_results$Contrast <- contrast_name
  
  # Сохраняем результаты в список
  results_list[[contrast_name]] <- all_results
}

# Объединяем все результаты в один датафрейм
final_results <- do.call(rbind, results_list)

final_results$Significance <- -log10(as.numeric(final_results$classicFisher))
final_results$GO.ID <- factor(final_results$GO.ID, levels = rev(unique(final_results$GO.ID)))

# Преобразование данных для общего dotplot
final_results_long <- final_results
final_results_long$Significance <- -log10(as.numeric(final_results_long$classicFisher))

# Сортировка GO.ID для удобства визуализации
final_results_long$GO.ID <- factor(final_results_long$GO.ID, levels = rev(unique(final_results_long$GO.ID)))

# Вычисление BgRatio
total_genes_in_background <- length(exp_g_go_background)

final_results_long <- final_results_long %>%
  mutate(BgRatio = Annotated / total_genes_in_background)

# Определяем порядок для оси X
contrast_order <- c(
  'LMF3_vs_AtF3',
  'up_LMF3_vs_AtF3',
  'down_LMF3_vs_AtF3',
  'LMF5_vs_AtF5',
  'up_LMF5_vs_AtF5',
  'down_LMF5_vs_AtF5',
  'LMK3_vs_AtK3',
  'up_LMK3_vs_AtK3',
  'down_LMK3_vs_AtK3',
  'LMF5_vs_AtK5',
  'up_LMF5_vs_AtK5',
  'down_LMF5_vs_AtK5'
)
# Преобразуем переменную Contrast в фактор с заданным порядком уровней
final_results_long$Contrast <- factor(final_results_long$Contrast, levels = contrast_order)

# Построение dotplot
png("At_LM_topGO_dotplot.png", width = 16, height = 12, units = "in", res = 150, type = "cairo")
ggplot(final_results_long, aes(x = Contrast, y = Term)) +
  geom_point(aes(size = BgRatio, color = Significance)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Contrast", y = "GO Terms", size = "BgRatio", color = "-log10(p-value)", title = "GO Enrichment Analysis Across Contrasts Atalante and LM98") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
