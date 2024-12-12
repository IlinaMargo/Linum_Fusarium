BiocManager::install("topGO")
library(topGO)
library(ggplot2)

# загрузка генов бэкграунда
load("go_background.RData")

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

################### для всех контрастов

# Подготовка всех контрастов
# Список значимых генов для каждого контраста
contrast_genes <- list(
  "AtK3_vs_AtF3" = sig_genes_lfc,
  "AtK5_vs_AtF5" = sig_genes_lfc_5,
  "AtF3_vs_AtF5" = sig_genes_lfc3_5,
  "Up_AtK3_vs_AtF3" = up_genes_lfc,
  "Up_AtK5_vs_AtF5" = up_genes_lfc_5,
  "Up_AtF3_vs_AtF5" = up_genes_lfc3_5,
  "Down_AtK3_vs_AtF3" = down_genes_lfc,
  "Down_AtK5_vs_AtF5" = down_genes_lfc_5,
  "Down_AtF3_vs_AtF5" = down_genes_lfc3_5
)

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


# Создаем список для хранения результатов
results_list <- list()

# Цикл по всем контрастам для atalante
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

# Сохраняем результаты в CSV
write.csv(final_results, "GO_enrichment_results_all_contrasts.csv", row.names = FALSE)

# Визуализация для всех контрастов

final_results$Significance <- -log10(as.numeric(final_results$classicFisher))
final_results$GO.ID <- factor(final_results$GO.ID, levels = rev(unique(final_results$GO.ID)))

# Преобразование данных для общего dotplot
final_results_long <- final_results
final_results_long$Significance <- -log10(as.numeric(final_results_long$classicFisher))

# Сортировка GO.ID для удобства визуализации
final_results_long$GO.ID <- factor(final_results_long$GO.ID, levels = rev(unique(final_results_long$GO.ID)))

# Добавляем вычисление BgRatio
total_genes_in_background <- length(exp_g_go_background)

final_results_long <- final_results_long %>%
  mutate(BgRatio = Annotated / total_genes_in_background)

# Определяем порядок для оси X
contrast_order <- c(
  'AtK3_vs_AtF3',
  'up_AtK3_vs_AtF3',
  'down_AtK3_vs_AtF3',
  'AtF3_vs_AtF5',
  'up_AtF3_vs_AtF5',
  'down_AtF3_vs_AtF5',
  'AtK5_vs_AtF5',
  'up_AtK5_vs_AtF5',
  'down_AtK5_vs_AtF5',
  'AtK3_vs_AtK5',
  'up_AtK3_vs_AtK5',
  'down_AtK3_vs_AtK5'
)

# Преобразуем переменную Contrast в фактор с заданным порядком уровней
final_results_long_LM$Contrast <- factor(final_results_long_LM$Contrast, levels = contrast_order)

# Построение dotplot
png("topGO_dotplot.png", width = 16, height = 12, units = "in", res = 150, type = "cairo")
ggplot(final_results_long, aes(x = Contrast, y = Term)) +
  geom_point(aes(size = BgRatio, color = Significance)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Contrast", y = "GO Terms", size = "BgRatio", color = "-log10(p-value)", title = "GO Enrichment Analysis Across Contrasts Atalante") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




## загрузка данных для LM
load("go_background_LM.RData")
load("sig_genes_LMF3_LMK3.RData")
load("sig_genes_LMF5_LMK5.RData")
load("sig_genes_LMF3_LMF5.RData")
load("sig_genes_LMK3_LMK5.RData")
load("up_genes_LMF3_LMK3.RData")
load("up_genes_LMF5_LMK5.RData")
load("up_genes_LMF3_LMF5.RData")
load("up_genes_LMK3_LMK5.RData")
load("down_sig_genes_LMF3_LMK3.RData")
load("down_sig_genes_LMK5_LMF5.RData")
load("down_sig_genes_LMF3_LMF5.RData")
load("down_sig_genes_LMK3_LMF5.RData")

# Список наборов генов
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
  'down_LMK3_LMK5' = down_genes_lfc_cont_LM
)

#### Для LM
results_list_LM <- list()

# Цикл по всем контрастам для LM
for (contrast_name_LM in names(contrast_genes_LM)) {
  message("Processing contrast: ", contrast_name_LM)
  
  # Текущий список генов
  sig_genes_LM <- contrast_genes_LM[[contrast_name_LM]]
  
  # Создаем geneList для данного контраста
  geneList_LM <- ifelse(exp_g_go_background_LM %in% sig_genes_LM, 1, 0)
  names(geneList_LM) <- exp_g_go_background_LM
  
  # Создаем объект topGOdata
  GOdata_LM <- new("topGOdata",
                description = paste("GO Analysis for", contrast_name_LM),
                ontology = "BP",  # Можно заменить на "MF" или "CC"
                allGenes = geneList_LM,
                geneSel = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = go_map,
                nodeSize = 10)
  
  # Выполняем тесты
  result_fisher_LM <- runTest(GOdata_LM, algorithm = "classic", statistic = "fisher")
  result_ks_LM <- runTest(GOdata_LM, algorithm = "elim", statistic = "ks")
  
  # Сбор результатов
  all_results_LM <- GenTable(GOdata_LM,
                          classicFisher = result_fisher_LM,
                          elimKS = result_ks_LM,
                          orderBy = "classicFisher",
                          topNodes = 20)
  
  # # Добавляем поправку на множественные сравнения
  all_results$adj_pval <- p.adjust(all_results$classicFisher, method = "BH")
  
  # Добавляем название контраста в результаты
  all_results_LM$Contrast <- contrast_name_LM
  
  # Сохраняем результаты в список
  results_list_LM[[contrast_name_LM]] <- all_results_LM
}

# Объединяем все результаты в один датафрейм
final_results_LM <- do.call(rbind, results_list_LM)

# Сохраняем результаты в CSV
write.csv(final_results_LM, "GO_enrichment_results_all_contrasts_LM.csv", row.names = FALSE)

# Визуализация для всех контрастов

final_results_LM$Significance <- -log10(as.numeric(final_results_LM$classicFisher))
final_results_LM$GO.ID <- factor(final_results_LM$GO.ID, levels = rev(unique(final_results_LM$GO.ID)))

# Преобразование данных для общего dotplot
final_results_long_LM <- final_results_LM
final_results_long_LM$Significance <- -log10(as.numeric(final_results_long_LM$classicFisher))

# Сортировка GO.ID для удобства визуализации
final_results_long_LM$GO.ID <- factor(final_results_long_LM$GO.ID, levels = rev(unique(final_results_long_LM$GO.ID)))

# Добавляем вычисление BgRatio
total_genes_in_background_LM <- length(exp_g_go_background_LM)

final_results_long_LM <- final_results_long_LM %>%
  mutate(BgRatio = Annotated / total_genes_in_background_LM)

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
  'up_LMK3_LMK5',
  'down_LMK3_LMK5'
)

# Преобразуем переменную Contrast в фактор с заданным порядком уровней
final_results_long_LM$Contrast <- factor(final_results_long_LM$Contrast, levels = contrast_order)

filtered_results <- final_results_long %>%
  filter(as.numeric(classicFisher) <= 0.05)
filtered_results_uo <- final_results_long %>%
  filter(as.numeric(adj_pval) <= 0.05)

# Преобразуем classicFisher в числовой
filtered_results$classicFisher <- as.numeric(filtered_results$classicFisher)

# Построение dotplot с заданным порядком
png("topGO_dotplot_LM.png", width = 18, height = 12, units = "in", res = 150, type = "cairo")
ggplot(final_results_long_LM, aes(x = Contrast, y = Term)) +
  geom_point(aes(size = BgRatio, color = Significance)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Contrast", 
    y = "GO Terms", 
    size = "BgRatio", 
    color = "-log10(p-value)", 
    title = "GO Enrichment Analysis Across Contrasts LM98"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
