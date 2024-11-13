# install.packages("dnet",repos="http://cran.r-project.org",type="source")
# BiocManager::install("hfang-bristol/XGR", dependencies=T)
# install.packages("/home/milina/R/x86_64-pc-linux-gnu-library/4.3/XGR_1.1.1.tar.gz", repos = NULL, type = "source")
library(XGR)
library(ggplot2)

## база данных reactome
file_path <- "/media/eternus1/projects/milina/R/RScripts/reactome_joint.xlsx"

# Загрузка данных из Excel-файла
reactome_data <- read_excel(file_path)
head(reactome_data)

load("go_background.RData")

## набор генов для анализа
load("sig_genes_AtK3_AtF3.RData")
load("sig_genes_AtK5_AtF5.RData")
load("sig_genes_AtF3_AtF5.RData")
load("up_sig_genes_AtK3_AtF3.RData")
load("up_sig_genes_AtK5_AtF5.RData")
load("up_sig_genes_AtF3_AtF5.RData")
load("down_sig_genes_AtK3_AtF3.RData")
load("down_sig_genes_AtK5_AtF5.RData")
load("down_sig_genes_AtF3_AtF5.RData")

# Список наборов генов
gene_sets <- list(
  AtK3_AtF3 = sig_genes_lfc,
  AtK5_AtF5 = sig_genes_lfc_5,
  AtF3_AtF5 = sig_genes_lfc3_5,
  up_AtK3_AtF3 = up_genes_lfc,
  up_AtK5_AtF5 = up_genes_lfc_5,
  up_AtF3_AtF5 = up_genes_lfc3_5,
  down_AtK3_AtF3 = down_genes_lfc,
  down_AtK5_AtF5 = down_genes_lfc_5,
  down_AtF3_AtF5 = down_genes_lfc3_5
)

# Проведение анализа обогащения для каждого набора генов
enrichment_results <- lapply(gene_sets, function(genes) {
  xEnricherYours(
    data.file = genes,
    annotation.file = reactome_data,
    background.file = exp_g_go_background,
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


## датафрейм для всех контрастов с построением столбца BgRatio
# Функция для вычисления BgRatio для одного объекта eTerm
calculate_BgRatio <- function(eTerm_obj) {
  # Общее количество генов в бэкграунде
  total_genes_in_background <- length(eTerm_obj$background)
  
  # Проверяем, что есть перекрывающиеся гены для расчета
  if (length(eTerm_obj$overlap) > 0) {
    # Вычисляем BgRatio для каждого термина в eTerm_obj
    bg_ratios <- sapply(eTerm_obj$overlap, function(overlap_genes) {
      # Количество генов, аннотированных для этого термина в бэкграунде
      num_genes_in_term_background <- length(overlap_genes)
      # Вычисляем BgRatio
      num_genes_in_term_background / total_genes_in_background
    })
    return(bg_ratios)
  } else {
    # Если данных нет, возвращаем NULL
    return(NULL)
  }
}

# Вектор с названиями контрастов
contrast_names <- c(
  "AtK3_AtF3",
  "AtK5_AtF5",
  "AtF3_AtF5",
  "up_AtK3_AtF3",
  "up_AtK5_AtF5",
  "up_AtF3_AtF5",
  "down_AtK3_AtF3",
  "down_AtK5_AtF5",
  "down_AtF3_AtF5"
)

# Применение функции ко всем объектам в enrichment_results и сбор результатов в единый датафрейм
all_results <- do.call(rbind, lapply(seq_along(enrichment_results), function(i) {
  eTerm_obj <- enrichment_results[[i]]
  bg_ratios <- calculate_BgRatio(eTerm_obj)
  
  # Проверка, что bg_ratios содержит данные
  if (!is.null(bg_ratios) && length(bg_ratios) > 0) {
    # Создание датафрейма для текущего контраста
    dotplot_data <- data.frame(
      TermID = eTerm_obj$term_info[, "id"],
      TermName = eTerm_obj$term_info[, "name"],
      adjp = eTerm_obj$adjp,
      BgRatio = bg_ratios,
      GeneSet = contrast_names[i],  # Используем название контраста
      stringsAsFactors = FALSE
    )
    return(dotplot_data)
  } else {
    # Если данных нет, возвращаем NULL, чтобы пропустить этот контраст
    return(NULL)
  }
}))

# Просмотр итогового датафрейма
head(all_results)


### это старый код без построение столбца BgRatio
# ## строим датафрейм со всеми контрастами для дальнейшей визуализации
# # Создаем список, чтобы хранить данные для всех контрастов
# all_results <- do.call(rbind, lapply(seq_along(enrichment_results), function(i) {
#   eTerm_obj <- enrichment_results[[i]]
#   
#   # Проверяем, есть ли данные в eTerm_obj перед извлечением
#   if (!is.null(eTerm_obj$term_info) && nrow(eTerm_obj$term_info) > 0) {
#     # Извлекаем данные для текущего контраста
#     dotplot_data <- data.frame(
#       TermID = eTerm_obj$term_info[, "id"],
#       TermName = eTerm_obj$term_info[, "name"],
#       Namespace = eTerm_obj$term_info[, "namespace"],
#       Distance = eTerm_obj$term_info[, "distance"],
#       adjp = eTerm_obj$adjp,
#       pvalue = eTerm_obj$pvalue,
#       odds_ratio = eTerm_obj$or,
#       fold_change = eTerm_obj$fc,
#       stringsAsFactors = FALSE
#     )
#     
#     # Добавляем имя контраста, чтобы различать их в итоговом датафрейме
#     dotplot_data$GeneSet <- paste0("Contrast_", i)
#     
#     return(dotplot_data)
#   } else {
#     # Если данных нет, возвращаем NULL, чтобы пропустить этот контраст
#     return(NULL)
#   }
# }))
# 
# 
# # Вектор с новыми именами для каждого контраста
# contrast_names <- c(
#   "AtK3_AtF3",
#   "AtK5_AtF5",
#   "AtF3_AtF5",
#   "up_AtK3_AtF3",
#   "up_AtK5_AtF5",
#   "up_AtF3_AtF5",
#   "down_AtK3_AtF3",
#   "down_AtK5_AtF5",
#   "down_AtF3_AtF5"
# )
# 
# # Обновляем значения в колонке GeneSet, используя новые имена контрастов
# all_results$GeneSet <- contrast_names[as.integer(gsub("Contrast_", "", all_results$GeneSet))]


## Рисуем dotplot
png("reactome_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(all_results, aes(x = GeneSet, y = TermName, size = BgRatio, color = adjp)) +
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
