
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


# Фильтрация результирующего датафрейма по padj < 0.01
filtered_results <- all_results %>%
  filter(adjp < 0.01)

# Определяем порядок для оси X
contrast_order <- c(
  'AtK3_AtF3',
  'up_AtK3_AtF3',
  'down_AtK3_AtF3',
  'up_AtF3_AtF5',
  'down_AtF3_AtF5',
  'AtK5_AtF5',
  'up_AtK5_AtF5',
  'down_AtK5_AtF5'
)

# Преобразуем переменную Contrast в фактор с заданным порядком уровней
filtered_results$GeneSet <- factor(filtered_results$GeneSet, levels = contrast_order)


## Рисуем dotplot
png("reactome_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(filtered_results, aes(x = GeneSet, y = TermName, size = BgRatio, color = adjp)) +
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
