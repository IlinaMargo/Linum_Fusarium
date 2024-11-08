# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("YuLab-SMU/ggtree")
# BiocManager::install("clusterProfiler")
BiocManager::install("DOSE", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("ReactomePA", force = TRUE)

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(tibble)
library(DBI)
library(ReactomePA)

# Загрузка бэкграунда для GO
load("go_background.RData")
print(exp_g_go_background)

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
filt_go_annotations <- go_annotations %>%
  filter(GeneID %in% exp_g_go_background)

# go_list <- go_annotations %>%
#   select(GeneID, GO_ID) %>%
#   distinct() %>%
#   group_by(GeneID) %>%
#   summarize(GO_terms = list(GO_ID)) %>%
#   deframe()
# 
# connection <- dbConnect(RSQLite::SQLite(), "/media/eternus1/projects/milina/db/Lu_2_go")
# DBI::dbDisconnect(connection)

# Выполнение анализа обогащения GO с использованием базы данных для льна
# Проведение анализа обогащения GO
# enricher чтобы установить свой собственный сет background genes
# гены 
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

# Подготовка TERM2GENE и TERM2NAME
TERM2GENE <- go_annotations %>%
  select(GO_ID, GeneID) %>%
  distinct()

TERM2NAME <- go_annotations %>%
  select(GO_ID, GO_term_description) %>%
  distinct()

# Создаем список генов для каждого контраста
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


# Объединение результатов анализа обогащения для всех контрастов
compare_results <- compareCluster(
  geneClusters = contrast_genes,
  fun = "enricher",
  pAdjustMethod = "BH",
  universe = filt_go_annotations$GeneID,
  minGSSize = 10,
  maxGSSize = 1000,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)


# Извлекаем результаты в виде data.frame для каждого контраста
compare_results_df <- as.data.frame(compare_results)
# Фильтрация строк, где qvalue >= 0.01
compare_results_df <- compare_results_df %>%
  dplyr::filter(qvalue <= 0.01)


# Построение dot plot с помощью ggplot2
library(ggplot2)

png("gg_GO_dotplot.png", width = 16, height = 10, units = "in", res = 150, type = "cairo")
ggplot(compare_results_df, aes(x = Cluster, y = Description, size = -log10(qvalue), color = qvalue)) +
  geom_point() +
  scale_size_continuous(range = c(3, 10)) +  # Устанавливаем диапазон размера точек
  scale_color_gradient(low = "red", high = "blue") +  # Градиент цвета для qvalue
  ggtitle("GO Enrichment Analysis Across Contrasts") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms", size = "-log10(qvalue)", color = "qvalue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Поворот подписей на оси X
dev.off()


png("GO_dotplot.png", width = 16, height = 12, units = "in", res = 150, type = "cairo")
# Построение dot plot для объединенного результата
dotplot(compare_results, showCategory = 20, color = 'qvalue') +
  ggtitle("GO Enrichment Analysis Across Contrasts") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## если гнать со всеми генами, а не только дифэкспрессирующимися

# Объединение результатов анализа обогащения для всех контрастов
compare_results_all <- compareCluster(
  geneClusters = contrast_genes,
  fun = "enricher",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

png("GO_dotplot_all.png", width = 12, height = 8, units = "in", res = 150, type = "cairo")
# Построение dot plot для объединенного результата
dotplot(compare_results_all, showCategory = 30, color = 'qvalue') +
  ggtitle("GO Enrichment Analysis Across Contrasts") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###только для всех контрастов без up and down
contrast_genes_onlyall <- list(
  "AtK3_vs_AtF3" = sig_genes_lfc,
  "AtK5_vs_AtF5" = sig_genes_lfc_5,
  "AtF3_vs_AtF5" = sig_genes_lfc3_5
)


# Объединение результатов анализа обогащения для всех контрастов
compare_results_onlyall <- compareCluster(
  geneClusters = contrast_genes_onlyall,
  fun = "enricher",
  pAdjustMethod = "bonferroni",
  universe = filt_go_annotations$GeneID,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.05,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

png("GO_dotplot_only_all.png", width = 12, height = 8, units = "in", res = 150, type = "cairo")
# Построение dot plot для объединенного результата
dotplot(compare_results_onlyall, showCategory = 30, color = 'qvalue') +
  ggtitle("GO Enrichment Analysis Across Contrasts") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### загружаем контрасты res_
res_AtF3_AtK3 <- readRDS("res_AtF3_AtK3.rds")
res_AtF3_AtF5 <- readRDS("res_AtF3_AtF5.rds")
res_AtK5_AtF5 <- readRDS("res_AtF5_AtK5.rds")


### попробуем использовать другую функцию, которая учитывает lfc при подсчете
# Конвертируем результаты для каждого контраста в data.frame и создаем gene list
# Фильтрация данных и создание отсортированных gene list для контраста AtF3 vs AtK3
df_AtF3_AtK3 <- as.data.frame(res_AtF3_AtK3)
df_AtF3_AtK3 <- df_AtF3_AtK3[!is.na(df_AtF3_AtK3$log2FoldChange) & abs(df_AtF3_AtK3$log2FoldChange) >= 1.5, ]
gene_list_AtF3_AtK3 <- setNames(df_AtF3_AtK3$log2FoldChange, rownames(df_AtF3_AtK3))
gene_list_AtF3_AtK3 <- sort(gene_list_AtF3_AtK3, decreasing = TRUE)

# Повторяем для контраста AtF3 vs AtF5
df_AtF3_AtF5 <- as.data.frame(res_AtF3_AtF5)
df_AtF3_AtF5 <- df_AtF3_AtF5[!is.na(df_AtF3_AtF5$log2FoldChange) & abs(df_AtF3_AtF5$log2FoldChange) >= 1.5, ]
gene_list_AtF3_AtF5 <- setNames(df_AtF3_AtF5$log2FoldChange, rownames(df_AtF3_AtF5))
gene_list_AtF3_AtF5 <- sort(gene_list_AtF3_AtF5, decreasing = TRUE)

# Повторяем для контраста AtK5 vs AtF5
df_AtK5_AtF5 <- as.data.frame(res_AtK5_AtF5)
df_AtK5_AtF5 <- df_AtK5_AtF5[!is.na(df_AtK5_AtF5$log2FoldChange) & abs(df_AtK5_AtF5$log2FoldChange) >= 1.5, ]
gene_list_AtK5_AtF5 <- setNames(df_AtK5_AtF5$log2FoldChange, rownames(df_AtK5_AtF5))
gene_list_AtK5_AtF5 <- sort(gene_list_AtK5_AtF5, decreasing = TRUE)

set.seed(42)

# Запуск GSEA для каждого контраста
gsea_result_AtF3_AtK3 <- GSEA(
  geneList = gene_list_AtF3_AtK3,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  minGSSize = 10,
  maxGSSize = 1000,
  pvalueCutoff = 0.05
)

gsea_result_AtF3_AtF5 <- GSEA(
  geneList = gene_list_AtF3_AtF5,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  minGSSize = 10,
  maxGSSize = 1000,
  pvalueCutoff = 0.05
)

gsea_result_AtK5_AtF5 <- GSEA(
  geneList = gene_list_AtK5_AtF5,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  minGSSize = 10,
  maxGSSize = 1000,
  pvalueCutoff = 0.05
)

# Преобразование результатов GSEA в data.frame и добавление меток контрастов
gsea_df_AtF3_AtK3 <- as.data.frame(gsea_result_AtF3_AtK3@result)
gsea_df_AtF3_AtK3$Contrast <- "AtF3 vs AtK3"

gsea_df_AtF3_AtF5 <- as.data.frame(gsea_result_AtF3_AtF5@result)
gsea_df_AtF3_AtF5$Contrast <- "AtF3 vs AtF5"

gsea_df_AtK5_AtF5 <- as.data.frame(gsea_result_AtK5_AtF5@result)
gsea_df_AtK5_AtF5$Contrast <- "AtK5 vs AtF5"

# Объединение результатов в один data.frame
gsea_all_results <- rbind(gsea_df_AtF3_AtK3, gsea_df_AtF3_AtF5, gsea_df_AtK5_AtF5)

png("GSEA_dotplot_contrasts.png", width = 12, height = 8, units = "in", res = 150, type = "cairo")

ggplot(gsea_all_results, aes(x = Contrast, y = Description, color = p.adjust, size = NES)) +
  geom_point() +
  scale_color_continuous(low = "blue", high = "red", name = "p.adjust") +
  ggtitle("GSEA Dot Plot Across Contrasts") +
  theme_minimal() +
  labs(x = "Contrast", y = "GO Terms") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
