## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(123)

## -----------------------------------------------------------------------------
library(tidypopgen)
vcf_path <-
  system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
anole_gt <-
  gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))

## -----------------------------------------------------------------------------
anole_gt

## -----------------------------------------------------------------------------
pops_path <- system.file("/extdata/anolis/punctatus_n46_meta.csv",
  package = "tidypopgen"
)
pops <- read.csv(pops_path)
pops

## ----error=TRUE---------------------------------------------------------------
try({
anole_gt <- anole_gt %>% left_join(pops, by = "id")
})

## -----------------------------------------------------------------------------
anole_gt %>% glimpse()

## -----------------------------------------------------------------------------
anole_gt <- gt_add_sf(anole_gt, c("longitude", "latitude"))
anole_gt

## -----------------------------------------------------------------------------
library(rnaturalearth)
library(ggplot2)

map <- ne_countries(
  continent = "South America",
  type = "map_units", scale = "medium"
)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = anole_gt$geometry) +
  coord_sf(
    xlim = c(-85, -30),
    ylim = c(-30, 15)
  ) +
  theme_minimal()

## ----error=TRUE---------------------------------------------------------------
try({
anole_pca <- anole_gt %>% gt_pca_partialSVD(k = 30)
})

## -----------------------------------------------------------------------------
anole_gt <- gt_impute_simple(anole_gt, method = "mode")

## -----------------------------------------------------------------------------
anole_pca <- anole_gt %>% gt_pca_partialSVD(k = 30)

## -----------------------------------------------------------------------------
anole_pca

## -----------------------------------------------------------------------------
tidy(anole_pca, matrix = "eigenvalues")

## ----fig.alt = "Scree Plot of eigenvalues for each Principal Component"-------
autoplot(anole_pca, type = "screeplot")

## ----fig.alt = "Score plot of individuals across the first and second Principal Components"----
autoplot(anole_pca, type = "scores")

## ----fig.alt = "Score plot of individuals across the first and second Principal Components, with individual samples coloured by population"----
library(ggplot2)
autoplot(anole_pca, type = "scores") +
  aes(color = anole_gt$population) +
  labs(color = "population")

## -----------------------------------------------------------------------------
anole_gt <- augment(anole_pca, data = anole_gt)

## ----fig.alt = "Another score plot of individuals across the first and second Principal Components, with individual samples coloured by population, using ggplot2"----
anole_gt %>% ggplot(aes(.fittedPC1, .fittedPC2, color = population)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", color = "Population")

## ----fig.alt = "Plot of the loadings of Principal Component 1 for each loci"----
autoplot(anole_pca, type = "loadings")

## -----------------------------------------------------------------------------
anole_gt_load <- augment_loci(anole_pca, data = anole_gt)

## ----fig.alt = "Plot of the cumulative loadings for each Principal Component"----
library(ggplot2)
tidy(anole_pca, matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC, y = cumulative)) +
  geom_point()

## -----------------------------------------------------------------------------
anole_clusters <- gt_cluster_pca(anole_pca, n_pca = 10)

## ----fig.alt = "Plot of BIC (Bayesian Information Criteria) for each value of k"----
autoplot(anole_clusters)

## -----------------------------------------------------------------------------
anole_clusters <- gt_cluster_pca_best_k(anole_clusters)

## -----------------------------------------------------------------------------
anole_clusters$best_k

## -----------------------------------------------------------------------------
anole_dapc <- gt_dapc(anole_clusters)

## -----------------------------------------------------------------------------
anole_dapc

## -----------------------------------------------------------------------------
tidy(anole_dapc, matrix = "eigenvalues")

## ----fig.alt = "Scree plot of the eigenvalues on the two discriminant axes (defined by `ld`)"----
autoplot(anole_dapc, type = "screeplot")

## ----fig.alt = "Bar plot of eigenvalues against the two discriminant axes"----
tidy(anole_dapc, matrix = "eigenvalues") %>%
  ggplot(aes(x = LD, y = eigenvalue)) +
  geom_col()

## ----fig.alt = "Scatterplot of the scores of each individual on the two discriminant axes (defined by `ld`)"----
autoplot(anole_dapc, type = "scores")

## ----fig.alt = "Bar plot showing the probability of assignment to each cluster"----
autoplot(anole_dapc, type = "components", group = anole_gt$population)

## ----fig.alt = "Plot of the loadings of Principal Component 1 for each loci"----
autoplot(anole_dapc, "loadings")

## ----fig.alt = "Scatterplot of individuals with ellipses surrounding each cluster of individuals identified by the DAPC analysis. Eigenvalues are displayed in an inset in the bottom right corner of the plot."----
library(adegenet)
scatter(anole_dapc, posi.da = "bottomright")

## -----------------------------------------------------------------------------
anole_gt <- anole_gt %>% mutate(dapc = anole_dapc$grp)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = anole_gt$geometry, aes(color = anole_gt$dapc)) +
  coord_sf(
    xlim = c(-85, -30),
    ylim = c(-30, 15)
  ) +
  labs(color = "DAPC cluster") +
  theme_minimal()

## ----results='hide', echo=TRUE, eval=FALSE------------------------------------
# anole_gt <- anole_gt %>% group_by(population)
# 
# anole_adm_original <- gt_admixture(
#   x = anole_gt,
#   k = 2:6,
#   n_runs = 1,
#   crossval = TRUE,
#   seed = 1
# )

## ----results="hide", echo=FALSE, eval=TRUE, warning=FALSE---------------------
anole_gt <- anole_gt %>% group_by(population)

anole_adm <- readRDS("a03_anole_adm.rds")

## ----fig.alt = "Plot of cross-entropy for each value of K, showing K = 3 has the lowest cross-entropy"----
autoplot(anole_adm, type = "cv")

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources"----
autoplot(anole_adm,
  k = 3, run = 1, data = anole_gt, annotate_group = FALSE,
  type = "barplot"
)

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources"----
autoplot(anole_adm,
  type = "barplot", k = 3, run = 1, data = anole_gt,
  annotate_group = TRUE, arrange_by_group = TRUE
)

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources, with individuals arranged according to their dominant ancestry propotion"----
autoplot(anole_adm,
  type = "barplot", k = 3, run = 1,
  data = anole_gt, annotate_group = TRUE, arrange_by_group = TRUE,
  arrange_by_indiv = TRUE, reorder_within_groups = TRUE
)

## -----------------------------------------------------------------------------
q_mat <- get_q_matrix(anole_adm, k = 3, run = 1)

## -----------------------------------------------------------------------------
anole_gt_adm <- augment(q_mat, data = anole_gt)
head(anole_gt_adm)

## -----------------------------------------------------------------------------
tidy_q <- tidy(q_mat, anole_gt)
head(tidy_q)

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources"----
tidy_q <- tidy_q %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(dominant_q = max(percentage)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(dominant_q)) %>%
  dplyr::mutate(
    plot_order = dplyr::row_number(),
    id = factor(id, levels = unique(id))
  )

plt <- ggplot2::ggplot(tidy_q, ggplot2::aes(x = id, y = percentage, fill = q)) +
  ggplot2::geom_col(
    width = 1,
    position = ggplot2::position_stack(reverse = TRUE)
  ) +
  ggplot2::labs(
    y = "Population Structure for K = 3",
    title = "ADMIXTURE algorithm on A. punctatus"
  ) +
  theme_distruct() +
  scale_fill_distruct()

plt

## -----------------------------------------------------------------------------
anole_adm <- gt_admix_reorder_q(anole_adm, group = anole_gt$pop)

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources, labelled Atlantic forest or Amazonian forest"----
autoplot(anole_adm,
  type = "barplot", k = 3, run = 1,
  annotate_group = TRUE, arrange_by_group = TRUE,
  arrange_by_indiv = TRUE, reorder_within_groups = TRUE
)

## ----echo=FALSE, results='hide'-----------------------------------------------
# reload the admxiture results (they were reordered above)
anole_adm <- readRDS("a03_anole_adm.rds")
runs <- c(1)
k_values <- c(2:6)
q_matrices <- list()
adm_dir <- file.path(tempdir(), "anolis_adm")
if (!dir.exists(adm_dir)) {
  dir.create(adm_dir)
}

for (x in k_values) {
  q_matrices[[x]] <- list()

  for (i in runs) {
    q_matrices[[x]][[i]] <- get_q_matrix(anole_adm, k = x, run = i)

    write.table(q_matrices[[x]][[i]],
      file = paste0(adm_dir, "/K", x, "run", i, ".Q"),
      col.names = FALSE, row.names = FALSE, quote = FALSE
    )
  }
}

## -----------------------------------------------------------------------------
adm_dir <- file.path(tempdir(), "anolis_adm")
list.files(adm_dir)

## -----------------------------------------------------------------------------
q_list <- read_q_files(adm_dir)
summary(q_list)

## -----------------------------------------------------------------------------
head(get_q_matrix(q_list, k = 3, run = 1))

## ----fig.alt = "Barplot of individuals coloured by predicted ancestry proportion (Q) from each of K ancestral sources"----
autoplot(get_q_matrix(q_list, k = 3, run = 1),
  data = anole_gt,
  annotate_group = TRUE, arrange_by_group = TRUE,
  arrange_by_indiv = TRUE, reorder_within_groups = TRUE
)

