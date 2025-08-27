## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("tidypopgen",
#   repos = c(
#     "https://evolecolgroup.r-universe.dev",
#     "https://cloud.r-project.org"
#   )
# )

## -----------------------------------------------------------------------------
library(tidypopgen)
library(ggplot2)

## -----------------------------------------------------------------------------
lobsters <- gen_tibble(
  x = system.file("extdata/lobster/lobster.bed", package = "tidypopgen"),
  quiet = TRUE, backingfile = tempfile()
)
head(lobsters)

## -----------------------------------------------------------------------------
lobsters %>% show_genotypes(indiv_indices = 1:5, loci_indices = 1:10)

## -----------------------------------------------------------------------------
head(lobsters %>% show_loci())

## -----------------------------------------------------------------------------
indiv_qc_lobsters <- lobsters %>% qc_report_indiv()

## ----fig.alt = "Scatter plot of missingness proportion and observed heterozygosity for each individual"----
autoplot(indiv_qc_lobsters, type = "scatter")

## -----------------------------------------------------------------------------
lobsters <- lobsters %>% filter(indiv_missingness(genotypes) < 0.2)

## -----------------------------------------------------------------------------
loci_qc_lobsters <- lobsters %>% qc_report_loci()

## -----------------------------------------------------------------------------
lobsters <- lobsters %>% group_by(population)
loci_qc_lobsters <- lobsters %>% qc_report_loci()

## ----fig.alt = "Histogram of minor allele frequency"--------------------------
autoplot(loci_qc_lobsters, type = "maf")

## ----fig.alt = "Histogram of the proportion of missing data"------------------
autoplot(loci_qc_lobsters, type = "missing")

## -----------------------------------------------------------------------------
lobsters <- lobsters %>% select_loci_if(loci_missingness(genotypes) < 0.05)

## -----------------------------------------------------------------------------
lobsters <- gt_update_backingfile(lobsters, backingfile = tempfile())

## -----------------------------------------------------------------------------
lobsters <- gt_impute_simple(lobsters, method = "random")

## -----------------------------------------------------------------------------
partial_pca <- gt_pca_partialSVD(lobsters)

## ----fig.alt = "Score plot of individuals across the first and second Principal Components"----
autoplot(partial_pca, type = "scores")

## ----fig.alt = "Score plot of individuals across the first and second Principal Components, with individual samples coloured by population"----
autoplot(partial_pca, type = "scores") +
  aes(color = lobsters$population) +
  labs(color = "population")

## -----------------------------------------------------------------------------
pcs <- augment(x = partial_pca, data = lobsters)

## -----------------------------------------------------------------------------
eigenvalues <- tidy(partial_pca, "eigenvalues")

xlab <- paste("Axis 1 (", round(eigenvalues[1, 3], 1), " %)",
  sep = ""
)
ylab <- paste("Axis 2 (", round(eigenvalues[2, 3], 1), " %)",
  sep = ""
)

## ----fig.alt = "Score plot of individuals across the first and second Principal Components, with individual samples coloured by population"----
ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = population),
    shape = 21, size = 3, show.legend = FALSE
  ) +
  scale_fill_distruct() +
  labs(x = xlab, y = ylab) +
  ggtitle("Lobster PCA")

## -----------------------------------------------------------------------------
# Calculate centre for each population
centroid <- aggregate(cbind(.fittedPC1, .fittedPC2, .fittedPC3) ~ population,
  data = pcs, FUN = mean
)

# Add these coordinates to our augmented pca object
pcs <- left_join(pcs, centroid, by = "population", suffix = c("", ".cen"))

## ----fig.alt = "Score plot of individuals across the first and second Principal Components. Individual samples are coloured by population, the centroid of points for each group is labelled with the population name, and lines are drawn from the centroid to each point."----
ggplot(data = pcs, aes(x = .fittedPC1, y = .fittedPC2)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(xend = .fittedPC1.cen, yend = .fittedPC2.cen),
    show.legend = FALSE
  ) +
  geom_point(aes(fill = population),
    shape = 21, size = 3, show.legend = FALSE
  ) +
  geom_label(
    data = centroid,
    aes(label = population, fill = population),
    size = 4, show.legend = FALSE
  ) +
  scale_fill_distruct() +
  labs(x = xlab, y = ylab) +
  ggtitle("Lobster PCA")

