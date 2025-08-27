## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)

## -----------------------------------------------------------------------------
library(tidypopgen)
data <- gen_tibble(
  system.file("extdata/related/families.bed",
    package = "tidypopgen"
  ),
  quiet = TRUE, backingfile = tempfile(),
  valid_alleles = c("1", "2")
)

## -----------------------------------------------------------------------------
individual_report <- qc_report_indiv(data)
summary(individual_report)

## ----fig.alt = "Scatter plot of missingness proportion and observed heterozygosity for each individual"----
autoplot(individual_report)

## -----------------------------------------------------------------------------
data <- data %>% filter(indiv_missingness(genotypes) < 0.045)
nrow(data)

## -----------------------------------------------------------------------------
mean_val <- mean(individual_report$het_obs)
sd_val <- stats::sd(individual_report$het_obs)

lower <- mean_val - 3 * (sd_val)
upper <- mean_val + 3 * (sd_val)

data <- data %>% filter(indiv_het_obs(genotypes) > lower)
data <- data %>% filter(indiv_het_obs(genotypes) < upper)
nrow(data)

## -----------------------------------------------------------------------------
individual_report <- qc_report_indiv(data, kings_threshold = 0.177)
summary(individual_report)

## -----------------------------------------------------------------------------
data <- data %>%
  filter(id %in% individual_report$id & individual_report$to_keep == TRUE)

## -----------------------------------------------------------------------------
summary(data)

## -----------------------------------------------------------------------------
loci_report <- qc_report_loci(data)
summary(loci_report)

## ----fig.alt = "UpSet plot giving counts of snps over the threshold for: missingness, minor allele frequency, and Hardy-Weinberg equilibrium P-value"----
autoplot(loci_report, type = "overview")

## ----fig.alt = "Upset plot as above, with adjusted thresholds"----------------
autoplot(loci_report,
  type = "overview",
  miss_threshold = 0.03,
  maf_threshold = 0.02,
  hwe_p = 0.01
)

## ----fig.alt = "Four panel plot, containing: a histogram of the proportion of missing data for snps with minor allele frequency above the threshold, a histogram of the proportion of missing data for snps with minor allele freqency below the threshold, a histogram of HWE exact test p-values, and a histogram of significant HWE exact test p-values"----
autoplot(loci_report,
  type = "all",
  miss_threshold = 0.03,
  maf_threshold = 0.02,
  hwe_p = 0.01
)

## ----fig.alt = "Histogram of minor allele frequency"--------------------------
autoplot(loci_report, type = "maf")

## -----------------------------------------------------------------------------
data <- data %>% select_loci_if(loci_maf(genotypes) > 0.02)
count_loci(data)

## ----fig.alt = "Histogram of the proportion of missing data"------------------
autoplot(loci_report, type = "missing", miss_threshold = 0.05)

## -----------------------------------------------------------------------------
data <- data %>% select_loci_if(loci_missingness(genotypes) < 0.05)
count_loci(data)

## ----fig.alt = "Histogram of significant HWE exact test p-values"-------------
autoplot(loci_report, type = "significant hwe", hwe_p = 0.01)

## -----------------------------------------------------------------------------
data <- data %>% select_loci_if(loci_hwe(genotypes) > 0.01)
count_loci(data)

## -----------------------------------------------------------------------------
data <- gt_update_backingfile(data)

## -----------------------------------------------------------------------------
imputed_data <- gt_impute_simple(data, method = "random")

## -----------------------------------------------------------------------------
to_keep_ld <- loci_ld_clump(imputed_data, thr_r2 = 0.2, size = 10)
head(to_keep_ld)

## -----------------------------------------------------------------------------
ld_data <- imputed_data %>%
  select_loci_if(loci_ld_clump(genotypes, thr_r2 = 0.2, size = 10))

## -----------------------------------------------------------------------------
gt_save(ld_data, file_name = tempfile())

## -----------------------------------------------------------------------------
data <- data %>% mutate(population = c(rep("A", 4), rep("B", 5)))

## -----------------------------------------------------------------------------
grouped_loci_report <- data %>%
  group_by(population) %>%
  qc_report_loci()
head(grouped_loci_report)

## -----------------------------------------------------------------------------
grouped_individual_report <- data %>%
  group_by(population) %>%
  qc_report_indiv(kings_threshold = 0.177)
head(grouped_individual_report)

## -----------------------------------------------------------------------------
loci_maf_grouped <- data %>%
  group_by(population) %>%
  loci_maf()

