## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(tidypopgen)
example_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e"),
  population = c("pop1", "pop1", "pop2", "pop2", "pop2")
)
example_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 0, 0, 0, NA, 0),
  c(1, 2, 0, 0, 1, 1),
  c(0, 2, 0, 1, 2, 1),
  c(1, 1, NA, 2, 1, 0)
)
example_loci <- data.frame(
  name = c("rs1", "rs2", "rs3", "rs4", "x1", "x2"),
  chromosome = c(1, 1, 1, 1, 2, 2),
  position = c(3, 5, 65, 343, 23, 456),
  genetic_dist = c(0, 0, 0, 0, 0, 0),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

## -----------------------------------------------------------------------------
example_gt <- gen_tibble(example_genotypes,
  indiv_meta = example_indiv_meta,
  loci = example_loci,
  backingfile = tempfile()
)

## -----------------------------------------------------------------------------
example_gt

## -----------------------------------------------------------------------------
example_gt %>% show_genotypes()

## -----------------------------------------------------------------------------
example_gt %>% show_loci()

## -----------------------------------------------------------------------------
example_gt %>% indiv_het_obs()

## -----------------------------------------------------------------------------
example_gt %>% mutate(het_obs = indiv_het_obs(.data$genotypes))

## -----------------------------------------------------------------------------
example_gt %>% mutate(het_obs = indiv_het_obs(genotypes))

## -----------------------------------------------------------------------------
example_pop2 <- example_gt %>% filter(population == "pop2")
example_pop2

## -----------------------------------------------------------------------------
example_gt %>% indiv_het_obs()

## -----------------------------------------------------------------------------
example_gt %>%
  filter(population == "pop2") %>%
  indiv_het_obs()

## -----------------------------------------------------------------------------
example_gt %>% mutate(het_obs = indiv_het_obs(genotypes))

## -----------------------------------------------------------------------------
loci_names(example_gt)

## -----------------------------------------------------------------------------
example_sub <- example_gt %>% select_loci(starts_with("rs"))
example_sub

## -----------------------------------------------------------------------------
loci_names(example_sub)

## -----------------------------------------------------------------------------
example_sub %>% indiv_het_obs()

## -----------------------------------------------------------------------------
example_gt %>%
  select_loci(c(2, 6, 1)) %>%
  show_loci()

## -----------------------------------------------------------------------------
example_gt %>%
  select_loci(c(2, 6, 1)) %>%
  show_genotypes()

## -----------------------------------------------------------------------------
example_gt %>% loci_maf()

## -----------------------------------------------------------------------------
sel_indices <- which((example_gt %>% loci_maf()) > 0.2)
example_gt %>%
  select_loci(all_of(sel_indices)) %>%
  show_loci()

## -----------------------------------------------------------------------------
example_gt_sub <- example_gt %>% select_loci_if(loci_maf(genotypes) > 0.2)
example_gt_sub %>% show_genotypes()

## -----------------------------------------------------------------------------
example_gt %>%
  select_loci_if(
    loci_chromosomes(genotypes) == 2 &
      loci_maf(genotypes) > 0.2
  ) %>%
  show_loci()

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  tally()

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  summarise(n = n(), mean_het = mean(indiv_het_obs(genotypes)))

## -----------------------------------------------------------------------------
example_gt %>%
  mutate(het_obs = indiv_het_obs(genotypes)) %>%
  group_by(population) %>%
  summarise(n = n(), mean_het = mean(het_obs))

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  reframe(loci_hwe = loci_hwe(genotypes))

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  loci_maf()

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  loci_maf(type = "list")

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  loci_maf(type = "matrix")

## -----------------------------------------------------------------------------
example_gt %>%
  group_by(population) %>%
  pop_fst()

## -----------------------------------------------------------------------------
example_gt %>%
  pairwise_ibs()

## -----------------------------------------------------------------------------
bed_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
pop_a_gt <- gen_tibble(bed_path_pop_a, backingfile = tempfile("pop_a_"))

## -----------------------------------------------------------------------------
gt_as_plink(example_gt, file = tempfile("new_bed_"))

## -----------------------------------------------------------------------------
bed_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
pop_a_gt <- gen_tibble(bed_path_pop_a, backingfile = tempfile("pop_a_"))
bed_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
pop_b_gt <- gen_tibble(bed_path_pop_b, backingfile = tempfile("pop_b_"))

## -----------------------------------------------------------------------------
pop_a_gt

## -----------------------------------------------------------------------------
pop_b_gt

## -----------------------------------------------------------------------------
report <- rbind_dry_run(pop_a_gt, pop_b_gt, flip_strand = TRUE)

## -----------------------------------------------------------------------------
# #create merge
merged_gt <- rbind(pop_a_gt, pop_b_gt,
  flip_strand = TRUE,
  backingfile = file.path(tempdir(), "gt_merged")
)

## -----------------------------------------------------------------------------
merged_gt

## -----------------------------------------------------------------------------
merged_gt %>% show_loci()

## -----------------------------------------------------------------------------
bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(bed_file, backingfile = tempfile("missing_"))
missing_gt

## ----fig.alt = "Histogram of loci missingness, showing that most loci have no missing data, while some have a small proportion of missing data"----
missing_gt %>%
  loci_missingness() %>%
  hist()

## ----error = TRUE-------------------------------------------------------------
try({
missing_pca <- missing_gt %>% gt_pca_autoSVD()
})

## -----------------------------------------------------------------------------
missing_gt <- gt_impute_simple(missing_gt, method = "mode")

## -----------------------------------------------------------------------------
gt_has_imputed(missing_gt)

## -----------------------------------------------------------------------------
gt_uses_imputed(missing_gt)

## ----fig.alt = "Histogram of loci missingness as above, showing that the use of imputed data is not automatic"----
missing_gt %>%
  loci_missingness() %>%
  hist()

## ----fig.alt = "Histogram of loci missingness after setting use of imputed data to true, showing that there is no missingness"----
gt_set_imputed(missing_gt, set = TRUE)
missing_gt %>%
  loci_missingness() %>%
  hist()

## ----fig.alt = "Histogram of loci missingness after setting use of imputed data to false again"----
gt_set_imputed(missing_gt, set = FALSE)
missing_gt %>%
  loci_missingness() %>%
  hist()

## -----------------------------------------------------------------------------
missing_pca <- missing_gt %>% gt_pca_partialSVD()
missing_pca

## ----fig.alt = "Histogram of loci missingness, showing that use of imputed data is not automatically set to true, after using a PCA function"----
gt_uses_imputed(missing_gt)
missing_gt %>%
  loci_missingness() %>%
  hist()

## -----------------------------------------------------------------------------
gt_file_name <- gt_save(example_gt)
gt_file_name

## -----------------------------------------------------------------------------
gt_get_file_names(example_gt)

## -----------------------------------------------------------------------------
new_example_gt <- gt_load(gt_file_name[1])
new_example_gt %>% show_genotypes()

## -----------------------------------------------------------------------------
new_example_gt <- new_example_gt %>% filter(!population == "pop1")

## ----error = TRUE-------------------------------------------------------------
try({
gt_impute_simple(new_example_gt)
})

## -----------------------------------------------------------------------------
new_example_gt <- gt_update_backingfile(new_example_gt,
  backingfile = tempfile()
)

## -----------------------------------------------------------------------------
gt_impute_simple(new_example_gt)

## -----------------------------------------------------------------------------
is_loci_table_ordered(new_example_gt)

## -----------------------------------------------------------------------------
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 2, 1, 1, 1, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(c(0.01, 0.01, 0.03, 0.03, 0.02, 0.015)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)

## -----------------------------------------------------------------------------
is_loci_table_ordered(test_gt)

## -----------------------------------------------------------------------------
reorder_test_gt <- gt_order_loci(test_gt)

## -----------------------------------------------------------------------------
is_loci_table_ordered(reorder_test_gt)

## -----------------------------------------------------------------------------
reorder_test_gt_again <- gt_order_loci(reorder_test_gt,
  ignore_genetic_dist = FALSE
)

## -----------------------------------------------------------------------------
is_loci_table_ordered(reorder_test_gt_again, ignore_genetic_dist = FALSE)

