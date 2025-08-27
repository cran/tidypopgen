## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include = FALSE----------------------------------------------------------
library(tidypopgen)
# Create gen_tibble for examples
bed_path <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path <- bigsnpr::snp_readBed(bed_path, backingfile = tempfile())
data <- gen_tibble(bigsnp_path, quiet = TRUE)

## -----------------------------------------------------------------------------
data %>% select_loci_if(loci_chromosomes(genotypes) %in% c(1:22))

## -----------------------------------------------------------------------------
my_snps <- c("rs4477212", "rs3094315", "rs3131972", "rs12124819", "rs11240777")

data %>%
  select_loci_if(loci_names(genotypes) %in% my_snps) %>%
  show_loci()

## -----------------------------------------------------------------------------
my_individuals <- c("GRC14300079", "GRC14300142", "GRC14300159")

data %>% filter(id %in% my_individuals)

