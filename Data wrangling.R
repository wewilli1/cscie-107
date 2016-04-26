## define folder
setwd("C:/Users/guy12/Documents/Programming/cscie107/Project/Annotated Cluster Multidimensional Enrichment (ACME) analysis")

## load excel file into r
library(xlsx)
dt <- read.xlsx("CTRPv2.2._INFORMER_SET.xlsx", 1)
## remove columns not needed
library(dplyr)
dt = select(dt, -top_test_conc_umol, -target_or_activity_of_compound, -cpd_status, -index_cpd)

## remove rows with NA
dt <- dt[complete.cases(dt), ]


## split commma separated gene symbol into rows
library(tidyr)
dt <- dt %>% 
mutate(gene_symbol_of_protein_target= strsplit(as.character(gene_symbol_of_protein_target), ";")) %>%
  unnest(gene_symbol_of_protein_target)

## load data of cancer cell line: one gene is associated with many cell lines. Use the ccle_primary_site instead
library(data.table)
ccl <- fread("v22.anno.ccl_mut_features.txt")
ccl <- select(ccl, hugo_gene_symbol, index_ccl)
cell_primary_site <- fread("v22.meta.per_cell_line.txt")
cell_primary_site <- select(cell_primary_site, index_ccl, ccle_primary_site, ccle_primary_hist)