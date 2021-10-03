

# Make a tse object from the avilable human gastrointestinal microbiota data

setwd("")

library("ggplot2"); packageVersion("ggplot2")
library(tidyr)
library(mia)
library(miaViz)
library(tidyverse)
library(scater)
library(ape)

#Read the data 

samples_df <- read.csv(file ="HGMA.web.metadata.csv", 
                       header = TRUE)

otu_mat <- read.csv(file ="HGMA.web.MSP.abundance.matrix.csv", 
                    header = TRUE)

tax_mat <- read.table(file='IGC2.1989MSPs.taxo.tsv',sep = '\t', header = TRUE)

tree <- read.tree("IGC2.1990MSPs.nwk",text = NULL, tree.names = NULL, skip = 0,comment.char = "", keep.multi = FALSE)


otu_mat <- otu_mat %>%coordinates
  tibble::column_to_rownames("X")
otu_mat <- as.matrix(otu_mat)
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("msp_name")
samples_df <- samples_df %>%
  tibble::column_to_rownames("sample.ID")


dim(otu_mat)
dim(samples_df)
dim(tax_mat)



tax_mat <- tax_mat[,c("phylum", "class", "order", "family", "genus","species" )]
se <- TreeSummarizedExperiment(assays = list(counts = otu_mat[,rownames(samples_df)]),
                               colData = samples_df,
                               rowData = tax_mat[rownames(otu_mat),])
se

#how many taxa and samples the data contains

dim(se)


#TreeSummarizedExperiment object which includes also a rowTree slot

tse <- as(se, "TreeSummarizedExperiment")
rowTree(tse) <- tree

tse


common.nodes <- intersect(rownames(tse), rowTree(tse)$tip.label)
tse <- TreeSummarizedExperiment::subsetByLeaf(x = tse, rowLeaf = common.nodes, updateTree = TRUE)


