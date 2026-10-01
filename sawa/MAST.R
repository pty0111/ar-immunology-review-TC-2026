setwd("/data1/lesliec/tyler/CCR7_DC/GSE278161-sawa-2026-tcp/")
Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")

suppressPackageStartupMessages({
    library(curl)
    library(Seurat)
    library(Rmagic)
    library(Matrix)
    library(ggplot2)
    library(cowplot)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(ggrepel)
    library(anndata)
})

set.seed(1)
options(future.globals.maxSize = Inf)
plan(strategy = "multicore", workers = 8)

source("../utils.R")
pal <- readRDS("plots/palette.rds")

# ############################################################################ #
# MAST ####
# ############################################################################ #
pref.sro <- 'Seurat/res1_ex_c1c5c13c14c15/'; pref.p.sro <- "plots/res1_ex_c1c5c13c14c15/"
dir.create(paste0(pref.sro, "markers"))
dir.create(paste0(pref.p.sro, "markers"))

sro <- readRDS(paste0(pref.sro, "SRO.rds"))
imp.expr <- readRDS("Seurat/imputed-expr.rds")

# ############################################################################ #
# 1. one-vs-rest ####
## 0.8 ####
group.name <- "RNA_snn_res.1"
fn <- paste0(group.name, "-all")
n <- 20

sro.input <- sro
DefaultAssay(sro.input) <- "RNA"

res <- MAST.DE.multiple(fn, sro.input, idents = group.name, output.dir = paste0(pref.sro, 'markers/'))

markers <- select.markers(fn, output.dir = paste0(pref.sro, 'markers/'), n=n)
write.csv(markers, paste0(pref.sro, 'markers/', fn, '-top-', n, '-markers.csv'), row.names = T, quote = F)

all.markers <- select.markers(fn, output.dir = paste0(pref.sro, 'markers/'))
write.csv(all.markers, paste0(pref.sro, 'markers/', fn, '-sig-markers.csv'), row.names = T, quote = F)

### heatmap ####
pal$named.Clusters <- pal$Clusters[1:length(unique(sro.input@meta.data[[group.name]]))]
names(pal$named.Clusters) <- sort(unique(sro.input@meta.data[[group.name]]))

### heatmap ####
pal$sample <- pal$Clusters[1:length(unique(sro.input@meta.data[['sample']]))]
names(pal$sample) <- sort(unique(sro.input@meta.data[['sample']]))



ca.col.list <- list(pal$sample, pal$named.Clusters)
names(ca.col.list) <- c('sample', group.name)

h <- nrow(markers)*0.2+2
pdf(width = 50, height = h, file = paste0(pref.p.sro, "markers/hm-", fn, "-top-", n, ".pdf"))
DE.heatmap(sro.input, expr = imp.expr, legend.name = "scaled\nimputed\nexpression",
           cells = colnames(sro.input), markers = markers$gene, md = sro.input@meta.data,
           split = group.name, cluster_column_slices = T, column_title_rot = 0,
           rowsplit = markers$cluster, cluster_rows = T, cluster_row_slices = T,
           split.columns = names(ca.col.list),
           ca.col = ca.col.list)
dev.off()
