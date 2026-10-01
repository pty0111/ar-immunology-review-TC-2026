setwd("~/0-workspace/CCR7_DC/CB-3141_pool/")

Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")
suppressPackageStartupMessages({
    library(curl)
    library(Seurat)
    library(Rmagic)
    library(Matrix)
    library(future)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(ComplexHeatmap)
    library(circlize)
    library(ggrepel)
    library(readxl)
    library(rasterpdf)
    library(readxl)
    library(pals)
    library(plotly)
    library(anndata)
})
source("../utils.R")
pal <- readRDS("plots/palette.rds")

# ############################################################################ #
sro1 <- subset(sro.subset, subset.annotation2 %in% c("ILCp", "ILC Ki67+", "ILC early", "ILC1", "ILC2", "LTi - 1", "LTi - 2", "LTi - 3", "LTip", "TLp", "TCp", "TC I-IV"))
gene.info <- sro1@assays$RNA@meta.features
sro1$source <- "CB3141"
sro1 <- RenameGenesSeurat(sro1, sro1@assays$RNA@meta.features$ID)

sro1 <- NormalizeData(sro1, verbose = FALSE)
sro1 <- FindVariableFeatures(sro1, selection.method = "vst", nfeatures = 5000, verbose = FALSE)

sro2 <- readRDS("../GSE312479-colonna-2026-tcp-fig2/Seurat/Fig2B/SRO.rds")
sro2$source <- "ColonnaFig2B"
sro2 <- RenameGenesSeurat(sro2, sro2@assays$RNA@meta.features$ID)

sro3 <- readRDS("../GSE278161-sawa-2026-tcp/Seurat/res1_ex_c1c5c13c14c15/SRO.rds")
sro3$source <- "SawaFL"
sro3$tissue <- "FL"
sro3[["RNA"]]$scale.data <- NULL
sro3 <- RenameGenesSeurat(sro3, sro3@assays$RNA@meta.features$ID)

sro4 <- readRDS("../TC_all_LN/results/SRO.rds") %>% subset(Annotation %ni% c("DC contaminant", "low QC"))
features <- read.table(
  '../TC_all_LN/filtered_feature_bc_matrix/features.tsv.gz', 
  header = F, sep = "\t", stringsAsFactors = F, col.names = c("ID", "symbol", "type")
)
features$symbol.unique <- toupper(make.unique(features$symbol))
rownames(features) <- features$symbol.unique
sro4@assays$RNA <- AddMetaData(sro4@assays$RNA, features[rownames(sro4), ])

sro4 <- RenameGenesSeurat(sro4, sro4@assays$RNA@meta.features$ID)
sro4$source <- "Cabric"
sro4$MtFrac_RNA <- sro4$percent.MT

genes <- intersect(intersect(rownames(sro1), rownames(sro2)), intersect(rownames(sro3), rownames(sro4)))
gene.info <- sro1@assays$RNA@meta.features[c('ID', "Symbol", "Symbol.unique")]
gene.info <- gene.info[match(genes, gene.info$ID),]
dim(gene.info)

# ############################################################################ #
pref.i <- "results/res-2-subset-progenitors-option2-v2exC18/integrate-exDCexProg-with-ColonnaFig2B-SawaRes1_ex_c1c5c13c14c15-Cabric/"
dir.create(pref.i)
# write.csv(gene.info, file = paste0(pref.i, "gene-info.csv"))

reduction <- "cca"; pref <- paste0(pref.i, reduction, "/")
dir.create(pref)

## Run integration ####
intg.anchors <- FindIntegrationAnchors(
    object.list = list(sro1, sro2, sro3, sro4),
    normalization.method = "LogNormalize",
    anchor.features = 5000,
    reduction = reduction
)
sro.i <- IntegrateData(
    anchorset = intg.anchors,
    features.to.integrate = rownames(gene.info),
    normalization.method = "LogNormalize"
) 
sro.i <- ScaleData(sro.i, features = rownames(sro.i)) %>%
    RunPCA(features = rownames(sro.i), npcs = 50) %>%
    FindNeighbors(dims = 1:30, k.param = 30) %>%
    RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

###
source.curr <- sources[1]
group.name.1 <- 'subset.annotation2'
p.1 <- plot.clusters.highlight.one(SRO = sro.i, 
                                   idx1 = colnames(sro.i), idx2 = sro.i$source == source.curr, 
                                   groups = get.named.vector.sro(sro.i, group.name.1), clusters.col = paste0(source.curr, '-', group.name.1), labels = T, label.size = 4,
                                   col = pal$subset.annotation2, label.pad = 1, pref.C = F, raster=T) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))

pal$fine.annotations <- pal$Clusters[1:length(unique(sro2$fine.annotations))]
names(pal$fine.annotations) <- unique(sro2$fine.annotations)
###
source.curr <- sources[2]
group.name.2 <- 'fine.annotations'
p.2 <- plot.clusters.highlight.one(SRO = sro.i, 
                                   idx1 = colnames(sro.i), idx2 = sro.i$source == source.curr, 
                                   groups = get.named.vector.sro(sro.i, group.name.2), clusters.col = paste0(source.curr, '-', group.name.2), labels = T, label.size = 4,
                                   col = pal$fine.annotations, label.pad = 1, pref.C = F, raster=T) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))

pal$SawaClusters <- pal$Clusters[1:length(unique(sro3$RNA_snn_res.0.4))]
names(pal$SawaClusters) <- unique(sro3$RNA_snn_res.0.4)
###
source.curr <- sources[3]
group.name.2 <- 'RNA_snn_res.0.4'
p.3 <- plot.clusters.highlight.one(SRO = sro.i, 
                                   idx1 = colnames(sro.i), idx2 = sro.i$source == source.curr, 
                                   groups = get.named.vector.sro(sro.i, group.name.2), clusters.col = paste0(source.curr, '-', group.name.2), labels = T, label.size = 4,
                                   col = pal$SawaClusters, label.pad = 1, pref.C = F, raster=T) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))

pal$Annotation <- pal$Clusters[1:length(unique(sro4$Annotation))]
names(pal$Annotation) <- unique(sro4$Annotation)
###
source.curr <- sources[4]
group.name.2 <- 'Annotation'
p.4 <- plot.clusters.highlight.one(SRO = sro.i, 
                                   idx1 = colnames(sro.i), idx2 = sro.i$source == source.curr, 
                                   groups = get.named.vector.sro(sro.i, group.name.2), clusters.col = paste0(source.curr, '-', group.name.2), labels = T, label.size = 4,
                                   col = pal$Annotation, label.pad = 1, pref.C = F, raster=T) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))

raster_pdf(file = paste0(pref.p, "UMAP_regCC-subset.annotation.pdf"), width = 12*2, height = 8*2, res = 300)
plot_grid(plotlist = list(p.1, p.2, p.3, p.4), ncol = 2)
dev.off()
