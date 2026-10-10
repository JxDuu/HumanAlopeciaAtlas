library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
library(tibble)
library(stringr)
library(RColorBrewer)
library(ComplexHeatmap)
library(parallel)
library(future)
options(future.global.maxSize = 500 * 1024^3)
set.seed(2026)

source('/data/work/00.script/colorlist.R')

color_djx = readxl::read_excel('/data/work/10.AGA_AA/01.RNA/Colorlist_DJX.xlsx')

brode_col = color_djx[color_djx$Group == 'BrodeCluster', ]$Color
names(brode_col) = color_djx[color_djx$Group == 'BrodeCluster', ]$Type

ct_col = color_djx[color_djx$Group == 'Celltype', ]$Color
names(ct_col) = color_djx[color_djx$Group == 'Celltype', ]$Type

sample_col = color_djx[color_djx$Group == 'Sample_type', ]$Color
names(sample_col) = color_djx[color_djx$Group == 'Sample_type', ]$Type

obj = readRDS('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca.rds')
meta = read.csv('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca260702.meta.csv.gz', row.names=1)
obj@meta.data = meta[match(rownames(obj@meta.data), rownames(meta)),]

head(obj,3)

atac = readRDS('/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi.rds')
meta = read.csv('/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi.meta260825.csv.gz',row.names=1)
atac@meta.data = meta[match(rownames(atac@meta.data), rownames(meta)),]

get_points_within_radius_grid <- function(points, r) {
  # points: matrix/data.frame, n x 2
  points <- as.matrix(points)
  stopifnot(ncol(points) >= 2)
  x <- points[, 1]
  y <- points[, 2]
  n <- nrow(points)
    
  gx <- floor(x / r)
  gy <- floor(y / r)
  
  cell_key <- paste(gx, gy, sep = "_")
  cell_index <- split(seq_len(n), cell_key)
  
  counts <- integer(n)
  
  for (i in seq_len(n)) {
    cx <- gx[i]
    cy <- gy[i]
    
    neighbor_keys <- as.vector(outer((cx-1):(cx+1), (cy-1):(cy+1), paste, sep = "_"))
    candidate_idx <- unlist(cell_index[neighbor_keys], use.names = FALSE)
    
    if (length(candidate_idx) == 0) {
      counts[i] <- 0L
      next
    }
    
    dx <- x[candidate_idx] - x[i]
    dy <- y[candidate_idx] - y[i]
    counts[i] <- sum(dx*dx + dy*dy <= r*r)
  }
  
  counts
}

# Figure 1b

dmt_mat <- Embeddings(obj, reduction = "dmt")[, 1:2]
obj$r_points <- get_points_within_radius_grid(dmt_mat, r = 0.5)
threshold <- quantile(obj$r_points, 0.005, na.rm = TRUE)
subdata <- subset(obj, subset = r_points > threshold)

options(repr.plot.width = 8.5, repr.plot.height = 8)
p = DimPlot(subdata, group.by = c('CellType'), reduction = 'dmt', pt.size=0.01, cols = ct_col, label = F, label.size = 4, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = '0.figure/RNA_Celltype_dimplot.pdf', width = 8.5, height = 8)

options(repr.plot.width = 8.5, repr.plot.height = 8)
p = DimPlot(subdata, group.by = c('Sample'), reduction = 'dmt', pt.size=0.1, cols = sample_col, shuffle = T, label = F, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = '0.figure/RNA_Sample_dimplot.pdf', width = 6.5, height = 6)

options(repr.plot.width = 8.5, repr.plot.height = 8)
p = DimPlot(subdata, group.by = c('Dataset'), reduction = 'dmt', pt.size=0.1, cols = best_color, shuffle = T, label = F, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = '0.figure/RNA_Dataset_dimplot.pdf', width = 6.5, height = 6)

# Figure 1c

dmt_mat <- Embeddings(atac, reduction = "dmt")[, 1:2]
atac$r_points <- get_points_within_radius_grid(dmt_mat, r = 0.6)
threshold <- quantile(atac$r_points, 0.005, na.rm = TRUE)
subdata <- subset(atac, subset = r_points > threshold)

p = DimPlot(subdata, group.by = c('CellType'), reduction = 'dmt', pt.size=0.1, cols = best_color, shuffle = T, label = F, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = 'ATAC_CellType_dimplot.pdf', width = 6.5, height = 6)

p = DimPlot(subdata, group.by = c('Dataset'), reduction = 'dmt', pt.size=0.1, cols = best_color, shuffle = T, label = F, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = 'ATAC_Dataset_dimplot.pdf', width = 6.5, height = 6)

options(repr.plot.width = 8.5, repr.plot.height = 8)
p = DimPlot(subdata, group.by = c('Sample'), reduction = 'dmt', pt.size=0.1, cols = sample_col, shuffle = T, label = F, raster = F)+NoLegend()+
    # theme_white+#coord_fixed()+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
# p
ggsave(p, file = 'ATAC_Sample_dimplot.pdf', width = 6.5, height = 6)





# Figure 1d

ct_marker = c(
    "GJA5", "ITGB8", "PDGFRB", "IGHG1", "KRT15", "KRTDAP", "CD3D", "IL7R", 
    "CLEC9A", "CA4", "RSPO4", "ITGB5", "COL11A1", "KRT19", "SLC12A2", "DCD", 
    "HBB", "GRHL3", "SPINK5", "KIT", "S100A8", "KRT27", "KRT82", "KRT79", 
    "FLT4", "CD86", "MSR1", "CPA3", "DSG4", "KRT8", "NRXN1", "KRT75", 
    "KRT6A", "RGS5", "EREG", "COLEC12", "MKI67", "SOX9", "SFRP1", "COL14A1", 
    "ACSBG1", "CDH1", "MITF", "ANK3", "LEF1", "WNT10B", "ACKR1"
)

DefaultAssay(obj) = 'RNA'
Idents(obj) = 'CellType'
p = DotPlot(obj, features = ct_marker, dot.scale = 6, scale.by = 'size')+ 
    scale_y_discrete(limits = rev(sort(unique(obj$CellType))))

data = p$data %>% select(features.plot,id,avg.exp.scaled) %>%
    pivot_wider(., names_from = features.plot, values_from = avg.exp.scaled) %>%
    column_to_rownames(.,'id')
data = data[sort(unique(obj$CellType)), ct_marker] %>% as.matrix
head(data,3)

mode(data) <- "numeric"
data[is.na(data)] <- 0

data_plot <- data
data_plot[data_plot > 2] <- 2
data_plot[data_plot < -2] <- -2

options(repr.plot.width = 7, repr.plot.height = 7)
col_fun = col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#67A9CF", "#F7F7F7", "#D73027"))
heatmap_legend_param <- list(
  title = "Exp.",
  at = c(-2, -1, 0, 1, 2), 
  labels = c("-2", "-1", "0", "1", "2"), 
  title_gp = gpar(fontsize = 8, col = "black"), 
  labels_gp = gpar(fontsize = 6, col = "black"), 
  border = "black"
)
ht <- ComplexHeatmap::Heatmap(
    data_plot, 
    name = "Exp.",
    col = col_fun,
    row_names_gp = gpar(fontsize = 6, col = "black"),
    column_names_gp = gpar(fontsize = 6, col = "black"),
    row_names_side = "left",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = "black",
    width = unit(14.1, "cm"),    
    height = unit(14.1, "cm"),
    heatmap_legend_param = heatmap_legend_param)
ht
# pdf('0.figure/RNA_CellMarker_heatmap_260608.pdf', width = 7, height = 11)
# ht
# dev.off()

# Figure 1e

atac = readRDS('/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi_dmt_ctpeak.rds')
meta = read.csv('/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi.meta260825.csv.gz',row.names=1)
atac@meta.data = meta[match(rownames(atac@meta.data), rownames(meta)),]
atac

Idents(atac) = 'CellType'
atac$CellType = factor(as.character(atac$CellType),ordered = T, levels = sort(unique(atac$CellType)))
# levels(atac$CellType)

peakmarker1_select = read.csv('Figure1e_markerpeak.csv', row.names=1)
peakmarker1_select

DefaultAssay(atac) = 'ATAC'
lapply(1:nrow(peakmarker1_select), function(x){
    Idents(atac) = 'CellType'
    ct = peakmarker1_select$celltype[x]
    peak = peakmarker1_select$peak[x]
    gene = peakmarker1_select$gene_name[x]
    p = CoveragePlot(
        object = atac,
        region = peak,
        peaks = F,
        # annotation=F,
        # region.highlight = regions_highlight,
        extend.upstream = 0,
        extend.downstream = 0
    ) & scale_fill_manual(values = ct_col)
    ggsave(p, file = paste0('MarkerPeak_plot_ctpeak/',ct,'_',gene,'_',peak,'.pdf'), width = 1.75, height = 20)
})



# Supplementary Figure 1a

meta = obj@meta.data

cellnum = meta %>% 
    group_by(Dataset, Sample, batch1, Age) %>%
    summarise(log1p_n = log1p(n()))
cellnum$Sample = factor(cellnum$Sample, ordered = T, levels = c('HC','AO','AF','AA'))
cellnum = cellnum %>% 
    arrange(Sample,Age)
cellnum

options(repr.plot.width = 4, repr.plot.height = 8)
p = ggplot(cellnum, aes(x = log1p_n, y = batch1, fill = Sample))+
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_y_discrete(limits = cellnum$batch1 %>% rev)+
    scale_fill_manual(values = sample_col) +
    theme_classic()+
    labs(title = 'scRNA')+
    theme_white+
    theme(legend.position = "right",
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6)
         )
p
# ggsave(p, file = 'RNA_batch_cellnumber.boxplot.pdf', width = 4, height = 8)



# Supplementary Figure 1b

meta = atac@meta.data

cellnum = meta %>% 
    group_by(Dataset, Sample, batch1, Age) %>%
    summarise(log1p_n = log1p(n()))
cellnum$Sample = factor(cellnum$Sample, ordered = T, levels = c('HC','AO','AF','AA'))
cellnum = cellnum %>% 
    arrange(Sample,Age)
cellnum

options(repr.plot.width = 4, repr.plot.height = 6.5)
p = ggplot(cellnum, aes(x = log1p_n, y = batch1, fill = Sample))+
    # geom_boxplot(width = 0.75, outliers = FALSE)+
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_y_discrete(limits = cellnum$batch1 %>% rev)+
    scale_fill_manual(values = sample_col) +
    theme_classic()+
    labs(title = 'scATAC-seq cell number')+
    theme_white+
    theme(legend.position = "right",
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
         )
    # facet_wrap(~Dataset)
p
# ggsave(p, file = '0.figure/ATAC_batch_cellnumber.boxplot.pdf', width = 4, height = 6.5)





# Supplementary Figure 1c

set.seed(1234)

rna = readRDS('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca.rds')
meta = read.csv('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca260702.meta.csv.gz', row.names=1)
rna@meta.data <- meta[rownames(rna@meta.data),]
rna

DefaultAssay(rna) <- "RNA"
rna <- DietSeurat(rna, assays = "RNA")

rna <- NormalizeData(
  rna,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

set.seed(1234)
rna_meta <- rna@meta.data
rna_meta$cellid = rownames(rna_meta)
prop = rna_meta %>% 
    group_by(CellType) %>%
    summarize(n=n()) %>%
    mutate(prop=n/sum(n))
# prop

total_cell = 85000
meta_sub = lapply(prop$CellType, function(ct){
    rna_meta[rna_meta$CellType == ct,] %>% slice_sample(n = round(total_cell*prop[prop$CellType == ct, ]$prop))
}) %>% Reduce(rbind, .)

rna1 <- subset(rna, cells = meta_sub$cellid)
rna1

rna1.list <- SplitObject(rna1, split.by = "batch1")

for (i in 1:length(rna1.list)) {
    rna1.list[[i]] <- NormalizeData(rna1.list[[i]], verbose = FALSE)
    rna1.list[[i]] <- FindVariableFeatures(rna1.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)%>%ScaleData(.,verbose = F) %>% RunPCA(., npcs = 30)
}
features <-SelectIntegrationFeatures(object.list = rna1.list)
rna1.anchors <- FindIntegrationAnchors(object.list = rna1.list, normalization.method = "LogNormalize",anchor.features = features,
                                      reduction = "rpca",  dims = 1:30)
rna1<- IntegrateData(anchorset = rna1.anchors, dims = 1:30)

rna1 <- ScaleData(rna1, verbose = FALSE)
rna1 <- RunPCA(rna1, npcs = 30, verbose = FALSE)
rna1 <- RunUMAP(rna1, reduction = "pca", dims = 1:30)

DefaultAssay(rna1) <- 'RNA'
rna1 = rna1 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
rna = rna1
DefaultAssay(rna) <- 'RNA'
rna$Tech = 'RNA'

atac = readRDS("/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi.rds")
meta2 = read.csv("/data/work/10.AGA_AA/02.ATAC/02.Signac/AGA_AA_ATAC_postQC_rlsi.meta260825.csv.gz",row.names = 1)
atac@meta.data <- meta2[rownames(atac@meta.data),]
atac

DefaultAssay(atac) <- "RNA"
atac$Tech = 'ATAC'
##### 03. CCA integrated #####
DefaultAssay(atac) <- "RNA"
features <- SelectIntegrationFeatures(object.list = list(rna, atac))
anchors <- FindIntegrationAnchors(
  object.list = list(rna, atac),
  anchor.features = features,
  assay = c('RNA', 'RNA')
)
# integrate data and create a new merged object
integrated <- IntegrateData(
  anchorset = anchors,
  weight.reduction = c('pca', 'integrated_lsi'),
  dims = 2:30,
  preserve.order = TRUE
)
integrated = ScaleData(integrated)
integrated <- RunPCA(
  object = integrated,
  npcs = 5,
  reduction.name = 'integr_PCA'
)
integrated <- RunUMAP(
  object = integrated,
  dims = 1:30,
  reduction = 'integr_PCA',
  min.dist = 0.3
)

# integrated <- ScaleData(obj, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE, reduction.name = 'integr_PCA')
integrated <- RunUMAP(integrated, reduction = "integr_PCA", dims = 1:30)
library(reticulate)
np <- import("numpy")
pca <- Embeddings(object = integrated, reduction = "integr_PCA") %>% as.matrix()
np$save("/data/work/03.ATAC/AA_AGA_RNA_ATAC_CCA_0531.npy", pca)

obj = integrated
# 导入dmt结果到Seurat对象
celltrain="40000"
nu="0.6"
np <- import("numpy")
dmt_highdim = np$load(paste0('/data/work/03.ATAC/01.dmt/DMT_celltrain_',celltrain,'_nu_',nu,'/X_dmt_highdim.npy'))
dmt = np$load(paste0("/data/work/03.ATAC/01.dmt/DMT_celltrain_",celltrain,'_nu_',nu,'/X_dmt.npy'))

rownames(dmt_highdim)= Cells(obj)
colnames(dmt_highdim)= paste0('dmt_highdim_',1:50)
rownames(dmt)= Cells(obj)
colnames(dmt)= paste0('dmt_',1:2)
obj[["dmt_highdim"]] <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(dmt_highdim), key = "dmt_highdim_ ", assay ="integrated")
obj[["dmt"]] <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(dmt), key = "dmt_", assay ="integrated")
obj = obj %>% FindNeighbors(., reduction = 'dmt_highdim', dims = 1:30) %>% 
    FindClusters(., resolution = 0.5)

p = DimPlot(obj, group.by = 'Tech',
             cols = best_color, pt.size = 0.01, label = FALSE, 
             label.size = 3.5, raster = FALSE, shuffle = FALSE, reduction = 'dmt')
p



# Supplementary Figure 1d

ratio = meta %>% 
    group_by(Dataset, Sample, batch1, Age, CellType) %>%
    summarise(n = n()) %>%
    mutate(Ratio = n/sum(n) * 100) %>% 
    arrange(desc(Sample),Age, toupper(CellType))
ratio

options(repr.plot.width = 12, repr.plot.height = 3)
p = ggplot(ratio, aes(x = batch1, y = Ratio, fill = CellType))+
    geom_col(position = "stack", width = 0.7) +
    scale_x_discrete(limits = cellnum$batch1)+
    scale_fill_manual(values = ct_col) +
    theme_classic()+
    labs(title = 'scRNA-seq cell type', y = 'Percentage (%)')+
    theme_white+
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
p
# ggsave(p, file = '0.figure/RNA_batch_CelltypeRaio.boxplot.pdf', width = 12, height = 3)



# Supplementary Figure 1e

meta = atac@meta.data

ratio = meta %>% 
    group_by(Dataset, Sample, batch1, Age, CellType) %>%
    summarise(n = n()) %>%
    mutate(Ratio = n/sum(n) * 100) %>% 
    arrange(desc(Sample),Age, toupper(CellType))
ratio

options(repr.plot.width = 12, repr.plot.height = 3)
p = ggplot(ratio, aes(x = batch1, y = Ratio, fill = CellType))+
    geom_col(position = "stack", width = 0.7) +
    scale_x_discrete(limits = cellnum$batch1)+
    scale_fill_manual(values = ct_col) +
    theme_classic()+
    labs(title = 'scATAC-seq cell type', y = 'Percentage (%)')+
    theme_white+
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
         )
p
# ggsave(p, file = 'ATAC_batch_CelltypeRaio.boxplot.pdf', width = 12, height = 3)

