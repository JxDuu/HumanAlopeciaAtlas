library(Seurat)
library(data.table)
library(ggplot2)
library(ggpubr)
library(viridis)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
library(tibble)
library(stringr)
library(purrr)
library(RColorBrewer)
library(ComplexHeatmap)
library(clusterProfiler)
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

# Figure 2a - RNA cell type proportion

meta = read.csv('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca260702.meta.csv.gz', row.names=1)
colnames(meta,3)

ratio = meta %>% 
    group_by(Sample, CellType) %>%
    summarise(n = n()) %>%
    mutate(Ratio = n/sum(n) * 100)
ratio$BrodeCluster = meta[match(ratio$CellType, meta$CellType), ]$BrodeCluster
ratio$BrodeCluster = factor(ratio$BrodeCluster, ordered = T, levels = c('Immune','ENDO','Stromal','FB','NC','EG','SG','SkinKC','HFKC'))
ratio = ratio %>% arrange(BrodeCluster, toupper(CellType))
ratio$CellType = factor(ratio$CellType, ordered = T, levels = unique(ratio$CellType))
head(ratio,3)

rect_data <- ratio %>%
    group_by(Sample) %>%
    arrange(desc(CellType)) %>% 
    mutate(xmax = cumsum(Ratio), xmin = xmax - Ratio) %>%
    ungroup()

target_boxes <- rect_data %>%
    filter(BrodeCluster %in% c('Immune', 'HFKC')) %>%
    group_by(Sample, BrodeCluster) %>%
    summarise(xmin = min(xmin), xmax = max(xmax), .groups = "drop")

y_levels <- rev(c('HC','AO','AF','AA'))

target_boxes <- target_boxes %>%
    mutate(y_center = as.numeric(factor(Sample, levels = y_levels)),
           ymin = y_center - 0.4,  
           ymax = y_center + 0.4)

options(repr.plot.width = 8, repr.plot.height = 2.5)
p = ggplot(ratio, aes(x = Ratio, y = Sample, fill = CellType))+
    geom_col(position = "stack", width = 0.8) +
    scale_y_discrete(limits = c('HC','AO','AF','AA') %>% rev)+
    # scale_x_discrete(limits = levels(meta$CellType) %>% rev)+
    scale_fill_manual(values = ct_col) +
    theme_classic()+
    labs(title = 'scRNA-seq cell type', y = 'Percentage (%)')+
    theme_white+
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6))+ 
    geom_rect(data = target_boxes, inherit.aes = FALSE, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "red", linewidth = 0.5)
p
# ggsave(p, file = '00.figure/RNA_CelltypeRaio.boxplot_260901.pdf', width = 8, height = 2.5)



# Figure 2b - IMM cell proportion

ratio1 = meta %>% 
    group_by(batch1, CellType) %>%
    summarise(n = n()) %>%
    pivot_wider(., names_from = batch1, values_from = n, values_fill = 0) %>%
    pivot_longer(., cols = 2:ncol(.), names_to = 'batch1', values_to = 'n') %>%
    group_by(batch1) %>%
    mutate(Ratio = n/sum(n) * 100) %>% 
    arrange(toupper(CellType))
ratio1$Sample = meta[match(ratio1$batch1, meta$batch1), ]$Sample
ratio1$Dataset = meta[match(ratio1$batch1, meta$batch1), ]$Dataset
ratio1

# sfig2_ct = meta %>% select(BrodeCluster, CellType) %>% distinct %>% filter(BrodeCluster %in% c('HFKC','Immune','FB')) %>% arrange(BrodeCluster,toupper(CellType)) %>% pull(CellType)
sfig2_ct = c('DP','DS','DSC','aHFSC','IRS','IRSC','MCC','ORSB','ORSS','qHFSC','TAC1','TAC2','BC','CD4_TC','CD8_TC','cDC','LGS','MACRO','MAST','pDC')
length(sfig2_ct)

p = lapply(sfig2_ct, function(ct){
    df1 = ratio1 %>% filter(CellType == ct)
    p = ggplot(df1, aes(x = Sample, y = Ratio, fill = Sample)) +
    stat_summary(fun = "mean", geom = "col", width = 0.8, color = "black", linewidth = 0.5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, linewidth = 0.5) +
    scale_fill_manual(values = sample_col) + # 替换为你的 sample 颜色
    scale_x_discrete(limits = c('HC', #'AO', 
                                'AF', 'AA'))+
    theme_white +
    labs(title = ct, y = "Ratio (%)") +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.text = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 8),
          legend.position = "none")
    return(p)
}) %>% wrap_plots(., nrow = 4)
# options(repr.plot.width = 45, repr.plot.height = 25)
# p
ggsave(p, file = 'CellType_Ratio.pdf', width = 10, height = 8)

# Figure 2c - IMM GSEA

gsea = read.csv('/data/work/10.AGA_AA/01.RNA/04.subcluster/IMM/GSEA_immune_v1/results/GSEA_significant_padj0.05.csv')
head(gsea,3)

OUT_DIR <- "/data/work/10.AGA_AA/01.RNA/04.subcluster/IMM/GSEA_immune_v1"

RES_DIR  <- file.path(OUT_DIR, "results")
PLOT_DIR <- file.path(OUT_DIR, "plots")

allres_path <- file.path(RES_DIR, "GSEA_all_results_fgsea.csv")
allres <- fread(allres_path)

head(allres,3)

selected <- data.table(
  pathway = c(
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "GOBP_MAST_CELL_ACTIVATION",
    "GOBP_TISSUE_REMODELING"),
  short = c(
    "IFN-g signaling",
    "IFN-a response",
    "MHC-I antigen presentation",
    "TNFb signaling via NF-kB",
    "Mast cell activation",
    "Tissue remodeling"),
  axis = c("AA", "AA", "AA", "AF", "AF", "AF")
)

dat <- allres[pathway %in% selected$pathway]
dat <- merge(dat, selected[, .(pathway, short, axis)], by = "pathway")
# head(dat,3)

main_cells <- c("MACRO", "cDC", "MAST", "CD4_TC")
all_cells  <- c("BC", "LGS", "MAST", "MACRO", "cDC", "CD8_TC", "CD4_TC")

dat[, short := factor(short, levels = selected$short)]
dat[, contrast_lab := fifelse(contrast == "AF-HC", "AF vs HC", "AA vs HC")]
dat[, contrast_lab := factor(contrast_lab, levels = c("AA vs HC", "AF vs HC"))]
dat[, sig := !is.na(padj) & padj < 0.05]
dat[, star := fifelse(sig, "*", "")]
head(dat,3)

plot_dat <- dat[, .(celltype, contrast, contrast_lab, collection, pathway, short,
                    axis, NES, pval, padj, sig)]
plot_dat <- plot_dat[order(celltype, contrast, short)]
# fwrite(plot_dat, file.path(RES_DIR, "GSEA_selected_pathways_plotting_data.csv"))

nes_wide <- dcast(dat, celltype + contrast ~ short, value.var = "NES")
padj_wide <- dcast(dat, celltype + contrast ~ short, value.var = "padj")
nes_mat <- merge(nes_wide, padj_wide,
                 by = c("celltype", "contrast"),
                 suffixes = c("_NES", "_padj"))
head(nes_mat,3)
# fwrite(nes_mat, file.path(RES_DIR, "GSEA_selected_pathways_NES_matrix.csv"))

make_plot <- function(d, cells, title, subtitle, cell_levels) {
  d <- d[celltype %in% cells]
  d[, celltype := factor(celltype, levels = cell_levels)]

  # rectangle highlight behind CD4_TC row (CellChat outgoing hub)
  n <- length(cell_levels)             # CD4_TC is the last level -> top -> numeric pos n
  cd4_y <- as.numeric(factor("CD4_TC", levels = cell_levels))

  p <- ggplot(d, aes(x = celltype, y = short)) +
    # annotate("rect", xmin = 0.5, xmax = length(selected$short) + 0.5,
    #          ymin = cd4_y - 0.5, ymax = cd4_y + 0.5,
    #          fill = "grey30", alpha = 0.08) +
    geom_point(aes(size = abs(NES), color = NES), shape = 16) +
    geom_text(aes(label = star), size = 4.5, vjust = 0.72, hjust = 0.5, color = "black") +
    scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                          midpoint = 0, name = "NES",
                          limits = c(-max(abs(range(d$NES, na.rm = TRUE))),
                                      max(abs(range(d$NES, na.rm = TRUE))))) +
    scale_size_area(max_size = 9, name = "|NES|", limits = c(0, 3)) +
    scale_y_discrete(limits = rev(selected$short))+
    scale_x_discrete(limits = sort(cell_levels))+
    facet_wrap(~ contrast_lab, nrow = 1) +
    labs(x = NULL, y = NULL, title = title) +
    # theme_minimal(base_size = 11) +
    theme_white+
    theme(
      base_size = 8,
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text = element_text(size = 8, face = "bold"),
      # plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 8, hjust = 0),
      panel.spacing = unit(0.5, "lines"),
      # panel.grid.major.x = element_blank(),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.4, "cm"))
  p
}

p_main <- make_plot(
  dat, main_cells,
  title = "Functional pathways",
  # subtitle = "Main figure (4 immune cells \u00d7 6 selected pathways); CD4_TC row shaded (CellChat outgoing hub)",
  cell_levels = main_cells)

options(repr.plot.width = 5, repr.plot.height = 2.5)
p_main
# ggsave(p_main, file = '00.figure/IMM_GSEA_260829.pdf', width = 5, height = 2.5)



# Figure 2d - cellchat result

samples = c('HC','AF','AA')

cco.list = lapply(samples, function(x){
    readRDS(paste0(x, "_cellchat_260916.rds"))
})
names(cco.list) = samples

cellchat = readRDS('immune_DSC_HF_sample_cellchat.rds')

Subtype = sort(unique(cellchat@meta$CellType))
# Subtype

hfsc = c('aHFSC','IRS','IRSC','MCC','ORSB','ORSS','qHFSC','TAC1','TAC2') %>% sort
dsc = c('DP','DS','DSC') %>% sort
imm = setdiff(Subtype,c(dsc,hfsc)) %>% sort

get_data_strength_dotplot = function(cco.list, select.celltype = NULL, color.use, ymax = 0.6, xmax = 0.6){
    object = cco.list
    title = names(cco.list)
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]][[x.measure]]
        incoming[, i] <- centr[[i]][[y.measure]]
    }

    outgoing.cells <- rowSums(outgoing)
    incoming.cells <- rowSums(incoming)
    num.link <- aggregateNet(object, return.object = FALSE, 
                             remove.isolate = FALSE)$count
    num.link <- rowSums(num.link) + colSums(num.link) - diag(num.link)
    df <- data.frame(x = outgoing.cells, y = incoming.cells, 
                     labels = names(incoming.cells), Count = num.link)

    if (!is.null(select.celltype)) df = df[select.celltype, ]
    df$labels <- factor(df$labels, levels = names(incoming.cells))
    return(df)
}

data = lapply(names(cco.list),function(x){
    df = cco.list[[x]]@net$weight %>% as.data.frame %>% tibble::rownames_to_column(., 'Outgoing') %>%
        pivot_longer(., cols = 2:ncol(.), names_to = 'Incoming', values_to = 'weight') %>%
        mutate(Sample = x)
    return(df)
}) %>% Reduce(rbind,.) %>% 
    select(Sample,Outgoing, everything()) %>%
    mutate(Sample = factor(Sample, ordered = T, levels = c('HC','AF','AA')))
data

order = data %>% filter(Sample == 'HC') %>% 
    filter(Outgoing %in% imm) %>% 
    arrange(desc(weight)) %>% pull(Outgoing) %>% unique

options(repr.plot.width = 6, repr.plot.height = 3)
p = ggplot(data[data$Outgoing %in% imm, ], aes(x = Outgoing, y = weight, fill = Sample)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    theme_classic(base_size = 10) +
    scale_fill_manual(values = sample_col) +
    labs(x = "Cell Type", 
         y = "Outgoing Interaction Strength", 
         fill = "Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limits = order)+
    theme_white
p
ggsave(p, file='0.figure/Outgoing.pdf', width = 6, height = 3)

order = data %>% filter(Sample == 'AF') %>% 
    filter(Incoming %in% c(dsc,hfsc)) %>% 
    arrange(desc(weight)) %>% pull(Incoming) %>% unique

options(repr.plot.width = 6.5, repr.plot.height = 3)
p = ggplot(data[data$Incoming %in% c(dsc,hfsc), ], aes(x = Incoming, y = weight, fill = Sample)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    theme_classic(base_size = 10) +
    scale_fill_manual(values = sample_col) +
    labs(x = "Cell Type", 
         y = "Incoming Interaction Strength", 
         fill = "Sample") +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limits = order)+
    theme_white
p
ggsave(p, file='0.figure/Incoming.pdf', width = 6.5, height = 3)



# Figure 2e - IMM-dermal cell chat

cco.list = lapply(samples, function(x){
    readRDS(paste0(x, "_cellchat_260703.rds"))
})
names(cco.list) = samples

cellchat = readRDS('immune_DSC_sample_cellchat.rds')

Subtype = sort(unique(cellchat@meta$CellType))
# Subtype

dsc = c('DP','DS','DSC') %>% sort
imm = setdiff(Subtype,dsc) %>% sort

# 数量与强度差异网络图 红色是case相当于control上调
options(repr.plot.width = 8, repr.plot.height = 8)
pdf('netVisual_diffInteraction_Sample.pdf', width = 20, height = 10)
par(mfrow = c(1,2))
p1 = netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",
                               comparison = c(3, 1), 
                               sources.use = imm, 
                               targets.use = dsc, # remove.isolate=T,
                               color.use=ct_col[Subtype], #c(color1[environment],color1[hair]),
                               margin=0, 
                               title.name = 'AA vs HC'
                              )
p2 = netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",
                               comparison = c(3, 2), 
                               sources.use = imm, 
                               targets.use = dsc,# remove.isolate=T,
                               color.use=ct_col[Subtype], #c(color1[environment],color1[hair]),
                               margin=0, 
                               title.name = 'AF vs HC'
                              )
dev.off()



# Figure 2f - CD4 TC-DP/DS1 LR expression

obj = readRDS('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca.rds')
meta = read.csv('/data/work/10.AGA_AA/01.RNA/AGA_AA_RNA_postQC_mt10_rpca260702.meta.csv.gz', row.names=1)
obj@meta.data = meta[match(rownames(obj@meta.data), rownames(meta)),]

DefaultAssay(obj) = 'RNA'
Idents(obj) = 'CellType'
obj

obj1 = subset(obj, subset = CellType == 'CD4_TC' & Sample != 'AO') %>% DietSeurat(., assay = 'RNA') %>% NormalizeData
obj2 = subset(obj, subset = CellType == 'DSC' & Sample != 'AO') %>% DietSeurat(., assay = 'RNA') %>% NormalizeData
obj3 = subset(obj, subset = CellType == 'DP' & Sample != 'AO') %>% DietSeurat(., assay = 'RNA') %>% NormalizeData

ligand <- c(
    # AA
    "IFNG", "CD40LG", #"MIF", "WNT10A", "WNT10B", "GDNF", "CLEC2D", "NMB", "TIGIT", "GZMA",'IFNA5',
    "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "HLA-DMA", "HLA-DPA1", # "TGFB1", "LIF", "LTA", "TNFSF14"
    # AF
    "MIF", "IL26","NPPC",'CCL8','CCL14','CXCL3','CXCL11','CCL17' #'CCL13','CXCL13',
)
receptor_DSC <- c(
    # "ITGB1", 
    "CD74", "IL10RA", "NPR2", "IFNGR2","ACKR1" #"CXCR4", "IL10RB", 
)
receptor_DP <- c(
    "IFNGR1", "CD40", #"ACKR3", "GFRA1", "RET", "AVPR1A", "NTRK3", 
    "IFNAR1"# , "CD4" "TGFBR1", "TGFBR2", "TNFRSF14", "LTBR", "LIFR", "IL6ST"
)

options(repr.plot.width = 6.25,repr.plot.height = 2)
Idents(obj1) = 'Sample'
p1 = DotPlot(obj1, features = ligand, dot.scale = 10, scale.by = 'size')+#coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(4,10)) +
    scale_y_discrete(limits = c('HC','AF','AA') %>% rev)+
    theme_white+labs(title = 'CD4 TC ligand')+
    theme(plot.title =element_text(size = 10,vjust=0.5, hjust=0.5),
          axis.text.x=element_text(vjust=0.5, hjust=1, size = 8, angle = 90, face = 'italic'), 
          axis.text.y=element_text(size = 8), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p1
ggsave(p1, file = '00.figure/AF_AA_CD4_ligand_dotplot.pdf', width = 6.25, height = 2)

options(repr.plot.width = 3, repr.plot.height = 2)
Idents(obj2) = 'Sample'
p2 = DotPlot(obj2, features = receptor_DSC, dot.scale = 10, scale.by = 'size')+#coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(4, 10)) +
    scale_y_discrete(limits = c('HC','AF','AA') %>% rev)+
    theme_white+labs(title = 'DSC receptor')+
    theme(plot.title =element_text(size = 10,vjust=0.5, hjust=0.5),
          axis.text.x=element_text(vjust=0.5, hjust=1, size = 8, angle = 90, face = 'italic'), 
          axis.text.y=element_text(size = 8), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p2
ggsave(p2, file = '00.figure/AF_AA_DSC_receptor_dotplot.pdf', width = 3, height = 2)

options(repr.plot.width = 2.5, repr.plot.height = 2)
Idents(obj3) = 'Sample'
p3 = DotPlot(obj3, features = receptor_DP, dot.scale = 10, scale.by = 'size')+#coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(4, 10)) +
    scale_y_discrete(limits = c('HC','AF','AA') %>% rev)+
    theme_white+labs(title = 'DP receptor')+
    theme(plot.title =element_text(size = 10,vjust=0.5, hjust=0.5),
          axis.text.x=element_text(vjust=0.5, hjust=1, size = 8, angle = 90, face = 'italic'), 
          axis.text.y=element_text(size = 8), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p3
ggsave(p3, file = '00.figure/AF_AA_DP_receptor_dotplot.pdf', width = 2.5, height = 2)

ls()

obj4 = subset(obj, subset = CellType == 'aHFSC' & Sample != 'AO') %>% DietSeurat(., assay = 'RNA') %>% NormalizeData

options(repr.plot.width = 2.5, repr.plot.height = 2)
Idents(obj4) = 'Sample'
p4 = DotPlot(obj4, features = c('NRP1','NRP2'), dot.scale = 10, scale.by = 'size')+#coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(2, 10)) +
    scale_y_discrete(limits = c('HC','AF','AA') %>% rev)+
    theme_white+labs(title = 'aHFSC receptor')+
    theme(plot.title =element_text(size = 10,vjust=0.5, hjust=0.5),
          axis.text.x=element_text(vjust=0.5, hjust=1, size = 8, angle = 90, face = 'italic'), 
          axis.text.y=element_text(size = 8), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p4
ggsave(p4, file = '00.figure/AF_AA_aHFSC_receptor_dotplot.pdf', width = 2.5, height = 2)

