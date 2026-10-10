library(Seurat)
library(monocle)
library(ggplot2)
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
library(circlize)
library(viridis)

library(parallel)
library(future)
options(future.global.maxSize = 500 * 1024^3)
set.seed(2026)

setwd('/data/work/10.AGA_AA/01.RNA/04.subcluster/DSC')
source('/data/work/00.script/colorlist.R')
source('/data/work/00.script/markerlist.r')

color_djx = readxl::read_excel('/data/work/10.AGA_AA/01.RNA/Colorlist_DJX.xlsx')

brode_col = color_djx[color_djx$Group == 'BrodeCluster', ]$Color
names(brode_col) = color_djx[color_djx$Group == 'BrodeCluster', ]$Type

ct_col = color_djx[color_djx$Group == 'Celltype', ]$Color
names(ct_col) = color_djx[color_djx$Group == 'Celltype', ]$Type

sample_col = color_djx[color_djx$Group == 'Sample_type', ]$Color
names(sample_col) = color_djx[color_djx$Group == 'Sample_type', ]$Type

# Figure 3b - DEG count

deg = read.csv('/data/work/10.AGA_AA/01.RNA/03.Ratio_DEG/AGA_AA_RNA_CellType_260702_Sample_DEG_mt10ds_logfc0.5_p0.05.csv', row.names=1)
head(deg,3)

deg_sum = deg %>%
    group_by(celltype,group1,group) %>%
    filter(group1 %in% c('AF-HC')) %>%
    summarise(n = n()) %>%
    arrange(toupper(celltype), desc(n))
deg_sum$n[deg_sum$group == 'Down'] = -deg_sum$n[deg_sum$group == 'Down']
# deg_sum

order.af = deg_sum %>% 
    filter(group == 'Down') %>%
    group_by(celltype) %>%
    summarise(n_total = sum(abs(n))) %>%
    arrange(desc(n_total)) %>%
    pull(celltype)
# order.af

p1 = ggplot(deg_sum, aes(x = celltype, y = n, fill = group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    # scale_fill_manual(values = c("Down" = '#6ad7fb', "Up" = '#e840ba'))+ 
    scale_x_discrete(limits = order.af) +
    theme_white + labs(y = 'AF-HC')+ 
    theme(base_size = 6,
          legend.position = 'right',
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    scale_fill_manual(values = c("Down" = '#B8D8AF', "Up" = '#B09FCB')) 
    # scale_x_continuous(limits = c(-2000,4700))

deg_sum = deg %>%
    group_by(celltype,group1,group) %>%
    filter(group1 %in% c('AA-HC')) %>%
    summarise(n = n()) %>%
    arrange(toupper(celltype), desc(n))
deg_sum$n[deg_sum$group == 'Down'] = -deg_sum$n[deg_sum$group == 'Down']
# deg_sum

order.aa = deg_sum %>% 
    filter(group == 'Down') %>%
    group_by(celltype) %>%
    summarise(n_total = sum(abs(n))) %>%
    arrange(desc(n_total)) %>%
    pull(celltype)
# order.aa

p2= ggplot(deg_sum, aes(x = celltype, y = n, fill = group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    # scale_fill_manual(values = c("Down" = '#6ad7fb', "Up" = '#e840ba'))+ 
    scale_x_discrete(limits = order.aa) +
    theme_white + labs(y = 'AA-HC')+
    theme(base_size = 6,
          legend.position = 'right',
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    scale_fill_manual(values = c("Down" = '#B8D8AF', "Up" = '#B09FCB'))
    # scale_x_continuous(limits = c(-2000,4700))

options(repr.plot.width = 7, repr.plot.height = 4)
p = wrap_plots(list(p1,p2), ncol=1)
p
# ggsave(p, file = '00.figure/DEG_mt10ds_all_barplot_AA-AF_split_260703.pdf', width = 7, height = 4)



# Figure 3c - dermal subcluster

options(repr.plot.width = 6, repr.plot.height = 5)
p = DimPlot(obj, group.by = c('Subtype'), pt.size=0.4, cols = ct_col, label = F, raster = F, reduction = 'dmt')+#NoLegend()+
    theme_white+coord_fixed(0.9)+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p
ggsave(p, file = '00.figure/DSC_subtype_dimplot.png', width = 6, height = 5)



# Figure 3d - dermal ratio

ratio = obj@meta.data %>% 
    group_by(Sample, Subtype) %>%
    summarise(n = n()) %>%
    mutate(Ratio = n/sum(n) * 100) %>% 
    arrange(Sample, Subtype)

options(repr.plot.width = 4, repr.plot.height = 3.5)
p = ggplot(ratio, aes(x = Sample, y = Ratio, fill = Subtype))+
    geom_col(position = "stack", width = 0.8) +
    scale_x_discrete(limits = c('HC','AO','AF','AA'))+
    # scale_x_discrete(limits = levels(meta$CellType) %>% rev)+
    scale_fill_manual(values = ct_col) +
    theme_classic()+
    labs(title = 'scRNA-seq cell type', y = 'Percentage (%)')+
    theme_white+
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
         )
p
ggsave(p, file = '00.figure/DSC_SubtypeRaio.boxplot.pdf', width = 4, height = 3.5)

ratio1 = ratio %>%
    filter(Subtype %in% c('DS1','dDS1'))
ratio1$Ratio[ratio1$Subtype == 'DS1'] = -ratio1$Ratio[ratio1$Subtype == 'DS1']

ratio1_plot <- ratio1 %>%
  mutate(Ratio_plot = ifelse(Ratio < 0, Ratio / 2, Ratio))
p1 = ggplot(ratio1_plot, aes(x = Sample, y = Ratio_plot, fill = Subtype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_discrete(limits = c('HC','AO','AF','AA')) +
    scale_y_continuous(
        breaks = c(-20, -10, 0, 10, 20), 
        labels = c("-40", "-20", "0", "10", "20")
    ) +
    theme_white + ggtitle('DSC state') + 
    scale_fill_manual(values = ct_col)
# p1

ratio1 = ratio %>%
    filter(Subtype %in% c('DS2','dDS2'))
ratio1$Ratio[ratio1$Subtype == 'DS2'] = -ratio1$Ratio[ratio1$Subtype == 'DS2']
p2 = ggplot(ratio1, aes(x = Sample, y = Ratio, fill = Subtype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_discrete(limits = c('HC','AO','AF','AA')) +
    theme_white + ggtitle('DS state')+ #theme(legend.position = 'none')+
    scale_fill_manual(values = ct_col)
    # scale_x_continuous(limits = c(-2000,4700))

ratio1 = ratio %>%
    filter(Subtype %in% c('DP','dDP'))
ratio1$Ratio[ratio1$Subtype == 'DP'] = -ratio1$Ratio[ratio1$Subtype == 'DP']
ratio1_plot <- ratio1 %>%
  mutate(Ratio_plot = ifelse(Ratio < 0, Ratio / 5, Ratio))
p3 = ggplot(ratio1_plot, aes(x = Sample, y = Ratio_plot, fill = Subtype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_discrete(limits = c('HC','AO','AF','AA')) +
    scale_y_continuous(
        breaks = c(-10, -5, 0, 5, 10), 
        labels = c("-50", "-25", "0", "5", "10")
    ) +
    theme_white + ggtitle('DP state')+ #theme(legend.position = 'none')+
    scale_fill_manual(values = ct_col)
    # scale_x_continuous(limits = c(-2000,4700))

options(repr.plot.width = 12, repr.plot.height = 3)
p = wrap_plots(list(p1,p2,p3), nrow=1)
p

ggsave(p,file = '00.figure/DSC_Subtype.ratio.pdf', width = 12, height=3)



# Figure 3e - dermal subcluster marker

obj$Subtype = factor(as.character(obj$Subtype), ordered = T, levels = c('DS1','dDS1','DS2','dDS2','DP','dDP'))

gene = c(
    "LGR5","MME","SFRP1", #DSC
    "MYL4","CTGF","CYR61", #dDSC
    "EGFR","SFRP4","PPARG", #DS "LGR4","THY1",
    "MMP2","CXCL12","ANGPTL5", #dDS
    "MOXD1","WNT10A","GREM1", #DP
    "SRD5A2","WIF1","BMP7" # dDP
)

options(repr.plot.width = 3, repr.plot.height = 4)
Idents(obj) = 'Subtype'
p = DotPlot(obj, features = rev(gene), dot.scale = 10, scale.by = 'size')+ coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(0, 6)) +
    scale_y_discrete(limits = levels(obj$Subtype))+
    theme_white+
    theme(axis.text.y=element_text(vjust=0.5, hjust=1, size = 10, face = 'italic'), 
          axis.text.x=element_text(vjust=0.5, hjust=1, size = 10, angle=90), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p
ggsave(p, file = '00.figure/DSC_Subtype_marker_dotplot.pdf', width = 3, height = 4)



# Figure 3f - DS1 vs dDS1 GO

data = data[c(
    'Developmental growth','Regulation of cell activation','Hair follicle development cytodifferentiation stage 3 of 3',
    'Transforming growth factor beta receptor superfamily signaling pathway','Epithelial cell differentiation','Negative regulation of apoptotic process','Cell chemotaxis',
    
    'Cellular response to starvation','Cellular responses to stress','Striated muscle contraction','Type i collagen synthesis in the context of osteogenesis imperfecta',
    'Response to oxidative stress','Regulation of intrinsic apoptotic signaling pathway','Response to reactive oxygen species'
),]

new_row_names <- str_wrap(rownames(data), width = 32)
new_row_names

options(repr.plot.width = 3, repr.plot.height = 6)
# col_fun = colorRamp2(c(0,1,2,2.5,3,3.5,4,4.5,5), magma(9))
col_fun = colorRamp2(c(0,-log10(0.05), -log10(0.01),-log10(1e-3),-log10(1e-5)), magma(5))
# col_fun = circlize::colorRamp2(quantile(data, c(0, 0.5, 0.9)), magma(11)[c(2,6,11)])
heatmap_legend_param = list(
    title = '-Log10P',
    title_gp = gpar(fontsize = 8, col = "black"),
    labels_gp = gpar(col = "black", fontsize = 6),
    border = 'black'
)
# cluster <- colnames(data)
# column_ha <- HeatmapAnnotation(
#     cluster = cluster,
#     col = list(cluster = setNames(c('#e840ba','#e840ba','#6ad7fb','#6ad7fb'), cluster)),
#     na_col = "grey"
# )
# row_ha = rowAnnotation(group = rownames(data), col = list(group = setNames(c('#e840ba','#e840ba','#6ad7fb','#6ad7fb'), rownames(data))))
ht = Heatmap(data,
        row_labels = new_row_names,
        col = col_fun,
        name = '-Log10P',
        heatmap_legend_param = heatmap_legend_param,
        # column_title = "Target", 
        # column_title_gp = gpar(fontsize = 10, col = "black"),
        row_title = "DSC vs dDSC",
        row_title_gp = gpar(fontsize = 8, col = "black"),
        # row_title_side = 'left',
        row_names_side = "right",
        border = "black",
        cluster_rows = FALSE, 
        cluster_columns = F,
        # top_annotation = column_ha,
        # left_annotation = row_ha,
        width = unit(2, "cm"),    
        height = unit(10, "cm"),
        row_names_gp = gpar(fontsize = 6, col = "black"),
        column_names_gp = gpar(fontsize = 6, col = "black")
)
# ht
pdf(paste0('00.figure/DSC_heatmap_260621.pdf'), width = 3, height = 5)
ht
dev.off()



# Figure 3g - DS1 vs dDS1

dsc_deg = FindMarkers(obj, ident.1 = 'DS1', ident.2 = 'dDS1', logfc.threshold=0, only.pos=F)
head(dsc_deg,3)

dsc_deg$gene = rownames(dsc_deg)
dsc_deg$Subtype = ifelse(dsc_deg$avg_log2FC > 0, 'DS1', 'dDS1')

dsc_deg = dsc_deg %>% arrange(desc(avg_log2FC))
write.csv(dsc_deg, file = 'DS1_vs_dDS1.deg.csv')

dsc_deg1 = dsc_deg %>% 
    filter(gene %in% c('APOD','FBLN1','LTBP2','TIMP3','ITGBL1','IGFBP7','FGF7', #DSC
                       'COL1A1','MYL4','CYR61','COL3A1','COL1A2','COL5A2','COX6C','COL4A1' #dDSC
                      ))

options(repr.plot.width = 3, repr.plot.height = 4)
p1 = ggplot(dsc_deg, aes(x = -log10(p_val), y = avg_log2FC)) +  
     geom_point(size = 0.8, color = 'grey') + 
     geom_point(data = dsc_deg1, aes(x = -log10(p_val), y = avg_log2FC, color = Subtype), size = 0.8) +  
     geom_hline(yintercept = c(log10(0.5), -log10(0.5)), linewidth = 0.5, color = "grey50", lty = 'dashed') +
     theme_white+
     scale_color_manual(values = c('#4e6691','#b8474d'))+
     ggrepel::geom_label_repel(data = dsc_deg1, aes(x = -log10(p_val), y = avg_log2FC, label = gene, color = Subtype),
                               fontface = 'italic', size = 3, seed = 233, min.segment.length = 0, force = 12, force_pull = 2, 
                               box.padding = 0.2,label.padding = 0.1,label.r = 0.1,label.size = 0.2,fill = "white",alpha = 0.8,
                               max.overlaps = Inf, segment.alpha = 0.5)
p1
# ggsave(p1, file = '00.figure/DSC_vs_dDSC_DEG_volcano_260825.pdf', width = 3, height = 4)



# Figure 3h - dDS1 DEG in samples

dsc = subset(obj, subset = Subtype == 'dDS1' & Sample != 'AO') %>% NormalizeData
options(repr.plot.width = 2.5, repr.plot.height = 2.5)

Idents(dsc) = 'Sample'
p = DotPlot(dsc, features = rev(c('SEMA3C','EREG','INHBA','AREG','TOP2A',"SFRP2","WIF1",'WNT5A','WNT11',"TGFB1") %>% unique), dot.scale = 8, scale.by = 'size')+ coord_flip() +
    scale_color_viridis(option = "D") +
    scale_size_continuous(range = c(2,6)) +
    scale_y_discrete(limits = c('HC','AF','AA'))+
    theme_white+
    theme(axis.text.x=element_text(vjust=0.5, hjust=0.5, size = 10), 
          axis.text.y=element_text(size = 10, face = 'italic'), 
          axis.title = element_blank())+
    guides(size = guide_legend(title = "Ratio (%)"))+
    guides(color = guide_colorbar(title = "Exp."))
p



# Figure 3i - SCENIC+

tf2g = read.table('/data/work/10.AGA_AA/02.ATAC/03.scenic/DSC_sub/out/tf_to_gene_adj.tsv', header=T)
head(tf2g,3)

dsc_tf = c(
    'SOX2','ALX4','MAF',
    'MSX1','MYC','NFIB','NFIA','AHR',
    'JUN','FOS','EGR1',
    'CEBPB','TWIST1',
    'LEF1','GLIS1','TFAP2A'
) %>% unique

eRegulon_name_select = c(
    'SOX2_direct_+/+','ALX4_direct_+/+','MAF_direct_+/+',
    'JUN_direct_+/+','FOS_direct_+/+','EGR1_direct_+/+',
    'MSX1_direct_+/-','MYC_direct_+/+','NFIB_direct_+/+','NFIA_direct_+/-','AHR_direct_+/+',
    'CEBPB_direct_+/+','TWIST1_direct_+/-',
    'LEF1_direct_+/+','GLIS1_direct_+/+','TFAP2A_direct_+/+'
)

unique(data_tf$eRegulon_name ) %>% sort

data_tf = read.csv('/data/work/10.AGA_AA/02.ATAC/03.scenic/DSC_sub/out/DSC_sub_All.plotting_df_all.csv', row.names=1) %>%
    dplyr::filter(repressor_activator == 'activator')
data_tf = data_tf %>%
    dplyr::filter(eRegulon_name %in% eRegulon_name_select)
data_tf$eRegulon = sapply(data_tf$eRegulon_name, function(x){str_split(x,'_')[[1]][1]})
# data_tf = data_tf %>%
#     dplyr::filter(eRegulon %in% dsc_tf)
data_tf$eRegulon = factor(as.character(data_tf$eRegulon), ordered = T, levels = dsc_tf)
data_tf$eRegulon1 = paste0(data_tf$eRegulon,'(+)')
data_tf = arrange(data_tf, eRegulon)
# data_tf

order = c('DS1','dDS1','DS2','dDS2','DP','dDP')

data_tf_mt = data_tf %>% dplyr::select(scRNA_counts.cell_type,direct_gene_based_AUC,eRegulon1) %>%
    pivot_wider(., names_from = eRegulon1, values_from = direct_gene_based_AUC) %>%
    tibble::column_to_rownames(.,'scRNA_counts.cell_type') %>% as.matrix
# data_tf_mt = data_tf %>% dplyr::select(scRNA_counts.cell_type,direct_gene_based_AUC,eRegulon_name) %>%
#     pivot_wider(., names_from = eRegulon_name, values_from = direct_gene_based_AUC) %>%
#     tibble::column_to_rownames(.,'scRNA_counts.cell_type') %>% as.matrix
data_tf_mt = data_tf_mt[order, ] %>% t
# data_tf_mt

options(repr.plot.width = 3, repr.plot.height = 3)
col_fun = circlize::colorRamp2(quantile(data_tf_mt, c(0, 0.35, 0.9)), brewer.pal(11, "RdYlBu")[c(9,6,2)])
heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, col = "black"),
    labels_gp = gpar(fontsize = 6, col = "black"),
    border = "black"
)
ht <- ComplexHeatmap::Heatmap(
    data_tf_mt, 
    name = "TF score",
    col = col_fun,
    row_names_gp = gpar(fontsize = 6, col = "black"),
    column_names_gp = gpar(fontsize = 6, col = "black"),
    row_names_side = "left",
    cluster_rows = F,
    cluster_columns = FALSE,
    border = "black",
    width = unit(ncol(data_tf_mt)*0.3, "cm"),    
    height = unit(nrow(data_tf_mt)*0.3, "cm"),
    heatmap_legend_param = heatmap_legend_param)
ht
pdf('DSC_sub_regulon_260625.pdf', width = 3, height = 3)
ht
dev.off()



# Figure 3j - GRN

dsc_tf = c(
    'SOX2','ALX4','MAF',
    'MSX1','MYC','NFIB','NFIA','AHR',
    'JUN','FOS','EGR1',
    'CEBPB','TWIST1',
    'LEF1','GLIS1','TFAP2A'
) %>% unique

tf2g1 = readxl::read_excel('DSC_sub_TF2G_select.xlsx',sheet=2)
tf2g1 = lapply(tf2g1$TF, function(x){
    df = data.frame(target = str_split(tf2g1[tf2g1$TF == x, ]$Select_target, ' / ')[[1]])
    df$TF = x
    return(df)
}) %>% Reduce(rbind,.)
tf2g1

tf2g2 = lapply(unique(tf2g1$TF), function(x){
    target = tf2g1[tf2g1$TF == x, ]$target
    df = tf2g[tf2g$TF == x & tf2g$target %in% target & tf2g$regulation == 1, ]
    return(df)
}) %>% Reduce(rbind,.)
tf2g2 = tf2g2 %>% mutate(TF = paste0(TF, '(+)'))
tf2g2

library(ggraph)
library(tidygraph)

edges <- tf2g2 %>% select(TF, target, importance)
graph_data <- as_tbl_graph(edges, directed = TRUE) %>%
  mutate(node_type = ifelse(name %in% tf2g2$TF, "TF", "Target"))

p <- ggraph(graph_data, layout = 'fr') +
    geom_edge_link(aes(color = importance), arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(4, 'mm'), alpha = 0.9,width = 0.8) +
    scale_edge_color_gradient(low = "#FFF59D", high = "#F57F17", name = "Importance") +
    geom_node_point(aes(color = node_type), size = 8) +
    scale_color_manual(values = c("TF" = "#D4C2D8", "Target" = "#B3DCE3"),
                       name = "Node Type") +
    geom_node_text(aes(label = name, size = node_type), repel = TRUE) +
    scale_size_manual(values = c("TF" = 4, "Target" = 3), guide = "none") +
    theme_graph(base_family = "")

options(repr.plot.width = 10, repr.plot.height = 10)
p
ggsave(p, file = '00.figure/DSC_TF2G_net_260903.pdf', width = 10, height = 10)



# Supplementary Figure 3e - dermal subcluster

options(repr.plot.width = 6, repr.plot.height = 5)
p = DimPlot(obj, group.by = c('Sample'), pt.size=0.4, cols = sample_col, label = F, raster = F, reduction = 'dmt')+#NoLegend()+
    theme_white+coord_fixed(0.9)+
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
p
ggsave(p, file = '00.figure/DSC_sample_dimplot.png', width = 6, height = 5)

