library(Seurat)
library(data.table)
library(sf)
library(hdf5r)
library(tidyverse)
library(scop)
library(ggplot2)
library(ggpubr)
library(STDistance)
library(scTenifoldKnk)
library(dplyr)
library(ggrepel)
library(grid)
library(reshape2)


seurat.res <- readRDS("seurat.NSCLC_Brip1.rds")
load("cols.Major_Spatial.Polygons.RData")   #cols.Major

##---------Fig. 6A
p <- FeatureStatPlot(
  srt = seurat.res, group.by = "celltype2", plot_type = "box",    
  stat.by = c("Ifnar1","Il2rb"), palcolor = cols.celltype, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
) 
ggsave('Fig.6A.pdf', p, width = 5, height = 5)



##---------Fig. 6B
fts <- c("Ifnar1","Ifnar2", "Il15","Il2rb")
p <- FeatureDimPlot(srt = seurat.res, features = fts, 
                     reduction = "UMAP", raster=T, pt.size = 2, theme_use = "theme_blank", ncol = 2)
ggsave("Fig.6B.pdf", p, width = 6, height = 6)




##---------Fig. 6C (Fig. 6E top)
ht <- GroupHeatmap(
  srt = seurat.res,
  features = c("Ifnb1","Ifnar1","Ifnar2"),
  group.by = c("group"),
  heatmap_palette = "YlOrRd",
  row_names_side = "left",
  add_dot = TRUE,
  add_reticle = TRUE
)
p1 <- ht$plot

obj <- subset(seurat.res, subset=celltype%in%c("CD8_Teff","CD8_Prol"))
p2 <- FeatureStatPlot(
  srt = obj, group.by = "group", plot_type = "violin",    
  stat.by = c("Ifnar1", "Il2rb"), palcolor = cols.group, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="CD8+ Teff cells", legend.position="none", ylab="Ifnar1 expression"
) 
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("Fig.6C.pdf", p, width = 6, height = 4)



##---------Fig. 6D (Fig. 6E bottom)
load("Pseudobulk_CD8Teff_group.RData")    
df <- reshape2::melt(exp.fpkm)
colnames(df) <- c("geneID", "sample", "fpkm")
exp.dt <- as.data.table(df)
exp.dt$group <- sub("_.*","",exp.dt$sample)
cols.group <- c("#3361A5","#B24745CC")
selected <- c("Ifnar1", "Ifnar2","Il2rb")
dt <- exp.dt[geneID%in%selected,]
dt$value <- log2(dt$fpkm+1)
dt$group <- factor(dt$group, levels = c("WT","KO"))
ps <- list()
for(i in selected){
  df_i <- dt[dt$geneID==i,]
  pi <- ggplot(df_i, aes(x=group, y=value, fill=group)) +
    geom_violin(trim=FALSE,color="white") + 
    geom_boxplot(width=0.1,position=position_dodge(0.1),outlier.shape = NA)+ 
    scale_fill_manual(values = cols.group)+ 
    stat_compare_means(method="t.test",label="p.signif")+
    xlab("") + ylab(paste0(i," (log2FPKM)")) + theme_bw() + 
    theme(legend.position="none",panel.background = element_rect(colour = "black"))
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], ps[[3]], ncol = 3, nrow = 1)   
ggsave("Fig.6D.pdf", p, width = 5, height = 4)




##---------Fig. 6F
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")
p <- SpatialFeaturePlot(seurat.Spatial, features = c("Il2rb"),
                        pt.size.factor = 0.5,      
                        image.alpha = 0.5,         
                        alpha = c(0, 1),         
                        images = 'slice1.polygons') 
ggsave("Fig.6F.png", p,  width = 5, height = 5, dpi = 300)





##---------Fig. 6G
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")
obj <- subset(seurat.res, subset=Major%in%"Tcells")
load("Pseudobulk_Spatial_Tcells.Polygons.RData")   #T cells
p1 <- FeatureStatPlot(
  srt = obj, group.by = "sample", plot_type = "violin",    
  stat.by = "Il2rb", palcolor = c("#eeca40","#fd763f","#23bac5"), 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="expression"
) 

cols.group <- c("#eeca40","#fd763f","#23bac5")
dt <- exp.dt[geneID%in%"Il2rb",]
dt$value <- log2(dt$fpkm+1)
dt$group <- factor(dt$group, levels = c("WT + IgG", "KO + IgG", "KO + PD-L1"))
my_comparisons <- list(c("WT + IgG", "KO + IgG"), c("WT + IgG","KO + PD-L1"))
p2 <- ggplot(dt, aes(x=group, y=value, fill=group)) +
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.1, position=position_dodge(0.1), outlier.shape = NA)+ 
  scale_fill_manual(values = cols.group)+ 
  stat_compare_means(comparisons = my_comparisons, label="p.signif") +
  xlab("") + ylab("Il2rb expression") + theme_bw() + 
  theme(legend.position="none",panel.background = element_rect(colour = "black"))
p <- ggarrange(p1, p2, ncol = 1, nrow = 2)   
ggsave("Fig.6G.pdf", p, width = 4, height = 6)






##---------Fig. 6L
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")

cells_in_umap <- rownames(seurat.Spatial@meta.data)    
sample <- fread("sample_from_LoupeBrowser.csv")    ##cloupe.cloupe
colnames(sample) <- c("barcode", "Sample")
square_008um <- fread("tissue_positions_square_008um.csv")    #tissue_positions_square_008um.csv
bc_map <- fread("barcode_mappings.csv")    #barcode_mappings
bc_map <- unique(bc_map[,c("square_008um","cell_id")])
bc_map_2 <- bc_map[cell_id%in%cells_in_umap,]
tissue_positions <- merge(bc_map_2, square_008um, by.x="square_008um", by.y="barcode")
tissue_positions$square_008um <- NULL
colnames(tissue_positions)[1] <- "barcode"
fac <- factor(tissue_positions$barcode)
x <- tapply(tissue_positions$"pxl_row_in_fullres", fac, mean)
y <- tapply(tissue_positions$"pxl_col_in_fullres", fac, mean)
spatial.locs_mean <- data.table(barcode=names(x), pxl_row_in_fullres=x, pxl_col_in_fullres=y[names(x)])
tissue_positions <- merge(spatial.locs_mean, sample, by="barcode")
tissue_positions$Sampleid <- tissue_positions$"Sample"
tissue_positions$Newbarcode <- paste0(tissue_positions$barcode, "_", tissue_positions$Sampleid)
tissue_positions$barcode <- NULL
metadata <- data.table(barcode=rownames(seurat.Spatial@meta.data), seurat.Spatial@meta.data)
metadata$Newbarcode <- paste0(metadata$barcode, "_", metadata$"sample")
tissue_posi_normalized <- normalize_spatial(tissue_positions)    
posi <- merge(x = tissue_posi_normalized, y = metadata, by = "Newbarcode", all.y = TRUE)
distance_results <- calculate_nearest_distances(
  posi,
  reference_type = "DCs",   
  target_types = c("Fibro","Macro", "Neut","Mono","Tcells", "Tumor"),   
  x_col = "pxl_row_in_fullres",   # x 
  y_col = "pxl_col_in_fullres",   # y 
  id_col = "Newbarcode",          # ID 
  type_col = "Major"   
)

target_types <- c("Fibro","Macro", "Neut","Mono","Tcells", "Tumor")   
dis_all <- data.table()
for(i in c("WT + IgG", "KO + IgG")){
  dis <- distance_results[grep(i, distance_results$Newbarcode, fixed = TRUE), target_types]
  dis_mean <- apply(dis[,target_types], 2, mean)
  dis.dt <- data.table(celltype=names(dis_mean), dis_mean)
  dis.dt$"sample" <- i
  dis_all <- rbind(dis_all, dis.dt)
}

celltype <- c("Tumor", "Tcells")   
ps <- list()
for(i in celltype){
  df <- dis_all[celltype%in%i,]
  df$sample <- factor(df$sample, levels = c("WT + IgG", "KO + IgG"))
  pi <- ggplot(df, aes(x = sample, y = dis_mean)) +
    geom_segment(aes(x = sample, xend = sample, y = 0, yend = dis_mean - 0.5), color = "gray") +
    geom_point(aes(x = sample, y = dis_mean, color = sample), size = 5) +
    geom_point(aes(x = sample, y = dis_mean, color = sample), shape = 1, size = 6) +
    geom_text(aes(label = round(dis_mean,1), x = sample, y = dis_mean + 1), hjust = -0.2, size = 3) + 
    scale_color_manual(values = c("#eeca40","#fd763f","#23bac5")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2*max(df$dis_mean))) +
    coord_flip() + ylab("Distance to the nearest DC cell (mean)") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.7),
          axis.ticks = element_line(size = 0.7),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 11, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"), 
          axis.title.y = element_blank(),  
          legend.position = "none")
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], labels=celltype, ncol = 2, nrow = 1)   
ggsave("Fig.6L.pdf", p, width = 6, height = 4)




##---------Fig. S12J
samples <- c("WT + IgG", "KO + IgG")
ps <- list()
for(i in samples){
  dis <- distance_results[grep(i, distance_results$Newbarcode, fixed = TRUE), ]
  
  #plot
  pi <- plot_radial_distance(
    dis,
    id_col = "Newbarcode",
    reference_type = "Tcells",  
    label_padding = 0.1,            
    show_labels = TRUE,             
    palette = "Dark2"               
  ) + scale_fill_manual(values = cols.Major)
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], labels=samples, ncol = 2, nrow = 1)  
ggsave("Fig.S12J.pdf", p, width = 10, height = 6)




##---------Fig. 6M
load("gKO_Il2rb_Lymphocytes.RData")
res <- gKO_res$"Il2rb"
selected.genes <- c("Runx1","Runx3","Ifngr1","Ifnar1","Mki67","Nfkb1","Itga1","Stat3","Stat4","Itgav","Ccr7")
labeled_genes <- res[gene%in%selected.genes,]
p <- ggplot(labeled_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') + coord_flip() + 
  labs(title="Representative regulated genes", x="", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ggsave("Fig.6M.pdf", p, width=4, height=3)








