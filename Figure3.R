library(Seurat)
library(data.table)
library(ggsci)
library(scales)
library(ggpubr)
library(scop)
library(pheatmap)
library(TCMNP)
library(ggplot2)
library(stringr)
library(parallel)
library(patchwork)
library(scales)
library(RColorBrewer)
library(scRNAtoolVis)  
library(readxl)       
library(corrplot)      
library(clusterProfiler)
library(org.Mm.eg.db)  
library(DOSE)
library(ggraph)




seurat.res <- readRDS("seurat.NSCLC_Brip1.rds")
load("cols.Major.RData")   #cols.Major
load("celltype_Immunecells.RData")    #cols.celltype
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")   #cols.group


# Fig. 3D
p <- DimPlot(seurat.res, reduction = "umap", group.by = "celltype", cols=cols.celltype, raster=TRUE,
             label = T,label.box = F, pt.size = 2)
ggsave("Fig.3D.pdf", p, width = 8, height = 5)


# Fig. 3E
p <- CellStatPlot(seurat.res,  stat.by = "celltype", group.by = "group", plot_type = "trend", palcolor = cols.celltype)
ggsave("Fig.3E.pdf", p, width = 8, height = 6)


# Fig. 3F
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")

p <- SpatialDimPlot(seurat.Spatial, group.by = 'Major', pt.size.factor = 0.5, 
                    images = 'slice1.polygons') + scale_fill_manual(values = cols.Major)  + NoLegend()
ggsave("Fig.3F.pdf", p,  width = 8, height = 8)


# Fig. 3G
p1 <- CellStatPlot(seurat.Spatial, stat.by = "Major", group.by = "Major", stat_type = "count", palcolor = cols.Major, position = "dodge", label = TRUE) + NoLegend()
p2 <- CellStatPlot(seurat.Spatial,  stat.by = "Major", group.by = "sample", plot_type = "trend", palcolor = cols.Major)
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)   
ggsave("Fig.3G.pdf", p, width = 6, height = 5)


# Fig. 3I
p1 <- FeatureStatPlot(
  srt = seurat.res, group.by = "group", plot_type = "violin",    
  stat.by = c("Cd274"), palcolor = cols.group, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
)

p2 <- FeatureStatPlot(
  srt = seurat.res, group.by = "group", plot_type = "violin",    
  stat.by = c("Pdcd1"), palcolor = cols.group, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
)
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)   
ggsave('Fig.3I.pdf', p, width = 4, height = 5)


# Fig. 3J
p <- SpatialFeaturePlot(seurat.Spatial, features = c("Cd274"),
                        pt.size.factor = 0.5,      
                        image.alpha = 0.5,         
                        alpha = c(0, 1),         
                        images = 'slice1.polygons') 
ggsave("Fig.3J.png", p,  width = 5, height = 5, dpi = 300)


# Fig. 3K (Fig. S3O)
cols.group_spatial <- c("#eeca40","#fd763f","#23bac5")
selected <- c("Cd274","Pdcd1","Cxcl10","Ccl5")
dt <- exp.dt[geneID%in%selected,]
dt$value <- log2(dt$fpkm+1)
dt$group <- factor(dt$group, levels = c("WT + IgG", "KO + IgG", "KO + PD-L1"))

ps <- list()
for(i in selected){
  df_i <- dt[dt$geneID==i,]
  df_i <- df_i[fpkm > 0,]   
  my_comparisons <- list(c("WT + IgG", "KO + IgG"), c("WT + IgG","KO + PD-L1"))
  pi <- ggplot(df_i, aes(x=group, y=value, fill=group)) +
    geom_violin(trim=FALSE,color="white") + 
    geom_boxplot(width=0.1,position=position_dodge(0.1),outlier.shape = NA)+ 
    scale_fill_manual(values = cols.group_spatial)+ 
    stat_compare_means(comparisons = my_comparisons, label="p.signif") +
    xlab("") + ylab(paste0(i," Expression")) + theme_bw() + 
    theme(legend.position="none",panel.background = element_rect(colour = "black"))
  
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 2, nrow = 2)   
ggsave("Fig.3K.pdf", p, width = 4, height = 4)


# Fig. 3L
dt <- data.table(cellName=rownames(seurat.res@meta.data), seurat.res@meta.data)
celltype_sorted <- sort(unique(seurat.res$celltype))
normalized_data_matrix <- GetAssayData(seurat.res, assay = "RNA", layer = "data")
num_filtered <- c()
for(i in celltype_sorted){
  mat <- normalized_data_matrix[,dt[celltype==i, cellName]]
  gene_expression_sums <- Matrix::rowSums(mat)
  threshold_sum <- ncol(mat) * 0.03
  genes_over_threshold <- gene_expression_sums > threshold_sum
  genes_expressed_in_cells <- Matrix::rowSums(mat > 0) > (ncol(mat) * 0.1)
  filtered_genes <- genes_over_threshold & genes_expressed_in_cells
  num_filtered <- c(num_filtered, sum(filtered_genes))
}
names(num_filtered) <- celltype_sorted

total_counts_df <- data.table(celltype=celltype_sorted, Genes_Tested=num_filtered, Freq=num_filtered)
total_counts_df$celltype <- factor(total_counts_df$celltype,levels =celltype_sorted )

p1 <- ggplot(total_counts_df, aes(x=celltype, y=Freq, fill = celltype)) +
  geom_bar(stat="identity",width=0.8)+theme_bw(base_size=6)+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ 
  scale_fill_manual(values=color_vector)+ylab('Genes\nTested')+
  theme(plot.margin=grid::unit(c(0,0,-2.5,0), "mm"))


up <- table(dt_deg[avg_log2FC>0,celltype])
dt_up <- data.table(celltype=names(up), count=as.numeric(up), Expression=1)   
down <- table(dt_deg[avg_log2FC<0,celltype])
dt_down <- data.table(celltype=names(down), count=as.numeric(down)*-1, Expression=0.7)   
deg_counts_df <- rbind(dt_up, dt_down)
deg_counts_df$celltype <- factor(deg_counts_df$celltype,levels = celltype_sorted)

p2 <- ggplot(deg_counts_df, aes(x = celltype, y = count, fill = celltype)) +
  geom_bar(stat = "identity",width=0.8,aes(alpha=factor(Expression)) ) + 
  scale_fill_manual(values=color_vector)+
  labs(y = "DEG\nCounts", x = "Cell Type") +
  scale_alpha_manual(values = c(0.7, 1)) +
  theme_bw(base_size=6)+
  theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.margin=grid::unit(c(0,0,0,-.25), "mm"))+
  scale_x_discrete(labels = function(x) gsub(" cell", "", x))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


res_all2 <- rbind(res_mye, res_lym)
res_all_4groups <- rbind(res_all, res_all2)

deg_counts_df_total <- data.table()
for(i in c("KO-WT&IgG", "KO-WT&PDL1", "PDL1-IgG&WT", "PDL1-IgG&KO")){
  dt_deg <- res_all_4groups[group==i & p_val_adj<=0.01 & abs(avg_log2FC)>=0.58,] 
  deg <- table(dt_deg$celltype)
  counts_df <- data.table(celltype=names(deg), Contrast=i, total_count_clipped=as.numeric(deg))   
  deg_counts_df_total <- rbind(deg_counts_df_total, counts_df)
}
deg_counts_df_total$celltype <- factor(deg_counts_df_total$celltype,levels =celltype_sorted )
deg_counts_df_total$Contrast <- factor(deg_counts_df_total$Contrast,levels =c("KO-WT&IgG", "KO-WT&PDL1", "PDL1-IgG&WT", "PDL1-IgG&KO") )


max_deg_number <- max(deg_counts_df_total$total_count_clipped)
p3 <- ggplot(deg_counts_df_total, aes(celltype, Contrast, fill = total_count_clipped)) +
  geom_tile() +
  theme_bw(base_size = 6) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.margin = grid::unit(c(0, 0, 0, -0.25), "mm")) +
  scale_x_discrete(labels = function(x) gsub(" cell", "", x)) +
  scale_fill_gradientn(
    colours = c("white", brewer.pal(9, "Blues")[0:9]),
    values = scales::rescale(c(0, seq(1, max_deg_number, 
                                      length.out = 9))),
    limits = c(0, max_deg_number)) +theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
p <- p1 / p2 / p3 + plot_layout(heights = c(1.5,3,1.3))
ggsave("Fig.3L.pdf", plot = p, width = 3, height = 3)




# Fig. S3A
p1 <- DimPlot(seurat.res, reduction = "umap", group.by = "seurat_clusters", raster=TRUE,  
              cols=cols, label = T, label.box = F)

##sample
p2 <- DimPlot(seurat.res, reduction = "umap", group.by = "sample", raster=TRUE, 
             cols=cols, label = F,label.box = F)
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("Fig.S3A.pdf", p, width = 14, height = 6)



# Fig. S3B
seurat.res[["Major"]] <- factor(seurat.res@meta.data$Major, 
                                levels = c("Macrophages","DC cells","Monocytes","Neutrophils","T cells","B cells","NK cells"))
p <- DimPlot(seurat.res, reduction = "umap", group.by = "Major", cols=cols.Major, raster=TRUE, order = T,
             label = T,label.box = F, pt.size = 1)
ggsave("Fig.S3B.pdf", p, width = 6, height = 5)



# Fig. S3C
p <- CellStatPlot(seurat.res,  stat.by = "Major", group.by = "group", plot_type = "trend", palcolor = cols.Major)
ggsave("Fig.S3C.pdf", p, width = 8, height = 5)




# Fig. S3E
final <- list(
  B=c("Cd79a","Cd19","Ighm","Ms4a1"), # B lymphocytes
  DC=c("Itgax","Cd83","H2-DMb1", "H2-DMa"), #DC
  Macro=c("Cd68","Lgmn","Ctsb","Trem2"), # Macrophage
  cMon=c("Cd14","Vcan","Clec4e","Fcn1"), # Monocyte 
  Neut=c("Ncf2","S100a8","S100a9","Fpr1","G0s2"), #Neutrophil
  NK=c("Nkg7","Gnly","Klrd1"), # Natural Killer cells
  T=c("Cd3d","Cd3e","Cd3g","Trac")  #T lymphocytes
)

p <- DotPlot(seurat.res, features = final, assay='RNA', group.by = "Major") + 
  RotatedAxis() + scale_color_gradientn(colours = viridis::viridis(20))
ggsave("Fig.S3E.pdf", p, width =10, height = 3)




# Fig. S3H
fts <- c("Cd86","Arg1","C1qc","Vegfa","Clec9a","Ccr7","G0s2","Gzmb","Mki67")
p <- FeatureDimPlot(srt = seurat.res, features = fts, 
                    reduction = "UMAP", raster=T, pt.size = 2, theme_use = "theme_blank", ncol = 3)
ggsave("Fig.S3H.pdf", p, width = 6, height = 6)




# Fig. S4A
seurat.res[["celltype"]] <- factor(seurat.res@meta.data$celltype, levels = sort(names(table(seurat.res@meta.data$celltype))))
p <- FeatureStatPlot(
  srt = seurat.res, group.by = "celltype", plot_type = "box",    
  stat.by = fts, add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
) 
ggsave('Fig.S4A.pdf', p, width = 20, height = 17)




# Fig. S3I
obj_WT <- subset(seurat.res, subset=grepl("WT", group))
obj_KO <- subset(seurat.res, subset=grepl("KO", group))
p1 <- DimPlot(obj_WT, reduction = "umap", group.by = "celltype", cols=cols.celltype, raster=T,
              label = F,label.box = F, pt.size = 1)
p2 <- DimPlot(obj_KO, reduction = "umap", group.by = "celltype", cols=cols.celltype, raster=T,
              label = F,label.box = F, pt.size = 1)
p <- ggarrange(p1, p2,ncol = 2, nrow = 1)   
ggsave("Fig.S3I.pdf", p, width = 16, height = 5)





# Fig. S3J
dt <- data.table(cellName=rownames(seurat.res@meta.data), seurat.res@meta.data)
dt.lst <- split(dt, dt$sample)
freq.mat <- matrix(0, nrow=length(unique(dt$celltype)), ncol=length(unique(dt$sample)))
rownames(freq.mat) <- sort(unique(dt$celltype))
colnames(freq.mat) <- sort(unique(dt$sample))
for(sample in names(dt.lst)){
  x <- dt.lst[[sample]]
  tb <- table(x$"celltype")
  f <- tb/sum(tb)
  freq.mat[names(f), sample] <- f
}
y <- unique(dt[,c("sample","group"), with = FALSE])
anno_col <- as.data.frame(y)
rownames(anno_col) <- anno_col$sample
anno_col$sample <- NULL
freq.mat <- na.omit(freq.mat[, rownames(anno_col)])

annColors <- list()
color <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
names(color) <- c("WT&IgG","KO&IgG", "WT&PDL1", "KO&PDL1")  
annColors[['group']] <- color
cols_h <- c("royalblue1", "royalblue2", "royalblue3", "royalblue4", "black", "yellow4","yellow3","yellow2", "yellow1")
pdf("Fig.S3J.pdf", width= 5, height=4)
pheatmap(freq.mat, cluster_cols = T, cluster_rows = T, 
         annotation_col=anno_col, 
         annotation_colors = annColors,
         show_rownames=T,show_colnames=T, scale='row',  
         treeheight_row=0, treeheight_col=0, 
         color = colorRampPalette(cols_h)(50)
)
dev.off()







# Fig. S3K
dt <- data.table(cellName=rownames(seurat.res@meta.data), seurat.res@meta.data)
celltype_sorted <- sort(unique(dt$celltype))
color_vector <- cols.celltype[celltype_sorted]

initial.num <- rep(0,length(celltype_sorted))
names(initial.num) <- celltype_sorted

tr_sample <- split(dt, dt$sample)
tr_num <- sapply(tr_sample, function(x){
  tb <- table(x$celltype)
  res <- initial.num
  res[names(tb)] <- tb
  return(res)
})
freq <- apply(tr_num,2,function(x){x/sum(x)})

mean_1 <- apply(freq[,unique(dt[group%in%c("WT&IgG"),sample])],1,mean)
mean_2 <- apply(freq[,unique(dt[group%in%c("KO&IgG"),sample])],1,mean)
mean_3 <- apply(freq[,unique(dt[group%in%c("WT&PDL1"),sample])],1,mean)
mean_4 <- apply(freq[,unique(dt[group%in%c("KO&PDL1"),sample])],1,mean)


fc1 <- mean_2[names(initial.num)]/mean_1[names(initial.num)]   
fc2 <- mean_4[names(initial.num)]/mean_3[names(initial.num)]   
res_all <- rbind(data.table(celltype=names(fc1), FC=fc1, Expression=1),  
                 data.table(celltype=names(fc2), FC=fc2*-1, Expression=0.8))
res_all$celltype <- factor(res_all$celltype,levels = celltype_sorted)

p <- ggplot(res_all, aes(x = celltype, y = FC, fill = celltype)) +
  geom_bar(stat = "identity",width=0.8,aes(alpha=factor(Expression)) ) + 
  scale_fill_manual(values=color_vector)+
  labs(y = "FC\n(Brip1-KO vs. -WT)", x = "") +
  scale_alpha_manual(values = c(0.8, 1)) +
  theme_bw() + theme(legend.position="none")+ 
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.S3K.pdf", plot = p, width = 6, height = 3)





# Fig. S4B
load("sce.markers_group_Myeloid.RData")
res$cluster <- res$celltype

data1 <- res[group=="KO-WT&IgG",]     
data_clean <- data1 %>% 
  filter(abs(avg_log2FC) > 0.58)

p <- jjVolcano(
  diffData = data_clean,                
  log2FC.cutoff = 0.58,               
  tile.col = cols.celltype,
  size = 3.5,                           
  fontface = 'italic',                  
  topGeneN = 3,                         
  cluster.order = sort(unique(data_clean$cluster)),  
  legend.position = c(0.08, 0.02),      
  aesCol = c('#5A8FC8','#E85A4F')      
)
ggsave('Fig.S4B.pdf', p, width = 8, height = 5)




# Fig. S5C
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")
infercnv_obj <- readRDS("run.final.infercnv_obj")
expr <- infercnv_obj@expr.data  

dt <- data.table(cellName=rownames(seurat.Spatial@meta.data), seurat.Spatial@meta.data)
dt <- dt[cellName%in%colnames(expr),]
dat <- sapply(sort(unique((dt$celltype))),function(i){
  cells <- dt[celltype==i, cellName]
  tmp <- apply(expr[,cells],1,mean)
  return(tmp)
})
colnames(dat) <- sort(unique((dt$celltype)))
cols5 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))

dat_2 <- t(dat)
genes_pos <- fread("GRCm39_genes_pos.txt", header = F)
colnames(genes_pos) <- c("gene","chr","start","end")
genes_pos <- genes_pos[gene%in%colnames(dat_2),]
genes_pos$"num" <- as.numeric(sub("chr","",genes_pos$chr))
genes_pos <- genes_pos[order(num),]
annotation_col <- data.frame(chr=factor(genes_pos$"chr", levels = unique(genes_pos$"chr")))
rownames(annotation_col) <- genes_pos$gene

ann_colors <- list(chr = rep(c("#D2D2D2", "#A8A8A8"),10))
names(ann_colors$chr) <- unique(annotation_col$chr)
dat_2 <- dat_2[,rownames(annotation_col)]
pdf("Fig.S5C.pdf", width= 8, height=3)
pheatmap(dat_2, cluster_cols = F, cluster_rows = F, 
         show_rownames=T,show_colnames=F, scale='column', 
         treeheight_row=0, treeheight_col=0,
         annotation_col=annotation_col, annotation_colors=ann_colors,
         color = colorRampPalette(cols5)(50)
)
dev.off()





# Fig. S5E
entrezID2geneSymbol <- function(entrezID){
  entrezID <-  strsplit(entrezID, split="/")[[1]]
  SYMBOL <- bitr(entrezID,"ENTREZID","SYMBOL",org.Mm.eg.db)[,"SYMBOL"]
  SYMBOL_chr <- paste0(SYMBOL,"",collapse ="/")
  return(SYMBOL_chr)
}

load("sce.markers_group.incelltype_Myeloid.RData")
dt_all <- res[p_val_adj<=0.01 & abs(avg_log2FC)>=0.58,]
dt <- dt_all[group=="KO-WT&IgG" & avg_log2FC>=0.58,]    ## KO vs. WT (IgG)
ids <- bitr(dt$gene,'SYMBOL','ENTREZID','org.Mm.eg.db') 
merged <- merge(dt, ids, by.x='gene', by.y='SYMBOL')
degs.list <- split(merged$ENTREZID, merged$celltype)
degs.list <- lapply(degs.list, function(x){as.numeric(unique(x))})

bp <- compareCluster(degs.list,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)
cols <- c("#118ab2", "#fdffb6", "#e63946")   
p <- dotplot(bp, label_format=100, font.size=14, showCategory = 14) + scale_fill_gradientn(colours=cols)+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1,size = 9))
ggsave("Fig. S5E.pdf", p, width= 9, height=10)






# Fig. S5F
load("DEG.LLC1.KO_WT.RData")  
res$obj <- rownames(res)
dt1 <- na.omit(as.data.table(res))
dt <- dt1[(pvalue <= 0.01 & abs(log2FoldChange) >= 0.58),]  
up <- head(dt$obj,200)   #Sorted in descending order based on log2FoldChange
eg <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
BP <- enrichGO(
  gene = eg$ENTREZID,
  "org.Mm.eg.db",  
  ont = "BP",  ##Biological Process
  pvalueCutoff = 0.01,
  readable = TRUE
)
p <- pathway_ccplot(BP, top = 15, root = "LLC1 Brip1-KO")    
ggsave("Fig.S5F.pdf", p, width= 6, height=6)




# Fig. S6C
load("DESeq2_DEG.LLC1-KO_mouse.KO_WT.IgG.RData")   #KO vs. WT (in IgG)
res$obj <- rownames(res)
dt1 <- na.omit(as.data.table(res))
dt <- dt1[(pvalue <= 0.01 & abs(log2FoldChange) >= 1),]  
up <- head(dt$obj,200)   #Sorted in descending order based on log2FoldChange
eg <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
BP <- enrichGO(
  gene = eg$ENTREZID,
  "org.Mm.eg.db",  
  ont = "BP",  ##Biological Process
  pvalueCutoff = 0.01,
  readable = TRUE
)
p <- dot_sankey(BP, top = 15, dot.y = 0.01, dot.height = 0.95)  
ggsave("Fig.S6C.pdf", p, width= 6, height=6)

















