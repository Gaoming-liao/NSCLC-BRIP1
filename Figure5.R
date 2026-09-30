library(Seurat)
library(data.table)
library(ggsci)
library(scales)
library(ggpubr)
library(scop)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(UCell)
library(magrittr)
library(enrichR)
library(GeneOverlap)
library(fgsea)
library(clusterProfiler)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Mmusculus.v79)    
library(BSgenome.Mmusculus.UCSC.mm10)    
library(GenomicRanges)
library(xgboost)
library(JASPAR2024)
library(RSQLite)
library(EnsDb.Mmusculus.v79) 
library(ggplot2)
library(ggrepel)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)
options(future.globals.maxSize = 40 * 1024^3)  



seurat.res <- readRDS("seurat.NSCLC_Brip1.rds")
seurat_obj <- subset(x = seurat.res, subset = Major%in%c("T cells","NK cells","B cells"))

load("celltype_Immunecells.RData")    #cols.celltype
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")   #cols.group


##---------Fig. 5A
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", cols=cols.celltype, raster=T,
             label = T,label.box = F, pt.size = 2, alpha = 0.8)
ggsave("Fig.5A.pdf", p, width = 5, height = 4)




##---------Fig. S10B-S10E
seurat_obj[['group2']] <- sub("&.*","",seurat_obj@meta.data$group)   #Brip1 KO and WT
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select ="fraction", fraction = 0.05, wgcna_name = "Brip1")
seurat_obj <- MetacellsByGroups( seurat_obj = seurat_obj,
                                 group.by = c("sample","batch","Major","celltype","group","group2"),   
                                 reduction ='harmony', # select the dimensionality reduction to perform KNN on
                                 k = 25,  #nearest-neighbors parameter
                                 max_shared = 10, # maximum number of shared cells between two metacells
                                 ident.group ='group2' # set the Idents of the metacell seurat object
)
# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
# Co-expression network analysis
obj <- SetDatExpr(seurat_obj, group_name = "KO", group.by='group2', assay = 'RNA', layer = 'data')


##---------Fig. S10B
obj <- TestSoftPowers(obj, networkType = 'signed')
plot_list <- PlotSoftPowers(obj)
pdf("Fig.S10B.pdf", width=8, height=7)
wrap_plots(plot_list, ncol=2)
dev.off()


##---------Fig. S10C
power_table <- GetPowerTable(obj)
obj <- ConstructNetwork(obj, tom_name = 'KO', overwrite_tom =T)
pdf("Fig.S10C.pdf", width=7, height=4)
PlotDendrogram(obj, main='Brip1-KO hdWGCNA Dendrogram')   
dev.off()

##-- Module Eigengenes and Connectivity
obj <- ScaleData(obj, features=VariableFeatures(obj))         # need to run ScaleData first or else harmony throws an error:
obj <- ModuleEigengenes(obj, group.by.vars="batch")          # compute all MEs in the full single-cell dataset. group.by.vars (groups to harmonize by)
hMEs <- GetMEs(obj)          # harmonized module eigengenes:
MEs <- GetMEs(obj, harmonized=FALSE)          # module eigengenes:

# compute eigengene-based connectivity (kME):  
obj <- ModuleConnectivity(obj, group.by = 'group2', group_name = 'KO')


##---------Fig. S10D
# rename the modules
obj <- ResetModuleNames(obj, new_name = "KO-M")
# plot genes ranked by kME for each module
p <- PlotKMEs(obj, ncol=7)
ggsave("Fig.S10D.pdf", p, width=18, height=8)


# get the module assignment table:
modules <- GetModules(obj) %>% subset(module != 'grey')
# get hub genes
hub_df <- GetHubGenes(obj, n_hubs = 10)
head(hub_df)
saveRDS(obj, file='hdWGCNA_obj_Lymphocytes.rds')



##---------Fig. 5B
pdf("Fig.5B.pdf", width=5, height=5)
HubGeneNetworkPlot(
  obj,
  n_hubs = 6, n_other=6,  
  edge_prop = 0.2,
  mods = 'all'
)
dev.off()



##---------Fig. S10E
# compute gene scoring for the top 25 hub genes by kME for each module with UCell method
obj <- ModuleExprScore(obj, n_genes = 25, method='UCell')   
plot_list <- ModuleFeaturePlot(obj, features='hMEs', raster = T, point_size = 0.1, order=TRUE)
p <- wrap_plots(plot_list, ncol=4)
ggsave("Fig.S10E.pdf", p, width=10, height=6)



##---------Fig. S10F (Fig. S10G)
gmt <- read.gmt("m5.all.v2025.1.Mm.symbols.gmt")    #MSigDB GMT file
gmt$term <- Hmisc::capitalize(tolower(gsub("_"," ",sub("[A-Z]+_","",gmt$term))))  
pathways <- lapply(unique(gmt$term), function(x){
  res <- unique(gmt[gmt$term==x,"gene"])
  return(res)
})
names(pathways) <- unique(gmt$term)

# get the modules table and remove grey genes
modules <- GetModules(obj) %>% subset(module != 'grey')

# rank all genes in this kME (eigengene-based connectivity (kME))
cur_mod <- 'KO-M1'    
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME 
names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]    #They are ordered (kME) from low to high

# run fgsea to compute enrichments  
gsea_df <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500)
gsea_df <- gsea_df[pval <= 0.05,]

select.keys1 <- c("Antigen","Mhc class", "nf kappab","Interferon", "immune",
                  "chemotaxis","Chemokine","cytokine", "Cytokinesis", "Cell killing", 
                  "Lymphocyte","tumor necrosis factor"
)   #chemotaxis
select.keys2 <- c("Myeloid", "Macrophage", "Neutrophil","mononuclear","T cell","B cell","Natural killer","t helper"
)   #cells

index <- c()
for(i in c(select.keys1, select.keys2)){   #c(select.keys1,select.keys2)
  pos <- grep(i, gsea_df$pathway, ignore.case = T)
  index <- c(index, pos)
}
gsea_df_deal <- gsea_df[unique(index),]
gsea_df_deal <- gsea_df_deal[order(NES, decreasing = F),]   
top_pathways <- gsea_df_deal[1:25,pathway]
p <- plotGseaTable(pathways[top_pathways], ranks, gsea_df_deal, gseaParam=1, colwidths = c(10, 4, 1, 1.5, 1.5))
ggsave("Fig.S10F.pdf", p, width=10, height=6)





##---------Fig. 5C
##Transcription factor regulatory network analysis
seurat_obj <- readRDS('hdWGCNA_obj_Lymphocytes.rds')
JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan
seurat_obj <- MotifScan(seurat_obj, species_genome = 'mm10', pfm = pfm_core, EnsDb = EnsDb.Mmusculus.v79)

# Construct TF Regulatory Network
motif_df <- GetMotifs(seurat_obj)

# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj, wgcna_name = "Brip1")
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
seurat_obj <- SetWGCNAGenes(seurat_obj, genes_use, wgcna_name = "Brip1")
seurat_obj <- SetDatExpr(seurat_obj, group_name = "KO", group.by='group2', wgcna_name = "Brip1")

#Now we are ready to run ConstructTFNetwork. Since this function models each gene, the runtime will scale wth the number of genes selected from the 
#previous step and it will scale with the number of metacells / metaspots that are used for this analysis. 
# define model params:
model_params <- list(
  objective = 'reg:squarederror',
  max_depth = 1,
  eta = 0.1,
  nthread=16,
  alpha=0.5
)

# construct the TF network. This function results in a table showing information about each potential TF-gene pair
seurat_obj <- ConstructTFNetwork(seurat_obj, model_params,wgcna_name = "Brip1")    
results <- GetTFNetwork(seurat_obj, wgcna_name = "Brip1")

res <- as.data.table(results)
#save(res, file="TFnetwork.RData")
#saveRDS(seurat_obj, file='hdWGCNA_obj_TFnetwork.rds')


# Define TF Regulons
seurat_obj <- AssignTFRegulons(
  seurat_obj,
  strategy = "B",
  reg_thresh = 0.01,
  n_genes = 50
)

# Calculate regulon expression signatures
# positive regulons
seurat_obj <- RegulonScores(
  seurat_obj,
  target_type = 'positive',
  cor_thresh = 0.05,
  ncores=8
)
# negative regulons
seurat_obj <- RegulonScores(
  seurat_obj,
  target_type = 'negative',
  cor_thresh = -0.05,
  ncores=8
)

##Differential regulon analysis
group1_cells <- seurat_obj@meta.data %>% subset(group2 == "KO") %>% rownames
group2_cells <- seurat_obj@meta.data %>% subset(group2 != "KO") %>% rownames
# calculate differential regulons
dregs <- FindDifferentialRegulons(
  seurat_obj,
  barcodes1 = group1_cells,
  barcodes2 = group2_cells, 
  wgcna_name = "Brip1"
)
p <- PlotDifferentialRegulons(seurat_obj, dregs)
ggsave("Fig.5C.pdf", p, width=4, height=4)



##---------Fig. 5D (Fig. 5L, Fig. S11C-S11D)
obj_teff <- subset(seurat.res, subset=celltype%in%c("CD8_Teff","CD8_Prol"))
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
fts <- c("Runx1","Ifngr1","Gzma","Gzmb","Prf1","Gzmk","Ifng")
p1 <- FeatureStatPlot(
  srt = obj_teff, group.by = "group", plot_type = "violin",   
  stat.by = fts, palcolor = cols.group, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
) 
ggsave('Fig.5D.pdf', p1, width = 7, height = 8)



##---------Fig. 5G
normalized_data_matrix <- GetAssayData(obj_teff, assay = "RNA", layer = "data")   
mat <- normalized_data_matrix["Runx1",]
cells_high <- names(mat)[mat>=quantile(mat,2/3)]   #top 1/3
cells_low <- names(mat)[mat<=quantile(mat,1/3)]   #bottom 1/3
obj_high <- subset(obj_teff, cells = cells_high) 
obj_high[['Runx1']] <- "High"
obj_low <- subset(obj_teff, cells = cells_low) 
obj_low[['Runx1']] <- "Low"

##merge two types data objects
obj_Runx1 <- merge(obj_high, obj_low)
obj_Runx1 <- JoinLayers(obj_Runx1, overwrite = TRUE)  
dt <- data.table(cellName=rownames(obj_Runx1@meta.data), obj_Runx1@meta.data)
tb <- table(dt[,c("group","Runx1")])

obj_Runx1[["Runx1"]] <- factor(obj_Runx1@meta.data$Runx1, levels = c("High","Low"))
obj_Runx1[["group"]] <- factor(obj_Runx1@meta.data$group, levels = c("WT&IgG","KO&IgG", "WT&PDL1", "KO&PDL1"))

p <- CellStatPlot(obj_Runx1,  stat.by = "Runx1", group.by = "group", plot_type = "trend", palcolor = c("#AAD0AC","#BFA6C9"))
ggsave("Fig.5G.pdf", p, width = 3, height = 4.5)



##---------Fig. 5H
load("sce.markers.CD8Teff.Runx1_group.RData")
df <- as.data.frame(res)
rownames(df) <- df$gene
df$Significant <- ifelse(df$p_val_adj < 0.01 & abs(df$avg_log2FC) > 0.4,
                         ifelse(df$avg_log2FC >= 0.4, "Up", "Down"), "Ns")
dt <- res[(p_val_adj <= 0.01 & abs(avg_log2FC) >= 0.4),]
select.genes <- c("Cxcr4","Nr4a1","Rora","Tgfb1","Gzmk","Gzma","Stat5a","Stat5b","Jak2","Cd6","Ap3b1",
                  "Cd226","Zeb2", "Nfkb1","Traf3","Satb1","Itgb3","Mki67","Mxi1","Smad3","Eomes","Ccr2","Cxcr3",
                  "Trerf1","Ptprk","Lef1","Rbm47","Myb")
select.genes <- c(select.genes, tail(dt[,gene],16))
df_sel <- df[select.genes, ]
p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("#acd6ec","#d3d5d4","#f5a889")) +
  geom_text_repel( 
    data = df[select.genes,],
    aes(x = df_sel$avg_log2FC, y = -log10(df_sel$"p_val_adj"),label = rownames(df_sel)), 
    size = 3, box.padding = unit(0.2, "lines"), point.padding = unit(0.2, "lines")) + 
  theme_test(base_size = 16)+ 
  geom_vline(xintercept=c(-0.4,0.4),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 0,lty=4,col="black",lwd=0.8) +
  xlab("Average log2FC") + ylab("-log10 (Padj)") +ggtitle('Runx1 high vs. low')+
  theme(legend.position = "top") 
ggsave("Fig.5H.pdf", p, width=6,height=6)





##---------Fig. 5M
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")
p <- SpatialFeaturePlot(seurat.Spatial, features = c("Gzma"),
                        pt.size.factor = 0.5,      
                        image.alpha = 0.5,         
                        alpha = c(0, 1),         
                        images = 'slice1.polygons') 
ggsave("Fig.5M.png", p,  width = 5, height = 5, dpi = 300)


##---------Fig. 5N (Fig. S11F)
selected <- c("Gzmb","Prf1","Gzma","Ifng")
dt <- exp.dt[geneID%in%selected,]
dt$value <- log2(dt$fpkm+1)
dt$group <- factor(dt$group, levels = c("WT + IgG", "KO + IgG", "KO + PD-L1"))
ps <- list()
for(i in selected){
  df_i <- dt[dt$geneID==i,]
  my_comparisons <- list(c("WT + IgG", "KO + IgG"), c("WT + IgG","KO + PD-L1"))
  pi <- ggplot(df_i, aes(x=group, y=value, fill=group)) +
    geom_violin(trim=FALSE,color="white") + 
    geom_boxplot(width=0.1,position=position_dodge(0.1),outlier.shape = NA)+ 
    scale_fill_manual(values = cols.group)+ 
    stat_compare_means(comparisons = my_comparisons, label="p.signif") +
    xlab("") + ylab(paste0(i," Expression")) + theme_bw() + 
    theme(legend.position="none",panel.background = element_rect(colour = "black"))
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 2, nrow = 2)   
ggsave("Fig.5N.pdf", p, width = 4, height = 4)




##---------Fig. 5R
load("clinical.data_TCGA-NSCLC.RData")
clinical.data[Immune_Subtype=="C1", Immune_Subtype:="C1 (wound healing)"]
clinical.data[Immune_Subtype=="C2", Immune_Subtype:="C2 (IFN-g dominant)"]
clinical.data[Immune_Subtype=="C3", Immune_Subtype:="C3 (inflammatory)"]
clinical.data[Immune_Subtype=="C4", Immune_Subtype:="C4 (lymphocyte depleted)"]
clinical.data[Immune_Subtype=="C6", Immune_Subtype:="C6 (TGF-b dominant)"]

tb <- table(clinical.data[,c("BRIP1_status","Immune_Subtype")])
Ratio <- tb/apply(tb,1,sum)
cols <- RColorBrewer::brewer.pal(n = 10, name = "Paired")[1:5]
pdf("Fig.5R.pdf", width = 3,height = 6)
par(mar=c(8, 5, 4, 2))
x <- barplot(height=t(Ratio), ylim = c(0, 1), beside=FALSE, las=2, cex.lab=1.3,
             border = 0.1, col = cols, xlab=" ",ylab="Percentage")
abline(h=seq(0.2, 1, by=0.2),col="gray",lty=2)
legend("right", legend=colnames(Ratio), pch=16, col=cols)
dev.off()




##---------Fig. S11H
seurat.res <- readRDS("seurat.NSCLC_Brip1.rds")
p <- FeatureStatPlot(
  srt = seurat.res, group.by = "celltype2", plot_type = "box",    
  stat.by = c("Ifng","Ifngr1","Ifngr2", "Ifnar2", "Il15","Il15ra"), palcolor = cols.celltype, 
  add_box = TRUE, multiplegroup_comparisons=F, 
  subtitle="", legend.position="none", ylab="Expression"
) 
ggsave('Fig.S11H.pdf', p, width = 15, height = 7)







