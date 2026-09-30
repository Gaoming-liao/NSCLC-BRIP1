library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(data.table)
library(GSVA)
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(scTenifoldKnk)
library(ActivePathways)
library(ggrepel)
library(tidygraph)
library(ggraph)
library(scop)
library(grid)
options(stringsAsFactors = FALSE)

cellchat <- readRDS("cellchat.rds")
seurat.res <- readRDS("seurat.NSCLC_Brip1.rds")

load("celltype_Immunecells.RData")    #cols.celltype
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")   #cols.group


##---------Fig. 4A (Fig. S7A)
h1 <- netVisual_heatmap(cellchat, title.name="KO vs. WT (IgG)", comparison=c(1,2), color.use=cell_colors, font.size = 7) 
h2 <- netVisual_heatmap(cellchat, title.name="KO vs. WT (PD-L1)", comparison=c(3,4), color.use=cell_colors, font.size = 7)
pdf("Fig.4A.pdf", width = 7, height = 4)
print(h1+h2)
dev.off()



##---------Fig.4B (Fig. S7B)
ct <- levels(cellchat@idents$joint)
sources.use <- c("M1-Macro","DC1","DC2","DC3","pDCs","B")   
targets.use <- c("CD8_Teff","CD8_Prol","CD4_Tfh", "CD4_Prol")  
p1 <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, 
                       comparison = c(1,2), color.text=cols_group[1:2], 
                       max.dataset = 2, title.name = "Increased signaling in Brip1-KO (IgG)", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, 
                       comparison = c(3,4), color.text=cols_group[3:4], 
                       max.dataset = 4, title.name = "Increased signaling in Brip1-KO (PD1)", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Fig. 4B.pdf", pc, width = 18, height = 12)



##---------Fig. 4C (Fig. S7C, Fig. S7H)
ps <- list()
for (i in 1:length(cco.list)) {
  ps[[i]] <- netVisual_heatmap(cco.list[[i]], measure="weight", signaling = "TNF", color.heatmap = "Reds", color.use=cell_colors,
                               title.name = paste0("TNF signaling, (", names(cco.list)[i], ")"), font.size = 10, width = 4, height = 5)
}
pdf("Fig.4C.pdf", width = 15, height = 5)
ComplexHeatmap::draw(ps[[1]]+ps[[2]]+ps[[3]]+ps[[4]], ht_gap = unit(0.5, "cm"))
dev.off()




##---------Fig. 4D
dt <- data.table(cellName=rownames(seurat.res@meta.data), seurat.res@meta.data)
raw_count_mat <- GetAssayData(seurat.res, assay = "RNA", layer = "counts")
celltypes <- levels(cellchat@idents$joint)
cell_cols <- cols.celltype[celltypes]
names(cell_cols) <- sub("c[12]\\.\\d+[-|_]","",names(cell_cols))
cell_levels <- sub("c[12]\\.\\d+[-|_]","",celltypes)
reg_levels <- c("Down", "Up", "None")
reg_cols <- c(
  "Down" = "#3A45A4",
  "Up"   = "#8E2D2D",
  "None" = "#9A9A9A"
)
bg_col <- "#F3F3F3"

left <- data.table(id=cell_levels, label=cell_levels, group=cell_levels, 
                   side="left", node_class="circle", x=0.17)
left$y <- seq(1,0, -0.05)

###ligand and receptor
ligand <- c("Tnf","Tnfsf9","Tnfsf4")
receptor <- c("Tnfrsf1b","Tnfrsf9","Tnfrsf4")
mid <- data.table(id=c(ligand, receptor), label=c(ligand, receptor), group=NA, side="mid", node_class="square")
mid$x <- c(0.28, 0.28, 0.28, 0.39, 0.39, 0.39)
mid$y <- c(0.2, 0.5, 0.8, 0.2, 0.5, 0.8)

right <- data.table(id=cell_levels, label=cell_levels, group=cell_levels, 
                    side="right", node_class="circle", x=0.5)
right$y <- seq(1,0, -0.05)

nodes <- rbind(left, rbind(mid, right))
nodes <- nodes %>%
  mutate(
    group = factor(group, levels = cell_levels),
    side = factor(side, levels = c("left", "mid", "right")),
    node_class = factor(node_class, levels = c("circle", "square"))
  )

res_diff <- rbind(res_mye, res_lym) %>% filter(gene %in% c(ligand,receptor))
res_up <- res_diff[avg_log2FC>0,]
res_down <- res_diff[avg_log2FC<0,]

mat <- sapply(celltypes, function(i){
  cells <- dt[celltype2==i, cellName]
  tmp <- apply(raw_count_mat[c(ligand,receptor), cells],1,mean)    #raw_count_mat, normalized_mat
  return(tmp)
})
mat.dt <- as.data.table(reshape2::melt(t(mat)))
colnames(mat.dt) <- c("celltype","gene","expr")

mat.dt$ID <- paste0(mat.dt$celltype, ":", mat.dt$gene)
mat.dt$reg <- "None"
mat.dt[ID%in%paste0(res_up$celltype, ":", res_up$gene), reg:="Up"]
mat.dt[ID%in%paste0(res_down$celltype, ":", res_down$gene), reg:="Down"]
mat.dt$celltype <- sub("c[12]\\.\\d+[-|_]","",mat.dt$celltype)
mat.dt$ID <- NULL

##-left_edges
left_edges <- mat.dt[gene%in%ligand,c("celltype","gene","reg","expr")]
colnames(left_edges) <- c("from","to","reg","expr")
left_edges <- left_edges %>% mutate(edge_class = "left")

##--mid_edges
mid_edges <- tribble(
  ~from,   ~to,      ~reg,   ~expr,
  "Tnf",  "Tnfrsf1b",  "None", 1,
  "Tnfsf9",  "Tnfrsf9",   "None", 1,
  "Tnfsf4",  "Tnfrsf4",  "None", 1
) %>%
  mutate(edge_class = "mid")

##--right_edges
right_edges <- mat.dt[gene%in%receptor,c("gene","celltype","reg","expr")]
colnames(right_edges) <- c("from","to","reg","expr")
right_edges <- right_edges %>% mutate(edge_class = "right")

edges <- bind_rows(left_edges, mid_edges, right_edges) %>%
  mutate(
    reg = factor(reg, levels = reg_levels),
    edge_class = factor(edge_class, levels = c("left", "mid", "right"))
  )
edges[expr>=3, expr:=3]

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
p <- ggraph(graph, layout = "manual", x = x, y = y) +
  
  # left edge
  geom_edge_link(
    aes(
      filter = edge_class == "left",
      edge_colour = reg,
      edge_width = 0.5*expr
    ),
    lineend = "round",
    alpha = 1
  ) +
  
  # mid dashed line
  geom_edge_link(
    aes(filter = edge_class == "mid"),
    colour = "black",
    linewidth = 0.7,
    linetype = "22",
    show.legend = FALSE
  ) +
  
  # right edge
  geom_edge_link(
    aes(
      filter = edge_class == "right",
      edge_colour = reg,
      edge_width = 0.5*expr
    ),
    lineend = "round",
    alpha = 1
  ) +
  
  # node_class
  geom_node_point(
    aes(
      filter = node_class == "circle",
      colour = group
    ),
    shape = 21,
    fill = "white",
    size = 4.4,
    stroke = 2.3
  ) +
  
  # mid
  geom_node_point(
    aes(filter = node_class == "square"),
    shape = 22,
    fill = "white",
    colour = "black",
    size = 3.5,
    stroke = 0.9,
    show.legend = FALSE
  ) +
  
  # left label
  geom_node_text(
    aes(
      filter = side == "left",
      label = label
    ),
    hjust = 1,
    nudge_x = -0.018,
    size = 4.8
  ) +
  
  # right label
  geom_node_text(
    aes(
      filter = side == "right",
      label = label
    ),
    hjust = 0,
    nudge_x = 0.020,
    size = 4.8
  ) +
  
  # ligand 
  geom_node_text(
    aes(
      filter = id %in% ligand,
      label = label
    ),
    fontface = "italic",
    hjust = 0,
    nudge_x = 0.005,
    nudge_y = 0.055,
    size = 5.5
  ) +
  
  # receptor
  geom_node_text(
    aes(
      filter = id %in% receptor,
      label = label
    ),
    fontface = "italic",
    hjust = 0,
    nudge_y = 0.055,
    size = 5.5
  ) +
  
  # title
  annotate(
    "text", x = 0.25, y = 1, label = "TNF/CD137/OX40 interactions",
    hjust = 0, vjust = 0.5, size = 6
  ) +
  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  
  scale_colour_manual(
    values = cell_cols,
    breaks = cell_levels,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  scale_edge_colour_manual(
    values = reg_cols,
    breaks = reg_levels,
    drop = FALSE,
    name = "Regulation with Brip1-KO"
  ) +
  
  scale_edge_width_continuous(   
    range = c(0.01, 3),
    breaks = 1:3,
    limits = c(0.1, 3),
    name = "Mean expression"
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        size = 5.4,
        stroke = 2.3
      )
    ),
    edge_colour = guide_legend(
      order = 2,
      override.aes = list(edge_width = 1.3)
    ),
    edge_width = guide_legend(
      order = 3,
      override.aes = list(edge_colour = "black")
    )
  ) +
  
  theme_void() +
  theme(
    plot.background = element_rect(fill = bg_col, colour = bg_col),
    panel.background = element_rect(fill = bg_col, colour = bg_col),
    legend.background = element_rect(fill = bg_col, colour = bg_col),
    legend.key = element_rect(fill = bg_col, colour = bg_col),
    plot.margin = margin(12, 20, 12, 12),
    
    legend.position = c(0.9,0.4),
    legend.box = "vertical",
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 13)
  )
ggsave("Fig.4D.pdf", p, width = 8, height = 5)





##---------Fig. 4E-4F
p1 <- plotGeneExpression(cellchat, signaling = c("TNF","CD137","OX40"), type = "dot", color.use=cell_colors, enriched.only = TRUE)
ggsave("Fig.4E.pdf", p1, width = 6, height = 3)
p2 <- plotGeneExpression(cellchat, signaling = c("TNF","CD137","OX40"), type = "dot", group.by = "group", color.use=cols_group, enriched.only = TRUE)
ggsave("Fig.4F.pdf", p2, width = 3, height = 3)



##---------Fig. 4G
seurat.Spatial <- readRDS("Spatial.Polygons.NSCLC_Brip1_celltype.rds")
p <- SpatialFeaturePlot(seurat.Spatial, features = c("Tnfrsf9"),
                        pt.size.factor = 0.3,      
                        image.alpha = 0.5,         
                        alpha = c(0, 1),         
                        images = 'slice1.polygons') 
ggsave("Fig.4G.png", p,  width = 5, height = 5, dpi = 300)






##---------Fig. 4M
load("dexp.data_OAK.RData")      #log2(TPM+1)
load("clinical.data_OAK.RData")
load("cGAS-STING.signature.RData")
clinical <- clinical.data[ACTARM=="MPDL3280A" & STUDYNAME=="OAK",]
con.sample <- intersect(clinical$SampleID, colnames(exp.data))
exp.data <- exp.data[,con.sample]
Param <- ssgseaParam(as.matrix(exp.data), gene.list) 
score <- gsva(Param)
score_CGAS.STING <- score["cGAS_STING",]
act <- data.table(SampleID=colnames(score_CGAS.STING), t(score_CGAS.STING))

act$"cgas.sting.sta" <- "Low"
act[cGAS_STING>=median(cGAS_STING), cgas.sting.sta:="High"]

ia.gene <- read.table("immunostimulators.txt")
ia.gene <- as.character(ia.gene$V1)
ic.gene <- read.table("immunoinhibitors.txt", header = FALSE)
ic.gene <- as.character(ic.gene$V1)
im.gene <- unique(c(ia.gene, ic.gene))

exp_im <- exp.data[intersect(im.gene,rownames(exp.data)),con.sample]
exp_long <- as.data.table(reshape2::melt(t(exp_im)))
colnames(exp_long) <- c("SampleID","im.gene","log2Exp")
df <- merge(act, exp_long, by="SampleID")
df$cgas.sting.sta <- factor(df$cgas.sting.sta, levels = c("High","Low"))
df$im.gene <- factor(df$im.gene, levels = sort(unique(as.character(df$im.gene))))

cols <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")
p <- ggplot(df, aes(x = im.gene, y = log2Exp)) +
  geom_boxplot(notch = FALSE, aes(color = cgas.sting.sta), width = 0.8,  
               position = position_dodge(0.7)) +  
  scale_color_manual(values = cols[c(4,3)]) + xlab("") + ylab("log2(TPM)") +
  stat_compare_means(aes(group = cgas.sting.sta), label = "p.signif") +  
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.4M.pdf", p, width = 10,height = 4)





##---------Fig. 4N (Fig. S8D-8F)
df <- reshape2::melt(score_CGAS.STING)
dt <- data.table(sample=df$Var2, cGAS_STING=df$value)
selected <- c("TNF","TNFRSF1B","TNFSF9","TNFRSF9","TNFSF4","TNFRSF4")
df2 <- reshape2::melt(as.matrix(exp.data[selected,]))
colnames(df2) <- c("variable", "sample", "value")
data2 <- merge(dt, df2, by="sample")
colorlist <- c("#e88f18","#a199be","#abd0a7","#ca6a6b","#58a6d6","#e1c548")
p <- ggplot(data2, aes(x=cGAS_STING, y=value, group=variable, color=variable)) + geom_point(size=2,shape=16, alpha = 0.3, aes(color=variable)) + 
  geom_smooth(method = lm) + stat_cor(method = "pearson") +
  scale_color_manual(values=colorlist) +
  xlab("cGAS-STING activity") + ylab("Expression (log2TPM), OAK cohort")  + theme(legend.position="right") + 
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='transparent', color='black')) 
ggsave("Fig.4N.pdf", p, width = 5, height = 6)





##---------Fig. S7D
ligand <- c("Tnf","Tnfsf9","Tnfsf4")
receptor <- c("Tnfrsf1b","Tnfrsf9","Tnfrsf4")
fts <- c(ligand, receptor)
p1 <- FeatureDimPlot(srt = seurat.res, features = fts, 
                     reduction = "UMAP", raster=T, pt.size = 2, theme_use = "theme_blank", ncol = 3)
ggsave("Fig.S7D.pdf", p1, width = 6, height = 5)




##---------Fig. S7G
#Myeloid
obj <- subset(x = seurat.res, subset = Major%in%c("Neutrophils","Monocytes", "Macrophages", "DC cells"))
mat <- GetAssayData(obj, layer = "counts")    
ligand <- c("Tnf","Tnfsf9","Tnfsf4")  

gKO_res <- list()
for(ko_gene in ligand){
  res <- scTenifoldKnk(
    countMatrix = mat,    
    gKO = ko_gene,                                   
    qc_mtThreshold = 0.1,       
    qc_minLSize = 1000,          
    nc_nNet = 10,                 
    nc_nCells = 500,              
    nc_nComp = 3                
  )
  diff_reg <- res$diffRegulation     
  tmp <- as.data.table(diff_reg)
  tmp <- tmp[(order(FC, decreasing = T)),]
  tmp.lst <- list(tmp)
  names(tmp.lst) <- ko_gene
  gKO_res <- c(gKO_res, tmp.lst)
}
save(gKO_res, file="gKO_TNF.ligand_Myeloid.RData")

#Lymphocytes
obj <- subset(x = seurat.res, subset = Major%in%c("T cells","NK cells","B cells"))
mat <- GetAssayData(obj, layer = "counts")    # extract expression data from single-cell object
receptor <- c("Tnfrsf1b","Tnfrsf9","Tnfrsf4")
gKO_res <- list()
for(ko_gene in receptor){
  res <- scTenifoldKnk(
    countMatrix = mat,    
    gKO = ko_gene,                                    
    qc_mtThreshold = 0.1,         
    qc_minLSize = 1000,           
    nc_nNet = 10,               
    nc_nCells = 500,             
    nc_nComp = 3                 
  )
  diff_reg <- res$diffRegulation     # extract differential regulation results after gene knockout
  tmp <- as.data.table(diff_reg)
  tmp <- tmp[(order(FC, decreasing = T)),]
  tmp.lst <- list(tmp)
  names(tmp.lst) <- ko_gene
  gKO_res <- c(gKO_res, tmp.lst)
}
save(gKO_res, file="gKO_TNF.receptor_Lymphocytes.RData")

load("gKO_TNF.ligand_Myeloid.RData")
gKO_res1 <- gKO_res
load("gKO_TNF.receptor_Lymphocytes.RData")
gKO_res2 <- gKO_res
gKO_res <- c(gKO_res1, gKO_res2)
res <- do.call(rbind, gKO_res)
library(dplyr)
scores <- matrix(1,nrow=length(unique(res$gene)),ncol=length(gKO_res))  
rownames(scores) <- unique(res$gene)
colnames(scores) <- names(gKO_res)
for(i in 1:nrow(scores)){
  for(j in 1:ncol(scores)){
    y <- gKO_res[[j]]
    pos <- which(y$"gene" == rownames(scores)[i])
    if(length(pos)>0){
      scores[i,j] <- y[gene == rownames(scores)[i],p.value]   #p.value, p.adj
    }
  }
}

fname_GMT <- "m5.all.v2025.1.Mm.symbols.gmt"    #MSigDB GMT file
gmt <- read.GMT(fname_GMT)
for(i in 1:length(gmt)){
  gmt[[i]]$name <- Hmisc::capitalize(tolower(gsub("_"," ",sub(".*/[A-Z]+_","",gmt[[i]]$name))))  
}
enrich <- ActivePathways(scores, gmt, significant = 0.01, geneset_filter = c(10, 500), correction_method="BH")
colnames(enrich) <- gsub("_",".",colnames(enrich))

select.keys1 <- c("Antigen","Mhc class", "nf kappab","Interferon", "immune",
                  "chemotaxis","Chemokine","cytokine", "Cytokinesis", "Cell killing", 
                  "Lymphocyte","tumor necrosis factor"
)   #chemotaxis
select.keys2 <- c("Myeloid", "Macrophage", "Neutrophil","mononuclear","T cell","B cell","Natural killer","t helper"
)   #immune cells

index <- c()
for(i in c(select.keys1)){   #select.keys1 or select.keys2
  pos <- grep(i, enrich$term.name, ignore.case = T)
  index <- c(index, pos)
}
res1 <- enrich[unique(index),]
unique(res1$term.name)
export_as_CSV(res1, "results_ActivePathways.csv")    
###This results then used as input for the EnrichmentMap plugin in Cytoscape (v3.9.1) to generate a network visualization of the enriched pathways




##---------Fig. S8B
genes.gmt <- unique(gmt[term%in%unique(res1$term.name), gene])
load("gKO_TNF.receptor_Lymphocytes.RData")
res <- gKO_res$"Tnfrsf1b"
res$log_Padj <- -log10(res$p.adj)
res$significant <- ifelse(res$p.adj < 0.05, "Significant", "Not significant")
selected.genes <- res[FC>=2 & p.adj<=0.05,gene]
labeled_genes <- res[gene%in%intersect(c("Tnfrsf1b", genes.gmt), selected.genes),]
y_axis_upper_limit <- quantile(res$log_Padj, 0.999, na.rm = TRUE)
p2 <- ggplot(res, aes(x=Z, y=log_Padj, color=significant)) +
  geom_point(alpha=0.7, size=1) +  
  scale_color_manual(values = c("Significant" = "black", "Not significant" = "gray50")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_text_repel(data=labeled_genes, aes(label=gene), size=3, max.overlaps=50) +
  labs(title="", x="Z-score", y="-log10(Padj)") +
  theme_classic() + 
  coord_cartesian(ylim = c(0, y_axis_upper_limit)) + theme(legend.position = "none")
ggsave("Fig.S8B.pdf", p, width=5, height=5)





##---------Fig. S9D-9E
fts <- c("H2-T23","H2-T22","H2-Q7","H2-Q6","H2-Q4", "H2-M3", "H2-K1", "H2-Eb1",
         "H2-DMb2", "H2-DMb1", "H2-DMa", "H2-D1", "H2-Ab1", "H2-Aa", "Klrd1", "Klrc1")
ht <- GroupHeatmap(
  srt = seurat.res,
  features = fts,
  group.by = c("celltype"),
  heatmap_palette = "YlOrRd",
  row_names_side = "left",
  add_dot = TRUE,
  add_reticle = TRUE
)
p <- ht$plot
ggsave("Fig.S9D.pdf", p, width = 8, height = 4)


ht <- GroupHeatmap(
  srt = seurat.res,
  features = fts,
  group.by = c("group"),
  heatmap_palette = "YlOrRd",
  row_names_side = "left",
  add_dot = TRUE,
  add_reticle = TRUE
)
p <- ht$plot
ggsave("Fig.S9E.pdf", p, width = 4, height = 4)






