library(data.table)
library(GSVA)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(homologene)
library(estimate)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(reshape2)


cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")   #cols.group



##---------Fig. 7F
dt <- fread("grouped_tumorweights_Brip1_LNP.csv")
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
dt$group <- factor(dt$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
my_comparisons <- list(c("PBS + IgG", "siRNA + IgG"),c("PBS + IgG", "PBS + PDL1"), c("PBS + PDL1", "siRNA + PDL1"),c("siRNA + IgG", "siRNA + PDL1"))
p <- ggboxplot(dt, x = "group", y = "Weight", color = "group", add = "jitter", palette = cols.group) + 
  stat_compare_means(comparisons=my_comparisons, method = "t.test") + xlab("") + ylab("Tumor weight (g)") +
  theme_pubclean() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.7F.pdf", p, width = 2.5, height = 4)






##---------Fig. 7G
load("exp.data_Brip1_LNP.FPKM.RData")    
dt.gene <- data.table(gene=rownames(exp.data))
df <- mouse2human(dt.gene$gene, db = homologene::homologeneData)
pos <- match(dt.gene$gene, df$mouseGene)
dt.gene$"humanGene" <- df$"humanGene"[pos]
dt.gene <- na.omit(dt.gene[,c("gene","humanGene")])   
dt.gene <- dt.gene[!duplicated(humanGene),]
int.gene <- intersect(rownames(exp.data), dt.gene$gene)
exp.data <- exp.data[int.gene,]
pos2 <- match(rownames(exp.data), dt.gene$gene)
rownames(exp.data) <- dt.gene$"humanGene"[pos2]
write.table(na.omit(exp.data), file="exp.data_Brip1_LNP.FPKM.txt", quote = FALSE, sep="\t", row.names = TRUE)

# Calculate immune score using ESTIMATE
expr <- "exp.data_Brip1_LNP.FPKM.txt"
filterCommonGenes(input.f=expr, output.f="Brip1_LNP_genes.gct", id="GeneSymbol")
estimateScore(input.ds = "Brip1_LNP_genes.gct", output.ds="Brip1_LNP_estimate_score.gct")
estimate.score <- read.table("Brip1_LNP_estimate_score.gct",skip = 2,header = T)
rownames(estimate.score) <- estimate.score[,1]
estimate.score <- t(estimate.score[,3:ncol(estimate.score)])
rownames(estimate.score) <- gsub("\\.","-", rownames(estimate.score))
estimate.score <- data.table(sample=rownames(estimate.score),estimate.score)
save(estimate.score, file="Brip1_LNP_estimate.score.RData")

dt <- fread("Brip1_LNP.mouse_groups.csv")   #group
load("Brip1_LNP_estimate.score.RData")  
df <- merge(dt, estimate.score, by="sample")
df$group <- factor(df$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
my_comparisons <- list(c("PBS + IgG", "siRNA + IgG"),c("PBS + PDL1", "siRNA + PDL1"))
p <- ggboxplot(df, x = "group", y = "Immune score",
                color = "group", add = "jitter", palette = cols.group) + 
  stat_compare_means(comparisons=my_comparisons, method = "t.test")+ 
  theme_pubclean()+theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.7G.pdf", p, width = 3, height = 4)




##---------Fig. 7H
df <- data.table(pathway="cGAS_STING", sample=colnames(score), value=as.numeric(score["cGAS_STING", ]))
dt <- fread("Brip1_LNP.mouse_groups.csv")
df <- merge(df, dt, by="sample")
df$group <- factor(df$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
my_comparisons <- list(c("PBS + IgG", "siRNA + IgG"),c("PBS + PDL1", "siRNA + PDL1"))
p <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(),alpha=0.7) +
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge( .9)) + ylab("cGAS-STING activities") + xlab("")+
  stat_compare_means(comparisons=my_comparisons, 
                     method="t.test", label="p.signif", 
                     bracket.size = 0.6,size=4) + 
  scale_fill_manual(values = cols.group) + 
  theme_test(base_size = 13)+ 
  theme(legend.position = "none",  
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        plot.title = element_text(hjust = 0.5,size = 13),
        axis.title = element_text(size = 13,color = 'black'))+
  geom_jitter(size=1, alpha=0.5, shape =21)  
ggsave("Fig.7H.pdf", p, width = 6, height = 4)




##---------Fig. S13B
# hallmark
gmt <- read.gmt("mh.all.v2026.1.Mm.symbols.gmt")  #MSigDB file
gmt$term <- Hmisc::capitalize(tolower(gsub("_"," ",sub("[A-Z]+_","",gmt$term))))  
dt.hallmark <- as.data.table(gmt)

# ImmPort
pathway <- fread("GeneList_from_ImmPort.txt")
unique(pathway$Category)
pathway <- pathway[Category!="Antimicrobials",]
genes.select <- pathway[,c("Category","Symbol")]
colnames(genes.select) <- c("term","gene")
dt.ImmPort <- genes.select

# Immune Signature
pathway <- read.delim("Immune Signature Gene.txt",sep="\t", header = TRUE)
pathway.list <- lapply(1:nrow(pathway),function(i){  
  x <- as.character(pathway$Genes[i])
  genes <- unique(unlist(strsplit(x,split=" ")))
  tmp <- data.table(term=pathway$Signatures[i], gene=genes)
  return(tmp)
})
genes.select <- do.call(rbind, pathway.list)
dt.ImmuneSig <- genes.select

# ImmuneCellTypeGeneSets
pathway <- read.delim("immuneCellTypeGeneSets.txt",sep="\t", header = TRUE)
pathway.list <- lapply(1:nrow(pathway),function(i){  
  x <- as.character(pathway[i,-1])
  genes <- unique(x[x!="NA" & x!=""])
  tmp <- data.table(term=pathway$Cell.type[i], gene=genes)
  return(tmp)
})
genes.select <- do.call(rbind, pathway.list)
dt.ImmuneCell <- genes.select
dt.all <- rbind(dt.ImmPort,
                rbind(dt.ImmuneSig,
                      dt.ImmuneCell))
df <- human2mouse(dt.all$gene, db = homologene::homologeneData)
pos <- match(dt.all$gene, df$humanGene)
dt.all$"mouseGene" <- df$"mouseGene"[pos]
dt.all2 <- na.omit(dt.all[,c("term","mouseGene")])   
colnames(dt.all2) <- c("term","gene")
dt.all <- rbind(dt.all2, dt.hallmark)   

load("DESeq2_DEG.Brip1_LNP_IgG.RData")   #siRNA vs. PBS (in IgG)
res$obj <- rownames(res)
dt1 <- na.omit(as.data.table(res))
temp <- dt1[order(log2FoldChange,decreasing = T),]  
geneList <- temp$log2FoldChange
names(geneList) <- temp$obj    
##GSEA
gsea1 <- GSEA(geneList, TERM2GENE = dt.all, minGSSize = 5,
              maxGSSize = 1000, pvalueCutoff=1) 
temp1 <- as.data.table(gsea1)[,1:9]
temp1 <- temp1[pvalue<=0.05,]

load("DESeq2_DEG.Brip1_LNP_PDL1.RData")     #siRNA vs. PBS (PDL1)
res$obj <- rownames(res)
dt1 <- na.omit(as.data.table(res))
temp <- dt1[order(log2FoldChange,decreasing = T),]  
geneList <- temp$log2FoldChange
names(geneList) <- temp$obj    
##GSEA
gsea2 <- GSEA(geneList, TERM2GENE = dt.all, minGSSize = 5,
              maxGSSize = 1000, pvalueCutoff=1) 
temp2 <- as.data.table(gsea2)[,1:9]
temp2 <- temp2[pvalue<=0.05,]

id <- unique(c(temp1$ID, temp2$ID))
res.gsea1 <- temp1[ID%in%id,c("ID","NES","pvalue"),]
res.gsea1$"group" <- "siRNA-PBS(IgG)"
res.gsea2 <- temp2[ID%in%id,c("ID","NES","pvalue"),]
res.gsea2$"group" <- "siRNA-PBS(PDL1)"
res.gsea <- rbind(res.gsea1, res.gsea2)
notIn <- c("Xenobiotic metabolism","Spermatogenesis","Uv response dn", "Uv response up", "Pi3k akt mtor signaling",
           "Allograft rejection","Androgen response","Adipogenesis","Apical junction", "Myogenesis",
           "Estrogen response early","Estrogen response late","Bile acid metabolism","Kras signaling dn",
           "Cholesterol homeostasis","Mitotic spindle","GEP","B_cells","Coagulation","Hedgehog signaling")
data.fr <- res.gsea[!ID%in%notIn, c("ID","group","NES","pvalue")]
colnames(data.fr) <- c("yname", "xname", "NES", "pval")
data.fr$pval <- -log(data.fr$pval,10)
cols5 <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
library(ggplot2)
p <- ggplot(data.fr, aes(x = factor(xname),  y = factor(yname))) + 
  geom_point(aes(colour = NES,  size =pval)) + 
  scale_color_gradientn(colours = cols5) +
  scale_size(range = c(0,6)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("") + ylab("") + theme_bw() +   
  theme(axis.text.x=element_text(size=12, face="plain", angle=30, hjust=1)) + 
  theme(axis.text.y=element_text(size=12, face="plain"))
ggsave("Fig.S13B.pdf", p, width=5, height=8)



##---------Fig. 7I
int.ters <- c("Co_inhibition_T_cell","Cytotoxic_cells","Co_stimulation_T_cell","NaturalKiller_Cell_Cytotoxicity","CD8+ T-cells","CD_8_T_effector","Neutrophils")
genes.select <- dt.all[term%in%int.ters,]
load("DESeq2_DEG.Brip1_LNP_PDL1.RData")     #siRNA vs. PBS (in PDL1)
res$obj <- rownames(res)
dt1 <- na.omit(as.data.table(res))
temp <- dt1[order(log2FoldChange,decreasing = T),]  
geneList <- temp$log2FoldChange
names(geneList) <- temp$obj
gsea2 <- GSEA(geneList, TERM2GENE = genes.select, minGSSize = 5,
              maxGSSize = 1000, pvalueCutoff=1) 
temp <- as.data.table(gsea2)[,1:9]
temp2 <- temp[pvalue<=0.05,]
selected.sets <- temp2$ID   
cols1 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
pdf("Fig.7I.pdf", width = 6, height = 6)
gseaplot2(gsea2, geneSetID = selected.sets, pvalue_table = T, 
          base_size=12, color=cols1[c(1:length(selected.sets))]) 
dev.off()






##---------Fig. 7J
load("exp.data_Brip1_LNP.FPKM.RData")
pathway <- read.delim("Immune Signature Gene.txt",sep="\t", header = TRUE)
pathway.list <- lapply(1:nrow(pathway),function(i){  
  x <- as.character(pathway$Genes[i])
  genes <- unique(unlist(strsplit(x,split=" ")))
  tmp <- data.table(term=pathway$Signatures[i], gene=genes)
  return(tmp)
})
names(pathway.list) <- pathway$Signatures
gene.list <- lapply(pathway.list, function(x){unique(x$gene)})
gene.list.mm <- lapply(gene.list, function(x){
  y <- homologene::homologene(x, inTax = 9606, outTax = 10090)$"10090"  
  return(y)
})

Param <- gsvaParam(as.matrix(exp.data), gene.list.mm) 
score <- gsva(Param)
df <- reshape2::melt(score)
colnames(df) <- c("pathway", "sample", "value")
df <- as.data.table(df)
dt <- fread("Brip1_LNP.mouse_groups.csv")
df <- merge(df, dt, by="sample")
df$group <- factor(df$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
y <- unique(df[,c("sample","group")])
anno_col <- as.data.frame(y)
rownames(anno_col) <- anno_col$sample
anno_col$sample <- NULL
annColors <- list()
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
names(cols.group) <- c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1")  
annColors[['group']] <- cols.group
color <- c("royalblue1", "royalblue2", "royalblue3", "royalblue4", "black", "yellow4","yellow3","yellow2", "yellow1")
pdf("Fig.7J.pdf", width= 5.2, height=3.6)
pheatmap(score, cluster_cols = T, cluster_rows = T, 
         annotation_col=anno_col, 
         annotation_colors = annColors,
         show_rownames=T,show_colnames=T, scale='row',  
         treeheight_row=0, treeheight_col=0,
         color = colorRampPalette(color)(50)
)
dev.off()



##---------Fig. 7K
selected <- c("Macrophages", "Monocytes", "Neutrophils", "aDC", 
              "CD8+ T-cells", "CD8+ Tem", "CD4+ T-cells","Th2 cells","NK cells","B-cells")
df2 <- df[pathway%in%selected,]
df2$pathway <- factor(df2$pathway, levels = selected)
df3 <- data.table()
for(i in selected){
  y <- as.data.frame(df)[df$pathway==i,]
  fac <- factor(y$group)
  value <- tapply(y$value, fac, mean_ci)
  tmp <- do.call(rbind, value)
  res <- data.table(group=rownames(tmp),cell=i,tmp)
  df3 <- rbind(df3,res)
}
df3$group <- factor(df3$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
df3$cell <- factor(df3$cell, levels = selected)
p <- ggplot(df3, aes(x = cell, y = y, fill = group))+
  geom_col(position = 'dodge', width = 0.8) +
  geom_errorbar(aes(x = cell, ymin = ymin, ymax = ymax),
                width = 0.2, linewidth = 0.1, color='black', position = position_dodge(0.8))+
  scale_fill_manual(values = cols.group) + 
  xlab("") + ylab("Activities") +  theme(legend.position = "top") +
  theme_pubclean()+theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.7K.pdf", p, width = 5, height = 4)







##---------Fig. S13D
load("exp.data_Brip1_LNP.FPKM.RData")
pathway <- read.delim("Immune Signature Gene.txt",sep="\t", header = TRUE)
pathway.list <- lapply(1:nrow(pathway),function(i){  
  x <- as.character(pathway$Genes[i])
  genes <- unique(unlist(strsplit(x,split=" ")))
  tmp <- data.table(term=pathway$Signatures[i], gene=genes)
  return(tmp)
})
names(pathway.list) <- pathway$Signatures
gene.list <- lapply(pathway.list, function(x){unique(x$gene)})
gene.list.mm <- lapply(gene.list, function(x){
  y <- homologene::homologene(x, inTax = 9606, outTax = 10090)$"10090"  
  return(y)
})

Param <- gsvaParam(as.matrix(exp.data), gene.list.mm) 
score <- gsva(Param)
df <- reshape2::melt(score)
colnames(df) <- c("pathway", "sample", "value")
df <- as.data.table(df)
dt <- fread("Brip1_LNP.mouse_groups.csv")
df <- merge(df, dt, by="sample")
df$group <- factor(df$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
y <- unique(df[,c("sample","group")])
anno_col <- as.data.frame(y)
rownames(anno_col) <- anno_col$sample
anno_col$sample <- NULL
annColors <- list()
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
names(cols.group) <- c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1")  
annColors[['group']] <- cols.group
color <- c("royalblue1", "royalblue2", "royalblue3", "royalblue4", "black", "yellow4","yellow3","yellow2", "yellow1")
pdf("Fig.S13D.pdf", width= 5.6, height=4)
pheatmap(score, cluster_cols = T, cluster_rows = T, 
         annotation_col=anno_col, 
         annotation_colors = annColors,
         show_rownames=T,show_colnames=T, scale='row',  
         treeheight_row=0, treeheight_col=0,
         color = colorRampPalette(color)(50)
)
dev.off()



##---------Fig. S13E
selected <- c("Co_stimulation_APC", "Co_stimulation_T_cell",
              "MHC_Class_I", "IFNG_sigture", "Cytotoxic_cells",
              "CD8","CD_8_T_effector","NK_cells","Neutrophils","aDC")
df2 <- df[pathway%in%selected,]
df2$pathway <- factor(df2$pathway, levels = selected)
df3 <- data.table()
for(i in selected){
  y <- as.data.frame(df)[df$pathway==i,]
  fac <- factor(y$group)
  value <- tapply(y$value, fac, mean_ci)
  tmp <- do.call(rbind, value)
  res <- data.table(group=rownames(tmp),cell=i,tmp)
  df3 <- rbind(df3,res)
}
df3$group <- factor(df3$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
df3$cell <- factor(df3$cell, levels = selected)
p <- ggplot(df3, aes(x = cell, y = y, fill = group))+
  geom_col(position = 'dodge', width = 0.8) +
  geom_errorbar(aes(x = cell, ymin = ymin, ymax = ymax),
                width = 0.2, linewidth = 0.1, color='black', position = position_dodge(0.8))+
  scale_fill_manual(values = cols.group) + 
  xlab("") + ylab("Activities") +  theme(legend.position = "top") +
  theme_pubclean()+theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.S13E.pdf", p, width = 5, height = 4)







##---------Fig. 7L
load("exp.data_Brip1_LNP.FPKM.RData")
exp.data <- log(exp.data+1, 2)
load("DESeq2_DEG.Brip1_LNP.siRNA_PBS.IgG.RData")   #siRNA vs. PBS (in IgG)
res$obj <- rownames(res)
res$"group" <- "siRNA_PBS(IgG)"
res1 <- na.omit(as.data.table(res))
load("DESeq2_DEG.Brip1_LNP.siRNA_PBS.PDL1.RData")     #siRNA vs. PBS (in PDL1)
res$obj <- rownames(res)
res$"group" <- "siRNA_PBS(PDL1)"
res2 <- na.omit(as.data.table(res))
res <- rbind(res1, res2)
res <- res[(pvalue <= 0.01 & abs(log2FoldChange) >= 0.58),]  

ia.gene <- read.table("immunostimulators.txt")
ia.gene <- as.character(ia.gene$V1)
ic.gene <- read.table("immunoinhibitors.txt", header = FALSE)
ic.gene <- as.character(ic.gene$V1)
immunogene <- c(ia.gene, ic.gene)
df.match <- human2mouse(immunogene, db = homologene::homologeneData)
selected_genes <- c("Cxcl8", "Cxcl9", "Cxcl10", "Ccl5", "Ccl8","Ccl4", "Ccl3", "Stat1", "Ccr2", "Irf3")
inte.gene <- intersect(c(df.match$mouseGene, selected_genes), res$obj)
exp.deal <- as.matrix(exp.data[inte.gene, , drop = FALSE])
df <- reshape2::melt(exp.deal)
colnames(df) <- c("gene", "sample", "value")
df <- as.data.table(df)
dt <- fread("Brip1_LNP.mouse_groups.csv")
df <- merge(df, dt, by="sample")
df$group <- factor(df$group, levels = c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1"))
y <- unique(df[,c("sample","group")])
anno_col <- as.data.frame(y)
rownames(anno_col) <- anno_col$sample
anno_col$sample <- NULL
annColors <- list()
cols.group <- c("#AECDE1", "#3C77AF","#FDAF91CC","#BC3C29CC")
names(cols.group) <- c("PBS + IgG", "siRNA + IgG", "PBS + PDL1", "siRNA + PDL1")  
annColors[['group']] <- cols.group
color <- c("#118ab2", "#fdffb6", "#e63946")
pdf("Fig.7L.pdf", width= 5.6, height=4)
pheatmap(exp.deal, cluster_cols = T, cluster_rows = T, 
         annotation_col=anno_col, 
         annotation_colors = annColors,
         show_rownames=T,show_colnames=T, scale='row', 
         treeheight_row=0, treeheight_col=0,
         color = colorRampPalette(color)(50)
)
dev.off()





##---------Fig. 7M
selected_genes <- c("Cd274", "Pdcd1")
my_comparisons <- list(c("PBS + IgG", "siRNA + IgG"), c("PBS + PDL1", "siRNA + PDL1"))
ps <- list()
for(i in selected_genes){
  tmp <- df[gene==i,]
  pi <- ggboxplot(tmp, x = "group", y = "value", color = "group", add = "jitter", palette = cols.group) + 
    stat_compare_means(comparisons=my_comparisons, method = "t.test") + xlab("") + ylab(i) +
    theme_pubclean() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
  ps <- c(ps, list(pi))
}
p <- ggarrange(ps[[1]], ps[[2]], ncol = 2, nrow = 1)   
ggsave("Fig.7M.pdf", p, width = 3.5, height = 3)




##---------Fig. S13F
selected_genes <- c("Cd86","Cd80", "Cd28", "Ccl5", "Cxcl10", "Cxcl9", "Ccl4", "Cxcl8", "Ccl8", "Cd274", "Pdcd1", "Lag3")
df2 <- df[gene%in%selected_genes,]
df2$gene <- factor(df2$gene, levels = selected_genes)
p <- ggplot(df2, aes(gene, value, fill = group)) +
  geom_violin(aes(fill = group),width = 0.7, scale = "width", size = 0.1) +
  geom_boxplot(width = 0.2, position = position_dodge(0.7),  
               outlier.shape = NA, cex = 0.1) + xlab("")+ylab("log2(FPKM+1)")+
  scale_fill_manual(values = cols.group) +
  stat_compare_means(method="anova", aes(group = group), label = "p.signif",
                     symnum.args = list(cutpoint=c(0,0.0001,0.001,0.01,0.05,1),
                                        symbols=c("****","***","**","*", NA))) +  
  theme_bw()+theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.S13F.pdf", p, width = 8, height = 3)














