library(data.table)
library(GSVA)
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(Hmisc)
library(reshape2)


cols.group <- c("#3361A5","#B24745CC")

##---------Fig. 2I
#siBRIP1
load("exp.data_A549.siRNA.FPKM.RData")
load("cGAS-STING.signature.RData")
gene.list <- lapply(sig_list, function(x){unique(x$gene)})
Param <- ssgseaParam(as.matrix(exp.data), gene.list) 
score <- gsva(Param)
df <- data.table(pathway="cGAS_STING", sample=colnames(score), value=as.numeric(score["cGAS_STING", ]))
df$group <- sub("-\\d+","",rownames(df))
df$group <- factor(df$group, levels = c("siNC","siBRIP1"))
p1 <- ggplot(df, aes(x=group, y=value, fill=group)) +
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.1,position=position_dodge(0.1))+ 
  scale_fill_manual(values = cols.group)+ stat_compare_means()+
  ggtitle("A549") + xlab("") + ylab("cGAS-STING activities") + theme_bw() + 
  theme(legend.position="none",panel.background = element_rect(colour = "black"))

#BRIP1 KO
load("A549.KO-RNAseq_expression.RData")
load("cGAS-STING.signature.RData")
gene.list <- lapply(sig_list, function(x){unique(x$gene)})
Param <- ssgseaParam(as.matrix(exp.data), gene.list) 
score <- gsva(Param)
df <- data.table(pathway="cGAS_STING", sample=colnames(score), value=as.numeric(score["cGAS_STING", ]))
df$group <- sub("-\\d+","", df$sample)
df$group <- factor(df$group, levels = c("WT","BRIP1 KO"))
p2 <- ggplot(df, aes(x=group, y=value, fill=group)) +
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.1,position=position_dodge(0.1))+ 
  scale_fill_manual(values = cols.group)+ stat_compare_means()+
  ggtitle("A549") + xlab("") + ylab("cGAS-STING activities") + theme_bw() + 
  theme(legend.position="none",panel.background = element_rect(colour = "black"))
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("Fig.2I.pdf", p, width = 4, height = 3)





##---------Fig. 2J
load("exp.FPKM_TCGA-NSCLC.RData")
load("clinical.data_TCGA-NSCLC.RData")
load("cGAS-STING.signature.RData")
gene.list <- lapply(sig_list, function(x){unique(x$gene)})
Param <- ssgseaParam(as.matrix(exp.data), gene.list) 
score <- gsva(Param)
df <- data.table(pathway="cGAS_STING", sample=colnames(score), value=as.numeric(score["cGAS_STING", ]))
df <- merge(df, clinical.data, by="PATIENT_ID")
df$BRIP1 <- factor(df$BRIP1, levels = c("WT","BRIP1 Mut"))
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(2,1)]
p <- ggboxplot(df, x = "BRIP1", y = "value",  
               color = "BRIP1", palette = cols, add = "boxplot") + ylab("cGAS-STING activities") +
  stat_compare_means() 
ggsave("Fig.2J.pdf", p, width = 2, height = 3)






##---------Fig. 2K
gmt <- read.gmt("h.all.v7.4.symbols.gmt")    #MSigDB file
genes.select <- gmt[gmt$term%in%"HALLMARK_TNFA_SIGNALING_VIA_NFKB",]
genes.select$term <- Hmisc::capitalize(tolower(gsub("_"," ",sub("[A-Z]+_","",genes.select$term))))    
load("DESeq2_DEG.A549-KO.RData")
dt1 <- na.omit(as.data.table(DESeq2_DEG))
temp <- dt1[order(log2FoldChange,decreasing = T),]  
geneList <- temp$log2FoldChange
names(geneList) <- temp$gene
gsea2 <- GSEA(geneList, TERM2GENE = genes.select, minGSSize = 5,
              maxGSSize = 1000, pvalueCutoff=1) 
temp2 <- as.data.frame(gsea2)[,1:9]
selected.sets <- temp2$ID  
cols2 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
pdf("Fig.2K.pdf", width = 5, height = 4)
gseaplot2(gsea2, geneSetID = selected.sets, pvalue_table = T, subplots=1:2,  
          base_size=12, color=cols2[c(1:length(selected.sets))], 
          title = "TNFa signaling via NF-kB") 
dev.off()




##---------Fig. S2C
#barplot
#siBRIP1
load("exp.data_A549.siRNA.FPKM.RData")
dt <- data.table(sample=colnames(exp.data), value=exp.data["BRIP1",])
dt$group <- sub("-\\d+","", dt$sample)
dt$group <- factor(dt$group, levels = c("siNC","siBRIP1"))

p <- ggplot(dt, aes(x = group, y = value, fill = group)) +  
  geom_bar(stat = "summary", fun ="mean", alpha=0.7, width = 1,   
           position = position_dodge(width = 0.8)) +   
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge( .9)) + 
  stat_compare_means() + scale_fill_manual(values = cols.group) + 
  theme_test(base_size = 12)+  
  xlab("") + ylab("BRIP1 FPKM (A549)") + ylim(0,8) + theme_bw() +  
  theme(legend.position="none",panel.background = element_rect(colour = "black"))
ggsave("Fig.S2C.pdf", p, width = 2, height = 3)











