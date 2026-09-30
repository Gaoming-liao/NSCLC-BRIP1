library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggpubr)

source("src/KM.curve.and.logrank.R")


# Fig. 1A
load("clinical.data_SU2C-MARK.RData")
load("mut.M_SU2C-MARK.RData")
load("DDRgene.RData")
hr.gene <- unique(unlist(DDRgene.core))
con.gene <- intersect(hr.gene, rownames(mut.M))
dt.status <- data.table(SampleID=colnames(mut.M),t(mut.M[con.gene,]))
clinical_su2c <- merge(clinical.data, dt.status, by="SampleID")

univ_formulas <- sapply(con.gene, function(x){as.formula(paste('Surv(OS, OS_Event)~', x))}) 
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = clinical_su2c)})
univ_results <- lapply(univ_models, function(x){
  tmp <-summary(x)
  pval <- round(tmp$coefficients[ ,5], digits = 4)
  HR <- round(tmp$coefficients[ ,2], digits = 3)    
  obj <- rownames(tmp$coefficients)
  all.data <- as.data.frame(cbind(obj, HR, pval))
})
univ_results <- do.call(rbind,univ_results)
univ_results$HR <- as.numeric(as.character(univ_results$HR))
univ_results$pval <- as.numeric(as.character(univ_results$pval))
univ_results$log2HR <- log(univ_results$HR,2)

dt <- as.data.table(univ_results)[pval<=0.05,]
dt <- dt[order(HR),]
dt$freq <- apply(mut.M[dt$obj,],1,sum)/ncol(mut.M)
dt2 <- dt[freq>0.04,]
df <- as.data.frame(dt2)

res <- as.data.frame(univ_results)
res <- res[!is.infinite(res$log2HR),]
res$Significant <- ifelse(res$pval < 0.05 & abs(res$log2HR) != 0,
                          ifelse(res$log2HR > 0, "Risk factors", "Favorable factors"), "Unsignificant")
select.genes <- intersect(con.gene,c("MLH3",DDRgene.core$"HR (Homologous Recombination)"))

p <- ggplot(res, aes(x = log2HR, y = -log10(pval))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("darkblue", "darkred", "grey")) +

  geom_text_repel( 
    data = res[select.genes,],
    aes(x = res[select.genes,"log2HR"], y = -log10(res[select.genes,"pval"]), label = select.genes), 
    box.padding = unit(1, "lines"), point.padding = unit(0.8, "lines")) + 
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8) +   
  geom_hline(yintercept = 1.3,lty=4,col="black",lwd=0.8) +    
  xlab("log2HR")+ylab("-log10 (pval)") + 
  theme(legend.position = "bottom") + theme(title=element_blank()) 
ggsave("Fig.1A.pdf", p, width=5, height=6)



# Fig. 1B
tmp <- data.table(SampleID=colnames(mut.M),t(mut.M[con.gene,]))
clinical.dt <- merge(clinical_su2c, tmp, by="SampleID")
clinical.dt$"int" <- clinical.dt$"BRIP1"
clinical.dt[int==1, int.status:="BRIP1 mut"]
clinical.dt[int==0, int.status:="WT"]
tb <- table(clinical.dt[,c("int.status","Response")])
Ratio <- tb/apply(tb,1,sum)
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(2,1)]
pdf("Fig.1B.pdf", width = 3,height = 6)
par(mar=c(8, 5, 4, 2))
x <- barplot(height=t(Ratio), ylim = c(0, 1), beside=FALSE, las=1, cex.lab=1.3,
             border = 0.1, col = cols, xlab=" ",ylab="Percentage")
abline(h=seq(0.1,1,by=0.1),col="gray",lty=2)
legend("right", legend=colnames(Ratio), pch=16, col=cols)
dev.off()



# Fig. 1C
clinical_sur <- clinical.dt[,c("SampleID","OS_Event", "OS", "int.status")]
colnames(clinical_sur) <- c("Patient_ID", "event", "time", "sample.label")
clinical_sur$"sample.label" <- as.factor(clinical_sur$"sample.label")  
pdf("Fig.1C.pdf",width = 5, height = 6, onefile = FALSE)
plot.surv(clinical_sur,
          upper.time=NULL,
          xscale = "d_m", 
          unit.xlabel = "month", 
          surv.median.line = "none", 
          risk.table = TRUE, 
          pval = TRUE, 
          conf = FALSE, 
          main = "BRIP1 status (SU2C-MARK cohort)", 
          survival.event = "Overall Survival")
dev.off()




# Fig. 1F
load("exp.FPKM_TCGA-NSCLC.RData")
load("clinical.data_TCGA-NSCLC.RData")
BRIP1.exp <- data.table(PATIENT_ID=colnames(exp.FPKM), log2Exp=log(exp.FPKM["BRIP1",]+1,2))
df <- merge(clinical.data, BRIP1.exp, by="PATIENT_ID")
df$int.status <- factor(df$int.status, levels = c("BRIP1 Mut","WT"))
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(2,1)]
p <- ggplot(df, aes(x = int.status, y = log2Exp)) +
  geom_boxplot(notch = FALSE, aes(color = int.status), width = 0.8,  
               position = position_dodge(0.7)) + stat_compare_means() +
  scale_color_manual(values = cols[c(1,2)]) + 
  ggtitle("NSCLC Patients") + xlab("") + ylab("BRIP1 log(FPKM)")
ggsave("Fig.1F.pdf", p, width = 4, height = 4)




# Fig. 1G
load("exp.TPM_cell.RData")
load("cell.data.RData")
BRIP1.exp <- data.table(CELL_ID=colnames(exp.TPM), log2Exp=log(exp.TPM["BRIP1",]+1,2))
df <- merge(cell.data, BRIP1.exp, by="CELL_ID")
df$int.status <- factor(df$int.status, levels = c("BRIP1 Mut","WT"))
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(2,1)]
p <- ggplot(df, aes(x = int.status, y = log2Exp)) +
  geom_boxplot(notch = FALSE, aes(color = int.status), width = 0.8,  
               position = position_dodge(0.7)) + stat_compare_means() +
  scale_color_manual(values = cols[c(1,2)]) + 
  ggtitle("NSCLC Celllines") + xlab("") + ylab("BRIP1 log(TPM)")
ggsave("Fig.1G.pdf", p, width = 4, height = 4)




# Fig. 1H
load("exp.TPM_OAK.RData")
load("clinical.data_OAK.RData")
BRIP1.exp <- data.table(Sample=colnames(exp.data), t(exp.data["BRIP1",]))
dt <- merge(clinical.data, BRIP1.exp, by="Sample")
dt$BRIP1.exp <- "High"
dt[BRIP1<=median(BRIP1), BRIP1.exp:="Low"]    
clinical <- dt[,c("Sample","OS_CENSOR", "OS_MONTHS", "BRIP1.exp")]
colnames(clinical) <- c("Patient_ID", "event", "time", "sample.label")
clinical$"sample.label" <- as.factor(clinical$"sample.label")  
pdf("Fig.1H.pdf",width = 5, height = 5, onefile = FALSE)
plot.surv(clinical,
          upper.time=NULL,
          xscale = "d_m",
          unit.xlabel = "month",
          surv.median.line = "none",
          risk.table = TRUE,
          pval = TRUE,
          conf = FALSE,
          main = "BRIP1 expression (OAK cohort)",
          survival.event = "Overall Survival")
dev.off()




# Fig. 1I
load("exp.data_GSE161537.RData")
load("clinical.data_GSE161537.RData")
BRIP1.exp <- data.table(Sample=colnames(exp.data), t(exp.data["BRIP1",]))
dt <- merge(clinical.data, BRIP1.exp, by="Sample")
res.cut <- surv_cutpoint(dt, time = "OS_months", event = "OS_event",
                         variables = "BRIP1")
optimal_cutpoint <- res.cut$cutpoint$cutpoint
dt$BRIP1.exp <- "High"
dt[BRIP1<optimal_cutpoint, BRIP1.exp:="Low"]  
clinical_sur <- dt[,c("Sample","OS_event", "OS_months", "BRIP1.exp")]
colnames(clinical_sur) <- c("Patient_ID", "event", "time", "sample.label")
clinical_sur$"sample.label" <- as.factor(clinical_sur$"sample.label") 
pdf("Fig.1I.pdf",width = 5, height = 5, onefile = FALSE)
plot.surv(clinical_sur,
          upper.time=NULL,
          xscale = "d_m",
          unit.xlabel = "month",
          surv.median.line = "none",
          risk.table = TRUE,
          pval = TRUE,
          conf = FALSE,
          main = "BRIP1 expression (GSE161537)",
          survival.event = "Overall Survival")
dev.off()




# Fig. 1J
load("exp.TPM_GSE207422.RData")
load("clinical.data_GSE207422.RData")
BRIP1.exp <- data.table(Sample=colnames(exp.TPM), t(log(exp.TPM["BRIP1",]+1,2)))
dt <- merge(clinical.data, BRIP1.exp, by="Sample")
dt$"Pathologic_Response" <- dt$"Pathologic Response"
dt[Pathologic_Response%in%c("NMPR"), Response:="NR"]
dt[Pathologic_Response%in%c("MPR (pCR)","MPR"), Response:="R"]
dt$Response <- factor(dt$Response, levels = c("NR","R"))

cols2 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

dt$patient <- as.character(1:nrow(dt))
dt <- dt[order(BRIP1)]
p1 <- ggplot(data = dt) + 
  geom_bar(aes(x = patient, y = BRIP1, fill = Response), stat = "identity") + 
  scale_fill_manual(values = cols2[c(1,2)]) + 
  scale_x_discrete(limits = dt$patient, expand = c(0.05,0.05)) + 
  xlab("Sample (n=24)") + ylab("BRIP1 (log2TPM)") + 
  theme(legend.position = c(0.01,1.05),
        legend.justification = c(0,1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) + theme_pubclean()
p2 <- ggboxplot(dt, x = "Response", y = "BRIP1", order = c("NR","R"),  
                color = "Response", palette = cols2[c(1,2)], add = "jitter") + 
  ylab("BRIP1 (log2TPM)") + stat_compare_means() + 
  theme(legend.position = 'none')
p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("Fig.1J.pdf",p, width = 6, height = 4)




# Fig. 1K
load("clinical.data_MSK-IMPACT.RData")
clinical.data$log2TMB <- log(clinical.data$TMB,2)
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
p <- ggviolin(clinical.data, x = "int.status", y = "log2TMB", fill = "int.status",
              palette = cols[1:2],
              order=c("BRIP1 Mut","WT"),
              alpha=1, width=0.7, add="boxplot",
              add.params = list(fill="white"))+
  xlab('BRIP1 status') + ylab('log2(TMB)') + ggtitle("MSK-IMPACT cohort") +
  stat_compare_means() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
ggsave("Fig.1K.pdf", p, width = 3, height = 5)




# Fig. S1F
load("clinical.data_TCGA-NSCLC.RData")
df <- clinical.data[,c("PATIENT_ID","BRIP1","Indel.Neoantigens","SNV.Neoantigens")]

df$Indel.Neoantigens <- log(df$Indel.Neoantigens,2)
df$SNV.Neoantigens <- log(df$SNV.Neoantigens,2)

df$BRIP1 <- factor(df$BRIP1, levels = c("BRIP1 Mut", "WT"))
cols <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
p1 <- ggboxplot(df, x = "BRIP1", y = "Indel.Neoantigens", outlier.shape = NA, 
                fill = "BRIP1", palette = cols[c(1,2)], add = "boxplot") + stat_compare_means() + ylab("log2Neoantigens(Indel)")
p2 <- ggboxplot(df, x = "BRIP1", y = "SNV.Neoantigens", outlier.shape = NA, 
                fill = "BRIP1", palette = cols[c(1,2)], add = "boxplot") + stat_compare_means() + ylab("log2Neoantigens(SNV)")
p <- ggarrange(p1,p2, ncol = 2, nrow = 1)   
ggsave("Fig.S1F.pdf", p, width = 4, height = 4)















