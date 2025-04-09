#### KI proteome part ####
title: "KI_study_proteome"
author: "Jeppe Kjærgaard Larsen"
date: "2023-05-11"
output: html_document
#### loading libraries ####
library(biomaRt)
library(clusterProfiler)
library(dplyr)
library(factoextra)
library(fgsea)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(gplots)
library(GOSemSim)
library(limma)
library(msigdbr)
library(org.Hs.eg.db)
library(patchwork)
library(pheatmap)
library(PhosR)
library(scales)
library(stringr)
library(tidyverse)


#### loading files ####
#Grouping information
Sample_ID <- read.delim("Input_SA_info_R.txt")
Sample_ID_prot <- Sample_ID[-c(24,45,57,81,125),]
Sample_ID_prot$Gender[grepl("Male",Sample_ID_prot$Group)] <- "Male"
Sample_ID_prot$Gender[grepl("Female",Sample_ID_prot$Group)] <- "Female"
Sample_ID_prot$disease[grepl("T2D",Sample_ID_prot$Group)] <- "T2D"
Sample_ID_prot$disease[grepl("NGT",Sample_ID_prot$Group)] <- "NGT"
write_rds(Sample_ID_prot,"Sample_ID_prot.rds")

#proteome file
proteome_study_raw <- read.delim("20230704_095840_20221104_Proteome_directDIA_KI_v18_Report.tsv")

#### processing data ####
rownames(proteome_study_raw) <- paste0(proteome_study_raw$PG.Genes,"_",proteome_study_raw$PG.ProteinAccessions)
proteome_study_raw <- proteome_study_raw[,grepl("PG.Quantity",colnames(proteome_study_raw))]
proteome_study <- log2(proteome_study_raw)
new_colnames <- word(colnames(proteome_study),9, sep="_")
new_colnames <- as.integer(sub("S","",new_colnames))
colnames(proteome_study) <- new_colnames
colnames(proteome_study) <- as.integer(colnames(proteome_study))
proteome_study <- proteome_study[,str_sort(colnames(proteome_study), numeric = TRUE)]

####  PCA - batch effect ####
complete_proteome <- na.omit(proteome_study)
t_exprs <- as.data.frame(t(complete_proteome))
pca_res <- prcomp(t_exprs[,1:1428], scale. = TRUE)
df_out<- as.data.frame(pca_res$x)
df_out$batch <- "day1 batch"
df_out$batch[c(81:154)] <- "day2 batch"
df_out$contam <- "valid"
df_out$contam[c(24,45,57,81,125)] <- "outlier"
df_out$set <-factor(c(rep(1, 80),rep(2, 19),rep(3,55))) 
theme <- theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black"),panel.background = element_blank(),
               panel.border=element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),strip.background=element_blank(),
               axis.text.x=element_text(colour="black", size = 18),
               axis.text.y=element_text(colour="black", size = 18),
               axis.ticks=element_line(colour="black"), legend.title  = element_text(size=18), legend.text = element_text(size = 16))

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )


p<-ggplot(df_out,aes(x=PC1,y=PC2, colour = set))
p<-p+geom_point(size = 2) + theme + xlab(percentage[1]) + ylab(percentage[2]) #stat_ellipse(level=0.8)
p

##PCA on batch corrected data to identify outliers. The outliers identified here were also noted as "bloody" or contaminated during sample preparation
set_batch_effect <- factor(c(rep(1, 80),rep(2, 19),rep(3,55))) 
corrected_prot <- removeBatchEffect(complete_proteome, set_batch_effect)
t_exprs <- as.data.frame(t(corrected_prot))
pca_res <- prcomp(t_exprs[,1:1423], scale. = TRUE)
df_out<- as.data.frame(pca_res$x)
df_out$batch <- "day1 batch"
df_out$batch[c(81:154)] <- "day2 batch"
df_out$contam <- "valid"
df_out$contam[c(24,45,57,81,125)] <- "outlier"
df_out$set <-factor(c(rep(1, 80),rep(2, 19),rep(3,55))) 

theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=18),
             legend.text = element_text(size = 16))
percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
p<-ggplot(df_out,aes(x=PC1,y=PC2, colour = contam))
p<-p+geom_point(size = 2) + theme + xlab(percentage[1]) + ylab(percentage[2]) #stat_ellipse(level=0.8)
p

### Final PCA without outliers and batch-corrected
corrected_prot <- corrected_prot[,-c(24,45,57,81,125)]
write_rds(corrected_prot,"corrected_prot.rds")

t_exprs <- as.data.frame(t(corrected_prot))
pca_res <- prcomp(t_exprs[,1:1428], scale. = TRUE)
df_out<- as.data.frame(pca_res$x)
df_out$batch <- "day1 batch"
df_out$batch[c(82:ncol(corrected_prot))] <- "day2 batch"

write_rds(pca_res,"pca_res_prot.rds")

df_out$disease <- Sample_ID_prot$disease
df_out$Gender <- Sample_ID_prot$Gender
df_out$Clamp <- Sample_ID_prot$Clamp
df_out$Subject <- Sample_ID_prot$Subject_ID
df_out$group <- Sample_ID_prot$Group

Sample_ID_prot$Subject_ID <- sub("IRS","IRS1",Sample_ID_prot$Subject_ID)
Our_clin_data_prot <- Clinical_data[Clinical_data$Subject.ID %in% Sample_ID_prot$Subject_ID,]
Our_clin_data_prot <- Our_clin_data_prot[order(match(Our_clin_data_prot[,1],Sample_ID_prot[,5])),]

# Duplicate rows based on unique IDs while preserving the order
duplicated_df_prot <- Our_clin_data_prot %>%
  uncount(weights = 2, .remove = FALSE)
duplicated_df_prot <- duplicated_df_prot[-c(24,45,57,81,125),]

df_out$GIR <- duplicated_df_prot$M..µmol..kg.min..
df_out$HOMA1_IR <- duplicated_df_prot$HOMA1.IR
df_out$FFA <- duplicated_df_prot$FFA..mmol.L.
df_out$BMI <- duplicated_df_prot$BMI..kg.m2.
df_out$HDL <- duplicated_df_prot$HDL.Chol..mmol.L.
df_out$Age <- duplicated_df_prot$Age

theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=18),
             legend.text = element_text(size = 16))

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )


#### Figure 1H - PCA proteome colored by M-value #### 
p<-ggplot(df_out,aes(x=PC1,y=PC2, colour = GIR))
p<-p+geom_point(size = 2) + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  scale_colour_gradient(low = "lightblue",high = "darkred", space = "Lab", na.value = "grey50",guide = "colourbar", aesthetics = "colour")
p

#### Figure 1 S2C - discrete NGT/T2D comparison PCA ####
p<-ggplot(df_out,aes(x=PC1,y=PC2, colour = disease)) + geom_point(size = 3, alpha=0.8) + 
  theme + xlab(percentage[1]) + ylab(percentage[2]) +
  scale_colour_manual(values=c("darkgrey","darkred"))
p + stat_ellipse(level=0.8)


#### Differentially abundant proteins - LIMMA ####
Proteome_muscle2 <- removeBatchEffect(proteome_study, set_batch_effect)
Proteome_muscle2 <- Proteome_muscle2[,-c(24,45,57,81,125)]
Proteome_muscle2 <- medianScaling(Proteome_muscle2)
total_SA <- rep(1,149)
Proteome_muscle2 <- selectGrps(Proteome_muscle2, total_SA, 0.25, n=1)
write_rds(Proteome_muscle2,"Proteome_muscle2.rds")

Clamp <- factor(Sample_ID_prot$Clamp, levels=c("Pre","Post"))
Disease <- factor(Sample_ID_prot$disease, levels=c("NGT","T2D"))
Gender <- factor(Sample_ID_prot$Gender, levels=c("Male","Female"))
design <- model.matrix(~ Disease + Clamp + Gender)
corfit <- duplicateCorrelation(Proteome_muscle2, design, block=Sample_ID_prot$Subject_ID)
corfit$consensus
fit <- lmFit(Proteome_muscle2,design,block=Sample_ID_prot$Subject_ID,correlation=corfit$consensus)
fit2 <- eBayes(fit)

Disease_main_effect <- topTable(fit2, 
                                coef=2, 
                                number = Inf, sort.by = "none")

Disease_main_effect$ID <- rownames(Disease_main_effect)

Insulin_main_effect <- topTable(fit2, 
                                coef=3, 
                                number = Inf, sort.by = "none")

Insulin_main_effect$ID <- rownames(Insulin_main_effect)

Gender_main_effect <- topTable(fit2, 
                               coef=4, 
                               number = Inf, sort.by = "none")

Gender_main_effect$ID <- rownames(Gender_main_effect)

#### Extracting individual proteins ####
Protein_example <- Proteome_muscle2[rownames(Proteome_muscle2) == "FHL1_Q13642",]
plot_Protein_example <- cbind(Protein_example, Sample_ID_prot)
colnames(plot_Protein_example)[1] <- "Intensity"
plot_Protein_example$paired_ID <- factor(str_sub(plot_Protein_example$Subject_ID, start= -3))
plot_Protein_example$Clamp <- factor(plot_Protein_example$Clamp, levels = c("Pre","Post"))

p <-ggplot(plot_Protein_example, aes(x = factor(Clamp,levels = c("Pre","Post")), y = Intensity)) + 
  geom_boxplot(aes(fill = Clamp), alpha = .2) +
  geom_line(aes(group = paired_ID)) + 
  geom_point(size = 2) + 
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14),
        axis.line.y = element_line(size = .5)) + ylab("Log2 Intensity") +
  ggtitle("RCSD1")
p


#### Protein abundance & insulin sensitivity associations ####
Clinical_data <- read.delim("Combined_clinical_KI_data.txt", dec=",") #Clinical data won't be provided
Sample_ID_prot$Subject_ID <- sub("IRS","IRS1",Sample_ID_prot$Subject_ID)
pre_proteome_intensities <- Proteome_muscle2[,which(Sample_ID_prot$Clamp == "Pre")]
pre_proteome_Sample_ID <- Sample_ID_prot[which(Sample_ID_prot$Clamp == "Pre"),]

Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% pre_proteome_Sample_ID$Subject_ID,]
target_proteome_df <- Our_clin_data[match(unique(pre_proteome_Sample_ID$Subject_ID),Our_clin_data$Subject.ID),]

clin_var <- c("M..µmol..kg.min..","HbA1c..mmol.mol.","BMI..kg.m2.","Age",
              "FS.Insulin..mIE.L.","HOMA1.IR","Creatinine..µmol.L.","P.TG..mmol.L.",
              "P.Chol..mmol.L.","FFA..mmol.L.","FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.")

Prot_baseline_p.val_matrix <- matrix(nrow=nrow(pre_proteome_intensities),ncol=length(clin_var))
Prot_baseline_estimate_matrix <- matrix(nrow=nrow(pre_proteome_intensities),ncol=length(clin_var))

empty_list <- list()
for(k in 1:length(clin_var)){
  try({
    for(i in 1:nrow(pre_proteome_intensities)){
      empty_list <- cor.test(scale(pre_proteome_intensities[i,]), scale(log2(target_proteome_df[,clin_var[k]])), method="kendall", na.action = "na.exclude")
      Prot_baseline_p.val_matrix[i,k] <- empty_list$p.value
      Prot_baseline_estimate_matrix[i,k] <- empty_list$estimate
    }
  }, silent = TRUE)
}

colnames(Prot_baseline_p.val_matrix) <- paste0(clin_var,"_p.value")
colnames(Prot_baseline_estimate_matrix) <- paste0(clin_var,"_estimate")

alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

Prot_new_matrix_adjust <- matrix(nrow=nrow(Prot_baseline_p.val_matrix),ncol=ncol(Prot_baseline_p.val_matrix))
for(i in 1:12){
  Prot_new_matrix_adjust[,i] <- p.adjust(Prot_baseline_p.val_matrix[,i],method="BH")
}
colnames(Prot_new_matrix_adjust) <- paste0(clin_var,"_adj.p.value")
Prot_complete_p.val_matrix <- alternate.cols(Prot_baseline_p.val_matrix,Prot_new_matrix_adjust)

Prot_complete_cor_matrix <- cbind(Prot_complete_p.val_matrix,Prot_baseline_estimate_matrix)
Prot_complete_cor_frame <- data.frame(Prot_complete_cor_matrix)
Prot_complete_cor_frame <- cbind(rownames(pre_proteome_intensities), Prot_complete_cor_frame)
Prot_complete_cor_frame <- Prot_complete_cor_frame[,c(1:3,26,4:5,27,6:7,28,8:9,29,10:11,30,12:13,31,
                                                      14:15,32,16:17,33,18:19,34,20:21,35,22:23,36,
                                                      24:25,37)]
colnames(Prot_complete_cor_frame)[1] <- "ID"

Prot_PRE_GIR_Cor <- Prot_complete_cor_frame[,c(1,2:4)]
colnames(Prot_PRE_GIR_Cor) <- c("ID","p.val","adj.pval","estimate")
Prot_PRE_GIR_Cor$diff <- "NO"
Prot_PRE_GIR_Cor[Prot_PRE_GIR_Cor$adj.pval <= 0.05,]$diff <- "yes"
Prot_PRE_GIR_Cor$short_ID <- word(word(Prot_PRE_GIR_Cor$ID,1,sep="_"),1,sep=";")
write_rds(Prot_PRE_GIR_Cor,"Prot_PRE_GIR_Cor.rds")

#### Figure 2A - make volcano correlation plot ####
ggplot() +
  geom_point(data = Prot_PRE_GIR_Cor[Prot_PRE_GIR_Cor$adj.pval < 0.05, ],
             aes(x = estimate, y = -log10(p.val), color = diff),
             size = 2, alpha = 0.3, stroke = NA) +
  scale_color_manual(values = c("grey", "darkred")) +
  geom_point(data = Prot_PRE_GIR_Cor[Prot_PRE_GIR_Cor$adj.pval > 0.05, ],
             aes(x = estimate, y = -log10(p.val), color = diff),
             size = 0.75, alpha = 0.3, stroke = NA) +
  geom_text_repel(
    data = Prot_PRE_GIR_Cor[Prot_PRE_GIR_Cor$adj.pval < 0.05, ],
    aes(label = short_ID, x = estimate, y = -log10(p.val)),
    size = 3,
    nudge_x = 0.05,
    nudge_y = 0.05,
    colour = "darkred"  # Set label colour to dark red
  ) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=16, hjust = 1),
        axis.title = element_text(size = 18), aspect.ratio=1/1, legend.position="none") + xlim(-0.5,0.5)



#### Figure 2F LDHA/LDHB ratio corrleaiton with insulin sensitivity ####
LDHB <- as.data.frame(cbind(pre_proteome_intensities[which(rownames(pre_proteome_intensities) == "LDHB_P07195"),],target_proteome_df$M..µmol..kg.min..,target_proteome_df$Group))
colnames(LDHB) <- c("Intensity","GIR","Group")             

LDH_ratio <- cbind(LDHB,pre_proteome_intensities[which(rownames(pre_proteome_intensities) == "LDHA_P00338"),])
colnames(LDH_ratio) <- c("Intensity","GIR","Group","LDHA")             

LDH_ratio$stoc <- LDH_ratio$LDHA-as.numeric(LDH_ratio$Intensity)
cor.test(as.numeric(LDH_ratio$stoc),as.numeric(LDH_ratio$GIR),method="kendall")
cor.test(as.numeric(LDH_ratio$LDHA),as.numeric(LDH_ratio$GIR),method="kendall")

ggplot(LDH_ratio, aes(x=as.numeric(GIR), y=stoc)) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("LDHA/LDHB ratio") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=4, label="tau = -0.41", color = "black",size=4) +
  annotate(geom = "text",x=60,y=3.7,label="p = 5.5e-07",size=4)

#### Figure 2 S1G LDHA only ####
ggplot(LDH_ratio, aes(x=as.numeric(GIR), y=as.numeric(LDHA))) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("LDHA") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=14.2, label="tau = -0.31", color = "black",size=4) +
  annotate(geom = "text",x=60,y=14,label="p = 1.6e-04",size=4)


#### Figure 2 S1H LDHB only ####
ggplot(LDH_ratio, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("LDHB") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=10.5, label="tau = 0.37", color = "black",size=4) +
  annotate(geom = "text",x=60,y=10.2,label="p = 5.8e-06",size=4)


#### Figure 2 S2E - LDHA/B proportion ####
LDHA_prop <- 2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHA_P00338",Sample_ID_prot$Clamp == "Pre"])/(2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHA_P00338",Sample_ID_prot$Clamp == "Pre"])+2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHB_P07195",Sample_ID_prot$Clamp == "Pre"]))
LDHB_prop <- 2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHB_P07195",Sample_ID_prot$Clamp == "Pre"])/(2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHA_P00338",Sample_ID_prot$Clamp == "Pre"])+2^(Proteome_muscle2[rownames(Proteome_muscle2) == "LDHB_P07195",Sample_ID_prot$Clamp == "Pre"]))
LDH_data <- data.frame(LDHA_prop,LDHB_prop)
LDH_data$Individual <- 1:nrow(LDH_data)
LDH_data <- LDH_data[order(LDH_data$LDHB_prop),]

# Assuming your data frame is named df
# Reshape the data to long format
long_df <- pivot_longer(LDH_data, cols = c(LDHA_prop, LDHB_prop), names_to = "Enzyme", values_to = "Proportion")
long_df <- long_df[order(long_df$Proportion),]

LDH_data$Rank <- rank(-LDH_data$LDHB_prop)

long_df <- pivot_longer(LDH_data, cols = c(LDHA_prop, LDHB_prop), names_to = "Enzyme", values_to = "Proportion")

# Reorder Individuals based on LDHB_prop
long_df <- long_df %>%
  mutate(Individual = reorder(Individual, Proportion, FUN = function(x) -mean(x[Enzyme == "LDHB_prop"])))

long_df$Individual <- factor(long_df$Individual, levels = LDH_data$Individual[order(LDH_data$Rank)])

# Create the bar plot
ggplot(long_df, aes(x = Individual, y = Proportion, fill = Enzyme)) +
  geom_bar(stat = "identity", position = "stack", alpha=0.6) +
  theme_minimal() + theme(panel.grid = element_blank(),
                          axis.text.x = element_blank()) +
  labs(y = "LDHA/LDHB Proportion", x = "", fill = "Isoform") +
  scale_fill_manual(values = c("LDHA_prop" = "lightgrey", "LDHB_prop" = "darkred"))

#### Figure 2 S1A - HSPA2 insulin sensitivity ####
HSPA2 <- as.data.frame(cbind(pre_proteome_intensities[which(rownames(pre_proteome_intensities) == "HSPA2_P54652"),],target_proteome_df$M..µmol..kg.min..,target_proteome_df$Group))
colnames(HSPA2) <- c("Intensity","GIR","Group")             
cor.test(as.numeric(HSPA2$Intensity),as.numeric(HSPA2$GIR), method="kendall")

ggplot(HSPA2, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("HSPA2") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=10,y=10.7, label="tau = -0.46", color = "black",size=4) +
  annotate(geom = "text",x=10,y=10.5,label="p = 1.2e-08",size=4)

ggplot(HSPA2, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("HSPA2") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=10,y=10.7, label="tau = -0.46", color = "black",size=4) +
  annotate(geom = "text",x=10,y=10.5,label="p = 1.2e-08",size=4)

#### Figure 2 S1B - BDH1 insulin sensitivity ####
BDH1 <- as.data.frame(cbind(pre_proteome_intensities[which(rownames(pre_proteome_intensities) == "BDH1_Q02338"),],target_proteome_df$M..µmol..kg.min..,target_proteome_df$Group))
colnames(BDH1) <- c("Intensity","GIR","Group")             
cor.test(as.numeric(BDH1$Intensity),as.numeric(BDH1$GIR), method="kendall")

ggplot(BDH1, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("BDH1") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=50,y=6.1, label="tau = -0.46", color = "black",size=4) +
  annotate(geom = "text",x=50,y=5.9,label="p = 1.2e-08",size=4)

#### Figure 2B - KEGG pathways of correlated proteome ####
organism = "org.Hs.eg.db"

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = original_gene_list[names(original_gene_list) %in% dedup_ids$SYMBOL]

df2 <- as.data.frame(df2)
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$df2

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

selec_KEGG_prot_cor <- c("Ubiquitin mediated proteolysis","Motor proteins","Protein processing in endoplasmic reticulum",
                         "Adrenergic signaling in cardiomyocytes","Wnt signaling pathway",
                         "Citrate cycle (TCA cycle)","Fatty acid degradation","Proteasome",
                         "Lipoic acid metabolism","Oxidative phosphorylation",
                         "Retrograde endocannabinoid signaling","Propanoate metabolism")
require(viridis)

dotplot(kk2, showCategory = selec_KEGG_prot_cor, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign) + scale_color_gradient(low = "#56B1F7", high = "#132B43")
KEGG_Prot_cor <- kk2@result
prot_terms_to_plot <- KEGG_Prot_cor[which(KEGG_Prot_cor$Description %in% selec_KEGG_prot_cor),]
prot_terms_to_plot <- prot_terms_to_plot[order(prot_terms_to_plot$enrichmentScore),]
prot_terms_to_plot$rank <- 1:12

ggplot(prot_terms_to_plot, aes(x = 1, y = rank, color = p.adjust, label=Description)) +
  geom_point(aes(size=abs(enrichmentScore))) + geom_text_repel() +
  scale_color_gradient(low = "darkblue", high = "grey") +
  theme_minimal() +
  geom_hline(yintercept = 6.5) +
  theme(aspect.ratio=3/1,axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

#### Mitochondria summed protein abundance analysis ####
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- c(word(Disease_main_effect$ID,1,sep="_"))
annotations <- getBM(attributes = c('external_gene_name', 'go_id', 'name_1006'), 
                     filters = 'external_gene_name', 
                     values = genes, 
                     mart = mart)

# Filter for mitochondrial annotations (using GO:0005739 for mitochondrion)
mitochondrial_genes <- annotations$external_gene_name[annotations$go_id == "GO:0005739"]

Sum_int_data_frame_prot <- as.data.frame(Proteome_muscle2)
Sum_int_data_frame_prot$protein_name <- word(rownames(Proteome_muscle2),1,sep="_")

mito_df <- Sum_int_data_frame_prot[Sum_int_data_frame_prot$protein_name %in% mitochondrial_genes, ]


raw_mito_df <- 2^mito_df[,1:149]
summed_mito_intensities <- colSums(raw_mito_df, na.rm = TRUE)

Sum_int_data_frame_prot <- 2^Sum_int_data_frame_prot[,-ncol(Sum_int_data_frame_prot)]
summed_whole_intensities <- colSums(Sum_int_data_frame_prot, na.rm = TRUE)
proportions_mito_vs_whole <- summed_mito_intensities / summed_whole_intensities

Mito_plot <- long_data_sum[long_data_sum$Compartment == "Mitochondria",]
pre_Mito_plot <- Mito_plot[Mito_plot$Clamp == "Pre",]
pre_Mito_plot$Sum <- pre_Mito_plot$Sum*100

#### Figure 2C - Summed mitochondria abundance insulin sensitivity correlation ####
ggplot(pre_Mito_plot,aes(x=disease,y=Sum, color=disease)) + geom_violin(width=1) +
  geom_boxplot(width=0.1, alpha=0.2, position=position_dodge(width = 1.39)) +
  theme_bw() + scale_color_manual(values=c("darkgrey","darkred")) + 
  theme(aspect.ratio=1.25/1,legend.position = "none") + xlab("") +
  ylab("Sum mitochondrial \n protein abundance (%)")


#ploting mitochondrial protein abundance (proportion) correlation and wilcox.test
cor.test(new_test[new_test$Clamp == "Pre",]$Mitochondria,target_proteome_df$M..µmol..kg.min.., method="kendall")
pre_mito_proportion_data <- cbind(new_test[new_test$Clamp == "Pre",],target_proteome_df$M..µmol..kg.min..)

wilcox.test(new_test[new_test$disease == "NGT" & new_test$Clamp=="Pre",]$Mitochondria,
            new_test[new_test$disease == "T2D" & new_test$Clamp=="Pre",]$Mitochondria, alternative = "two.sided")

#### Figure 2D - mitochondria proportion insulin sensitivity ####
ggplot(pre_mito_proportion_data, aes(x=as.numeric(`target_proteome_df$M..µmol..kg.min..`), y=Mitochondria)) + geom_point(aes(color=disease,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  ggtitle("") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity")

#### Glycolysis summed abundance correlation ####
Sum_int_data_frame_prot <- as.data.frame(Proteome_muscle2)
Sum_int_data_frame_prot$protein_name <- word(word(rownames(Proteome_muscle2),1,sep="_"),1,sep=";")
Glycolysis_genes <- annotations$external_gene_name[annotations$go_id == "GO:0006096"] #TCA genes (27 genes)
Glycolysis_df <- Sum_int_data_frame_prot[Sum_int_data_frame_prot$protein_name %in% Glycolysis_genes, ]
non_glyc_proteins <- c("LDHA","PGM1","TIGAR","OGDH","PRKAG3")
Glycolysis_df <- Glycolysis_df[-c(which(Glycolysis_df$protein_name %in% non_glyc_proteins)),]
raw_Glycolysis_df <- 2^Glycolysis_df[,1:149]
summed_Glycolysis_intensities <- colSums(raw_Glycolysis_df, na.rm = TRUE)
Glycolysis_matrix <- cbind(summed_Glycolysis_intensities, Sample_ID_prot)
Glycolysis_matrix$proportion <- Glycolysis_matrix$summed_Glycolysis_intensities/summed_cytosol_cel_intensities
Glycolysis_matrix$proportion_whole <- Glycolysis_matrix$summed_Glycolysis_intensities/summed_whole_intensities
Glycolysis_matrix_data <- Glycolysis_matrix[Glycolysis_matrix$Clamp == "Pre",]
Glycolysis_matrix_data <- cbind(target_proteome_df$M..µmol..kg.min..,Glycolysis_matrix_data)
colnames(Glycolysis_matrix_data)[1] <- "GIR"
cor.test(log2(as.numeric(Glycolysis_matrix_data$GIR)),Glycolysis_matrix_data$proportion,method="kendall") # Proportion to cytosol proteome correlation
cor.test(log2(as.numeric(Glycolysis_matrix_data$GIR)),Glycolysis_matrix_data$proportion_whole,method="kendall") # Proportion to whole proteome correlation



#### Figure 2E heatmap/proportion of mito-complexes plot####
Sum_int_data_frame_prot <- as.data.frame(Proteome_muscle2)
Sum_int_data_frame_prot$protein_name <- word(rownames(Proteome_muscle2),1,sep="_")
complex_1_data <- read.delim("Complex_1_data.txt",header=T)
complex1_dataframe <- Sum_int_data_frame_prot[which(Sum_int_data_frame_prot$protein_name %in% complex_1_data$Approved.symbol),]
complex1_intensities <-2^complex1_dataframe[,-c(ncol(complex1_dataframe))]
summed_complex1_intensities <- colSums(complex1_intensities, na.rm = TRUE)
complex1_matrix <- cbind(summed_complex1_intensities, Sample_ID_prot)
complex1_matrix$proportion <- complex1_matrix$summed_complex1_intensities/summed_mito_intensities
complex1_matrix$proportion_whole <- complex1_matrix$summed_complex1_intensities/summed_whole_intensities
pre_complex1 <- complex1_matrix[complex1_matrix$Clamp == "Pre",]
pre_complex1 <- cbind(target_proteome_df$M..µmol..kg.min..,pre_complex1)
colnames(pre_complex1)[1] <- "GIR"

### complex2
complex_2_data <- read.delim("Complex_2_data.txt",header=T)
complex2_dataframe <- Sum_int_data_frame_prot[which(Sum_int_data_frame_prot$protein_name %in% complex_2_data$Approved.symbol),]
complex2_intensities <- 2^complex2_dataframe[,-c(ncol(complex2_dataframe))]
summed_complex2_intensities <- colSums(complex2_intensities, na.rm = TRUE)
complex2_matrix <- cbind(summed_complex2_intensities, Sample_ID_prot)
complex2_matrix$proportion <- complex2_matrix$summed_complex2_intensities/summed_mito_intensities
complex2_matrix$proportion_whole <- complex2_matrix$summed_complex2_intensities/summed_whole_intensities
pre_complex2 <- complex2_matrix[complex2_matrix$Clamp == "Pre",]
pre_complex2 <- cbind(target_proteome_df$M..µmol..kg.min..,pre_complex2)
colnames(pre_complex2)[1] <- "GIR"

### complex 3
complex_3_data <- read.delim("Complex_3_data.txt",header=T)
complex3_dataframe <- Sum_int_data_frame_prot[which(Sum_int_data_frame_prot$protein_name %in% complex_3_data$Approved.symbol),]
complex3_intensities <-2^complex3_dataframe[,-c(ncol(complex3_dataframe))]
summed_complex3_intensities <- colSums(complex3_intensities, na.rm = TRUE)
complex3_matrix <- cbind(summed_complex3_intensities, Sample_ID_prot)
complex3_matrix$proportion <- complex3_matrix$summed_complex3_intensities/summed_mito_intensities
complex3_matrix$proportion_whole <- complex3_matrix$summed_complex3_intensities/summed_whole_intensities
pre_complex3 <- complex3_matrix[complex3_matrix$Clamp == "Pre",]
pre_complex3 <- cbind(target_proteome_df$M..µmol..kg.min..,pre_complex3)
colnames(pre_complex3)[1] <- "GIR"

### complex 4
complex_4_data <- read.delim("Complex_4_data.txt",header=T)
complex4_dataframe <- Sum_int_data_frame_prot[which(Sum_int_data_frame_prot$protein_name %in% complex_4_data$Approved.symbol),]
complex4_intensities <- 2^complex4_dataframe[,-c(ncol(complex4_dataframe))]
summed_complex4_intensities <- colSums(complex4_intensities, na.rm = TRUE)
complex4_matrix <- cbind(summed_complex4_intensities, Sample_ID_prot)
complex4_matrix$proportion <- complex4_matrix$summed_complex4_intensities/summed_mito_intensities
complex4_matrix$proportion_whole <- complex4_matrix$summed_complex4_intensities/summed_whole_intensities
pre_complex4 <- complex4_matrix[complex4_matrix$Clamp == "Pre",]
pre_complex4 <- cbind(target_proteome_df$M..µmol..kg.min..,pre_complex4)
colnames(pre_complex4)[1] <- "GIR"

### complex 5
complex_5_data <- read.delim("Complex_5_data.txt",header=T)
complex5_dataframe <- Sum_int_data_frame_prot[which(Sum_int_data_frame_prot$protein_name %in% complex_5_data$Approved.symbol),]
complex5_intensities <-2^complex5_dataframe[,-c(ncol(complex5_dataframe))]
summed_complex5_intensities <- colSums(complex5_intensities, na.rm = TRUE)
complex5_matrix <- cbind(summed_complex5_intensities, Sample_ID_prot)
complex5_matrix$proportion <- complex5_matrix$summed_complex5_intensities/summed_mito_intensities
complex5_matrix$proportion_whole <- complex5_matrix$summed_complex5_intensities/summed_whole_intensities
pre_complex5 <- complex5_matrix[complex5_matrix$Clamp == "Pre",]
pre_complex5 <- cbind(target_proteome_df$M..µmol..kg.min..,pre_complex5)
colnames(pre_complex5)[1] <- "GIR"

est_complex <- matrix(ncol=2,nrow=5)
est_complex[,1] <- c(cor.test(log2(as.numeric(pre_complex1$GIR)),pre_complex1$proportion,method="kendall")$estimate[[1]],
                     cor.test(log2(as.numeric(pre_complex2$GIR)),pre_complex2$proportion,method="kendall")$estimate[[1]],
                     cor.test(log2(as.numeric(pre_complex3$GIR)),pre_complex3$proportion,method="kendall")$estimate[[1]],
                     cor.test(log2(as.numeric(pre_complex4$GIR)),pre_complex4$proportion,method="kendall")$estimate[[1]],
                     cor.test(log2(as.numeric(pre_complex5$GIR)),pre_complex5$proportion,method="kendall")$estimate[[1]])
est_complex[,2] <- c(cor.test(log2(as.numeric(pre_complex1$GIR)),pre_complex1$proportion,method="kendall")$p.value[[1]],
                     cor.test(log2(as.numeric(pre_complex2$GIR)),pre_complex2$proportion,method="kendall")$p.value[[1]],
                     cor.test(log2(as.numeric(pre_complex3$GIR)),pre_complex3$proportion,method="kendall")$p.value[[1]],
                     cor.test(log2(as.numeric(pre_complex4$GIR)),pre_complex4$proportion,method="kendall")$p.value[[1]],
                     cor.test(log2(as.numeric(pre_complex5$GIR)),pre_complex5$proportion,method="kendall")$p.value[[1]])

est_complex <- as.data.frame(est_complex)
colnames(est_complex) <- c("Estimate","p.value")
est_complex$complex <- c("I","II","III","IV","V")

Ratio_plot <- matrix(ncol=2,nrow=5)
Ratio_plot[,1] <- c(mean(pre_complex1$proportion)*100,
             mean(pre_complex2$proportion)*100,
             mean(pre_complex3$proportion)*100,
             mean(pre_complex4$proportion)*100,
             mean(pre_complex5$proportion)*100)
Ratio_plot[,2] <- c(sd(pre_complex1$proportion*100),
             sd(pre_complex2$proportion*100),
             sd(pre_complex3$proportion*100),
             sd(pre_complex4$proportion*100),
             sd(pre_complex5$proportion*100))

Ratio_plot <- as.data.frame(Ratio_plot)
colnames(Ratio_plot) <- c("mean","sd")
Ratio_plot$complex <- c("I","II","III","IV","V")

ggplot(Ratio_plot, aes(x=complex,y=mean)) +
  geom_crossbar(aes(ymin = mean, ymax = mean), width = 0.3,size=0.3,color="darkred") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,color="darkred") +
  theme_bw() + theme(aspect.ratio =0.3/1) + ylab("Proportion of mitochondria \n (%)")

# Plot
ggplot(est_complex, aes(x = complex, y = 1, fill = -log10(p.value))) +
  geom_tile(color = "black", size = 0.1) +
  geom_text(aes(label = sprintf("%.2f", Estimate))) +  # Format tau to 2 decimal places
  scale_fill_gradient(low = "white", high = "darkred", 
                      name = "-log10(p.value)",
                      breaks = pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(), aspect.ratio=0.1/1
  ) +
  labs(x = "Complex")


#### Figure 2G - plotting glycolysis, TCA, ECT together in binned GIR values ####
breaks_gir <- c(seq(0, 40, by = 10), Inf)
labels_gir <- c(seq(10, 40, by = 10), ">40")

# Bin the data
tca_summary <- TCA_matrix_data[,c(1,11)] %>% 
  mutate(GIR_bin = cut(GIR, breaks = breaks_gir, labels = labels_gir, include.lowest = TRUE)) %>% 
  group_by(GIR_bin) %>%
  summarise(
    mean_proportion = mean(proportion_whole),
    sem_proportion = sd(proportion_whole) / sqrt(n())
  )
tca_summary <- tca_summary[-c(nrow(tca_summary)),]

# Plotting
ggplot(tca_summary, aes(x = as.numeric(GIR_bin), y = mean_proportion)) +
  geom_line(aes(group = 1)) + 
  geom_point() +
  geom_errorbar(aes(ymin = mean_proportion - sem_proportion, ymax = mean_proportion + sem_proportion), width = 0.5) +
  labs(x = "Insulin Sensitivity (GIR)", y = "Proportion of Pathway")

# Repeat the binning and summarization for other dataframes
glycolysis_summary <- Glycolysis_matrix_data[,c(1,11)] %>% 
  mutate(GIR_bin = cut(GIR, breaks = breaks_gir, labels = labels_gir, include.lowest = TRUE)) %>% 
  group_by(GIR_bin) %>%
  summarise(
    mean_proportion = mean(proportion_whole),
    sem_proportion = sd(proportion_whole) / sqrt(n())
  ) %>% mutate(pathway = "Glycolysis")
glycolysis_summary <- glycolysis_summary[-c(nrow(glycolysis_summary)),]

electron_transport_summary <- all_ECT_matrix[,c(1,11)] %>% 
  mutate(GIR_bin = cut(GIR, breaks = breaks_gir, labels = labels_gir, include.lowest = TRUE)) %>% 
  group_by(GIR_bin) %>%
  summarise(
    mean_proportion = mean(proportion_whole),
    sem_proportion = sd(proportion_whole) / sqrt(n())
  ) %>% mutate(pathway = "Electron Transport")
electron_transport_summary <- electron_transport_summary[-c(nrow(electron_transport_summary)),]


all_data <- bind_rows(tca_summary %>% mutate(pathway = "TCA"), glycolysis_summary, electron_transport_summary)

# Plotting all together
ggplot(all_data, aes(x = as.numeric(GIR_bin), y = mean_proportion, color = pathway)) +
  geom_line(aes(group = pathway)) + 
  geom_point() +
  geom_errorbar(aes(ymin = mean_proportion - sem_proportion, ymax = mean_proportion + sem_proportion), width = 0.5) +
  labs(x = "Insulin Sensitivity (GIR)", y = "Proportion of Pathway") + theme_bw() + scale_color_manual(values=c("darkblue","darkred","darkgrey"))


#### Figure 2H - Glycolytic/oxidtive metabolism insulin sensitivity correlation ####
Oxidative_phosphorylaiton <- cbind(TCA_matrix_data$GIR,(TCA_matrix_data$summed_TCA_intensities+all_ECT_matrix$summed_all_electron_proteins_intensities)/summed_whole_intensities[which(TCA_matrix$Clamp == "Pre")])
Oxidative_phosphorylaiton <- as.data.frame(Oxidative_phosphorylaiton)
colnames(Oxidative_phosphorylaiton) <- c("GIR","proportion_whole")

cor.test(as.numeric(Oxidative_phosphorylaiton$GIR),Oxidative_phosphorylaiton$proportion_whole, method="kendall")

Oxidative_phosphorylaiton_summary <- Oxidative_phosphorylaiton %>% 
  mutate(GIR_bin = cut(GIR, breaks = breaks_gir, labels = labels_gir, include.lowest = TRUE)) %>% 
  group_by(GIR_bin) %>%
  summarise(
    mean_proportion = mean(proportion_whole),
    sem_proportion = sd(proportion_whole) / sqrt(n())
  ) %>% mutate(pathway = "Electron Transport")
Oxidative_phosphorylaiton_summary <- Oxidative_phosphorylaiton_summary[-c(nrow(Oxidative_phosphorylaiton_summary)),]

all_data <- bind_rows(glycolysis_summary %>% mutate(pathway = "TCA"), Oxidative_phosphorylaiton_summary)
ggplot(all_data, aes(x = as.numeric(GIR_bin), y = mean_proportion, color = pathway)) +
  geom_line(aes(group = pathway)) + 
  geom_point() +
  geom_errorbar(aes(ymin = mean_proportion - sem_proportion, ymax = mean_proportion + sem_proportion), width = 0.5) +
  labs(x = "Insulin Sensitivity (GIR)", y = "Proportion of Pathway") + theme_bw() + scale_color_manual(values=c("darkblue","darkred"))



cor.test(as.numeric(Oxidative_phosphorylaiton$GIR),Glycolysis_matrix_data$summed_Glycolysis_intensities/(TCA_matrix_data$summed_TCA_intensities+all_ECT_matrix$summed_all_electron_proteins_intensities), method="kendall")

Oxidative_phosphorylaiton$gly_ox_ratio <- Glycolysis_matrix_data$summed_Glycolysis_intensities/(TCA_matrix_data$summed_TCA_intensities+all_ECT_matrix$summed_all_electron_proteins_intensities)
Oxidative_phosphorylaiton$disease <- all_ECT_matrix$disease
cor.test(as.numeric(Oxidative_phosphorylaiton$GIR),Oxidative_phosphorylaiton$gly_ox_ratio, method="kendall")

ggplot(Oxidative_phosphorylaiton, aes(x=as.numeric(GIR), y=gly_ox_ratio)) + geom_point(aes(color=disease,stroke=NA),alpha=0.5, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgrey", "darkred")) +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Glycolytic/Oxidative metabolism") +
  annotate(geom="text",x=60,y=1, label="tau = -0.44", color = "black",size=4) +
  annotate(geom = "text",x=60,y=0.9,label="p = 5.8e-08",size=4) + scale_y_continuous(breaks = round(seq(0.3, 1.5, by = 0.2),1))








#### Figure 3J - glucose-into glycogen data ####
GU_glycogen_data <- read.delim("Glu_glycogen_data_HskM_for_R.txt", dec=",") # This clinical data won't be provided

GU_glycogen_summary <- GU_glycogen_data %>%
  group_by(Treatment, Condition) %>%
  summarise(
    mean_GU = mean(GU),
    SEM_GU = sd(GU) / sqrt(n())
  )

p <- ggplot(GU_glycogen_summary, aes(x = factor(Treatment, levels = c("Basal","Insulin")), y = mean_GU)) + 
  geom_col(aes(fill = Treatment), alpha = .2, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = mean_GU, ymax = mean_GU + SEM_GU),alpha=0.2, width = 0.2, position = position_dodge(width = 0.7)) +
  geom_line(data = GU_glycogen_data, aes(x = factor(Treatment, levels = c("Basal","Insulin")), y = GU, group = Donor, color=Donor),alpha=0.4, inherit.aes = FALSE,show.legend = FALSE) + 
  geom_point(data = GU_glycogen_data, aes(x = factor(Treatment, levels = c("Basal","Insulin")), y = GU, color = Donor), size = 4, alpha=0.4, inherit.aes = FALSE, show.legend = FALSE) + 
  facet_wrap(~ Condition, switch = "x", nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14, color="black"),
        axis.line.y = element_line(size = .5), aspect.ratio=2/1) + 
  ylab("nmol/mg protein/hr") +
  ggtitle("Glucose incorporation into glycogen") + 
  scale_fill_manual(values=c("lightgrey","black")) +
  scale_color_manual(values=c("blue","purple","black","darkgreen"))

p

#### Figure 3L - Hskm cells preprocessing and PCA ####
proteome_data <- read.delim("20240806_073639_Aurora25cm_KI_donor_deep_hskm_proteome_Report.tsv")
proteome_data <- proteome_data[,c(1:3,11:19,4:10)]
data_labels <- read.delim("Sample_label.txt")
colnames(proteome_data)[4:19] <- data_labels$Condition

log2_Hskm_proteome_data <- log2(proteome_data[,4:19]) #log2 transform data
rownames(log2_Hskm_proteome_data) <- proteome_data$PG.ProteinGroups
log2_Hskm_proteome_data_noNA <- na.omit(log2_Hskm_proteome_data) #remove all missing values


log2_Hskm_proteome_data_noNA <- medianScaling(log2_Hskm_proteome_data_noNA)

t_exprs <- as.data.frame(t(log2_Hskm_proteome_data_noNA))
pca_res <- prcomp(t_exprs, scale. = TRUE)

new_data_frame <- data.frame(pca_res$x)
new_data_frame$ID <- colnames(log2_Hskm_proteome_data_noNA)

new_data_frame$donor <- word(colnames(log2_Hskm_proteome_data_noNA),1,sep="_")
new_data_frame$kd <- word(colnames(log2_Hskm_proteome_data_noNA),2,sep="_")
new_data_frame$treatment <- word(colnames(log2_Hskm_proteome_data_noNA),3,sep="_")
new_data_frame$plate <- c("Plate_1","Plate_1","Plate_2","Plate_2","Plate_3","Plate_3","Plate_4","Plate_4",
                          "Plate_5","Plate_5","Plate_6","Plate_6","Plate_7","Plate_7","Plate_8","Plate_8")
new_data_frame$number <- 1:16
theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=18),
             legend.text = element_text(size = 16))

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(pca_res), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(new_data_frame,aes(x=PC1,y=PC2, color=donor)) +geom_point(size = 6, alpha = 0.4, aes(shape=kd), stroke=NA) + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values=c("blue","purple","black","darkgreen"))
p + stat_ellipse(level=0.8, alpha=0.2)


#### Differential regulation Hskm cells ####
log2_data_proteome_filtered <- selectGrps(log2_Hskm_proteome_data,
                                          word(colnames(log2_Hskm_proteome_data),2,3,sep="_"), 0.75,n=1)

# Median normalising data
log2_data_proteome_filtered <- medianScaling(log2_data_proteome_filtered)

Treatment <- factor(word(colnames(log2_Hskm_proteome_data),3,sep="_"), levels=c("basal","Insulin"))
kd <- factor(word(colnames(log2_Hskm_proteome_data),2,sep="_"), levels=c("Scr","AMPKy3"))
design <- model.matrix(~ kd * Treatment)
colnames(design) <- c("Intercept","KD","Treatment","KD_Treatment")
corfit <- duplicateCorrelation(log2_data_proteome_filtered, design, block=word(colnames(log2_Hskm_proteome_data),1,sep="_"))
fit <- lmFit(log2_data_proteome_filtered,design,block=word(colnames(log2_Hskm_proteome_data),1,sep="_"),correlation=corfit$consensus)
fit2 <- eBayes(fit)

proteome_KD_main_effect <- topTable(fit2, 
                                    coef=2, 
                                    number = Inf, sort.by = "none")
proteome_KD_main_effect$ID <- proteome_data[which(proteome_data$PG.ProteinGroups %in% rownames(proteome_KD_main_effect)),]$PG.Genes

proteome_Treatment_main_effect <- topTable(fit2, 
                                           coef=3, 
                                           number = Inf, sort.by = "none")

proteome_Treatment_main_effect$ID <- proteome_data[which(proteome_data$PG.ProteinGroups %in% rownames(proteome_KD_main_effect)),]$PG.Genes

proteome_Interaction_diet_treatment <- topTable(fit2, 
                                                coef=4, 
                                                number = Inf, sort.by = "none")

proteome_Interaction_diet_treatment$ID <- proteome_data[which(proteome_data$PG.ProteinGroups %in% rownames(proteome_KD_main_effect)),]$PG.Genes

For_export_hskm_prot <- cbind(proteome_KD_main_effect[,1:6],proteome_Treatment_main_effect[,1:7],rownames(proteome_Treatment_main_effect))

write.table(For_export_hskm_prot,"Hskm_cells_proteome.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### Figure 3K - AMPKy3 total abundance ####
Protein_plot <- log2_Hskm_proteome_data[rownames(log2_Hskm_proteome_data)=="Q9UGI9",]
sample_ID <- data_labels
Protein_plot_data <- data.frame(cbind(t(Protein_plot),sample_ID))
colnames(Protein_plot_data)[1] <- "Intensity"
Protein_plot_data$paired_ID <- factor(word(Protein_plot_data$Condition,1,sep="_"))
Protein_plot_data$treatment <- factor(word(Protein_plot_data$Condition,3,sep="_"), levels=c("basal","Insulin"))
Protein_plot_data$kd <- factor(word(Protein_plot_data$Condition,2,sep="_"), levels=c("Scr","AMPKy3"))

Protein_plot_data$Intensity <- as.numeric(Protein_plot_data$Intensity)

Protein_plot_data_summary <- Protein_plot_data %>%
  group_by(treatment, kd) %>%
  summarise(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    SEM_intensity = sd(Intensity, na.rm = TRUE) / sqrt(sum(!is.na(Intensity)))
  )

p <- ggplot(Protein_plot_data_summary, aes(x = factor(treatment, levels = c("basal","Insulin")), y = mean_intensity)) + 
  geom_col(aes(fill = treatment), alpha = .2, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = mean_intensity, ymax = mean_intensity + SEM_intensity),alpha=0.2, width = 0.2, position = position_dodge(width = 0.7)) +
  geom_line(data = Protein_plot_data, aes(x = factor(treatment, levels = c("basal","Insulin")), y = Intensity, group = paired_ID, color=paired_ID),alpha=0.4, inherit.aes = FALSE,show.legend = FALSE) + 
  geom_point(data = Protein_plot_data, aes(x = factor(treatment, levels = c("basal","Insulin")), y = Intensity, color = paired_ID), size = 4, alpha=0.4, inherit.aes = FALSE, show.legend = FALSE) + 
  facet_wrap(~ kd, switch = "x", nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14, color="black"),
        axis.line.y = element_line(size = .5), aspect.ratio=2/1) + 
  ylab("Log2 Intensity") +
  ggtitle("AMPKy3") + 
  scale_fill_manual(values=c("lightgrey","black")) +
  scale_color_manual(values=c("blue","purple","black","darkgreen"))

p

#### Figure 3M - KEGG pathway enrichment ####
genes <- proteome_KD_main_effect[proteome_KD_main_effect$adj.P.Val < 0.01,]$ID

ids4<-bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
ids4 = ids4[!duplicated(ids4[c("SYMBOL")]),]
universe_ID <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)

kegg_organism = "hsa"
kk <- enrichKEGG(gene=ids4$ENTREZID,
                 universe=universe_ID$ENTREZID,
                 organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

kk_results <- kk@result

Significant_KEGG_terms <- kk_results[kk_results$p.adjust < 0.05,]
Significant_KEGG_terms <- Significant_KEGG_terms[!(Significant_KEGG_terms$Description == "Longevity regulating pathway - multiple species"),]
Significant_KEGG_terms$Description <- factor(Significant_KEGG_terms$Description,
                                             levels=c("Insulin signaling pathway","mTOR signaling pathway",
                                                      "Carbohydrate digestion and absorption","Sphingolipid signaling pathway",
                                                      "Adipocytokine signaling pathway","T cell receptor signaling pathway",
                                                      "Insulin resistance","Autophagy - animal"))


ggplot(Significant_KEGG_terms, aes(y=Description, x=Count, fill=p.adjust)) + geom_col() +
  scale_fill_gradient(low="darkblue",high="lightgrey") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                            panel.grid.minor = element_blank())


#### Figure 5 S1A - sex-specific differences PCA ####
p<-ggplot(df_out,aes(x=PC1,y=PC5, colour = Gender)) + geom_point(shape=17, size = 4, alpha=0.8) + theme + xlab(percentage[1]) + ylab(percentage[5]) +
  scale_colour_manual(values=c("#3674b3","#bd626b")) + stat_ellipse()
p

#### Figure 5A - volcano proteome ####
ggplot() + geom_point(data=Gender_main_effect[which(Gender_main_effect$logFC > 0 & Gender_main_effect$adj.P.Val <= 0.05),], aes(x = logFC, y = -log10(P.Value), color = -log10(P.Value)),size=1.5) +
  scale_color_gradient(low = "#8ba0b5", high = "#3674b3") + 
  new_scale_color() + 
  geom_point(data = Gender_main_effect[which(Gender_main_effect$logFC < 0 & Gender_main_effect$adj.P.Val <= 0.05),], aes(x = logFC, y = -log10(P.Value), color = -log10(P.Value)),size=1.5) +
  scale_color_gradient(low = "#bda6a8", high = "#bd626b") +
  geom_point(data = Gender_main_effect[which(Gender_main_effect$adj.P.Val >0.05),], aes(x = logFC, y = -log10(P.Value)),color="grey", size=1, alpha=0.5, stroke=NA) +
  theme_bw() + theme(aspect.ratio=4/3,panel.grid = element_blank(),panel.border = element_blank()) + theme(legend.position = "none") + scale_x_continuous(limits = c(-0.75, 0.75), breaks = c(-0.75,-0.5, -0.25, 0, 0.25, 0.5,0.75)) +
  scale_y_continuous(limits = c(0, 13.5), breaks = c(2,4,6,8,10,12)) + geom_vline(xintercept = 0) 

#### Figure 5 S1K - Sex-specific glycolytiv/oxidative metabolism insulin sensitivity correlation ####
Oxidative_phosphorylaiton_gender <- cbind(Oxidative_phosphorylaiton,Sample_ID_prot[Sample_ID_prot$Clamp == "Pre",])

Oxidative_phosphorylaiton_gender$gly_ox_ratio <- Glycolysis_matrix_data$summed_Glycolysis_intensities/(TCA_matrix_data$summed_TCA_intensities+all_ECT_matrix$summed_all_electron_proteins_intensities)

ggplot(Oxidative_phosphorylaiton_gender, aes(x=as.numeric(GIR), y=gly_ox_ratio, group=Gender)) + geom_point(aes(color=Gender,stroke=NA),alpha=0.5, size=2.5) +
  geom_smooth(data = filter(Oxidative_phosphorylaiton_gender, Gender == "Male"), 
              method = "lm", se = TRUE, colour = "darkred", size = 0.75, alpha=0.05, fill="darkred") +
  geom_smooth(data = filter(Oxidative_phosphorylaiton_gender, Gender == "Female"), 
              method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.05, fill="darkblue") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),aspect.ratio=1/1.2, panel.grid = element_blank()) + scale_color_manual(values=c("darkblue", "darkred")) +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Glycolytic/Oxidative metabolism") +
  annotate(geom="text",x=60,y=1, label="tau = -0.46", color = "darkred",size=4) +
  annotate(geom = "text",x=60,y=0.9,label="p = 2.0e-05",size=4, color="darkred") + scale_y_continuous(breaks = round(seq(0.3, 1.5, by = 0.2),1)) +
  annotate(geom="text",x=20,y=0.5, label="tau = -0.46", color = "darkblue",size=4) +
  annotate(geom = "text",x=20,y=0.4,label="p = 3.4e-04",size=4, color = "darkblue") + scale_y_continuous(breaks = round(seq(0.3, 1.5, by = 0.2),1))
#### sex GSEA ####
original_gene_list <- Gender_main_effect$logFC
names(original_gene_list) <- sub(';.*$','', word(Gender_main_effect$ID,1,sep="_"))
original_gene_list = sort(original_gene_list, decreasing = TRUE)
original_gene_list <- original_gene_list[!duplicated(names(original_gene_list))]

set.seed(123)
gse_CC <- gseGO(geneList=original_gene_list,
                ont ="CC",
                keyType = "SYMBOL",
                nPerm = 10000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH",
                seed = T)

gse_CC_merged <- simplify(gse_CC)
dotplot(gse_CC, showCategory=20, split=".sign", font.size = 8, title = "T2D training (Post - Pre) - GOBP") + facet_grid(.~.sign) +
  scale_color_gradient(low = "#CC6666", high = "Dark grey")

gse_CC <- gseGO(geneList=original_gene_list,
                ont ="BP",
                keyType = "SYMBOL",
                nPerm = 10000,
                minGSSize = 3,
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                verbose = TRUE,
                OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH",
                seed = T)

dotplot(gse_CC, showCategory=20, split=".sign", font.size = 8, title = "T2D training (Post - Pre) - GOBP") + facet_grid(.~.sign) +
  scale_color_gradient(low = "#CC6666", high = "Dark grey")

gse_CC_merged <- simplify(gse_CC)

selected_GOBP <- c("neuron projection development","cytokine production","regulation of cytokine production involved in immune response","ventricular cardiac muscle tissue morphogenesis","protein localization to chromatin","nucleosome assembly",
                   "oxidative phosphorylation",
                   "glucose catabolic process",
                   "mitochondrial gene expression",
                   "proton motive force-driven ATP synthesis",
                   "NADH metabolic process",
                   "respiratory electron transport chain")

selected_GOCC <- c("NADH dehydrogenase complex",
                   "cytochrome complex", "mitochondrial inner membrane",
                   "mitochondrial ribosome","mitochondrial respiratory chain complex III","mitochondrial respiratory chain complex IV",
                   "complex of collagen trimers","nucleosome","axon terminus","ER to Golgi transport vesicle membrane","catalytic step 2 spliceosome","endocytic vesicle")

subset_GOBP$gender <- ""
subset_GOBP[which(subset_GOBP$enrichmentScore > 0),]$gender <- rep("Female",6)
subset_GOBP[which(subset_GOBP$enrichmentScore < 0),]$gender <- rep("Male",6)
subset_GOBP <- subset_GOBP[order(subset_GOBP$gender),]

subset_GOCC$gender <- ""
subset_GOCC[which(subset_GOCC$enrichmentScore > 0),]$gender <- rep("Female",6)
subset_GOCC[which(subset_GOCC$enrichmentScore < 0),]$gender <- rep("Male",6)
subset_GOCC <- subset_GOCC[order(subset_GOCC$gender),]

subset_GOBP$GO <- "GOBP"
subset_GOCC$GO <- "GOCC"
all_subset <- rbind(subset_GOBP,subset_GOCC)

ggplot(all_subset, aes(y=abs(enrichmentScore), x=reorder(Description, enrichmentScore), fill=p.adjust)) + geom_col() +
  theme_bw() + scale_fill_gradient(low="darkred",high="lightblue") +
  theme(aspect.ratio = 1/1.5) + coord_flip() + xlab("") + ylim(-1,1) +
  coord_polar(start = 0)


#### extracting lipid/glycogen proteins differentially expressed between sex ####
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

protein_info <- getBM(attributes = c('uniprotswissprot', 'external_gene_name', 'go_id'),
                      filters = 'uniprotswissprot',
                      values = Significant_proteins_gender$Prot_ID,
                      mart = ensembl)

lipid_go_terms <- c("GO:0006869", "GO:0008610","GO:0005811","GO:0019915")  # Add more relevant GO terms if necessary

# Check if the protein is involved in lipid-related processes
is_lipid_related <- protein_info[which(protein_info$go_id %in% lipid_go_terms),]

lipid_oxidation_terms <- c("GO:0003995","GO:0006635","GO:0009062","GO:0016042","GO:0044242","GO:0046355")  # Add more relevant GO terms if necessary
is_lipid_oxidation_related <- protein_info[which(protein_info$go_id %in% lipid_oxidation_terms),]

lipid_synthesis_terms <- c("GO:0006633","GO:0008610","GO:0072339","GO:0004312","GO:0030497","GO:0016053")
is_lipid_synthesis_related <- protein_info[which(protein_info$go_id %in% lipid_synthesis_terms),]

Glycogen_synthesis_terms <- c("GO:0005977")
is_Glycogen_synthesis_related <- protein_info[which(protein_info$go_id %in% Glycogen_synthesis_terms),]

Significant_proteins_gender$process <- "NO"
Significant_proteins_gender[which(Significant_proteins_gender$Prot_ID %in% unique(is_lipid_related$uniprotswissprot)),]$process <- "lipid_uptake"
Significant_proteins_gender[which(Significant_proteins_gender$Prot_ID %in% unique(is_lipid_oxidation_related$uniprotswissprot)[-1]),]$process <- "lipid_oxidation"
Significant_proteins_gender[which(Significant_proteins_gender$Prot_ID %in% unique(is_lipid_synthesis_related$uniprotswissprot)),]$process <- "lipid_synthesis"

Lipid_category <- Significant_proteins_gender[!Significant_proteins_gender$process == "NO",c(7,9)]
Lipid_category <- Lipid_category[-c(2,7),]
Lipid_names <- rownames(Lipid_category[order(Lipid_category$process),])
# Adding two additional lipogenesis genes LRG1, RTN3
Lipid_names <- c(Lipid_names,"LRG1_P02750","RTN3_O95197")
All_names <- c(Lipid_names, paste0(is_Glycogen_synthesis_related$external_gene_name,"_",is_Glycogen_synthesis_related$uniprotswissprot))

Single_gene <- list()
Multi_genes <- data.frame(1:149)
for(i in 1:length(All_names)){
  Single_gene <- Proteome_muscle2[rownames(Proteome_muscle2) == All_names[i],]
  Multi_genes <- cbind(Multi_genes,Single_gene)
}

Multi_genes <- Multi_genes[,-c(1)]
colnames(Multi_genes) <- word(All_names,1,sep="_")

plot_Multi_genes <- cbind(Multi_genes, Sample_ID_prot)
plot_Multi_genes <- plot_Multi_genes[,c(1:19,24:25)]
long_data <- gather(plot_Multi_genes, key = "Protein", value = "Intensity", -disease, -Gender, -Group)

summary_data <- long_data %>%
  group_by(Protein, Gender, disease) %>%
  summarise(
    Mean = mean(Intensity, na.rm = TRUE),  # Will calculate mean ignoring NAs
    SEM = ifelse(sum(!is.na(Intensity)) > 1, 
                 sd(Intensity, na.rm = TRUE) / sqrt(sum(!is.na(Intensity))),
                 NA)  # Will calculate SEM if there are more than one non-NA values, otherwise returns NA
  )

my_fill <- c("#a6cee3", "#c7b0af")
my_colors <- c("darkblue", "darkred")


ggplot(summary_data, aes(y = Protein, x = Mean, fill = Gender, color = disease)) +
  geom_bar(stat = "identity", position=position_dodge(0.9), width = 0.7, linewidth=0.25) +
  geom_errorbar(aes(xmin = Mean - SEM, xmax = Mean + SEM), width = 0.3, size = 0.1,
                position = position_dodge(width = 0.9)) +
  labs(x = "Log2 Intensity", y = NULL, fill = "Gender", pattern = "disease") +
  scale_fill_manual(values = my_fill) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        aspect.ratio=2/1) +
  scale_x_continuous(breaks = seq(5, 10, 1), minor_breaks = NULL) +
  coord_flip() +
  coord_cartesian(xlim = c(5.75, 14))

corrected_data <- summary_data %>%
  unite("Group", Gender, disease, sep = "_") %>%
  group_by(Protein, Group) %>%
  summarize(Mean = mean(Mean, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean)

new_All_names <- word(All_names,1,sep="_")
ordered_data <- corrected_data %>%
  arrange(match(Protein, new_All_names))

heatmap_matrix_corrected <- as.matrix(ordered_data[-1])
rownames(heatmap_matrix_corrected) <- ordered_data$Protein

heatmap_matrix_transposed <- t(heatmap_matrix_corrected)
heatmap(heatmap_matrix_transposed,
        Rowv = NA,
        Colv = NA,
        col = colorRampPalette(c("#bd626b", "white", "#3674b3"))(50),
        scale = "column",
        labRow = colnames(heatmap_matrix_corrected),
        labCol = rownames(heatmap_matrix_corrected),
        margins = c(47, 0),
        ColSideColors = c(rep("#bf999a",5),rep("#ab6163",2),rep("#5c080a",6),rep("darkblue",5))
)

legend(x="topleft")
legend(x="topleft",cex = 1.2, legend=c("High", "Medium", "Low"),fill=c("#3674b3", "white","#bd626b"))
?legend
legend(x="topleft",legend=c(-1:1),fill=colorRampPalette(c("#bd626b", "white", "#3674b3"))(50))
pheatmap(heatmap_matrix_transposed, scale = "column", cluster_rows="none", Colv=NA, de)


Significant_proteins_gender[Significant_proteins_gender$ID %in% Disease_main_effect[Disease_main_effect$adj.P.Val <= 0.05,]$ID,]
new_Significant_proteins_gender <- Significant_proteins_gender[which(Significant_proteins_gender$adj.P.Val <= 0.05),]

Sex_disease_overlap <- new_Significant_proteins_gender[new_Significant_proteins_gender[which(new_Significant_proteins_gender$adj.P.Val <= 0.05),]$ID %in% Prot_complete_cor_frame[Prot_complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05,]$ID,]
Sex_disease_overlap$col_group <- 1
Sex_disease_overlap[Sex_disease_overlap$logFC < 0 & Sex_disease_overlap$GIR_estimate < 0,]$col_group <- 2
Sex_disease_overlap[Sex_disease_overlap$logFC > 0 & Sex_disease_overlap$GIR_estimate < 0,]$col_group <- 3
Sex_disease_overlap$col_group <- factor(Sex_disease_overlap$col_group, levels=c(1,2,3))

Sex_disease_overlap$GIR_estimate <- Prot_complete_cor_frame[which(Prot_complete_cor_frame$ID %in% Sex_disease_overlap$ID),]$M..µmol..kg.min.._estimate
ggplot(Sex_disease_overlap,aes(x=logFC,y=GIR_estimate, label=ID, color=col_group)) + geom_point(size=4) + geom_text_repel(size =5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + theme_bw() + theme(panel.border = element_blank(), legend.position="none",aspect.ratio=1/1) +
  scale_color_manual(values=c("#ab6163","#7e7080","#3674b3"))
  
#### Validation cohort proteome processing ####
# Load proteomics and clinical data
KH_proteome_data <- read.delim("20231109_084853_20231031_KHvsl_SNv18.4_direct_Report.tsv")
KH_clinical_data <- read.delim("KH_validation_clinical_data_for_R.txt", dec = ",")
KH_clinical_data_additional <- read.delim("KH_validation_clinical_data_for_R_additional_data.txt", dec = ",")

# Add calculated clinical values
KH_clinical_data_additional <- KH_clinical_data_additional %>%
  mutate(Fast.insul_converted = round(Fast.insul / 6, 2),
         HOMA_IR = round(Fast.gluco * Fast.insul_converted / 22.5, 2))

# Merge clinical datasets
KH_clinical_data_combined <- KH_clinical_data %>%
  left_join(KH_clinical_data_additional, by = c("ID" = "PatientID")) %>%
  dplyr::select(-c(13:19))
colnames(KH_clinical_data_combined) <- sub("\\.x", "", colnames(KH_clinical_data_combined))

# Remove outliers and prepare sample IDs
KH_clinical_data_prot <- KH_clinical_data_combined[-c(34, 52), ]
KH_clinical_data_prot$SA_ID <- paste0(word(KH_clinical_data_prot$ID, 2, sep = "-"), "_", KH_clinical_data_prot$code_diab, "_", KH_clinical_data_prot$Clamp)

# Prepare proteomics matrix
KH_proteome_data <- KH_proteome_data[, -38]  # remove last metadata column
colnames(KH_proteome_data)[4:93] <- KH_clinical_data_prot$SA_ID
KH_val_log2_prot <- log2(KH_proteome_data[, 4:93])
complete_KH_val_prot <- na.omit(KH_val_log2_prot)

# PCA pre-correction
t_complete_KH_val_prot <- t(complete_KH_val_prot) %>% as.data.frame()
t_complete_KH_val_prot <- cbind(t_complete_KH_val_prot, KH_clinical_data_prot)
pca_res <- prcomp(t_complete_KH_val_prot[, 1:1139], scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>% mutate(batch = factor(KH_clinical_data_prot$batch))

# Variance explained
percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage_labels <- paste0("PC", 1:length(percentage), " (", percentage, "%)")

# PCA plot pre-correction
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2, alpha = 0.5) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = c("#9c4141", "#4c4794")) +
  xlab(percentage_labels[1]) +
  ylab(percentage_labels[2])
print(p)

# Batch correction and scaling
corrected_KH_val_prot <- removeBatchEffect(complete_KH_val_prot, batch = KH_clinical_data_prot$batch) %>%
  medianScaling()

# PCA post-correction
t_corrected_df <- t(corrected_KH_val_prot) %>% as.data.frame()
t_corrected_df <- cbind(t_corrected_df, KH_clinical_data_prot)
pca_res_corrected <- prcomp(t_corrected_df[, 1:1139], scale. = TRUE)
saveRDS(pca_res_corrected, "pca_res_prot_validation.rds")

pca_df <- as.data.frame(pca_res_corrected$x) %>%
  mutate(batch = factor(KH_clinical_data_prot$batch),
         disease = factor(KH_clinical_data_prot$code_diab),
         gender = factor(KH_clinical_data_prot$Gender),
         GIR = as.numeric(KH_clinical_data_prot$GIR),
         M.value = as.numeric(KH_clinical_data_prot$M.value),
         phenotype = factor(KH_clinical_data_prot$Phenotype))

percentage <- round(factoextra::get_eig(pca_res_corrected)$variance.percent, 2)
percentage_labels <- paste0("PC", 1:length(percentage), " (", percentage, "%)")

#### Figure 1 S2E Final PCA colored by M-value ####
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = M.value)) +
  geom_point(size = 4, alpha = 1) +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  xlab(percentage_labels[1]) +
  ylab(percentage_labels[2])
print(p2)


#### Differential Expression in Validation Cohort ####
# Set rownames and filter for proteins with at least 25% valid values
rownames(KH_val_log2_prot) <- KH_proteome_data$PG.ProteinGroups
Filt_KH_prot <- selectOverallPercent(KH_val_log2_prot, 0.25)

# Batch correction and normalization
Filt_KH_prot <- removeBatchEffect(Filt_KH_prot, batch = KH_clinical_data_prot$batch)
Filt_KH_prot <- medianScaling(Filt_KH_prot)

# Model setup
Clamp <- factor(KH_clinical_data_prot$Clamp, levels = c("Pre", "Post"))
Disease <- factor(KH_clinical_data_prot$code_diab, levels = c("con", "t2d"))
Gender <- factor(KH_clinical_data_prot$Gender, levels = c("Male", "Female"))
design <- model.matrix(~ Disease + Clamp + Gender)

# Linear model with duplicate correlation
corfit <- duplicateCorrelation(Filt_KH_prot, design, block = KH_clinical_data_prot$ID)
fit <- lmFit(Filt_KH_prot, design, block = KH_clinical_data_prot$ID, correlation = corfit$consensus)
fit2 <- eBayes(fit)

# Extract main effects
KH_Disease_main_effect <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
KH_Disease_main_effect$ID <- KH_proteome_data$PG.Genes[match(rownames(KH_Disease_main_effect), KH_proteome_data$PG.ProteinGroups)]

KH_Insulin_main_effect <- topTable(fit2, coef = 3, number = Inf, sort.by = "none")
KH_Insulin_main_effect$ID <- KH_proteome_data$PG.Genes[match(rownames(KH_Insulin_main_effect), KH_proteome_data$PG.ProteinGroups)]

KH_Gender_main_effect <- topTable(fit2, coef = 4, number = Inf, sort.by = "none")
KH_Gender_main_effect$ID <- KH_proteome_data$PG.Genes[match(rownames(KH_Gender_main_effect), KH_proteome_data$PG.ProteinGroups)]

#### Validation cohort Correlation with M-value at baseline (Pre) ####
KH_pre_prot_intensities <- Filt_KH_prot[, KH_clinical_data_prot$Clamp == "Pre"]
KH_pre_Sample_ID <- KH_clinical_data_prot[KH_clinical_data_prot$Clamp == "Pre", ]

baseline_p <- baseline_est <- numeric(nrow(KH_pre_prot_intensities))

for (i in seq_len(nrow(KH_pre_prot_intensities))) {
  res <- try(cor.test(KH_pre_prot_intensities[i, ], log2(KH_pre_Sample_ID$M.value), method = "kendall"), silent = TRUE)
  if (!inherits(res, "try-error")) {
    baseline_p[i] <- res$p.value
    baseline_est[i] <- res$estimate
  } else {
    baseline_p[i] <- NA
    baseline_est[i] <- NA
  }
}

adj_p <- p.adjust(baseline_p, method = "BH")

KH_prot_cor_frame <- data.frame(
  ID = rownames(KH_pre_prot_intensities),
  P.Value = baseline_p,
  Adj.P.Value = adj_p,
  Estimate = baseline_est,
  Gene_ID = KH_proteome_data$PG.Genes[match(rownames(KH_pre_prot_intensities), KH_proteome_data$PG.ProteinGroups)],
  Protein_name = KH_proteome_data$PG.ProteinDescriptions[match(rownames(KH_pre_prot_intensities), KH_proteome_data$PG.ProteinGroups)]
)

# Join with discovery cohort correlations (assuming already available as Prot_complete_cor_frame)
KH_prot_cor_frame$correct_ID <- word(KH_prot_cor_frame$ID, 1, sep = ";")
Prot_complete_cor_frame$correct_ID <- word(word(Prot_complete_cor_frame$ID, 2, sep = "_"), 1, sep = ";")

Both_cohort_correlation_IS <- inner_join(Prot_complete_cor_frame, KH_prot_cor_frame, by = "correct_ID")
Both_cohort_correlation_IS_sign <- Both_cohort_correlation_IS %>%
  filter(M..µmol..kg.min.._adj.p.value <= 0.05) %>%
  mutate(label = word(word(ID.x, 1, sep = "_"), 1, sep = ";"))

# Plot


#### Figure 2 S1C - Cohort Insulin Sensitivity correlation proteome  ####
ggplot(Both_cohort_correlation_IS_sign, aes(x = M..µmol..kg.min.._estimate, y = Estimate, label = label, color = abs(Estimate))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel() +
  theme_minimal() +
  scale_color_gradient(low = "grey", high = "darkred") +
  labs(
    title = "Protein insulin sensitivity associations",
    x = "Discovery cohort\n (coefficient)",
    y = "Validation cohort\n (coefficient)"
  ) +
  xlim(-0.5, 0.5) + ylim(-0.5, 0.5) +
  annotate("text", x = -0.25, y = 0.45, label = "tau = 0.47", size = 4) +
  annotate("text", x = -0.25, y = 0.38, label = "p = 2.1e-15", size = 4)


#### Figure 5 S1C - correlate FC between cohorts ####
Gender_overlap_datasets_proteome <-Gender_main_effect[which(word(KH_Gender_main_effect$ID,1,sep=";") %in% word(word(Gender_main_effect$ID,1,sep=";"),1,sep="_")),]
Gender_overlap_datasets_proteome$ID <- word(word(Gender_overlap_datasets_proteome$ID,1,sep=";"),1,sep="_")
KH_Gender_main_effect$ID <- word(KH_Gender_main_effect$ID,1,sep=";")

Gender_overlap_datasets_proteome <- inner_join(Gender_overlap_datasets_proteome,KH_Gender_main_effect,by="ID")
sign_Gender_overlap_datasets_proteome <- Gender_overlap_datasets_proteome[Gender_overlap_datasets_proteome$adj.P.Val.x <= 0.05 | Gender_overlap_datasets_proteome$adj.P.Val.y <= 0.05,]

ggplot(sign_Gender_overlap_datasets_proteome, aes(x=logFC.x, y=logFC.y, label=ID)) + geom_point(alpha=0.4, size=4, shape=1,color="darkred", stroke=0.4) +
  #geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkred", "darkgrey")) +
  ggtitle("Sex specific proteome") + geom_text_repel() +
  xlab("Discovery cohort \n Log2 FC (Female-Male") + ylab("Validation cohort \n Log2 FC (Female-Male) ") + xlim(-0.6,0.6) + ylim(-0.6,1) +
  annotate(geom="text",x=-0.3,y=0.6, label="tau = 0.60", color = "black",size=4) +
  annotate(geom = "text",x=-0.3,y=0.5,label="P = 4.0e-16",size=4)

cor.test(sign_Gender_overlap_datasets_proteome$logFC.x,sign_Gender_overlap_datasets_proteome$logFC.y, method="kendall")

nrow(sign_Gender_overlap_datasets_proteome[sign_Gender_overlap_datasets_proteome$logFC.x < 0 & sign_Gender_overlap_datasets_proteome$logFC.y < 0,]) #49
nrow(sign_Gender_overlap_datasets_proteome[sign_Gender_overlap_datasets_proteome$logFC.x > 0 & sign_Gender_overlap_datasets_proteome$logFC.y > 0,]) #49
(49+29)/85*100

#### Figure 1 S2H & 1 S2J - Clinical MDS plots ####
target_proteome_df2 <- target_proteome_df %>%
  dplyr::rename(
    BMI   = `BMI..kg.m2.`,
    Mval  = `M..µmol..kg.min..`,
    hba1c = `HbA1c..mmol.mol.`,
    FPG   = `FP.Glucose..mmol.L.`,
    TG    = `P.TG..mmol.L.`,
    HDL   = `HDL.Chol..mmol.L.`,
    Age   = `Age`
  )
# Define the clinical variables
clin_vars <- c("BMI", "Mval", "hba1c", "Age", "FPG", "TG", "HDL")

# Subset by group FROM THE RENAMED DATAFRAME
df_NGT <- target_proteome_df2 %>% filter(Group == "NGT")
df_T2D <- target_proteome_df2 %>% filter(Group == "T2D")

# Subset and scale clinical variables (column-wise Z-scaling)
df_out_NGT <- df_NGT[, clin_vars] %>% scale() %>% as.data.frame()
df_out_T2D <- df_T2D[, clin_vars] %>% scale() %>% as.data.frame()


# Compute distance and MDS separately
dist_NGT <- dist(df_out_NGT)
dist_T2D <- dist(df_out_T2D)

mds_df_NGT <- as.data.frame(cmdscale(dist_NGT, k = 2))
mds_df_T2D <- as.data.frame(cmdscale(dist_T2D, k = 2))

colnames(mds_df_NGT) <- c("Dim1", "Dim2")
colnames(mds_df_T2D) <- c("Dim1", "Dim2")

# Add clinical data back for correlation later
mds_df_NGT <- bind_cols(mds_df_NGT, df_out_NGT)
mds_df_T2D <- bind_cols(mds_df_T2D, df_out_T2D)


# Correlations
cor_df_NGT <- df_out_NGT %>%
  summarise(across(everything(), list(
    Dim1 = ~ cor(.x, mds_df_NGT$Dim1, method = "pearson"),
    Dim2 = ~ cor(.x, mds_df_NGT$Dim2, method = "pearson")
  ))) %>%
  pivot_longer(everything(), names_to = "key", values_to = "value") %>%
  separate(key, into = c("Variable", "Dim"), sep = "_Dim") %>%
  pivot_wider(names_from = Dim, values_from = value)

cor_df_T2D <- df_out_T2D %>%
  summarise(across(everything(), list(
    Dim1 = ~ cor(.x, mds_df_T2D$Dim1, method = "pearson"),
    Dim2 = ~ cor(.x, mds_df_T2D$Dim2, method = "pearson")
  ))) %>%
  pivot_longer(everything(), names_to = "key", values_to = "value") %>%
  separate(key, into = c("Variable", "Dim"), sep = "_Dim") %>%
  pivot_wider(names_from = Dim, values_from = value)

cor_df_NGT <- cor_df_NGT %>%
  dplyr::rename(Dim1 = `1`, Dim2 = `2`)
cor_df_T2D <- cor_df_T2D %>%
  dplyr::rename(Dim1 = `1`, Dim2 = `2`)


# Scale arrows
scaling_factor <- 2
cor_df_NGT <- cor_df_NGT %>% mutate(Dim1 = Dim1 * scaling_factor, Dim2 = Dim2 * scaling_factor)
cor_df_T2D <- cor_df_T2D %>% mutate(Dim1 = Dim1 * scaling_factor, Dim2 = Dim2 * scaling_factor)


# NGT Plot
plot_NGT <- ggplot(mds_df_NGT, aes(x = Dim1, y = Dim2)) +
  geom_point(size = 3, color = "darkblue", alpha = 0.2) +
  geom_segment(data = cor_df_NGT,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = cor_df_NGT, aes(x = Dim1, y = Dim2, label = Variable),
            size = 5, vjust = -0.5, hjust = 0.5) +
  labs(title = "MDS plot - NGT", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

mds_df_T2D$Dim2 <- -mds_df_T2D$Dim2
cor_df_T2D$Dim2 <- -cor_df_T2D$Dim2

# T2D Plot
plot_T2D <- ggplot(mds_df_T2D, aes(x = Dim1, y = Dim2)) +
  geom_point(size = 3, color = "darkred", alpha = 0.2) +
  geom_segment(data = cor_df_T2D,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = cor_df_T2D, aes(x = Dim1, y = Dim2, label = Variable),
            size = 5, vjust = -0.5, hjust = 0.5) +
  labs(title = "MDS plot - T2D", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")



#### Figure 1 S2I & 1 S2K - MDS clinical x proteome correlation GeneSet ####
# Function to compute spearman correlation to MDS dims
get_tau_correlations <- function(proteome_matrix, mds_df) {
  cor_df <- apply(proteome_matrix, 1, function(protein) {
    c(Tau_PC1 = cor(protein, mds_df$Dim1, method = "spearman", use = "pairwise.complete.obs"),
      Tau_PC2 = cor(protein, mds_df$Dim2, method = "spearman", use = "pairwise.complete.obs"))
  })
  as.data.frame(t(cor_df))
}

results_NGT <- get_tau_correlations(proteome_NGT, mds_df_NGT)
results_T2D <- get_tau_correlations(proteome_T2D, mds_df_T2D)

results_NGT$ID <- rownames(results_NGT)
results_T2D$ID <- rownames(results_T2D)

results_NGT$Gene <- word(word(results_NGT$ID,1,sep="_"),1,sep=";")
results_T2D$Gene <- word(word(results_T2D$ID,1,sep="_"),1,sep=";")

# Use SYMBOLs and Tau_PC1
gene_list_PC1 <- results_NGT %>%
  filter(!duplicated(Gene)) %>%
  drop_na(Tau_PC1) %>%
  arrange(desc(Tau_PC1)) %>%
  pull(Tau_PC1, name = Gene)

gene_list_PC2 <- results_NGT %>%
  filter(!duplicated(Gene)) %>%
  drop_na(Tau_PC2) %>%
  arrange(desc(Tau_PC2)) %>%
  pull(Tau_PC2, name = Gene)

set.seed(123)
gse_NGT_PC1 <- gseGO(geneList = gene_list_PC1,
                     ont = "BP",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 1,
                     verbose = FALSE,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = "BH")

gse_NGT_PC2 <- gseGO(geneList = gene_list_PC2,
                     ont = "BP",
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 1,
                     verbose = FALSE,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = "BH")

GOBP_PC1_NGT <- gse_NGT_PC1@result
GOBP_PC2_NGT <- gse_NGT_PC2@result

gene_list_PC1_T2D <- results_T2D %>%
  filter(!duplicated(Gene)) %>%
  drop_na(Tau_PC1) %>%
  arrange(desc(Tau_PC1)) %>%
  pull(Tau_PC1, name = Gene)

gene_list_PC2_T2D <- results_T2D %>%
  filter(!duplicated(Gene)) %>%
  drop_na(Tau_PC2) %>%
  arrange(desc(Tau_PC2)) %>%
  pull(Tau_PC2, name = Gene)

set.seed(123)
gse_T2D_PC1 <- gseGO(geneList = gene_list_PC1_T2D, ont = "BP", keyType = "SYMBOL",
                     nPerm = 10000, minGSSize = 3, maxGSSize = 800,
                     pvalueCutoff = 1, verbose = FALSE,
                     OrgDb = org.Hs.eg.db, pAdjustMethod = "BH")

gse_T2D_PC2 <- gseGO(geneList = gene_list_PC2_T2D, ont = "BP", keyType = "SYMBOL",
                     nPerm = 10000, minGSSize = 3, maxGSSize = 800,
                     pvalueCutoff = 1, verbose = FALSE,
                     OrgDb = org.Hs.eg.db, pAdjustMethod = "BH")

GOBP_PC1_T2D <- gse_T2D_PC1@result
GOBP_PC2_T2D <- gse_T2D_PC2@result


# Keep relevant GSEA columns
full_PC1_NGT <- GOBP_PC1_NGT %>%
  dplyr::select(ID, Description, ES_PC1 = enrichmentScore, p.adjust_PC1 = p.adjust)

full_PC2_NGT <- GOBP_PC2_NGT %>%
  dplyr::select(ID, Description, ES_PC2 = enrichmentScore, p.adjust_PC2 = p.adjust)

# Merge and annotate
merged_results_NGT <- full_join(full_PC1_NGT, full_PC2_NGT, by = c("ID", "Description")) %>%
  filter(p.adjust_PC1 < 0.05 | p.adjust_PC2 < 0.05) %>%
  mutate(Significance = case_when(
    p.adjust_PC1 < 0.05 & p.adjust_PC2 < 0.05 ~ "Both",
    p.adjust_PC1 < 0.05 ~ "PC1 Only",
    p.adjust_PC2 < 0.05 ~ "PC2 Only"
  ),
  highlight = ifelse(Description %in% c(
    "extracellular matrix organization", 
    "cellular respiration",
    "adaptive immune response",
    "ribosome biogenesis",
    "immunoglobulin mediated immune response"
  ), Description, NA))

plot_filtered_results_NGT <- ggplot(merged_results_NGT, aes(x = ES_PC1, y = ES_PC2)) +
  geom_point(size = 2, alpha = 0.3, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel(
    aes(label = highlight),
    size = 5, 
    color = "blue", 
    na.rm = TRUE, 
    segment.color = "black",
    segment.size = 0.5,
    box.padding = 0.3,
    point.padding = 0.5,
    nudge_y = 0.3
  ) +
  labs(
    title = "NGT",
    x = "Enrichment Score (PC1)",
    y = "Enrichment Score (PC2)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, color = "black"),
    aspect.ratio = 1
  )



# Prepare and merge
full_PC1_T2D <- GOBP_PC1_T2D %>%
  dplyr::select(ID, Description, ES_PC1 = enrichmentScore, p.adjust_PC1 = p.adjust)

full_PC2_T2D <- GOBP_PC2_T2D %>%
  dplyr::select(ID, Description, ES_PC2 = enrichmentScore, p.adjust_PC2 = p.adjust)

merged_results_T2D <- full_join(full_PC1_T2D, full_PC2_T2D, by = c("ID", "Description")) %>%
  filter(p.adjust_PC1 < 0.05 | p.adjust_PC2 < 0.05) %>%
  mutate(Significance = case_when(
    p.adjust_PC1 < 0.05 & p.adjust_PC2 < 0.05 ~ "Both",
    p.adjust_PC1 < 0.05 ~ "PC1 Only",
    p.adjust_PC2 < 0.05 ~ "PC2 Only"
  ),
  highlight = ifelse(Description %in% c(
    "branched-chain amino acid catabolic process", 
    "cellular respiration",
    "'de novo' protein folding",
    "telomere maintenance",
    "regulation of proteasomal ubiquitin-dependent protein catabolic process",
    "positive regulation of DNA replication",
    "macrophage derived foam cell differentiation",
    "response to unfolded protein"
  ), Description, NA))



plot_filtered_results_T2D <- ggplot(merged_results_T2D, aes(x = ES_PC1, y = ES_PC2)) +
  geom_point(size = 2, alpha = 0.3, color = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel(
    aes(label = highlight),
    size = 5, 
    color = "darkred", 
    na.rm = TRUE, 
    segment.color = "black",
    segment.size = 0.5,
    box.padding = 0.3,
    point.padding = 0.5,
    nudge_y = 0.3
  ) +
  labs(
    title = "T2D",
    x = "Enrichment Score (PC1)",
    y = "Enrichment Score (PC2)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, color = "black"),
    aspect.ratio = 1
  )

plot_filtered_results_NGT + plot_filtered_results_T2D
