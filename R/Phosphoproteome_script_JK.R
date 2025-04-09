#### KI phosphoproteome part ####
---
  title: "KI Phospho project"
author: "Jeppe Kjærgaard"
date: "2023-03-10"
output: html_document
---
  
#### Load libraries ####
library(biomaRt)
library(caret)
library(circlize)
library(cluster)
library(ComplexHeatmap)
library(dplyr)
library(factoextra)
library(fmsb)
library(GenomicRanges)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(karyoploteR)
library(limma)
library(mgcv)
library(pheatmap)
library(PhosR)
library(pROC)
library(protr)
library(RColorBrewer)
library(randomForest)
library(stats)
library(stringr)
library(tidyverse)
library(VennDiagram)

#### Loading files ####
#Clinical data loading - IS NOT PROVIDED
Sample_ID <- read.delim("Input_SA_info_R.txt")
Sample_ID$Gender[grepl("Male",Sample_ID$Group)] <- "Male"
Sample_ID$Gender[grepl("Female",Sample_ID$Group)] <- "Female"
Sample_ID$disease[grepl("T2D",Sample_ID$Group)] <- "T2D"
Sample_ID$disease[grepl("NGT",Sample_ID$Group)] <- "NGT"
Sample_ID$Group <- factor(Sample_ID$Group)

Clinical_data <- read.delim("Clinical data/Combined_clinical_KI_data.txt", dec=",")
Sample_ID$Subject_ID <- sub("IRS","IRS1",Sample_ID$Subject_ID)
Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% Sample_ID$Subject_ID,]
colnames(Our_clin_data)[1] <- "Subject_ID"
Our_clin_data <- Our_clin_data[order(match(Our_clin_data[,1],Sample_ID[,5])),]
duplicated_Our_clin_data <- Our_clin_data %>%
  uncount(weights = 2, .remove = FALSE)

pre_Sample_ID <- Sample_ID[which(Sample_ID$Clamp == "Pre"),]
Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% pre_Sample_ID$Subject_ID,]
target_df <- Our_clin_data[match(unique(pre_Sample_ID$Subject_ID),Our_clin_data$Subject.ID),]

#Phosphoproteome - site table
phospho_study <-read.delim("20230612_070600_20230525_directDIA_hybrid_KI_Phospho_collapsed.txt")
# Phosphoproteome - peptide table
phospho_study2 <- read.delim("20230612_070600_20230525_directDIA_hybrid_KI_Phospho_Report.txt")

#### Figure 1B - Clinical characterizationn ####
long_target_df <- gather(data=target_df, key="number", value="Value",-c(Gender,Group))
GIR_data <- long_target_df[long_target_df$number == "M..µmol..kg.min..",]
GIR_data$Value <- as.numeric(GIR_data$Value)

ggplot(GIR_data, aes(y=Value,x=Group, fill=Group, color=Group, Group=Gender)) + geom_boxplot(alpha = .2) + geom_point(size=3,alpha=0.75, stroke=NA,position = position_dodge(width = .75)) +
  ylab("umol/kg/min)") + theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 14)) +  scale_fill_manual(values = c("darkgrey","darkred")) + scale_color_manual(values=c("#4d4d4d","darkred"))

#### Figure 1C - Clinical characterizationn ####

data_GIR_plot <- Our_clin_data[order(Our_clin_data$M..µmol..kg.min..),]
data_GIR_plot$rank <- 1:77
ggplot(data=data_GIR_plot,aes(y=M..µmol..kg.min..,x=rank, fill=Group)) +
  geom_col() + theme_bw() + theme(aspect.ratio=1/1,axis.text.x = element_blank()) + 
  ylab("umol/kg/min") + xlab("") + scale_fill_manual(values=c("grey","#CC6666"))

#### Figure 1E - Number of identified sites/peptides/proteins ####
phospho_study2 <- phospho_study2[grepl("Phospho",phospho_study2$EG.PrecursorId),] ### filter for phosphopeptides
### number of phosphopeptides identified
length(unique(phospho_study2$EG.PrecursorId))
length(unique(phospho_study$PG.ProteinGroups)) ## number of phosphoproteins

#phosphosites found in at least 5 samples
total_SA <- rep(1,154)
rownames(phospho_study) <- phospho_study$PTM_collapse_key
Identified_sites <- selectGrps(phospho_study[,c(1:154)], total_SA, 0.03, n=1)
tyrosine_sites <- length(which(grepl("_Y",rownames(Identified_sites)))) ## number of phosphoserine
threonine_sites <- length(which(grepl("_T",rownames(Identified_sites)))) ## number of phosphoserine
serine_sites <- length(which(grepl("_S",rownames(Identified_sites)))) ## number of phosphoserine

#phosphopeptides found in at least 5 samples
required_samples <- 5
total_samples <- 154
# Group the data by peptide and count the unique sample IDs for each peptide
peptide_counts <- aggregate(R.Replicate ~ EG.PrecursorId, data = phospho_study2, FUN = function(x) length(unique(x)))
num_peptides <- sum(peptide_counts$R.Replicate >= required_samples)
print(num_peptides) # number of phosphopeptides

#barplot showing distribution of sites on S/T/Y
site_data <- data.frame(numbers = c(tyrosine_sites, threonine_sites, serine_sites),
                        residue = c("tyrosine","threonine","serine"),
                        percentage = c(paste0(round(tyrosine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%"),
                                       paste0(round(threonine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%"),
                                       paste0(round(serine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%")))

ggplot(site_data, aes(x=residue,y=numbers)) + geom_col(fill="#e39aa2") +
  geom_text(aes(label = percentage), vjust = -0.5, size=6) + theme_minimal() +
  theme(aspect.ratio=1/3, panel.grid = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black",size=14)) + xlab("") + coord_flip()

#### Batch-effect control ####
# In this section we investigate the data by a scaled PCA plot. Clearly there is an expected day-batch effect
log2_phospho_study <- log2(phospho_study[,1:154]) #log2 transform data
rownames(log2_phospho_study) <- phospho_study$PTM_collapse_key
new_study <- na.omit(log2_phospho_study) #remove all missing values
new_colnames <- word(colnames(new_study),10, sep="_")
new_colnames <- as.integer(sub("S","",new_colnames))
colnames(new_study) <- new_colnames
colnames(new_study) <- as.integer(colnames(new_study))
new_study <- new_study[,str_sort(colnames(new_study), numeric = TRUE)]

#defining known pre-analytical batch-effects
day_batch_effect <- factor(c(rep(1, 80),rep(2, 74))) 
c18_batch_effect <- factor(c(rep(1, 40),rep(2, 40),rep(3,40),rep(4,34))) 
evotip_batch_effect <- rep(1,154)
evotip_batch_effect[c(104,105,116,117,121,123,129,132,133,134,136,137,139,140,1:12,45,46,81:93,126)] <- 2
evotip_batch_effect  <- factor(evotip_batch_effect)

#preparing data for PC-analysis
t_exprs <- as.data.frame(t(new_study))
t_exprs <- cbind(t_exprs, day_batch_effect)
pca_res <- prcomp(t_exprs[,1:1589], scale. = TRUE)

## visualising and defining three batch clusters
new_data_frame <- data.frame(pca_res$x)
new_data_frame$day <- day_batch_effect

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


p<-ggplot(new_data_frame,aes(x=PC1,y=PC2, color = day))
p<-p+geom_point(size = 2, alpha = 0.5) + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  stat_ellipse(level=0.8)
p + scale_color_manual(values=c("#9c4141", "#4c4794", "#56B4E9","pink"))

#Conclusion: obvious day batch effect can be seen


#We will now correct for this batch effect
#removing batch effect
corrected_phos <- removeBatchEffect(new_study, day_batch_effect)
t_exprs <- as.data.frame(t(corrected_phos))
t_exprs <- cbind(t_exprs, day_batch_effect)
pca_res <- prcomp(t_exprs[,1:1589], scale. = TRUE)
write_rds(pca_res,"pca_res_phos.rds")

#adding grouping data to explore distribution on PCA-plot
new_data_frame <- data.frame(pca_res$x)
new_data_frame$day <- day_batch_effect
new_data_frame$group <- Sample_ID$Group
new_data_frame$clamp <- Sample_ID$Clamp
new_data_frame$gender <- Sample_ID$Gender
new_data_frame$disease <- Sample_ID$disease

#Clinical data
new_data_frame$HOMA1_IR <- duplicated_Our_clin_data$HOMA1.IR
# for HOMA1.IR PCA there is one outlier subject which needs to be taken out.
new_data_frame$GIR <- duplicated_Our_clin_data$M..µmol..kg.min..

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


#### Figure 1I - PCA colored by M-value ####
p<-ggplot(new_data_frame,aes(x=PC1,y=PC2, color = GIR))
p<-p+geom_point(size = 2, alpha = 0.8) + theme + xlab(percentage[1]) + ylab(percentage[2]) + scale_colour_gradient(low = "lightblue",high = "darkred", space = "Lab", na.value = "grey50",guide = "colourbar", aesthetics = "colour")
p

#### Figure 1 S2D ####
p<-ggplot(new_data_frame,aes(x=PC1,y=PC2, colour = disease)) + geom_point(size = 3, alpha=0.8) + 
  theme + xlab(percentage[1]) + ylab(percentage[2]) +
  scale_colour_manual(values=c("darkgrey","darkred"))
p + stat_ellipse(level=0.8)


#### Figure 5 S1B ####
p<-ggplot(new_data_frame,aes(x=PC4,y=PC5, colour = gender, fill=gender)) +geom_point(shape=24,size = 4, alpha=0.8) + 
  theme + xlab(percentage[4]) + ylab(percentage[5]) +
  scale_colour_manual(values=c("#3674b3","#bd626b")) +
  scale_fill_manual(values=c("#3674b3","#bd626b"))
p + stat_ellipse(level=0.8)


#### Figure 1F - New CV figure ####
corrected_prot <- readRDS(file = "corrected_prot.rds")
Sample_ID_prot <- readRDS(file="Sample_ID_prot.rds")
DIA_targeted_data <- read.delim("R1326_directDIA_targeted/20240611_144716_DIAphospho_Targeted_spikein_KI_revision_SNv19_collapsed.txt")
DIA_targeted_data_intensity <- DIA_targeted_data[,1:15]
rownames(DIA_targeted_data_intensity) <- DIA_targeted_data$PTM_collapse_key
DIA_targeted_data_intensity <- log2(DIA_targeted_data_intensity)
DIA_targeted_data_intensity <- na.omit(DIA_targeted_data_intensity)
DIA_targeted_data_intensity <- DIA_targeted_data_intensity[which(rownames(DIA_targeted_data_intensity) %in% rownames(corrected_phos)),] # Extracting sites in the new data set also identified in original data for fair comparison. 85% of the sites with complete quantification identified in the new experiment
DIA_targeted_data_intensity <- removeBatchEffect(DIA_targeted_data_intensity, batch=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)) #Removing potential batch effect between different spike-in ratios.

# Plotting CV distribution for each protein
cv_values_1 <- apply(2^corrected_prot[,Sample_ID_prot$Clamp == "Pre"],1, function(x) sd(x) / mean(x) * 100)
cv_values_2 <- apply(2^corrected_phos[,Sample_ID$Clamp == "Pre"], 1, function(x) sd(x) / mean(x) * 100)
cv_values_3 <- apply(2^corrected_phos[,Sample_ID$Clamp == "Post"], 1, function(x) sd(x) / mean(x) * 100)
cv_values_4 <- apply(2^DIA_targeted_data_intensity, 1, function(x) sd(x) / mean(x) * 100)

# Create a data frame with protein names and CV values
df1 <- data.frame(Protein = rownames(corrected_prot), CV = cv_values_1)
df2 <- data.frame(Protein = rownames(corrected_phos), CV = cv_values_2)
df3 <- data.frame(Protein = rownames(corrected_phos), CV = cv_values_3)
df4 <- data.frame(Protein = rownames(DIA_targeted_data_intensity), CV = cv_values_4)

# Initialize the plot
df1 <- df1[order(df1$CV), ]
df2 <- df2[order(df2$CV), ]
df3 <- df3[order(df3$CV), ]
df4 <- df4[order(df4$CV), ]

df1$distri <- factor("Proteome")
df2$distri <- factor("Phosphoproteome_Basal")
df3$distri <- factor("Phosphoproteome_Insulin")
df4$distri <- factor("Phosphoproteome_technical_CV")
df <- rbind(df1,df2,df3,df4)
df$Protein <- reorder(df$Protein, df$CV)

# Plot the distributions
ggplot() +
  geom_density(data = subset(df, distri == "Proteome"), aes(x = CV, fill = "Proteome"), alpha = 0.2) +
  geom_density(data = subset(df, distri == "Phosphoproteome_Basal"), aes(x = CV, fill = "Baseline Phosphoproteome"), alpha = 0.1) +
  geom_density(data = subset(df, distri == "Phosphoproteome_Insulin"), aes(x = CV, fill = "Insulin Phosphoproteome"), alpha = 0.1) +
  geom_density(data = subset(df, distri == "Phosphoproteome_technical_CV"), aes(x = CV, fill = "Phospho technical CV"), alpha = 0.1) + 
  labs(x = "Inter-subject variation (%)", y = "Density", title = "") +
  theme_minimal() +
  xlim(0, 100) + theme + 
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),  # Move legend inside plot (adjust as needed)
        legend.background = element_rect(fill = "white", color = "black"),  # Add background to make it more readable
        legend.box.background = element_rect(color = "black")) +  # Optional: add a border around the legend box
  scale_fill_manual(values = c("Proteome" = "darkblue", "Baseline Phosphoproteome" = "darkred", "Insulin Phosphoproteome" = "darkred", "Phospho technical CV" = "grey")) +
  scale_y_continuous(breaks = seq(0, 0.06, 0.01))

#### Figure 1J & Figure 1 S2G - Total weigthed contribution GAM ####
# Load PCA results and clinical data for discovery and validation cohorts
pca_res_prot <- readRDS("pca_res_prot.rds")
pca_res_phos <- readRDS("pca_res_phos.rds")
pca_res_val_prot <- readRDS("pca_res_prot_validation.rds")
pca_res_val_phos <- readRDS("pca_res_phos_validation.rds")

# Calculate explained variance
explained_variance_discovery_prot <- pca_res_prot$sdev^2 / sum(pca_res_prot$sdev^2)
explained_variance_discovery_phos <- pca_res_phos$sdev^2 / sum(pca_res_phos$sdev^2)
explained_variance_validation_prot <- pca_res_val_prot$sdev^2 / sum(pca_res_val_prot$sdev^2)
explained_variance_validation_phos <- pca_res_val_phos$sdev^2 / sum(pca_res_val_phos$sdev^2)

# Extract PC scores
pc_scores_prot_new <- as.data.frame(pca_res_prot$x)
pc_scores_phos_new <- as.data.frame(pca_res_phos$x)
pc_scores_prot <- as.data.frame(pca_res_val_prot$x)
pc_scores_phos <- as.data.frame(pca_res_val_phos$x)

# Load and clean clinical data
Proteome_PC_clinical_discovery <- duplicated_df_prot
Phospho_PC_clinical_discovery <- duplicated_Our_clin_data
Proteome_PC_clinical_validation <- read.delim("Figures/Figure_1/Proteome_PC_clinical.txt", dec=",")
Phospho_PC_clinical_validation <- read.delim("Figures/Figure_1/Phospho_PC_clinical.txt", dec=",")

# Standardize clinical parameter names
rename_clinical_columns <- function(df) {
  df %>%
    rename_with(~ gsub("M\\.\\.µmol\\.\\.kg\\.min\\.\\.", "M.value", .)) %>%
    rename_with(~ gsub("FP\\.Glucose\\.\\.mmol\\.L\\.", "Fast.gluco", .)) %>%
    rename_with(~ gsub("HbA1c\\.\\.mmol\\.mol\\.", "HbA1c", .)) %>%
    rename_with(~ gsub("FS\\.Insulin\\.\\.mIE\\.L\\.", "Fast.insul_converted", .)) %>%
    rename_with(~ gsub("HOMA1\\.IR|HOMA\\.IR|HOMA_IR", "HOMA1.IR", .))
}



Proteome_PC_clinical_discovery <- rename_clinical_columns(Proteome_PC_clinical_discovery)
Phospho_PC_clinical_discovery <- rename_clinical_columns(Phospho_PC_clinical_discovery)
Proteome_PC_clinical_validation <- rename_clinical_columns(Proteome_PC_clinical_validation)
Phospho_PC_clinical_validation <- rename_clinical_columns(Phospho_PC_clinical_validation)

# Function to run GAM and calculate weighted contributions
analyze_gam <- function(pc_scores, clinical_data, explained_variance, cohort_name, modality) {
  gam_results <- list()
  for (i in 1:ncol(pc_scores)) {
    gam_results[[paste0("PC", i)]] <- gam(
      pc_scores[, i] ~ s(M.value) + s(Fast.gluco) + 
        s(HbA1c) + s(Fast.insul_converted) + s(HOMA1.IR),
      data = clinical_data,
      method = "REML"
    )
  }
  
  predictor_contributions <- data.frame()
  for (i in seq_along(gam_results)) {
    gam_summary <- summary(gam_results[[i]])
    contributions <- data.frame(
      Predictor = rownames(gam_summary$s.table),
      PC = paste0("PC", i),
      F_stat = gam_summary$s.table[, "F"]
    )
    predictor_contributions <- rbind(predictor_contributions, contributions)
  }
  
  predictor_contributions <- predictor_contributions %>%
    mutate(
      ExplainedVariance = explained_variance[as.numeric(sub("PC", "", PC))],
      WeightedContribution = F_stat * ExplainedVariance
    )
  
  total_contributions <- predictor_contributions %>%
    group_by(Predictor) %>%
    summarise(TotalWeightedContribution = sum(WeightedContribution, na.rm = TRUE)) %>%
    arrange(desc(TotalWeightedContribution))
  
  list(
    total_contributions = total_contributions,
    cohort_name = cohort_name,
    modality = modality
  )
}

# Run GAM analysis for discovery and validation cohorts
cohort1_prot_results <- analyze_gam(pc_scores_prot_new, Proteome_PC_clinical_discovery, explained_variance_discovery_prot, "Cohort 1", "Proteomics")
cohort1_phos_results <- analyze_gam(pc_scores_phos_new, Phospho_PC_clinical_discovery, explained_variance_discovery_phos, "Cohort 1", "Phosphoproteomics")
cohort2_prot_results <- analyze_gam(pc_scores_prot, Proteome_PC_clinical_validation, explained_variance_validation_prot, "Cohort 2", "Proteomics")
cohort2_phos_results <- analyze_gam(pc_scores_phos, Phospho_PC_clinical_validation, explained_variance_validation_phos, "Cohort 2", "Phosphoproteomics")

# Optional: Plot total contributions (one example)
plot_total_contributions <- function(total_contributions, cohort_name, modality) {
  ggplot(total_contributions, aes(x = reorder(Predictor, -TotalWeightedContribution), y = TotalWeightedContribution)) +
    geom_bar(stat = "identity", fill = "lightgrey", alpha = 0.5) +
    labs(
      x = "",
      y = "Total Weighted Contribution",
      title = paste(cohort_name, modality)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      aspect.ratio = 1.5 / 1
    )
}

# Plot for Discovery Cohort – Proteomics
plot_total_contributions(cohort1_prot_results$total_contributions, "Cohort 1", "Proteomics")

# Plot for Discovery Cohort – Phosphoproteomics
plot_total_contributions(cohort1_phos_results$total_contributions, "Cohort 1", "Phosphoproteomics")

# Plot for Validation Cohort – Proteomics
plot_total_contributions(cohort2_prot_results$total_contributions, "Cohort 2", "Proteomics")

# Plot for Validation Cohort – Phosphoproteomics
plot_total_contributions(cohort2_phos_results$total_contributions, "Cohort 2", "Phosphoproteomics")


#### Figure 1K - RandomForest classification ####

corrected_prot <- readRDS(file="corrected_prot.rds")

#proteome
Proteome_full <- as.data.frame(t(corrected_prot))
colnames(Proteome_full) <- rownames(corrected_prot)
rownames(Proteome_full) <- colnames(corrected_prot)
Proteome_full$disease <- duplicated_df_prot$Group
Proteome_full$GIR <- duplicated_df_prot$M..µmol..kg.min..
Proteome_full$IS <- "low"
Proteome_full[which(Proteome_full$GIR > median(Proteome_full$GIR, na.rm = T)),]$IS <- "high"
Proteome_full$IS <- as.factor(Proteome_full$IS)
Proteome_full$disease <- as.factor(Proteome_full$disease)
Proteome_full <- Proteome_full[!is.na(Proteome_full$GIR),]
colnames(Proteome_full) <- make.names(colnames(Proteome_full), unique = TRUE)

# Phospho
Phospho_full <- as.data.frame(t(corrected_phos[,Sample_ID_prot$New_SA_ID]))
colnames(Phospho_full) <- rownames(corrected_phos[,Sample_ID_prot$New_SA_ID])
rownames(Phospho_full) <- colnames(corrected_phos[,Sample_ID_prot$New_SA_ID])
Phospho_full$disease <- duplicated_df_prot$Group
Phospho_full$GIR <-duplicated_df_prot$M..µmol..kg.min..
Phospho_full$IS <- "low"
Phospho_full[which(Phospho_full$GIR > median(Phospho_full$GIR, na.rm = T)),]$IS <- "high"
Phospho_full$IS <- as.factor(Phospho_full$IS)
Phospho_full$disease <- as.factor(Phospho_full$disease)
Phospho_full <- Phospho_full[!is.na(Phospho_full$GIR),]
colnames(Phospho_full) <- make.names(colnames(Phospho_full), unique = TRUE)

train_rf_multiple <- function(data, outcome_col, n_repeats = 25, n_points = 100) {
  specificity_grid <- seq(0, 1, length.out = n_points)
  roc_list <- list()
  auc_list <- numeric(n_repeats)
  
  for (i in 1:n_repeats) {
    set.seed(i)
    train_index <- createDataPartition(data[[outcome_col]], p = 0.7, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    rf_model <- randomForest(as.formula(paste(outcome_col, "~ .")), 
                             data = train_data, 
                             ntree = 100)
    
    rf_probs <- randomForest:::predict.randomForest(rf_model, test_data, type = "prob")[, 2]
    test_labels <- factor(test_data[[outcome_col]], levels = levels(data[[outcome_col]]))
    
    roc_curve <- roc(test_labels, rf_probs, levels = levels(test_labels))
    auc_list[i] <- auc(roc_curve)
    
    interp_sens <- approx(1 - roc_curve$specificities, roc_curve$sensitivities, 
                          xout = specificity_grid, rule = 2)$y
    
    roc_list[[i]] <- data.frame(
      specificity = specificity_grid,
      sensitivity = interp_sens,
      run = i,
      outcome = outcome_col
    )
  }
  
  return(list(
    roc_data = bind_rows(roc_list),
    auc_values = auc_list
  ))
}

rf_result_IS <- train_rf_multiple(Proteome_full, "IS")
rf_result_disease <- train_rf_multiple(Proteome_full, "disease")
rf_result_IS_phospho <- train_rf_multiple(Phospho_full, "IS")
rf_result_disease_phospho <- train_rf_multiple(Phospho_full, "disease")

roc_data_IS <- rf_result_IS$roc_data
roc_data_disease <- rf_result_disease$roc_data
roc_data_IS_phospho <- rf_result_IS_phospho$roc_data
roc_data_disease_phospho <- rf_result_disease_phospho$roc_data

# Merge both datasets
roc_data <- bind_rows(roc_data_IS, roc_data_disease)

# Compute mean ROC curve with confidence intervals
roc_summary <- roc_data %>%
  group_by(specificity, outcome) %>%
  summarise(
    mean_sensitivity = mean(sensitivity),
    lower_ci = quantile(sensitivity, 0.025),  # 2.5% CI
    upper_ci = quantile(sensitivity, 0.975)   # 97.5% CI
  )

# Define colors
outcome_colors <- c("IS" = "#5a30cf", "disease" = "#bd342f")

ROC_proteome <- ggplot() +
  # Add confidence interval for IS
  geom_ribbon(data = roc_summary %>% filter(outcome == "IS"), 
              aes(x = specificity, ymin = lower_ci, ymax = upper_ci, fill = outcome), 
              alpha = 0.1) +
  # Add confidence interval for disease
  geom_ribbon(data = roc_summary %>% filter(outcome == "disease"), 
              aes(x = specificity, ymin = lower_ci, ymax = upper_ci, fill = outcome), 
              alpha = 0.1) +
  
  # Add mean ROC curve for IS
  geom_line(data = roc_summary %>% filter(outcome == "IS"), 
            aes(x = specificity, y = mean_sensitivity, color = outcome), 
            size = 1) +
  # Add mean ROC curve for disease
  geom_line(data = roc_summary %>% filter(outcome == "disease"), 
            aes(x = specificity, y = mean_sensitivity, color = outcome), 
            size = 1) +
  
  # Add diagonal reference line (random classifier)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.7) +
  
  # Labels and theme
  scale_color_manual(values = outcome_colors) +
  scale_fill_manual(values = outcome_colors) +
  labs(
    title = paste("Proteome"),
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Outcome",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio=1.5/1,
    legend.position = c(0.75, 0.25),  # Moves legend inside plot (adjust as needed)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5) 
  )


roc_data_phospho <- bind_rows(roc_data_IS_phospho, roc_data_disease_phospho)

# Compute mean ROC curve with confidence intervals
roc_summary_phospho <- roc_data_phospho %>%
  group_by(specificity, outcome) %>%
  summarise(
    mean_sensitivity = mean(sensitivity),
    lower_ci = quantile(sensitivity, 0.025),  # 2.5% CI
    upper_ci = quantile(sensitivity, 0.975)   # 97.5% CI
  )

# Define colors
outcome_colors <- c("IS" = "#5a30cf", "disease" = "#bd342f")

ROC_phospho <- ggplot() +
  # Add confidence interval for IS
  geom_ribbon(data = roc_summary_phospho %>% filter(outcome == "IS"), 
              aes(x = specificity, ymin = lower_ci, ymax = upper_ci, fill = outcome), 
              alpha = 0.1) +
  # Add confidence interval for disease
  geom_ribbon(data = roc_summary_phospho %>% filter(outcome == "disease"), 
              aes(x = specificity, ymin = lower_ci, ymax = upper_ci, fill = outcome), 
              alpha = 0.1) +
  
  # Add mean ROC curve for IS
  geom_line(data = roc_summary_phospho %>% filter(outcome == "IS"), 
            aes(x = specificity, y = mean_sensitivity, color = outcome), 
            size = 1) +
  # Add mean ROC curve for disease
  geom_line(data = roc_summary_phospho %>% filter(outcome == "disease"), 
            aes(x = specificity, y = mean_sensitivity, color = outcome), 
            size = 1) +
  
  # Add diagonal reference line (random classifier)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.7) +
  
  # Labels and theme
  scale_color_manual(values = outcome_colors) +
  scale_fill_manual(values = outcome_colors) +
  labs(
    title = paste("Phospho"),
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Outcome",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio=1.5/1,
    legend.position = c(0.75, 0.25),  # Moves legend inside plot (adjust as needed)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5) 
  )



ROC_proteome + ROC_phospho


get_mean_auc <- function(data, outcome_col, n_repeats = 25) {
  auc_values <- numeric(n_repeats)
  
  for (i in 1:n_repeats) {
    set.seed(i)
    train_index <- createDataPartition(data[[outcome_col]], p = 0.7, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    
    rf_model <- randomForest(as.formula(paste(outcome_col, "~ .")),
                             data = train_data, 
                             ntree = 100)
    
    rf_probs <- randomForest:::predict.randomForest(rf_model, test_data, type = "prob")[, 2]
    test_labels <- factor(test_data[[outcome_col]], levels = levels(data[[outcome_col]]))
    
    roc_curve <- roc(test_labels, rf_probs, levels = levels(test_labels))
    auc_values[i] <- auc(roc_curve)
  }
  
  mean(auc_values)
}

mean_auc_IS_proteome <- get_mean_auc(Proteome_full, "IS")
mean_auc_disease_proteome <- get_mean_auc(Proteome_full, "disease")
mean_auc_IS_phospho <- get_mean_auc(Phospho_full, "IS")
mean_auc_disease_phospho <- get_mean_auc(Phospho_full, "disease")

cat("Mean AUCs:\n")
cat("Proteome – IS:", round(mean_auc_IS_proteome, 3), "\n")
cat("Proteome – disease:", round(mean_auc_disease_proteome, 3), "\n")
cat("Phospho – IS:", round(mean_auc_IS_phospho, 3), "\n")
cat("Phospho – disease:", round(mean_auc_disease_phospho, 3), "\n")

#### Differentially abundant phosphosites - LIMMA ####
phospho_study2 <- log2(phospho_study[,1:154])
new_colnames <- word(colnames(phospho_study2),10, sep="_")
new_colnames <- as.integer(sub("S","",new_colnames))
colnames(phospho_study2) <- new_colnames
colnames(phospho_study2) <- as.integer(colnames(phospho_study2))
new_study <- phospho_study2[,str_sort(colnames(phospho_study2), numeric = TRUE)]
rownames(new_study) <- phospho_study$PTM_collapse_key

# Filtering for 25% quantified sites in the whole dataset
Phospho_muscle2 <- selectOverallPercent(new_study, 0.25, n=1)

# Removing batch effect
Phospho_muscle2 <- removeBatchEffect(Phospho_muscle2, day_batch_effect)
# Median normalising data
Phospho_muscle2 <- medianScaling(Phospho_muscle2)

Clamp <- factor(Sample_ID$Clamp, levels=c("Pre","Post"))
Disease <- factor(Sample_ID$disease, levels=c("NGT","T2D"))
Gender <- factor(Sample_ID$Gender, levels=c("Male","Female"))
design <- model.matrix(~ Disease * Clamp + Gender)

colnames(design) <- c("Intercept","Disease","Clamp","Gender","Disease_Clamp")
corfit <- duplicateCorrelation(Phospho_muscle2, design, block=Sample_ID$Subject_ID)
corfit$consensus

fit <- lmFit(Phospho_muscle2,design,block=Sample_ID$Subject_ID,correlation=corfit$consensus)
fit2 <- eBayes(fit)

Disease_main_effect_phos <- topTable(fit2, 
                                     coef=2, 
                                     number = Inf, sort.by = "none")
Disease_main_effect_phos$ID <- rownames(Disease_main_effect_phos)

Insulin_main_effect_phos <- topTable(fit2, 
                                     coef=3, 
                                     number = Inf, sort.by = "none")

Insulin_main_effect_phos$ID <- rownames(Insulin_main_effect_phos)

Gender_main_effect_phos <- topTable(fit2, 
                                    coef=4, 
                                    number = Inf, sort.by = "none")

Gender_main_effect_phos$ID <- rownames(Gender_main_effect_phos)

Interaction_main_effect_phos <- topTable(fit2, 
                                         coef=5, 
                                         number = Inf, sort.by = "none")

Interaction_main_effect_phos$ID <- rownames(Interaction_main_effect_phos)

g_sig <- Gender_main_effect_phos[Gender_main_effect_phos$adj.P.Val <=0.05,]
d_sig <- Disease_main_effect_phos[Disease_main_effect_phos$adj.P.Val <=0.05,]
i_sig <- Insulin_main_effect_phos[Insulin_main_effect_phos$adj.P.Val <=0.05,]

#### Targeted phospho processing ####
Targeted_phospho_plate_1 <- read.delim("R1405_Sky_Export_Plate1_2.csv",sep=",")
Targeted_phospho_plate_1 <- Targeted_phospho_plate_1[-grep("R1382_PRM_",Targeted_phospho_plate_1$Replicate),] # removing 3 technical standard spike in samples
Targeted_phospho_plate_1[Targeted_phospho_plate_1 == "#N/A"] <- NA # replacing #N/A with NA
Targeted_plate_1_and_2_Sample_ID <- read.delim("Targeted_plate_1_and_2_sample_ID.txt") # Loading in sample ID

#manual curration of sites#
phospho_peptides <- c("TBC1D4_S588_M1","TBC1D4_T642_M1","RAF1_S296_M1","TSC2_T1462_M1","HSPB1_S15_M1","HSPB1_S82_M1",
                      "MAPK3_T202_M1","MAPK3_Y204_M1","IRS1_S307_M1","MAPKAPK2_T334_M1","GSK3B_S9_M1","GRB10_S476_M1",
                      "AKT1S1_S183_M1","AKT1S1_T246_M1") 


Targeted_matrix <- matrix(nrow = 150, ncol = 15)
Targeted_matrix[,1] <- Targeted_plate_1_and_2_Sample_ID$ID
unique_peptides <- unique(Targeted_phospho_plate_1$Peptide.Modified.Sequence)
Targeted_matrix <- matrix(nrow = 150, ncol = length(unique_peptides) + 1)
Targeted_matrix[,1] <- Targeted_plate_1_and_2_Sample_ID$ID

# Loop through each unique peptide
for(i in 1:length(unique_peptides)){
  subset_data <- Targeted_phospho_plate_1[Targeted_phospho_plate_1$Peptide.Modified.Sequence == unique_peptides[i],]
  
  if(all(is.na(subset_data$light.Total.Area.Fragment)) | all(is.na(subset_data$heavy.Total.Area.Fragment))){
    # Calculate the ratio using MS1 areas if Fragment areas are not valid
    Targeted_matrix[,i+1] <- as.numeric(subset_data$light.Total.Area.MS1) /
      as.numeric(subset_data$heavy.Total.Area.MS1)
  } else {
    # Calculate the ratio using Fragment areas
    Targeted_matrix[,i+1] <- as.numeric(subset_data$light.Total.Area.Fragment) /
      as.numeric(subset_data$heavy.Total.Area.Fragment)
  }
}

# Convert the matrix to a data frame if needed
Targeted_matrix_df <- as.data.frame(Targeted_matrix)

# Assign column names, assuming you have the unique peptide sequences as columns
colnames(Targeted_matrix_df) <- c("Sample_ID", phospho_peptides)

Targeted_matrix_df[10:ncol(Targeted_matrix_df)] <- t(removeBatchEffect(t(Targeted_matrix_df[,10:ncol(Targeted_matrix_df)]), batch=c(rep(1,79),rep(2,71))))
Targeted_matrix_df <- cbind(Sample_ID[Targeted_matrix_df$Sample_ID,],duplicated_Our_clin_data[Targeted_matrix_df$Sample_ID,]$M..µmol..kg.min..,Targeted_matrix_df)
colnames(Targeted_matrix_df)[8] <- "Mvalue"
phospho_peptides_untargeted <- phospho_peptides[-7]

Untargeted_matrix <- matrix(nrow = 150, ncol = 13)
for(i in 1:length(phospho_peptides_untargeted)){
  Untargeted_matrix[,i] <- Phospho_muscle2[rownames(Phospho_muscle2) == phospho_peptides_untargeted[i],Targeted_matrix_df$Sample_ID]
}
colnames(Untargeted_matrix) <- paste0(phospho_peptides_untargeted,"_","DIA")
Targeted_and_untargeted_phos <- cbind(Targeted_matrix_df,Untargeted_matrix)

Plot_data <- Targeted_and_untargeted_phos
Plot_data <- Plot_data[!is.na(Plot_data$Mvalue),]
Plot_data$IS <- "> 40"
Plot_data[Plot_data$Mvalue > 30 & Plot_data$Mvalue < 40,]$IS <- "< 30-40"
Plot_data[Plot_data$Mvalue > 20 & Plot_data$Mvalue < 30,]$IS <- "< 20-30"
Plot_data[Plot_data$Mvalue > 10 & Plot_data$Mvalue < 20,]$IS <- "< 10-20"
Plot_data[Plot_data$Mvalue < 10,]$IS <- "< 10"

Plot_data$Clamp <- factor(Plot_data$Clamp, levels=c("Pre","Post"))
Plot_data$IS <- factor(Plot_data$IS, levels=c("< 10","< 10-20","< 20-30","< 30-40","> 40"))

#### Figure 1G & 1 S2B - Targeted ####
plot1 <- ggplot(Plot_data, aes(x=HSPB1_S15_M1_DIA,y=HSPB1_S15_M1)) + 
  geom_point(shape=21,size=4,alpha=0.25, color="black",fill="darkgrey") + theme_bw() +
  geom_smooth(method = "lm", se = TRUE, alpha=0.2, color="#4d4d4d") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="HSPB1 S15 abundance human skeletal muscle") + 
  ylab("Targeted PRM quantification (L/H AUC ratio)") + xlab("Untargeted DIA quantification (log2)") +
  annotate(geom="text",x=13,y=10, label="tau = 0.42", color = "black",size=5) + 
  annotate(geom = "text",x=13,y=9.2,label="P = 1.35e-14",size=5)


plot2 <- ggplot(Plot_data, aes(x=AKT1S1_S183_M1_DIA,y=AKT1S1_S183_M1)) + 
  geom_point(shape=21,size=4,alpha=0.25, color="black",fill="darkgrey") + theme_bw() +
  geom_smooth(method = "lm", se = TRUE, alpha=0.2, color="#4d4d4d") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="AKT1S1 S183 abundance human skeletal muscle") + 
  ylab("Targeted PRM quantification (L/H AUC ratio)") + xlab("Untargeted DIA quantification (log2)") +
  annotate(geom="text",x=11.5,y=0.6, label="tau = 0.63", color = "black",size=5) + 
  annotate(geom = "text",x=11.5,y=0.55,label="P < 2.2e-16",size=5)


grid.arrange(plot1, plot2, ncol = 2, nrow = 1)

#### Figure 3H - Targeted phospho validation ####
ggplot(Plot_data[Plot_data$Clamp == "Pre",], aes(x=Mvalue,y=HSPB1_S82_M1)) + 
  geom_point(shape=21,size=4,alpha=0.25, color="black",fill="darkgrey") + theme_bw() +
  geom_smooth(method = "lm", se = TRUE, alpha=0.2, color="#4d4d4d") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="HSPB1 S82 abundance human skeletal muscle") + 
  ylab("Targeted PRM quantification (L/H AUC ratio)") + xlab("M-value (umol/kg/min)") +
  annotate(geom="text",x=60,y=45, label="tau = -0.3", color = "black",size=5) + 
  annotate(geom = "text",x=60,y=40,label="P = 1.4e-04",size=5)

cor.test(Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Pre",]$HSPB1_S82_M1,Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Pre",]$Mvalue, method="kendall")

#### Figure 3 S1A - T2D vs NGT phosphosites ####
d_sig_ordered <- d_sig[order(d_sig$logFC),]
d_sig_ordered$rank <- 1:43
d_sig_ordered$ID <- sub("_"," ",word(rownames(d_sig),1,2,sep="_"))

ggplot(data=d_sig_ordered,aes(y=logFC, x=rank, color=adj.P.Val, label=ID)) + geom_point(size=5) +
  theme_bw() + theme(aspect.ratio=5/3) + coord_flip() + geom_text_repel() +
  scale_colour_gradient(low="#5fa35f",high="grey") +
  xlab("")

#### Figure 3 S1B - T2D vs NGT kinase enrichment ####
KinLib_Disease <- read.delim("Disease_enrichment-analysis-result-table.txt",
                             sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)


converted_gene_names <- read.delim("Converted_gene_names.txt",sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
colnames(converted_gene_names)[1] <- "kinase"

muscle_match_list <- inner_join(converted_gene_names,KinLib_Disease, by = "kinase")

protein_lib <- read.delim("Expression_atlas.tsv",
                          sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)

Found_muscle_match <- muscle_match_list[which(muscle_match_list$Approved.symbol %in% protein_lib$Gene),]
Found_muscle_match <- Found_muscle_match[,c(1,3,4,7,14,16,20,22)]
Found_muscle_match$UP_pAdj <- p.adjust(Found_muscle_match$upregulated_p_value, method="BH")
Found_muscle_match$DOWN_pAdj <- p.adjust(Found_muscle_match$downregulated_p_value, method="BH")

NGT_enrich_kinases <- Found_muscle_match[Found_muscle_match$DOWN_pAdj <= 0.05,c(1:4,7,8,10)]
NGT_enrich_kinases$Disease <- "NGT"
T2D_enrich_kinases <- Found_muscle_match[Found_muscle_match$UP_pAdj <= 0.05,c(1:6,9)]
T2D_enrich_kinases$Disease <- "T2D"
colnames(NGT_enrich_kinases)[5:7] <- c("Enrichment","p.val","adj.p.val")
colnames(T2D_enrich_kinases)[5:7] <- c("Enrichment","p.val","adj.p.val")

Disease_kinase_table <- rbind(T2D_enrich_kinases,NGT_enrich_kinases)
Disease_kinase_table[Disease_kinase_table$Disease == "NGT",]$Enrichment <- 
  -Disease_kinase_table[Disease_kinase_table$Disease == "NGT",]$Enrichment

ggplot(Disease_kinase_table, aes(x=rank, y=-log10(p.val), label=Approved.symbol, color=Disease)) + geom_point(aes(size=abs(Enrichment)), alpha=0.4,stroke=NA) +
  geom_text_repel(force_pull   = 0, nudge_y      = 0.05,
                  direction    = "y",
                  angle        = 90,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1,
                  max.overlaps = 12,
                  size=2) +
  xlab("") + ylab("-log10 p.value") +
  scale_color_manual(values=c("darkgreen","darkblue")) + theme(aspect.ratio=1/2,
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank())

#### Figure 4A - Volcano plot insulin regulation ####
Insulin_main_effect_phos$diffexpressed <- "NO"
Insulin_main_effect_phos$Short_ID <- sub("_"," ",word(Insulin_main_effect_phos$ID,1,2,sep="_"))
Insulin_main_effect_phos$diffexpressed[Insulin_main_effect_phos$logFC > 0 & Insulin_main_effect_phos$adj.P.Val < 0.05] <- "UP"
Insulin_main_effect_phos$diffexpressed[Insulin_main_effect_phos$logFC < 0 & Insulin_main_effect_phos$adj.P.Val < 0.05] <- "DOWN"
Insulin_main_effect_phos$delabel <- NA
Insulin_main_effect_phos$delabel[Insulin_main_effect_phos$diffexpressed != "NO"] <- Insulin_main_effect_phos$Short_ID[Insulin_main_effect_phos$diffexpressed != "NO"]
Insulin_main_effect_phos$point_size <- 1
Insulin_main_effect_phos[which(Insulin_main_effect_phos$diffexpressed != "NO"),]$point_size <- 2
plot_Insulin_main_effect_phos <- ggplot(data=Insulin_main_effect_phos, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.5,aes(size = point_size), stroke=NA) + scale_color_manual(values=c("#8494fa", "grey", "darkblue")) + ggtitle("Insulin signaling") + theme_minimal() + geom_text_repel() + scale_size_identity()
plot_Insulin_main_effect_phos + theme(legend.position = "none", panel.border=element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(size=.5))
table(Insulin_main_effect_phos$logFC > 0 & Insulin_main_effect_phos$adj.P.Val < 0.05)["TRUE"]
table(Insulin_main_effect_phos$logFC < 0 & Insulin_main_effect_phos$adj.P.Val < 0.05)["TRUE"]



#### Figure 4B - Insulin kinase enrichment ####
### import generated Gender enrichment results file
KinLib_Insulin <- read.delim("INSULIN_Kinase_enrichment-analysis-result-table.txt",
                             sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)


converted_gene_names <- read.delim("Converted_gene_names.txt",sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
colnames(converted_gene_names)[1] <- "kinase"

muscle_match_list <- inner_join(converted_gene_names,KinLib_Insulin, by = "kinase")

protein_lib <- read.delim("Expression_atlas.tsv",
                          sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)

Found_muscle_match <- muscle_match_list[which(muscle_match_list$Approved.symbol %in% protein_lib$Gene),]

Found_muscle_match <- Found_muscle_match[,c(1,3,4,7,14,16,20,22)]
Found_muscle_match$UP_pAdj <- p.adjust(Found_muscle_match$upregulated_p_value, method="BH")
Found_muscle_match$DOWN_pAdj <- p.adjust(Found_muscle_match$downregulated_p_value, method="BH")

colnames(Found_muscle_match)[c(5:6,9)]<- c("Enrichment","p.val","adj.p.val")
Found_muscle_match$diff <- "NO"
Found_muscle_match[Found_muscle_match$adj.p.val <= 0.05,]$diff <- "YES" 

ggplot(Found_muscle_match, aes(x=Enrichment, y=-log10(p.val), label=Approved.symbol, color=diff)) + geom_point(size=3, alpha=0.4,stroke=NA) +
  geom_text_repel(size=4) +
  theme + xlab("") + ylab("-log10 p.value") +
  scale_color_manual(values=c("grey","darkblue")) + theme(aspect.ratio=1/1, legend.position="none")


#### Extracting individual sites ####
Akt_site <- Phospho_muscle2[rownames(Phospho_muscle2)=="HBS1L_S67_M1",]
plot_Akt_site <- cbind(Akt_site, Sample_ID)
colnames(plot_Akt_site)[1] <- "Intensity"
plot_Akt_site$paired_ID <- factor(str_sub(plot_Akt_site$Subject_ID, start= -3))
plot_Akt_site$Clamp <- factor(plot_Akt_site$Clamp, levels = c("Pre","Post"))

p <-ggplot(plot_Akt_site, aes(x = factor(Clamp,levels = c("Pre","Post")), y = Intensity)) + 
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
  ggtitle("HBS1L_S67_M1")
p

#### Figure 4 S1A - Targeted phospho validation ####
Box_TSC2_T1462 <- ggplot(Plot_data,aes(x=Clamp, y=TSC2_T1462_M1, color=Clamp, fill = Clamp)) + 
  geom_boxplot(alpha=0.1) + 
  geom_point(stroke=NA, size=3, alpha=0.4) + 
  geom_line(aes(group=Subject_ID), alpha=0.4) + 
  #facet_wrap(~ Gender+disease, switch = "x",nrow=1) + 
  theme_minimal() +
  ylab("Targeted PRM quantification (L/H AUC ratio)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.position = "none") + 
  ggtitle(label="TSC2 T1462") + scale_color_manual(values = c("grey","darkblue")) + scale_fill_manual(values = c("grey","darkblue"))

Box_AKT1S1_T246 <- ggplot(Plot_data,aes(x=Clamp, y=AKT1S1_T246_M1, color=Clamp, fill = Clamp)) + 
  geom_boxplot(alpha=0.1) + 
  geom_point(stroke=NA, size=3, alpha=0.4) + 
  geom_line(aes(group=Subject_ID), alpha=0.4) + 
  #facet_wrap(~ Gender+disease, switch = "x",nrow=1) + 
  theme_minimal() +
  ylab("Targeted PRM quantification (L/H AUC ratio)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.position = "none") + 
  ggtitle(label="AKT1S1 T246") + scale_color_manual(values = c("grey","darkblue")) + scale_fill_manual(values = c("grey","darkblue"))

Box_AKT1S1_S183 <- ggplot(Plot_data,aes(x=Clamp, y=AKT1S1_S183_M1, color=Clamp, fill = Clamp)) + 
  geom_boxplot(alpha=0.1) + 
  geom_point(stroke=NA, size=3, alpha=0.4) + 
  geom_line(aes(group=Subject_ID), alpha=0.4) + 
  #facet_wrap(~ Gender+disease, switch = "x",nrow=1) + 
  theme_minimal() +
  ylab("Targeted PRM quantification (L/H AUC ratio)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.position = "none") + 
  ggtitle(label="AKT1S1 S183") + scale_color_manual(values = c("grey","darkblue")) + scale_fill_manual(values = c("grey","darkblue"))

grid.arrange(Box_TSC2_T1462, Box_AKT1S1_T246,Box_AKT1S1_S183, ncol = 3, nrow = 1)

# Arrange data by Subject_ID and Clamp, ensuring correct pairing
paired_data <- Plot_data %>%
  filter(Clamp %in% c("Pre", "Post")) %>%
  arrange(Subject_ID, Clamp) %>%
  group_by(Subject_ID) %>%
  filter(n() == 2)  # Ensure only subjects with both Pre and Post are included

# Split into Pre and Post
pre_data <- paired_data %>% filter(Clamp == "Pre") %>% pull(TSC2_T1462_M1)
post_data <- paired_data %>% filter(Clamp == "Post") %>% pull(TSC2_T1462_M1)
pre_data <- paired_data %>% filter(Clamp == "Pre") %>% pull(AKT1S1_T246_M1)
post_data <- paired_data %>% filter(Clamp == "Post") %>% pull(AKT1S1_T246_M1)
pre_data <- paired_data %>% filter(Clamp == "Pre") %>% pull(AKT1S1_S183_M1)
post_data <- paired_data %>% filter(Clamp == "Post") %>% pull(AKT1S1_S183_M1)

# Perform the paired t-test
paired_test <- t.test(post_data,pre_data, paired = TRUE)
# Output the results
print(paired_test)

#### Figure 4 S1B - Insulin category enrichment ####
# Load and process regulatory site annotations
regulatory_sites <- read.delim("Regulatory_sites.txt")
regulatory_sites <- regulatory_sites %>%
  filter(grepl("-p", MOD_RSD)) %>%
  mutate(MOD_RSD = sub("-p", "", MOD_RSD),
         comb_ID = paste0(ACC_ID, "_", MOD_RSD))

# Prepare phosphosite data
tmp_file_ins <- inner_join(Insulin_main_effect_phos, Gender_main_effect_phos, by = "ID")
Insulin_main_effect_phos$Protein <- tmp_file_ins$protein
Insulin_main_effect_phos$comb_ID <- paste0(Insulin_main_effect_phos$Protein, "_", word(Insulin_main_effect_phos$ID, 2, sep = "_"))

# Functional subsets
activated_phos <- regulatory_sites %>% filter(grepl("activity, induced", ON_FUNCTION))
inhibited_phos <- regulatory_sites %>% filter(grepl("activity, inhibited", ON_FUNCTION))
mol_assoc <- regulatory_sites %>% filter(grepl("molecular association", ON_FUNCTION))
intra_loc <- regulatory_sites %>% filter(grepl("intracellular localization", ON_FUNCTION))
disease_sites <- read.delim("Disease-associated_sites.txt") %>%
  filter(grepl("-p", MOD_RSD)) %>%
  mutate(MOD_RSD = sub("-p", "", MOD_RSD),
         comb_ID = paste0(ACC_ID, "_", MOD_RSD))

# Annotate phosphosites
add_flag <- function(df, flag_vector, column_name) {
  df[[column_name]] <- ifelse(df$comb_ID %in% flag_vector$comb_ID, "YES", "NO")
  return(df)
}

Insulin_main_effect_phos <- Insulin_main_effect_phos %>%
  add_flag(regulatory_sites, "regulatory") %>%
  add_flag(activated_phos, "activated") %>%
  add_flag(inhibited_phos, "inhibited") %>%
  add_flag(mol_assoc, "mol_assoc") %>%
  add_flag(intra_loc, "intra_loc") %>%
  add_flag(disease_sites, "disease")

# Fisher test function
run_fisher <- function(flag_column) {
  tbl <- table(
    significance = ifelse(Insulin_main_effect_phos$adj.P.Val < 0.05, "significant", "non-significant"),
    category = Insulin_main_effect_phos[[flag_column]]
  )
  test <- fisher.test(tbl)
  return(c(p.value = test$p.value, OR = unname(test$estimate)))
}

# Apply to all features
results <- data.frame(
  t(sapply(c("activated", "inhibited", "mol_assoc", "intra_loc", "disease"), run_fisher))
)
results$type <- c("Activation", "Inhibition", "Molecular Association", "Intracellular Localization", "Disease-associated")

# Reorder factor levels for desired top-to-bottom order in the plot
site_terms_enriched <- results %>%
  mutate(
    log10.pval = -log10(p.value),
    type = factor(type, levels = c(
      "Disease-associated",
      "Molecular Association",
      "Intracellular Localization",
      "Inhibition",
      "Activation"
    ))
  )


# Plot theme (same as before)
clean_theme <- theme(
  axis.title = element_text(size = 12),
  axis.line = element_line(colour = "black"),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12)
)

# Final plot
ggplot(site_terms_enriched, aes(x = type, y = OR, fill = log10.pval)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_bw() +
  clean_theme +
  ylab("Odds Ratio") +
  xlab("") +
  scale_fill_gradient(low = "grey", high = "#4078bd")

#### Figure 4 S1G-H - targeted boxplot/correlation ####

TBC1D4_T642_plot_1 <- ggplot(Plot_data,aes(x=Clamp, y=TBC1D4_T642_M1, color=Clamp, fill=Clamp)) + 
  geom_boxplot(alpha=0.15) + 
  geom_point(stroke=NA, size=3, alpha=0.4) + 
  geom_line(aes(group=Subject_ID), color="black", alpha=0.4) + 
  facet_wrap(~ IS, switch = "x",nrow=1) + 
  theme_bw() +
  ylab("Targeted PRM quantification (L/H AUC ratio)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_blank(),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="TBC1D4 T642 insulin responsiveness") + 
  scale_color_manual(values=c("darkgrey","darkblue")) + scale_fill_manual(values=c("darkgrey","darkblue"))

AKT1S1_T246_plot_1 <- ggplot(Plot_data,aes(x=Clamp, y=AKT1S1_T246_M1, color=Clamp, fill=Clamp)) + 
  geom_boxplot(alpha=0.15) + 
  geom_point(stroke=NA, size=3, alpha=0.4) + 
  geom_line(aes(group=Subject_ID), color="black", alpha=0.4) + 
  facet_wrap(~ IS, switch = "x",nrow=1) + 
  theme_bw() +
  ylab("Targeted PRM quantification (L/H AUC ratio)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_blank(),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="AKT1S1 T246 insulin responsiveness") + 
  scale_color_manual(values=c("darkgrey","darkblue")) + scale_fill_manual(values=c("darkgrey","darkblue"))


grid.arrange(TBC1D4_T642_plot_1,AKT1S1_T246_plot_1, ncol=1,nrow=2)



TBC1D4_T642_Mval <- ggplot(Plot_data[Plot_data$Clamp == "Post",], aes(x=Mvalue,y=TBC1D4_T642_M1)) + 
  geom_point(shape=21,size=4,alpha=0.4, color="black",fill="darkblue") + theme_bw() +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1, color="darkblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="TBC1D4 T642 & insulin sensitivity") + 
  ylab("Targeted PRM quantification (L/H AUC ratio)") + xlab("M-value (umol/kg/min)") +
  annotate(geom="text",x=60,y=0.4, label="tau = 0.07", color = "black",size=5) + 
  annotate(geom = "text",x=60,y=0.35,label="P = 0.40",size=5)

AKT1S1_T246_Mval <- ggplot(Plot_data[Plot_data$Clamp == "Post",], aes(x=Mvalue,y=AKT1S1_T246_M1)) + 
  geom_point(shape=21,size=4,alpha=0.4, color="black",fill="darkblue") + theme_bw() +
  geom_smooth(method = "lm", se = TRUE, alpha=0.1, color="darkblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),
        axis.title.x = element_text(size = 18),  # Increase x-axis label size
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 16),  # Reduce legend title size
        legend.text = element_text(size = 14)) + 
  ggtitle(label="AKT1S1 T246 & insulin sensitivity") + 
  ylab("Targeted PRM quantification (L/H AUC ratio)") + xlab("M-value (umol/kg/min)") +
  annotate(geom="text",x=60,y=0.3, label="tau = 0.09", color = "black",size=5) + 
  annotate(geom = "text",x=60,y=0.25,label="P = 0.28",size=5)

grid.arrange(TBC1D4_T642_Mval,AKT1S1_T246_Mval,ncol=1,nrow=2)


cor.test(Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Post",]$TBC1D4_T642_M1,Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Post",]$Mvalue, method="kendall")
cor.test(Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Post",]$AKT1S1_T246_M1,Targeted_and_untargeted_phos[Targeted_and_untargeted_phos$Clamp == "Post",]$Mvalue, method="kendall")

#### Figure 4D - stratified insulin heatmap ####
ins_heatmap_data <- cbind(target_df,t(Phospho_muscle2[rownames(Phospho_muscle2) %in% i_sig$ID,Sample_ID$Clamp == "Post"]))
colnames(ins_heatmap_data)[22] <- "GIR"
ins_heatmap_data <- ins_heatmap_data %>% filter(!is.na(GIR))
ins_heatmap_data$GIR_bins <- cut(ins_heatmap_data$GIR, breaks = c(0, 10, 20, 30, 40, Inf), right = FALSE, include.lowest = TRUE, labels = c("0-10", "10-20", "20-30", "30-40", ">40"))

# Calculate averages for each bin
average_intensities <- ins_heatmap_data %>% 
  group_by(GIR_bins) %>%
  summarise(across(contains("_"), mean, na.rm = TRUE))

# Z-score the calculated averages
z_scores <- as.data.frame(scale(average_intensities[,-1])) # Exclude the first column (GIR_bins)
# Convert to matrix, if necessary
z_scores_matrix <- as.matrix(z_scores)
rownames(z_scores_matrix) <- z_scores$GIR_bins
z_scores_matrix <- z_scores_matrix[,-1] # Remove the GIR_bins column

# Transpose the matrix if your rows are intervals and columns are phosphosites
transposed_matrix <- t(z_scores)

# Perform hierarchical clustering on the transposed matrix
phospho_dist <- dist(transposed_matrix)  # Calculate the distance matrix
phospho_clust <- hclust(phospho_dist, method = "complete")  # Perform clustering

# Silhouette analysis for different numbers of clusters
sil_widths <- sapply(2:10, function(k) {
  clusters <- cutree(phospho_clust, k)
  sil_obj <- silhouette(clusters, dist(transposed_matrix))
  mean(sil_obj[, "sil_width"])
})

# Plotting the average silhouette width for each number of clusters
plot(2:10, sil_widths, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters", ylab = "Average silhouette width")

color_palette <- brewer.pal(n = 9, name = "Blues") 
z_scores_matrix_transposed <- t(z_scores_matrix)

# Convert the matrix to a Heatmap object
ht <- Heatmap(z_scores_matrix_transposed, 
              name = "z-score", 
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = TRUE,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100))

# Draw the heatmap with increased spacing between clusters
# Define the gap sizes as numeric values
gap_sizes <- unit(c(3, 3, 3, 3), "mm")  # Adjust the number and size of gaps as needed

# Draw the heatmap with the specified gaps
draw(ht, row_split = 5, row_gap = gap_sizes)

#### Figure 4G - Insulin sensitivity & Insulin LogFC ####
LOGFC_INS_COR_plot <- data.frame(Post_complete_cor_frame[which(Post_complete_cor_frame$ID %in% Insulin_main_effect_phos_sign$ID),2:4], Insulin_main_effect_phos_sign)

Known_substrates <- Kinase_Substrate_Dataset[Kinase_Substrate_Dataset$UI %in% word(LOGFC_INS_COR_plot$ID,1,2,sep="_"),]
MTOR_sites <- Known_substrates[Known_substrates$GENE == "MTOR",]$UI
AKT_sites <- Known_substrates[Known_substrates$GENE == "AKT1"|Known_substrates$GENE == "AKT2",]$UI

LOGFC_INS_COR_plot$kinase <- ""
LOGFC_INS_COR_plot[LOGFC_INS_COR_plot$ID %in% c("TFEB_S122_M3","RPTOR_S863_M1","AKT1S1_S183_M1","GRB10_S476_M1","EIF4EBP1_T70_M1","PATL1_S179_M2","PATL1_S184_M2","EIF4EBP1_T37_M2","MAF1_S75_M1","PRKAA2_S377_M1"),]$kinase <- "mTOR"
LOGFC_INS_COR_plot[LOGFC_INS_COR_plot$ID %in% c("FLNC_S2233_M1","PANK4_T406_M1","PANK2_S189_M1","NIBAN1_S602_M2","GSK3A_S21_M2","TSC2_T1462_M1","NDRG2_T348_M2","GSK3B_S9_M1","RANBP3_S126_M1","AKT1S1_T246_M1","TBC1D4_S341_M1","TBC1D4_T642_M1","TBC1D4_S588_M1"),]$kinase <- "AKT"

ggplot(LOGFC_INS_COR_plot, aes(x=M..µmol..kg.min.._estimate, y= logFC, col =kinase, label=kinase)) + geom_point(aes(size=kinase, shape=kinase),alpha=0.5, stroke=0.1) + geom_vline(xintercept = -0.1, linetype = "dashed") + geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() + theme(panel.grid = element_blank()) + xlab("Insulin sensitivity association (Tau)") + ylab("Insulin response (LogFC)") + scale_color_manual(values=c("lightgrey","darkred","darkblue")) +
  scale_size_manual(values = c("mTOR" = 3, "AKT"=3, 3)) + scale_shape_manual(values=c(19,19,19)) + scale_alpha_manual(values=c(0.1,0.85,0.85))


# Wilcoxon geneset test
set.seed(123)
correlation <- LOGFC_INS_COR_plot$M..µmol..kg.min.._estimate

LOGFC_INS_COR_plot$MTOR <- 1
LOGFC_INS_COR_plot[LOGFC_INS_COR_plot$kinase == "mTOR",]$MTOR <- 2
LOGFC_INS_COR_plot$AKT <- 1
LOGFC_INS_COR_plot[LOGFC_INS_COR_plot$kinase == "AKT",]$AKT <- 2

mtorSubstrateIdx <- which(LOGFC_INS_COR_plot$MTOR == 2)
AktSubstrateIdx <- which(LOGFC_INS_COR_plot$AKT == 2)

# Perform the Wilcoxon Gene Set Test
result_mTOR <- wilcoxGST(mtorSubstrateIdx, correlation)
result_Akt <- wilcoxGST(AktSubstrateIdx, correlation)

print(result_mTOR)
print(result_Akt)


#### Figure 4H - Preserved Insulin signaling ####
Akt_site <- Phospho_muscle2[rownames(Phospho_muscle2)=="TBC1D4_T642_M1",]
plot_Akt_site <- cbind(Akt_site, Sample_ID)
colnames(plot_Akt_site)[1] <- "Intensity"
plot_Akt_site$paired_ID <- factor(str_sub(plot_Akt_site$Subject_ID, start= -3))
plot_Akt_site$Clamp <- factor(plot_Akt_site$Clamp, levels = c("Pre","Post"))

plot_Akt_site$Substrate <- "TBC1D4 T642"

Akt_site2 <- Phospho_muscle2[rownames(Phospho_muscle2)=="GSK3B_S9_M1",]
plot_Akt_site2 <- cbind(Akt_site2, Sample_ID)
colnames(plot_Akt_site2)[1] <- "Intensity"
plot_Akt_site2$paired_ID <- factor(str_sub(plot_Akt_site2$Subject_ID, start= -3))
plot_Akt_site2$Clamp <- factor(plot_Akt_site2$Clamp, levels = c("Pre","Post"))

plot_Akt_site2$Substrate <- "GSK3B S9"


plot_Akt_site_comb <- rbind(plot_Akt_site,plot_Akt_site2)


site <- cbind(plot_Akt_site_comb,duplicated_Our_clin_data $M..µmol..kg.min..)
site <- site %>% filter(!is.na(`duplicated_Our_clin_data$M..µmol..kg.min..`))

site$GIR_bins <- cut(site$`duplicated_Our_clin_data$M..µmol..kg.min..`, breaks = c(0, 10, 20, 30, 40, Inf), right = FALSE, include.lowest = TRUE, labels = c("0-10", "10-20", "20-30", "30-40", ">40"))


site_summary <- site %>%
  group_by(GIR_bins,Clamp,Substrate) %>%
  summarise(
    Mean = mean(Intensity,na.rm = TRUE),
    SD = sd(Intensity,na.rm = TRUE),
    N = n(),
    SEM = SD / sqrt(N),
    ymin = Mean - SEM,
    ymax = Mean + SEM
  )

my_labeller <- function(labels) {
  paste0(" ", labels, " ")  # Adding spaces to create a box-like appearance
}

p <- ggplot(site, aes(x = factor(Clamp, levels = c("Pre", "Post")), y = Intensity)) + 
  geom_boxplot(aes(fill = Clamp), alpha = .2,size = 0.1) +
  geom_line(aes(group = paired_ID),size=0.1) + 
  geom_point(size = 1.5, alpha = 0.5, aes(col = Clamp)) + 
  facet_grid(rows = vars(Substrate), cols = vars(GIR_bins), scales = "free_y", switch = "x", labeller = labeller(Substrate = my_labeller)) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 18),
    panel.grid = element_blank(),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    strip.text.x = element_text(size = 18, angle = 0, hjust = 0.5),  # Centered GIR_bins labels
    strip.text.y = element_text(size = 18, angle = 270, hjust = 0.5),  # Centered Substrate labels
    strip.placement = "outside",
    axis.title = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.line.x = element_line(size=.5),
    axis.text.y = element_text(size = 14),
    axis.line.y = element_line(size = .5),
    strip.background = element_rect(fill = "grey90", colour = "black")
  ) + 
  ylab("Log2 Intensity") +
  ggtitle("Preserved Insulin signaling") + 
  scale_fill_manual(values = c("grey", "darkblue")) + 
  scale_color_manual(values = c("grey", "darkblue"))

p

cor.test(Phospho_muscle2[rownames(Phospho_muscle2) == "TBC1D4_T642_M1",which(Sample_ID$Clamp == "Post")],target_df$M..µmol..kg.min.., method="kendall")
cor.test(Phospho_muscle2[rownames(Phospho_muscle2) == "GSK3B_S9_M1",which(Sample_ID$Clamp == "Post")],target_df$M..µmol..kg.min.., method="kendall")

#### Figure 4I - IS associated signaling ####
Akt_site <- Phospho_muscle2[rownames(Phospho_muscle2)=="RPTOR_S863_M1",]
plot_Akt_site <- cbind(Akt_site, Sample_ID)
colnames(plot_Akt_site)[1] <- "Intensity"
plot_Akt_site$paired_ID <- factor(str_sub(plot_Akt_site$Subject_ID, start= -3))
plot_Akt_site$Clamp <- factor(plot_Akt_site$Clamp, levels = c("Pre","Post"))

plot_Akt_site$Substrate <- "RPTOR S863"

Akt_site2 <- Phospho_muscle2[rownames(Phospho_muscle2)=="EIF4EBP1_T70_M1",]
plot_Akt_site2 <- cbind(Akt_site2, Sample_ID)
colnames(plot_Akt_site2)[1] <- "Intensity"
plot_Akt_site2$paired_ID <- factor(str_sub(plot_Akt_site2$Subject_ID, start= -3))
plot_Akt_site2$Clamp <- factor(plot_Akt_site2$Clamp, levels = c("Pre","Post"))

plot_Akt_site2$Substrate <- "EIF4EBP1 T70"

plot_Akt_site_comb <- rbind(plot_Akt_site,plot_Akt_site2)


site <- cbind(plot_Akt_site_comb,duplicated_Our_clin_data $M..µmol..kg.min..)
site <- site %>% filter(!is.na(`duplicated_Our_clin_data$M..µmol..kg.min..`))

site$GIR_bins <- cut(site$`duplicated_Our_clin_data$M..µmol..kg.min..`, breaks = c(0, 10, 20, 30, 40, Inf), right = FALSE, include.lowest = TRUE, labels = c("0-10", "10-20", "20-30", "30-40", ">40"))


site_summary <- site %>%
  group_by(GIR_bins,Clamp,Substrate) %>%
  summarise(
    Mean = mean(Intensity,na.rm = TRUE),
    SD = sd(Intensity,na.rm = TRUE),
    N = n(),
    SEM = SD / sqrt(N),
    ymin = Mean - SEM,
    ymax = Mean + SEM
  )

my_labeller <- function(labels) {
  paste0(" ", labels, " ")  # Adding spaces to create a box-like appearance
}

p <- ggplot(site, aes(x = factor(Clamp, levels = c("Pre", "Post")), y = Intensity)) + 
  geom_boxplot(aes(fill = Clamp), alpha = .2,size = 0.1) +
  geom_line(aes(group = paired_ID),size=0.1) + 
  geom_point(size = 1.5, alpha = 0.5, aes(col = Clamp)) + 
  facet_grid(rows = vars(Substrate), cols = vars(GIR_bins), scales = "free_y", switch = "x", labeller = labeller(Substrate = my_labeller)) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 18),
    panel.grid = element_blank(),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    strip.text.x = element_text(size = 18, angle = 0, hjust = 0.5),  # Centered GIR_bins labels
    strip.text.y = element_text(size = 18, angle = 270, hjust = 0.5),  # Centered Substrate labels
    strip.placement = "outside",
    axis.title = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.line.x = element_line(size=.5),
    axis.text.y = element_text(size = 14),
    axis.line.y = element_line(size = .5),
    strip.background = element_rect(fill = "grey90", colour = "black")
  ) + 
  ylab("Log2 Intensity") +
  ggtitle("Preserved Insulin signaling") + 
  scale_fill_manual(values = c("grey", "darkblue")) + 
  scale_color_manual(values = c("grey", "darkblue"))

p

cor.test(Phospho_muscle2[rownames(Phospho_muscle2) == "RPTOR_S863_M1",which(Sample_ID$Clamp == "Post")],target_df$M..µmol..kg.min.., method="kendall")
cor.test(Phospho_muscle2[rownames(Phospho_muscle2) == "EIF4EBP1_T70_M1",which(Sample_ID$Clamp == "Post")],target_df$M..µmol..kg.min.., method="kendall")



#### Baseline Phosphosite X insulin sensitivity associations ####
pre_phospho_intensities <- Phospho_muscle2[,which(Sample_ID$Clamp == "Pre")]
pre_Sample_ID <- Sample_ID[which(Sample_ID$Clamp == "Pre"),]
Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% pre_Sample_ID$Subject_ID,]
target_df <- Our_clin_data[match(unique(pre_Sample_ID$Subject_ID),Our_clin_data$Subject.ID),]

clin_var <- c("M..µmol..kg.min..","HbA1c..mmol.mol.","BMI..kg.m2.","Age",
              "FS.Insulin..mIE.L.","HOMA1.IR","Creatinine..µmol.L.","P.TG..mmol.L.",
              "P.Chol..mmol.L.","FFA..mmol.L.","FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.")

baseline_empty_matrix <- matrix(nrow=nrow(pre_phospho_intensities),ncol=4)
baseline_p.val_matrix <- matrix(nrow=nrow(pre_phospho_intensities),ncol=length(clin_var))
baseline_estimate_matrix <- matrix(nrow=nrow(pre_phospho_intensities),ncol=length(clin_var))

empty_list <- list()
for(k in 1:length(clin_var)){
  try({
    for(i in 1:nrow(pre_phospho_intensities)){
      empty_list <- cor.test(scale(pre_phospho_intensities[i,]), scale(log2(target_df[,clin_var[k]])), method="kendall", na.action = "na.exclude")
      baseline_p.val_matrix[i,k] <- empty_list$p.value
      baseline_estimate_matrix[i,k] <- empty_list$estimate
    }
  }, silent = TRUE)
}

colnames(baseline_p.val_matrix) <- paste0(clin_var,"_p.value")
colnames(baseline_estimate_matrix) <- paste0(clin_var,"_estimate")

alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

new_matrix_adjust <- matrix(nrow=nrow(baseline_p.val_matrix),ncol=ncol(baseline_p.val_matrix))
for(i in 1:12){
  new_matrix_adjust[,i] <- p.adjust(baseline_p.val_matrix[,i],method="BH")
}
colnames(new_matrix_adjust) <- paste0(clin_var,"_adj.p.value")
complete_p.val_matrix <- alternate.cols(baseline_p.val_matrix,new_matrix_adjust)

complete_cor_matrix <- cbind(complete_p.val_matrix,baseline_estimate_matrix)
complete_cor_frame <- data.frame(complete_cor_matrix)
complete_cor_frame <- cbind(rownames(pre_phospho_intensities), complete_cor_frame)

complete_cor_frame <- complete_cor_frame[,c(1:3,26,4:5,27,6:7,28,8:9,29,10:11,30,12:13,31,
                                            14:15,32,16:17,33,18:19,34,20:21,35,22:23,36,
                                            24:25,37)]
colnames(complete_cor_frame)[1] <- "ID"


#### Figure 3D - AMPKy3 S65 M-value correlation ####
test <- as.data.frame(cbind(pre_phospho_intensities[which(rownames(pre_phospho_intensities) == "PRKAG3_S65_M1"),],target_df$M..µmol..kg.min..,target_df$Group))
colnames(test) <- c("Intensity","GIR","Group")             

ggplot(test, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.4, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgreen", "black")) +
  ggtitle("PRKAG3 S65") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=7.3, label="tau = -0.42", color = "black",size=4) +
  annotate(geom = "text",x=60,y=7,label="p = 8.3e-08",size=4)

PRE_GIR_Cor <- complete_cor_frame[,1:4]
colnames(PRE_GIR_Cor) <- c("ID","p.val","adj.pval","estimate")
PRE_GIR_Cor$diff <- "NO"
PRE_GIR_Cor[PRE_GIR_Cor$adj.pval <= 0.05,]$diff <- "yes"

### make volcano correlation plot
ggplot(PRE_GIR_Cor, aes(x=estimate,y=-log10(p.val), color=diff, label=ID)) +
  geom_point(size=2) + xlim(-0.6,0.6) + scale_color_manual(values=c("grey","darkgreen")) +
  theme_bw() + geom_text_repel()

GIR_phos_complete_cor_frame <- complete_cor_frame[,1:4]
colnames(GIR_phos_complete_cor_frame) <- c("ID","p.val","adj.pval","estimate")
GIR_phos_complete_cor_frame$diff <- "NO"
GIR_phos_complete_cor_frame[GIR_phos_complete_cor_frame$adj.pval <= 0.05,]$diff <- "YES"

phos_sig_GIR <- GIR_phos_complete_cor_frame[GIR_phos_complete_cor_frame$diff == "YES",]

phos_sig_GIR <- phos_sig_GIR[order(phos_sig_GIR$estimate),]
phos_sig_GIR$rank <-1:118

ggplot(phos_sig_GIR, aes(x=estimate,y=-log10(p.val), label=ID)) +
  geom_point(aes(size=-log10(p.val))) + scale_color_manual(values=c("grey","darkgreen")) +
  theme_bw() + geom_text_repel() + theme(aspect.ratio=1/2)


#### Figure 3B - integrating proteome correlation with phospho correlation data ####
Prot_PRE_GIR_Cor <- readRDS(file="Prot_PRE_GIR_Cor.rds")
Sign_phos <- unique(word(phos_sig_GIR$ID,1,sep="_"))

PROT_phos_cor_overlap <- Prot_PRE_GIR_Cor[which(word(word(Prot_PRE_GIR_Cor$ID,1,sep="_"),1,sep=";") %in% Sign_phos),]
PROT_phos_cor_overlap$Gene <- word(PROT_phos_cor_overlap$ID,1,sep="_")
PHOS_prot_cor_overlap <- phos_sig_GIR[which(word(phos_sig_GIR$ID,1,sep="_") %in% word(PROT_phos_cor_overlap$ID,1,sep="_")),]
PHOS_prot_cor_overlap$Gene <- word(PHOS_prot_cor_overlap$ID,1,sep="_")

Comb_cor_overlap <- inner_join(PROT_phos_cor_overlap,PHOS_prot_cor_overlap,by="Gene",suffix=c("_Prot","_Phos"))

Comb_cor_overlap$alpha_pam <- ifelse(Comb_cor_overlap$diff_Prot == 'NO', 1, 0.2)

ggplot(Comb_cor_overlap,aes(x=estimate_Phos,y=-log10(p.val_Phos), color=diff_Prot, alpha=alpha_pam)) + geom_point(stroke=NA) +
  scale_x_break(c(-0.275, 0.275)) + xlim(-0.6,0.6) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values=c("darkgreen","black")) + theme_bw()



#### Figure 3I - Barplot Hskm cells phospho ####
Hskm_phospho_data <- read.delim("20240121_184924_DIAphospho_MK2i_quant_Report_collapsed.txt")

Hskm_conditions <- read.delim("Condition_names_Hskm.txt")
colnames(Hskm_phospho_data)[1:12] <- Hskm_conditions$Biological.Sample.ID 
log2_Hskm_phospho_data <- log2(Hskm_phospho_data[,1:12]) #log2 transform data
rownames(log2_Hskm_phospho_data) <- Hskm_phospho_data$PTM_collapse_key

Filtered_inhib_data <- medianScaling(log2_Hskm_phospho_data)
colnames(Filtered_inhib_data) <- word(colnames(Filtered_inhib_data),1,sep="_")

PRKAG3_data <- t(Filtered_inhib_data[rownames(Filtered_inhib_data) == "PRKAG3_S65_M1",])
ID_g3 <- rownames(PRKAG3_data)
PRKAG3_data <- as.data.frame(PRKAG3_data)
PRKAG3_data$condition <- factor(ID_g3, levels = c("DMSO","MK2i"))
colnames(PRKAG3_data)[1] <- "Intensity"

mean_df <- PRKAG3_data %>%
  group_by(condition) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE),
            SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

ggplot() + geom_bar(data = mean_df, aes(x = condition, y = MeanIntensity, fill=condition),
                    width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = condition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM),
                position = position_dodge(0.8), width = 0.25) + geom_jitter(data = PRKAG3_data, aes(x=condition,y=Intensity),size=2, width=0.05) +
  theme_minimal() + theme(aspect.ratio=1.5/1, axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black",size=14),
                          axis.title=element_text(size=20),plot.title = element_text(hjust = 0.5,size=18),legend.position="none", panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(6.25, 7.75)) + xlab("") + ylab("log2 intensity") + labs(title = "PRKAG3 S65") + scale_fill_manual(values=c("#ebebeb","#9c9c9c"))


t.test(Filtered_inhib_data[rownames(Filtered_inhib_data) == "PRKAG3_S65_M1",grepl("DMSO",colnames(Filtered_inhib_data))],
       Filtered_inhib_data[rownames(Filtered_inhib_data) == "PRKAG3_S65_M1",grepl("MK2i",colnames(Filtered_inhib_data))], var.equal = T)


HSPB1_data <- t(Filtered_inhib_data[rownames(Filtered_inhib_data) == "HSPB1_S15_M1",])
ID_g3 <- rownames(HSPB1_data)
HSPB1_data <- as.data.frame(HSPB1_data)
HSPB1_data$condition <- factor(ID_g3, levels = c("DMSO","MK2i"))
colnames(HSPB1_data)[1] <- "Intensity"

mean_df <- HSPB1_data %>%
  group_by(condition) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE),
            SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")


ggplot() + geom_bar(data = mean_df, aes(x = condition, y = MeanIntensity, fill=condition),
                    width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = condition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM),
                position = position_dodge(0.8), width = 0.25) + geom_jitter(data = HSPB1_data, aes(x=condition,y=Intensity),size=2, width=0.25) +
  theme_minimal() + theme(aspect.ratio=1.5/1, axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black",size=14),
                          axis.title=element_text(size=20),plot.title = element_text(hjust = 0.5,size=18),legend.position="none", panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(13, 15.5)) + xlab("") + ylab("log2 intensity") + labs(title = "HSPB1 S15") + scale_fill_manual(values=c("#ebebeb","#cccccc","#9c9c9c"))

t.test(Filtered_inhib_data[rownames(Filtered_inhib_data) == "HSPB1_S15_M1",grepl("DMSO",colnames(Filtered_inhib_data))],
       Filtered_inhib_data[rownames(Filtered_inhib_data) == "HSPB1_S15_M1",grepl("MK2i",colnames(Filtered_inhib_data))], var.equal = T)

#### Figure 1L - radarplot ####
Base_phos <- c(nrow(complete_cor_frame[which(complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),]),
nrow(complete_cor_frame[which(complete_cor_frame$HbA1c..mmol.mol._adj.p.value <= 0.05),]),
nrow(complete_cor_frame[which(complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]),
nrow(complete_cor_frame[which(complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),]),
nrow(complete_cor_frame[which(complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),]))

complete_cor_frame[which(complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]$ID %in% complete_cor_frame[which(complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),]$ID
complete_cor_frame[which(complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),]$ID %in% complete_cor_frame[which(complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]$ID
complete_cor_frame[which(complete_cor_frame$HbA1c..mmol.mol._adj.p.value <= 0.05),]$ID %in% complete_cor_frame[which(complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]$ID


ins_phos <- c(nrow(Post_complete_cor_frame[which(Post_complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),]),
nrow(Post_complete_cor_frame[which(Post_complete_cor_frame$HbA1c..mmol.mol._adj.p.value <= 0.05),]),
nrow(Post_complete_cor_frame[which(Post_complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]),
nrow(Post_complete_cor_frame[which(Post_complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),]),
nrow(Post_complete_cor_frame[which(Post_complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),]))

nrow(complete_cor_frame[which(complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),])
nrow(complete_cor_frame[which(complete_cor_frame$HbA1c..mmol.mol._adj.p.value <= 0.05),])
nrow(complete_cor_frame[which(complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),])
nrow(complete_cor_frame[which(complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),])
nrow(complete_cor_frame[which(complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),])

Base_prot <- c(nrow(Prot_complete_cor_frame[which(Prot_complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),]),
nrow(Prot_complete_cor_frame[which(Prot_complete_cor_frame$HbA1c..mmol.mol._adj.p.value <= 0.05),]),
nrow(Prot_complete_cor_frame[which(Prot_complete_cor_frame$FS.Insulin..mIE.L._adj.p.value <= 0.05),]),
nrow(Prot_complete_cor_frame[which(Prot_complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),]),
nrow(Prot_complete_cor_frame[which(Prot_complete_cor_frame$FP.Glucose..mmol.L._adj.p.value <= 0.05),]))


# Create a data frame for the plot
categories <- c("Category 1", "Category 2", "Category 3", "Category 4", "Category 5")

# Maximum and minimum values for each category
max_values <- c(160, 160, 160, 160, 160)
min_values <- c(0, 0, 0, 0, 0)

# Data for two lines (e.g., two different groups or observations)
group1 <- c(3, 4, 2, 5, 3)
group2 <- c(4, 3, 4, 3, 4)

# Combine data into a data frame
radar_data <- data.frame(rbind(max_values, min_values, Base_phos, Base_prot))
radar_data <- data.frame(rbind(max_values, min_values, Base_phos, ins_phos, Base_prot))

colnames(radar_data) <- c("M-value","HbA1c","Fasting Insulin","HOMA1-IR","Fasting Glucose")
ring_labels <- seq(0, max(max_values), by = 40)

# Create the radar chart
radarchart(radar_data, axistype = 1,
           # Customization options
           pcol = c("#c43333", "#4967c9"), plty = 1, plwd = 1,
           cglcol="lightgrey", cglty=1, axislabcol="black", caxislabels=ring_labels,
           cglwd=0.8)

radarchart(radar_data, axistype = 1,
           # Customization options
           pcol = c("#c43333", "purple","#4967c9"), plty = 1, plwd = 1,
           cglcol="lightgrey", cglty=1, axislabcol="black", caxislabels=ring_labels,
           cglwd=0.8)


# Add a legend
legend(x = 0.5, y = 1,  # Adjust the position as needed
       legend = c("Baseline Phosphoproteome", "Insulin Phosphoproteome","Proteome"), 
       col = c("#c43333", "purple","#4967c9"), 
       lty = 1, lwd = 1,
       bty = "n")  # 'bty = "n"' removes the box around the legend



#### Figure 1L - Venn diagram of overlapping phosphosites and proteins ####
data_comb <- list(
  HOMA1_IR = complete_cor_frame[which(complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),]$ID,
  M_value = complete_cor_frame[which(complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),]$ID)

data_comb_prot <- list(
  HOMA1_IR = Prot_complete_cor_frame[which(Prot_complete_cor_frame$HOMA1.IR_adj.p.value <= 0.05),]$ID,
  M_value = Prot_complete_cor_frame[which(Prot_complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05),]$ID)

# Create a Venn diagram
venn_data_IS <- list(
  HOMA1_IR = data_comb$HOMA1_IR,
  M_value = data_comb$M_value
)

venn.plot <- venn.diagram(
  x = venn_data_IS,
  filename = NULL,
  category.names = c("HOMA1_IR", "M_value"),
  output = TRUE,
  fill= c("#d9bdbd", "darkred")
)

# Display the plot
if (!is.null(venn.plot)) {
  grid.draw(venn.plot)
}

#prot
venn_data_IS_prot <- list(
  HOMA1_IR = data_comb_prot$HOMA1_IR,
  M_value = data_comb_prot$M_value
)

venn.plot <- venn.diagram(
  x = venn_data_IS_prot,
  filename = NULL,
  category.names = c("HOMA1_IR", "M_value"),
  output = TRUE,
  fill= c("#d4d1e0", "darkblue")
)
# Display the plot
if (!is.null(venn.plot)) {
  grid.draw(venn.plot)
}



#### 30min time point associations with insulin sensitivity ####
post_phospho_intensities <- Phospho_muscle2[,which(Sample_ID$Clamp == "Post")]
post_Sample_ID <- Sample_ID[which(Sample_ID$Clamp == "Post"),]
Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% post_Sample_ID$Subject_ID,]
target_df <- Our_clin_data[match(unique(post_Sample_ID$Subject_ID),Our_clin_data$Subject.ID),]

baseline_empty_matrix <- matrix(nrow=nrow(post_phospho_intensities),ncol=4)
baseline_p.val_matrix <- matrix(nrow=nrow(post_phospho_intensities),ncol=length(clin_var))
baseline_estimate_matrix <- matrix(nrow=nrow(post_phospho_intensities),ncol=length(clin_var))

clin_var <- c("M..µmol..kg.min..","HbA1c..mmol.mol.","BMI..kg.m2.","Age",
              "FS.Insulin..mIE.L.","HOMA1.IR","Creatinine..µmol.L.","P.TG..mmol.L.",
              "P.Chol..mmol.L.","FFA..mmol.L.","FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.")


empty_list <- list()
for(k in 1:length(clin_var)){
  try({
    for(i in 1:nrow(post_phospho_intensities)){
      empty_list <- cor.test(scale(post_phospho_intensities[i,]), scale(log2(target_df[,clin_var[k]])), method="kendall", na.action = "na.exclude")
      baseline_p.val_matrix[i,k] <- empty_list$p.value
      baseline_estimate_matrix[i,k] <- empty_list$estimate
    }
  }, silent = TRUE)
}


colnames(baseline_p.val_matrix) <- paste0(clin_var,"_p.value")
colnames(baseline_estimate_matrix) <- paste0(clin_var,"_estimate")



new_matrix_adjust <- matrix(nrow=nrow(baseline_p.val_matrix),ncol=ncol(baseline_p.val_matrix))
for(i in 1:12){
  new_matrix_adjust[,i] <- p.adjust(baseline_p.val_matrix[,i],method="BH")
}
colnames(new_matrix_adjust) <- paste0(clin_var,"_adj.p.value")
complete_p.val_matrix <- alternate.cols(baseline_p.val_matrix,new_matrix_adjust)

Post_complete_cor_matrix <- cbind(complete_p.val_matrix,baseline_estimate_matrix)
Post_complete_cor_frame <- data.frame(Post_complete_cor_matrix)
Post_complete_cor_frame <- cbind(rownames(post_phospho_intensities), Post_complete_cor_frame)

Post_complete_cor_frame <- Post_complete_cor_frame[,c(1:3,26,4:5,27,6:7,28,8:9,29,10:11,30,12:13,31,
                                                      14:15,32,16:17,33,18:19,34,20:21,35,22:23,36,
                                                      24:25,37)]
colnames(Post_complete_cor_frame)[1] <- "ID"


write.table(complete_cor_frame,"Updated_baseline_Cor_values.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")
write.table(Post_complete_cor_frame,"Updated_insulin_Cor_values.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")
write.table(Prot_complete_cor_frame,"Updated_Prot_Cor_values.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")


#### Figure 3 S1D & Figure 4 S1E - plotting individual correlations 30min time point ####
test <- as.data.frame(cbind(post_phospho_intensities[which(rownames(post_phospho_intensities) == "PRKAG3_S65_M1"),],target_df$M..µmol..kg.min..,target_df$Group))
colnames(test) <- c("Intensity","GIR","Group")             


ggplot(test, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.4, size=2) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkgreen", "black")) +
  ggtitle("PRKAG3 S65") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=7.3, label="tau = -0.37", color = "black",size=4) +
  annotate(geom = "text",x=60,y=7,label="p = 2.1e-06",size=4)

test <- as.data.frame(cbind(post_phospho_intensities[which(rownames(post_phospho_intensities) == "IRS1_Y612_M1"),],target_df$M..µmol..kg.min..,target_df$Group))
colnames(test) <- c("Intensity","GIR","Group")             

ggplot(test, aes(x=as.numeric(GIR), y=as.numeric(Intensity))) + geom_point(aes(color=Group,stroke=NA),alpha=0.4, size=3) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkblue", "black")) +
  ggtitle("IRS1 Y612") +
  xlab("Glucose infusion rate (umol/kg*min)") + ylab("Log2 intensity") +
  annotate(geom="text",x=60,y=7.3, label="tau = -0.37", color = "black",size=4) +
  annotate(geom = "text",x=60,y=7,label="p = 2.1e-06",size=4)

#### overlapping pre/post correlated phosphosites ####
sig_pre_cor <- complete_cor_frame[complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05,1:4]
sig_post_cor <- Post_complete_cor_frame[Post_complete_cor_frame$M..µmol..kg.min.._adj.p.value <= 0.05,1:4]

sig_post_cor <- cbind(sig_post_cor,complete_cor_frame[which(complete_cor_frame$ID %in% sig_post_cor$ID),2:4])
sig_post_cor <- sig_post_cor[,c(1,5:7,2:4)]
colnames(sig_post_cor) <- c("ID","pre.p.val","pre.adj.p.val","pre.est",
                            "post.p.val","post.adj.p.val","post.est")
sig_pre_cor <- cbind(sig_pre_cor,Post_complete_cor_frame[which(Post_complete_cor_frame$ID %in% sig_pre_cor$ID),2:4])
colnames(sig_pre_cor) <- c("ID","pre.p.val","pre.adj.p.val","pre.est",
                           "post.p.val","post.adj.p.val","post.est")

comb_pre_post_cor <- rbind(sig_pre_cor,sig_post_cor[which(!(sig_post_cor$ID %in% sig_pre_cor$ID)),])

comb_pre_post_cor$sign <- "Pre"
comb_pre_post_cor[comb_pre_post_cor$post.adj.p.val <= 0.05,]$sign <- "Post"
comb_pre_post_cor[comb_pre_post_cor$post.adj.p.val <= 0.05 & comb_pre_post_cor$pre.adj.p.val <= 0.05,]$sign <- "Both"


ggplot(comb_pre_post_cor,aes(x=pre.est,y=post.est, color=sign)) + geom_point(stroke=NA,alpha=0.5, size=3) +
  theme_bw() + theme(aspect.ratio=1/1) + scale_color_manual(values=c("grey","darkblue","darkgreen"))


Post_30min_correlation_only <- comb_pre_post_cor[comb_pre_post_cor$sign == "Post",]


nrow(comb_pre_post_cor[comb_pre_post_cor$sign == "Pre",]) #78
nrow(comb_pre_post_cor[comb_pre_post_cor$sign == "Post",]) # 66
nrow(comb_pre_post_cor[comb_pre_post_cor$sign == "Both",]) # 40

nrow(comb_pre_post_cor[comb_pre_post_cor$sign == "Post" & comb_pre_post_cor$post.est > 0,]) # 66
nrow(comb_pre_post_cor[comb_pre_post_cor$sign == "Post" & comb_pre_post_cor$post.est < 0,]) # 66


unique(comb_pre_post_cor$ID)


#### Figure 4E - U-shaped 30 minute time point ####
Insulin_dependent_matrix <- Post_complete_cor_frame[,1:4]
Insulin_dependent_matrix$dependent <- "NO" 
Insulin_dependent_matrix[Insulin_dependent_matrix$ID %in% Post_30min_correlation_only$ID,]$dependent <- "YES"
colnames(Insulin_dependent_matrix) <- c("ID","p.val","adj.p.val","estimate","dependent")
Insulin_dependent_matrix$ID <- sub("_"," ",word(Insulin_dependent_matrix$ID,1,2,sep="_"))

ggplot() + geom_point(data=Insulin_dependent_matrix[Insulin_dependent_matrix$adj.p.val < 0.05,],aes(x=estimate,y=-log10(p.val), color=dependent),size=2,alpha=0.3, stroke=NA) +
  geom_hline(yintercept = -log10(0.000339)) + scale_color_manual(values=c("grey","darkblue")) +
  geom_point(data=Insulin_dependent_matrix[Insulin_dependent_matrix$adj.p.val > 0.05,],aes(x=estimate,y=-log10(p.val), color=dependent),size=0.75, alpha=0.3, stroke=NA) +
  theme_bw() + theme(aspect.ratio=1/1, panel.border=element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(size=.5)) #+ geom_text_repel()


#### Figure 3G - Predicting upstream kinase of AMPKy3 S65 ####
Score_g3_site <- read.delim("Figures/Figure_3/score-site-result-table.tsv")

top_proteins <- Score_g3_site %>%
  arrange(desc(score_log2)) %>%
  head(8)

ggplot(Score_g3_site,aes(x=score_rank, y=score_log2, color=score_rank)) + geom_point(size=1,alpha=0.5) + geom_text_repel(data = top_proteins, aes(label = kinase), 
                                                                                                                          box.padding = 0.5,  # Further reduced box padding
                                                                                                                          point.padding = 0.5,
                                                                                                                          segment.color = 'grey50',
                                                                                                                          nudge_x = 1, hjust = 0,
                                                                                                                          size = 3)  +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") + scale_color_gradient(low= "darkblue", high = "darkgreen") + xlab("Kinase rank") + ylab("Score (log2)") +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) + geom_hline(yintercept = 0,linetype="dashed",size=0.25)

#### creating sequence window ####
Insulin_main_effect_phos$Protein <- phospho_study[which(Insulin_main_effect_phos$ID %in% phospho_study$PTM_collapse_key),]$PG.UniProtIds
Insulin_main_effect_phos$position <- phospho_study[which(Insulin_main_effect_phos$ID %in% phospho_study$PTM_collapse_key),]$PTM_collapse_key
Insulin_main_effect_phos$position <- word(Insulin_main_effect_phos$position,2,sep="_")
Insulin_main_effect_phos$position <- as.numeric(gsub("[STY]", "", Insulin_main_effect_phos$position))

#Downloading sequence from uniprot. Takes a while to download (2hours ish)
testir <- getUniProt(Insulin_main_effect_phos$Protein)
testir2 <- testir

#write.table(testir2,"sequence_window_nonprocessed.txt")
Insulin_main_effect_phos$new_uniprotID <- gsub("_.*","",testir2)
Insulin_main_effect_phos$down_site <- Insulin_main_effect_phos$position - 9
Insulin_main_effect_phos$up_site <- Insulin_main_effect_phos$position + 9

Insulin_main_effect_phos$dupl_down_site <- Insulin_main_effect_phos$position - 9
Insulin_main_effect_phos$dupl_down_site[Insulin_main_effect_phos$dupl_down_site < 0] <- 0
Insulin_main_effect_phos$window <- str_sub(testir, Insulin_main_effect_phos$dupl_down_site, Insulin_main_effect_phos$up_site)
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -7] <- paste0("________", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -7])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -6] <- paste0("_______", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -6])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -5] <- paste0("______", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -5])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -4] <- paste0("_____", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -4])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -3] <- paste0("____", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -3])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -2] <- paste0("___", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -2])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -1] <- paste0("__", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == -1])
Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == 0] <- paste0("_", Insulin_main_effect_phos$window[Insulin_main_effect_phos$down_site == 0])

Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 10)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 10)],"_________")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 11)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 11)],"________")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 12)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 12)],"_______")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 13)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 13)],"______")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 14)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 14)],"_____")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 15)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 15)],"____")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 16)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 16)],"___")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 17)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 17)],"__")
Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 18)] <- paste0(Insulin_main_effect_phos$window[which(str_length(Insulin_main_effect_phos$window) == 18)],"_")

#### KEGG pathway phosphoproteins ####
original_gene_list <- as.matrix(Insulin_main_effect_phos[,c(1)])
rownames(original_gene_list) <- word(Insulin_main_effect_phos[,7],1,sep="_")
rownames(original_gene_list) <- original_gene_list[,1]
original_gene_list[,2] <- as.numeric(original_gene_list[,2])
phosCollapse(original_gene_list, rownames(original_gene_list), stat, by='max')
abs(original_gene_list)

organism = "org.Hs.eg.db"

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = original_gene_list[names(original_gene_list) %in% dedup_ids$SYMBOL]
df2 <- as.data.frame(df2)
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$df2
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)

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


#### Figure 5 - Sex differences in the phosphoproteome ####
# making volcano plot #
Gender_main_effect_phos$diffexpressed <- "NO"
Gender_main_effect_phos$diffexpressed[Gender_main_effect_phos$logFC > 0 & Gender_main_effect_phos$adj.P.Val < 0.05] <- "UP"
Gender_main_effect_phos$diffexpressed[Gender_main_effect_phos$logFC < 0 & Gender_main_effect_phos$adj.P.Val < 0.05] <- "DOWN"
Gender_main_effect_phos$delabel <- NA
Gender_main_effect_phos$delabel[Gender_main_effect_phos$diffexpressed != "NO"] <- Gender_main_effect_phos$ID[Gender_main_effect_phos$diffexpressed != "NO"]
table(Gender_main_effect_phos$logFC > 0 & Gender_main_effect_phos$adj.P.Val < 0.05)["TRUE"]
table(Gender_main_effect_phos$logFC < 0 & Gender_main_effect_phos$adj.P.Val < 0.05)["TRUE"]

ggplot() + geom_point(data=Gender_main_effect_phos[which(Gender_main_effect_phos$logFC > 0 & Gender_main_effect_phos$adj.P.Val <= 0.05),], aes(x = logFC, y = -log10(P.Value), color = -log10(P.Value)),size=1.5) +
  scale_color_gradient(low = "#8ba0b5", high = "#3674b3") + 
  new_scale_color() + 
  geom_point(data = Gender_main_effect_phos[which(Gender_main_effect_phos$logFC < 0 & Gender_main_effect_phos$adj.P.Val <= 0.05),], aes(x = logFC, y = -log10(P.Value), color = -log10(P.Value)),size=1.5) +
  scale_color_gradient(low = "#bda6a8", high = "#bd626b") +
  geom_point(data = Gender_main_effect_phos[which(Gender_main_effect_phos$adj.P.Val >0.05),], aes(x = logFC, y = -log10(P.Value)),color="grey", size=1, alpha=0.5, stroke=NA) +
  theme_bw() + theme(aspect.ratio=4/3) + theme(legend.position = "none", panel.grid = element_blank(), panel.border = element_blank()) + geom_vline(xintercept = 0) 

#individual boxplots (export in 4x4)
Akt_site <- Phospho_muscle2[rownames(Phospho_muscle2)=="UBB_S57_M1",]
#Akt_site <- t(Akt_site[,1:154])
#colnames(Akt_site)
#Sample_ID2 <- t(Sample_ID)
plot_Akt_site <- cbind(Akt_site, Sample_ID)
colnames(plot_Akt_site)[1] <- "Intensity"
plot_Akt_site$paired_ID <- factor(str_sub(plot_Akt_site$Subject_ID, start= -3))
plot_Akt_site$Clamp <- factor(plot_Akt_site$Clamp, levels = c("Pre","Post"))
plot_Akt_site$Gender <- factor(word(plot_Akt_site$Group,1,sep="-"),levels=c("Male","Female"))


plot_Akt_site$Subject_ID <- as.factor(plot_Akt_site$Subject_ID)
plot_Akt_site$disease <- as.factor(plot_Akt_site$disease)
plot_Akt_site$Clamp <- factor(plot_Akt_site$Clamp, levels = c("Pre", "Post"))
colnames(plot_Akt_site)[1]<-"Intensity"
mean_df <- plot_Akt_site %>%
  group_by(disease, Clamp) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Create the plot
ggplot(plot_Akt_site, aes(x = factor(Gender,levels=c("Male","Female")), y = Intensity)) +
  # Add points
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(color = Gender), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, alpha=0.4, stroke=NA) +
  labs(x="",y = "Log2 Intensity") +
  scale_fill_manual(values = c("Male" = "#bda6a8", "Female" = "#8ba0b5")) +
  scale_color_manual(values = c("Male" = "#bd626b", "Female" = "#3674b3")) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),) +
  coord_cartesian(ylim = c(7, 9)) + theme(aspect.ratio=1.5/1)



#### Figure 5D - Sex enriched kinases ####
KinLib_Gender <- read.delim("GENDER_FDR_enrichment-analysis-result-table.txt",
                            sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)


converted_gene_names <- read.delim("Converted_gene_names.txt",sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
colnames(converted_gene_names)[1] <- "kinase"

muscle_match_list <- inner_join(converted_gene_names,KinLib_Gender, by = "kinase")

protein_lib <- read.delim("Expression_atlas.tsv",
                          sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)

Found_muscle_match <- muscle_match_list[which(muscle_match_list$Approved.symbol %in% protein_lib$Gene),]

Found_muscle_match <- Found_muscle_match[,c(1,3,4,7,14,16,20,22)]
Found_muscle_match$UP_pAdj <- p.adjust(Found_muscle_match$upregulated_p_value, method="BH")
Found_muscle_match$DOWN_pAdj <- p.adjust(Found_muscle_match$downregulated_p_value, method="BH")

Male_enrich_kinases <- Found_muscle_match[Found_muscle_match$DOWN_pAdj <= 0.05,c(1:4,7,8,10)]
Male_enrich_kinases$gender <- "Male"
Female_enrich_kinases <- Found_muscle_match[Found_muscle_match$UP_pAdj <= 0.05,c(1:6,9)]
Female_enrich_kinases$gender <- "Female"

colnames(Male_enrich_kinases)[5:7] <- c("Enrichment","p.val","adj.p.val")
colnames(Female_enrich_kinases)[5:7] <- c("Enrichment","p.val","adj.p.val")

Gender_kinase_table <- rbind(Female_enrich_kinases,Male_enrich_kinases)


Gender_kinase_table[Gender_kinase_table$gender == "Male",]$Enrichment <- 
  -Gender_kinase_table[Gender_kinase_table$gender == "Male",]$Enrichment


Gender_kinase_table <- Gender_kinase_table[order(Gender_kinase_table$Enrichment),]
Gender_kinase_table$rank <- as.numeric(1:nrow(Gender_kinase_table))

ggplot(Gender_kinase_table, aes(x=rank, y=-log10(p.val), label=Approved.symbol, color=gender)) + geom_point(aes(size=abs(Enrichment)), alpha=0.4, stroke=NA) +
  geom_text_repel() +
  theme + xlab("") + ylab("-log10 p.value") +
  scale_color_manual(values=c("darkblue","darkred")) + theme(aspect.ratio=1/2) +
  scale_y_continuous(limits=c(3,9),breaks=seq(2, 8, 2)) + xlim(0,13) + scale_size_continuous(range = c(8, 20))


#### Figure 5E x-chromosome regulated genes ####
# Create a list of your differentially expressed genes
diff_expressed_genes <- word(word(Gender_main_effect[which(Gender_main_effect$adj.P.Val <= 0.05),]$ID,1,sep = "_"),1,sep=";")
diff_expressed_genes <- Gender_main_effect_phos[which(Gender_main_effect_phos$adj.P.Val <= 0.05),]$ID
diff_expressed_genes <- unique(word(diff_expressed_genes,1,sep="_"))

# Use biomaRt to get the location data for your genes
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), 
               filters = 'hgnc_symbol', 
               values = diff_expressed_genes, 
               mart = mart)

# Remove genes not located in standard chromosomes
genes <- genes[genes$chromosome_name %in% c("X"), ]
genes_gr <- makeGRangesFromDataFrame(genes,
                                     seqnames.field = "chromosome_name",
                                     start.field = "start_position",
                                     end.field = "end_position",
                                     keep.extra.columns = TRUE)

# Change the seqnames to include 'chr' prefix
seqlevels(genes_gr) <- paste0("chr", seqlevels(genes_gr))

# Plot again
kp <- plotKaryotype(genome = "hg19", chromosomes = "chrX")
kpPoints(kp, data = genes_gr, y=0, col = "black", cex = 1.5)
kpText(kp, data = genes_gr, y = 0.1, labels = genes_gr$hgnc_symbol, cex = 1)
#### Figure 5J Sex-spexific insulin sensitivity correlation ####
pre_phospho_intensities_Male <- Phospho_muscle2[,which(Sample_ID$Clamp == "Pre" & Sample_ID$Gender == "Male")]
pre_Sample_ID_male <- Sample_ID[which(Sample_ID$Clamp == "Pre" & Sample_ID$Gender == "Male"),]
Our_clin_data_male_phospho <- Clinical_data[Clinical_data$Subject.ID %in% pre_Sample_ID$Subject_ID,]
target_df_male_phospho <- Our_clin_data_male_phospho[match(unique(pre_Sample_ID_male$Subject_ID),Our_clin_data_male_phospho$Subject.ID),]

clin_var <- c("M..µmol..kg.min..","HbA1c..mmol.mol.","BMI..kg.m2.","Age",
              "FS.Insulin..mIE.L.","HOMA1.IR","Creatinine..µmol.L.","P.TG..mmol.L.",
              "P.Chol..mmol.L.","FFA..mmol.L.","FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.")

baseline_empty_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Male),ncol=4)
baseline_p.val_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Male),ncol=length(clin_var))
baseline_estimate_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Male),ncol=length(clin_var))

empty_list <- list()
for(k in 1:length(clin_var)){
  try({
    if(shapiro.test(log2(target_df[,clin_var[k]]))[[2]] <0.05){
      for(i in 1:nrow(pre_phospho_intensities)){
        empty_list <- cor.test(scale(as.numeric(pre_phospho_intensities_Male[i,])),scale(log2(target_df_male_phospho[,clin_var[k]])),method="kendall",na.action = "na.exclude")
        baseline_p.val_matrix[i,k] <- empty_list$p.value
        baseline_estimate_matrix[i,k] <- empty_list$estimate
      }
    } else {
      for(i in 1:nrow(pre_phospho_intensities)){
        empty_list <- cor.test(scale(as.numeric(pre_phospho_intensities_Male[i,])),scale(log2(target_df_male_phospho[,clin_var[k]])),method="pearson",na.action = "na.exclude")
        baseline_p.val_matrix[i,k] <- empty_list$p.value
        baseline_estimate_matrix[i,k] <- empty_list$estimate
      }
    }
  }, silent = TRUE)
}

colnames(baseline_p.val_matrix) <- paste0(clin_var,"_p.value")
colnames(baseline_estimate_matrix) <- paste0(clin_var,"_estimate")

alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

new_matrix_adjust <- matrix(nrow=nrow(baseline_p.val_matrix),ncol=ncol(baseline_p.val_matrix))
for(i in 1:12){
  new_matrix_adjust[,i] <- p.adjust(baseline_p.val_matrix[,i],method="BH")
}
colnames(new_matrix_adjust) <- paste0(clin_var,"_adj.p.value")
complete_p.val_matrix <- alternate.cols(baseline_p.val_matrix,new_matrix_adjust)

complete_cor_matrix <- cbind(complete_p.val_matrix,baseline_estimate_matrix)
complete_cor_frame_male_phospho <- data.frame(complete_cor_matrix)
complete_cor_frame_male_phospho <- cbind(rownames(pre_phospho_intensities_Male), complete_cor_frame_male_phospho)

complete_cor_frame_male_phospho <- complete_cor_frame_male_phospho[,c(1:3,26,4:5,27,6:7,28,8:9,29,10:11,30,12:13,31,
                                                                      14:15,32,16:17,33,18:19,34,20:21,35,22:23,36,
                                                                      24:25,37)]
colnames(complete_cor_frame_male_phospho)[1] <- "ID"


complete_cor_frame_male_phospho <- complete_cor_frame_male_phospho[,c(1,2:4)]
colnames(complete_cor_frame_male_phospho) <- c("ID","p.val","adj.pval","estimate")


#female phospho
pre_phospho_intensities_Female <- Phospho_muscle2[,which(Sample_ID$Clamp == "Pre" & Sample_ID$Gender == "Female")]
pre_Sample_ID_Female <- Sample_ID[which(Sample_ID$Clamp == "Pre" & Sample_ID$Gender == "Female"),]
Our_clin_data_Female_phospho <- Clinical_data[Clinical_data$Subject.ID %in% pre_Sample_ID$Subject_ID,]
target_df_Female_phospho <- Our_clin_data_Female_phospho[match(unique(pre_Sample_ID_Female$Subject_ID),Our_clin_data_Female_phospho$Subject.ID),]

clin_var <- c("M..µmol..kg.min..","HbA1c..mmol.mol.","BMI..kg.m2.","Age",
              "FS.Insulin..mIE.L.","HOMA1.IR","Creatinine..µmol.L.","P.TG..mmol.L.",
              "P.Chol..mmol.L.","FFA..mmol.L.","FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.")

baseline_empty_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Female),ncol=4)
baseline_p.val_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Female),ncol=length(clin_var))
baseline_estimate_matrix <- matrix(nrow=nrow(pre_phospho_intensities_Female),ncol=length(clin_var))

empty_list <- list()
for(k in 1:length(clin_var)){
  try({
    if(shapiro.test(log2(target_df[,clin_var[k]]))[[2]] <0.05){
      for(i in 1:nrow(pre_phospho_intensities)){
        empty_list <- cor.test(scale(as.numeric(pre_phospho_intensities_Female[i,])),scale(log2(target_df_Female_phospho[,clin_var[k]])),method="kendall",na.action = "na.exclude")
        baseline_p.val_matrix[i,k] <- empty_list$p.value
        baseline_estimate_matrix[i,k] <- empty_list$estimate
      }
    } else {
      for(i in 1:nrow(pre_phospho_intensities)){
        empty_list <- cor.test(scale(as.numeric(pre_phospho_intensities_Female[i,])),scale(log2(target_df_Female_phospho[,clin_var[k]])),method="pearson",na.action = "na.exclude")
        baseline_p.val_matrix[i,k] <- empty_list$p.value
        baseline_estimate_matrix[i,k] <- empty_list$estimate
      }
    }
  }, silent = TRUE)
}

colnames(baseline_p.val_matrix) <- paste0(clin_var,"_p.value")
colnames(baseline_estimate_matrix) <- paste0(clin_var,"_estimate")

alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

new_matrix_adjust <- matrix(nrow=nrow(baseline_p.val_matrix),ncol=ncol(baseline_p.val_matrix))
for(i in 1:12){
  new_matrix_adjust[,i] <- p.adjust(baseline_p.val_matrix[,i],method="BH")
}
colnames(new_matrix_adjust) <- paste0(clin_var,"_adj.p.value")
complete_p.val_matrix <- alternate.cols(baseline_p.val_matrix,new_matrix_adjust)

complete_cor_matrix <- cbind(complete_p.val_matrix,baseline_estimate_matrix)
complete_cor_frame_female_phospho <- data.frame(complete_cor_matrix)
complete_cor_frame_female_phospho <- cbind(rownames(pre_phospho_intensities_Female), complete_cor_frame_female_phospho)

complete_cor_frame_female_phospho <- complete_cor_frame_female_phospho[,c(1:3,26,4:5,27,6:7,28,8:9,29,10:11,30,12:13,31,
                                                                          14:15,32,16:17,33,18:19,34,20:21,35,22:23,36,
                                                                          24:25,37)]
colnames(complete_cor_frame_female_phospho)[1] <- "ID"

complete_cor_frame_female_phospho <- complete_cor_frame_female_phospho[,c(1,2:4)]
colnames(complete_cor_frame_female_phospho) <- c("ID","p.val","adj.pval","estimate")


plot(complete_cor_frame_female_phospho$M..µmol..kg.min.._estimate,complete_cor_frame_male_phospho$M..µmol..kg.min.._estimate)
Phos_Comb_gender_cor <- cbind(complete_cor_frame_female_phospho,complete_cor_frame_male_phospho)

colnames(Phos_Comb_gender_cor)[5:8] <- paste0("males_",colnames(Phos_Comb_gender_cor)[5:8])
Phos_Comb_gender_cor$sign <- "NO"
Phos_Comb_gender_cor[which(Phos_Comb_gender_cor$ID %in% complete_cor_frame_sign$ID),]$sign <- "Yes"
Phos_Comb_gender_cor$distance <- sqrt(Phos_Comb_gender_cor$estimate^2 + Phos_Comb_gender_cor$males_estimate^2)

ggplot(Phos_Comb_gender_cor, aes(x=estimate, y=males_estimate, color=distance)) +
  geom_point(alpha=0.75) +
  scale_color_gradient(low="lightgrey", high="black") +  # Adjust colors as needed
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.75, alpha=0.15, fill="darkgrey") +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none") + 
  xlab("Females Insulin Sensitivity (Tau)") + 
  ylab("Males Insulin Sensitivity (Tau)") +
  annotate(geom="text",x=-0.55,y=0.2, label="tau = 0.14", color = "black",size=4) +
  annotate(geom = "text",x=-0.55,y=0.1,label="P < 2.2e-16",size=4)



cor.test(Phos_Comb_gender_cor$males_estimate, Phos_Comb_gender_cor$estimate, method="kendall")

Female_Prot_PRE_GIR_Cor <- Female_Prot_complete_cor_frame[,c(1,2:4)]
colnames(Female_Prot_PRE_GIR_Cor) <- c("ID","p.val","adj.pval","estimate")
#### Validation cohort phosphoproteome loading data ####
KH_validation_phospho <- read.delim("20231111_085108_DIAphosphoPRI_SNv18.4_KHval_cohort_plusSearchArchive_Collapsed.txt")
KH_validation_phosphopeptide <- read.delim("20231111_085108_DIAphosphoPRI_SNv18.4_KHval_cohort_plusSearchArchive_Report.txt")

KH_clinical_data <- read.delim("KH_validation_clinical_data_for_R.txt", dec=",")
KH_clinical_data_additional <- read.delim("KH_validation_clinical_data_for_R_additional_data.txt", dec=",")

KH_clinical_data_additional$Fast.insul_converted <- round(KH_clinical_data_additional$Fast.insul/6,2) #converting to mIU/L which is required for HOMA-IR calculation
KH_clinical_data_additional$HOMA_IR <- round(KH_clinical_data_additional$Fast.gluco*KH_clinical_data_additional$Fast.insul_converted/22.5,2)


# Join the additional data with the main data
KH_clinical_data_combined <- KH_clinical_data %>%
  left_join(KH_clinical_data_additional, 
            by = c("ID" = "PatientID"))

# Display the combined data
KH_clinical_data_combined <- KH_clinical_data_combined[,-c(13:19)]
colnames(KH_clinical_data_combined) <- sub(".x","",colnames(KH_clinical_data_combined))


KH_clinical_data_combined$run <- 1:92
KH_clinical_data_combined$SA_ID <- paste0(word(KH_clinical_data_combined$ID,2,sep="-"),"_",KH_clinical_data_combined$code_diab,"_",KH_clinical_data_combined$Clamp)
colnames(KH_validation_phospho)[1:92] <- KH_clinical_data_combined$SA_ID

#### Figure 1 S2A - Number of identified sites/peptides/proteins ####
length(grepl("Phospho",phospho_study2$EG.PrecursorId)) ### number of total identified peptides
phospho_study2 <- KH_validation_phosphopeptide[grepl("Phospho",KH_validation_phosphopeptide$EG.PrecursorId),] ### filter for phosphopeptides
### number of phosphopeptides identified
length(unique(phospho_study2$EG.PrecursorId))
length(unique(KH_validation_phospho$PG.ProteinGroups)) ## number of phosphoproteins

#phosphosites found in at least 5 samples
total_SA <- rep(1,92)
rownames(KH_validation_phospho) <- KH_validation_phospho$PTM_collapse_key
Identified_sites <- selectGrps(KH_validation_phospho[,c(1:92)], total_SA, 0.05, n=1)

tyrosine_sites <- length(which(grepl("_Y",rownames(Identified_sites)))) ## number of phosphoserine
threonine_sites <- length(which(grepl("_T",rownames(Identified_sites)))) ## number of phosphoserine
serine_sites <- length(which(grepl("_S",rownames(Identified_sites)))) ## number of phosphoserine

length(unique(word(KH_validation_phospho[which(KH_validation_phospho$PTM_collapse_key %in% rownames(Identified_sites)),]$PG.UniProtIds,1,sep=";"))) #number of phosphoproteins

#phosphopeptides found in at least 5 samples
required_samples <- 5
total_samples <- 92
# Group the data by peptide and count the unique sample IDs for each peptide
peptide_counts <- aggregate(R.Replicate ~ EG.PrecursorId, data = phospho_study2, FUN = function(x) length(unique(x)))
num_peptides <- sum(peptide_counts$R.Replicate >= required_samples)
print(num_peptides)


site_data <- data.frame(numbers = c(tyrosine_sites, threonine_sites, serine_sites),
                        residue = c("tyrosine","threonine","serine"),
                        percentage = c(paste0(round(tyrosine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%"),
                                       paste0(round(threonine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%"),
                                       paste0(round(serine_sites/(tyrosine_sites+threonine_sites+serine_sites)*100,1),"%")))

ggplot(site_data, aes(x=residue,y=numbers)) + geom_col(fill="#e39aa2") +
  geom_text(aes(label = percentage), vjust = -0.5, size=6) + theme_minimal() +
  theme(aspect.ratio=1/3, panel.grid = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black",size=14)) + xlab("") + coord_flip()

#### Validation phosphoproteome processing ####
KH_val_log2_phospho <- log2(KH_validation_phospho[,1:92]) #log2 transform data
rownames(KH_val_log2_phospho) <- KH_validation_phospho$PTM_collapse_key
complete_KH_val_phospho <- na.omit(KH_val_log2_phospho)

t_complete_KH_val_phospho <- as.data.frame(t(complete_KH_val_phospho))
t_complete_KH_val_phospho <- cbind(t_complete_KH_val_phospho, KH_clinical_data)
pca_res <- prcomp(t_complete_KH_val_phospho[,1:2085], scale. = TRUE)

t_complete_KH_val_phospho <- data.frame(pca_res$x)
t_complete_KH_val_phospho$batch <- factor(KH_clinical_data$batch)

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(pca_res), "(", paste( as.character(percentage), "%", ")", sep="") )


p<-ggplot(t_complete_KH_val_phospho,aes(x=PC1,y=PC2, color = batch))
p<-p+geom_point(size = 2, alpha = 0.5) + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  stat_ellipse(level=0.8)
p + scale_color_manual(values=c("#9c4141", "#4c4794"))


corrected_KH_val_phospho <- removeBatchEffect(complete_KH_val_phospho, batch = KH_clinical_data$batch)
corrected_KH_val_phospho <- medianScaling(corrected_KH_val_phospho)
t_complete_KH_val_phospho <- as.data.frame(t(corrected_KH_val_phospho))
t_complete_KH_val_phospho <- cbind(t_complete_KH_val_phospho, KH_clinical_data)
pca_res <- prcomp(t_complete_KH_val_phospho[,1:2085], scale. = TRUE)
write_rds(pca_res,"pca_res_phos_validation.rds")

t_complete_KH_val_phospho <- data.frame(pca_res$x)
t_complete_KH_val_phospho$batch <- factor(KH_clinical_data$batch)
t_complete_KH_val_phospho$disease <- factor(KH_clinical_data$code_diab)
t_complete_KH_val_phospho$gender <- factor(KH_clinical_data$Gender)
t_complete_KH_val_phospho$GIR <- as.numeric(KH_clinical_data$GIR)
t_complete_KH_val_phospho$phenotype <- factor(KH_clinical_data$Phenotype)
t_complete_KH_val_phospho$M.value <- as.numeric(KH_clinical_data$M.value)
t_complete_KH_val_phospho$M.value <- as.numeric(KH_clinical_data)

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(t_complete_KH_val_phospho), "(", paste( as.character(percentage), "%", ")", sep="") )


Phospho_PC_clinical <- cbind(t_complete_KH_val_phospho,KH_clinical_data_combined)
write.table(Phospho_PC_clinical,"Phospho_PC_clinical.txt",row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### Figure 1 S2F ####
p<-ggplot(t_complete_KH_val_phospho,aes(x=PC1,y=PC2, color = M.value))
p<-p+geom_point(size = 4, alpha = 1) + theme + xlab(percentage[1]) + ylab(percentage[2])
p + scale_color_gradient(low="lightblue", high="darkred")




#### Validation cohort differentially regulated phosphosites ####
rownames(KH_val_log2_phospho) <- KH_validation_phospho$PTM_collapse_key
Filt_KH_phospho <- selectOverallPercent(KH_val_log2_phospho, 0.25)

# Removing batch effect
Filt_KH_phospho <- removeBatchEffect(Filt_KH_phospho, batch = KH_clinical_data$batch)
# Median normalising data
Filt_KH_phospho <- medianScaling(Filt_KH_phospho)

Clamp <- factor(KH_clinical_data$Clamp, levels=c("Pre","Post"))
Disease <- factor(KH_clinical_data$code_diab, levels=c("con","t2d"))
Gender <- factor(KH_clinical_data$Gender, levels=c("Male","Female"))
design <- model.matrix(~ Disease + Clamp + Gender)
colnames(design) <- c("Intercept","Disease","Clamp","Gender")

corfit <- duplicateCorrelation(Filt_KH_phospho, design, block=KH_clinical_data$ID)
corfit$consensus

fit <- lmFit(Filt_KH_phospho,design,block=KH_clinical_data$ID,correlation=corfit$consensus)
fit2 <- eBayes(fit)

KH_Disease_main_effect_phos <- topTable(fit2, 
                                        coef=2, 
                                        number = Inf, sort.by = "none")
KH_Disease_main_effect_phos$ID <- rownames(KH_Disease_main_effect_phos)

KH_Insulin_main_effect_phos <- topTable(fit2, 
                                        coef=3, 
                                        number = Inf, sort.by = "none")

KH_Insulin_main_effect_phos$ID <- rownames(KH_Insulin_main_effect_phos)

KH_Gender_main_effect_phos <- topTable(fit2, 
                                       coef=4, 
                                       number = Inf, sort.by = "none")

KH_Gender_main_effect_phos$ID <- rownames(KH_Gender_main_effect_phos)



#### Figure 3 S1C - Baseline Phosphosite X insulin sensitivity associations across cohots ####
KH_pre_phospho_intensities <- Filt_KH_phospho[,which(KH_clinical_data$Clamp == "Pre")]
KH_pre_Sample_ID <- KH_clinical_data[which(KH_clinical_data$Clamp == "Pre"),]

baseline_p.val_matrix <- matrix(nrow=nrow(KH_pre_phospho_intensities),ncol=1)
baseline_estimate_matrix <- matrix(nrow=nrow(KH_pre_phospho_intensities),ncol=1)

empty_list <- list()
for(i in 1:length(baseline_p.val_matrix)){
  try({
    empty_list <- cor.test(as.numeric(KH_pre_phospho_intensities[i,]), log2(as.numeric(KH_pre_Sample_ID$M.value)),method="kendall",na.action = "na.exclude")
    baseline_p.val_matrix[i,1] <- empty_list$p.value
    baseline_estimate_matrix[i,1] <- empty_list$estimate
  }, silent = TRUE)
}

alternate.cols <- function(m1, m2) {
  cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

new_matrix_adjust <- matrix(nrow=nrow(baseline_p.val_matrix),ncol=ncol(baseline_p.val_matrix))
new_matrix_adjust <- p.adjust(baseline_p.val_matrix,method="BH")

complete_p.val_matrix <- alternate.cols(baseline_p.val_matrix,new_matrix_adjust)

KH_complete_cor_matrix <- cbind(complete_p.val_matrix,new_matrix_adjust,baseline_estimate_matrix)
KH_complete_cor_frame <- data.frame(KH_complete_cor_matrix)
KH_complete_cor_frame$ID <- rownames(KH_pre_phospho_intensities)

PRE_cohort_correlations <- inner_join(KH_complete_cor_frame,complete_cor_frame,by="ID")
PRE_cohort_correlations_sign <- PRE_cohort_correlations[PRE_cohort_correlations$M..µmol..kg.min.._adj.p.value <= 0.05,]
PRE_cohort_correlations_sign$label <- paste(word(PRE_cohort_correlations_sign$ID,1,sep="_"),word(PRE_cohort_correlations_sign$ID,2,sep="_"))

ggplot(PRE_cohort_correlations_sign, aes(x=M..µmol..kg.min.._estimate,y=V3, label=label, color=abs(V3))) + geom_point(size=3,alpha=0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.position="none", aspect.ratio=1/1.2) + scale_color_gradient(low="grey", high="darkred") +
  ggtitle("Phosphosite insulin sensitivity associations") + geom_text_repel() + xlim(-0.6,0.6) + ylim(-0.6,0.6) +
  xlab("Discovery cohort \n (coefficient)") + ylab("Validation cohort \n (coefficient) ") +
  annotate(geom="text",x=-0.25,y=0.45, label="tau = 0.24", color = "black",size=4) +
  annotate(geom = "text",x=-0.25,y=0.38,label="p = 3.1e-04",size=4)

cor.test(PRE_cohort_correlations_sign$M..µmol..kg.min.._estimate,PRE_cohort_correlations_sign$V3, method="kendall")
nrow(PRE_cohort_correlations_sign[PRE_cohort_correlations_sign$M..µmol..kg.min.._estimate > 0 & PRE_cohort_correlations_sign$V3 > 0,])
nrow(PRE_cohort_correlations_sign[PRE_cohort_correlations_sign$M..µmol..kg.min.._estimate < 0 & PRE_cohort_correlations_sign$V3 < 0,])
nrow(PRE_cohort_correlations_sign)
(19+75)/101*100


#### Figure 5 S1D - LogFc sex correlation across cohorts ####
Gender_overlap_datasets <- inner_join(Gender_main_effect_phos,KH_Gender_main_effect_phos,by="ID")
sign_Gender_overlap_datasets <- Gender_overlap_datasets[Gender_overlap_datasets$adj.P.Val.x <= 0.05 | Gender_overlap_datasets$adj.P.Val.y <= 0.05,]
ggplot(sign_Gender_overlap_datasets, aes(x=logFC.x, y=logFC.y)) + geom_point()

sign_Gender_overlap_datasets$short_ID <- paste(word(sign_Gender_overlap_datasets$ID,1,sep="_"),word(sign_Gender_overlap_datasets$ID,2,sep="_"))

ggplot(sign_Gender_overlap_datasets, aes(x=logFC.x, y=logFC.y, label=short_ID)) + geom_point(alpha=0.4, size=4, shape=1,color="darkred", stroke=0.4) +
  #geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(),axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1.2) + scale_color_manual(values=c("darkred", "darkgrey")) +
  ggtitle("Sex specific phosphoproteome") + geom_text_repel() +
  xlab("Discovery cohort \n Log2 FC (Female-Male") + ylab("Validation cohort \n Log2 FC (Female-Male) ") + xlim(-1.5,1.5) + ylim(-1.6,1.6) +
  annotate(geom="text",x=-1,y=1.2, label="tau = 0.42", color = "black",size=4) +
  annotate(geom = "text",x=-1,y=1,label="p < 2.2e-16",size=4)

cor.test(sign_Gender_overlap_datasets$logFC.x,sign_Gender_overlap_datasets$logFC.y, method="kendall")

nrow(sign_Gender_overlap_datasets[sign_Gender_overlap_datasets$logFC.x < 0 & sign_Gender_overlap_datasets$logFC.y < 0,])
nrow(sign_Gender_overlap_datasets[sign_Gender_overlap_datasets$logFC.x > 0 & sign_Gender_overlap_datasets$logFC.y > 0,])
125+114
239/294*100