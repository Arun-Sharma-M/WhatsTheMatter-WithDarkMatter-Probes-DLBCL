# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")

#BiocManager::install('TCGAbiolinks')

#library(TCGAbiolinks)
library(survminer)
library(survival)
library(tidyverse)
library(stringr)

# getting clinical data for TCGA-DLBCL cohort -------------------
#clinical_DLBCL <- GDCquery_clinic("TCGA-DLBCL")
clinical_DLBCL <- read.csv('TCGA_DLBCL_meta.csv', header = TRUE, sep = ",")
#clinical_DLBCL[] <- lapply(clinical_DLBCL, gsub, pattern = "'--", replacement = "", fixed = TRUE)

any(colnames(clinical_DLBCL) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_DLBCL) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
# [1] 13 16 21
clinical_DLBCL[,c(13,16,21)]

# looking at some variables associated with survival 
table(clinical_DLBCL$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_DLBCL$deceased <- ifelse(clinical_DLBCL$vital_status == "Alive", FALSE, TRUE)
table(clinical_DLBCL$deceased)

# Converting Days to Month by dividing Number of days / 30.417

clinical_DLBCL$Months_to_last_follow_up <- clinical_DLBCL$days_to_last_follow_up/30.417
clinical_DLBCL$Months_to_death <- clinical_DLBCL$days_to_death/30.417
#clinical_DLBCL$Months_to_death <- as.numeric(ifelse(clinical_DLBCL$days_to_death == '--', NA, clinical_DLBCL$days_to_death)) / 30.417
#summary(clinical_DLBCL$Months_to_death)


# create an "overall survival" variable that is equal to Months_to_death
# for dead patients, and to Months_to_last_follow_up for patients who are still alive
clinical_DLBCL$overall_survival <- ifelse(clinical_DLBCL$vital_status == "Alive",
                                         clinical_DLBCL$Months_to_last_follow_up,
                                         clinical_DLBCL$Months_to_death)

#Remove entire rows with NA's in important columns with data for KM plot

#clinical_DLBCL <- clinical_DLBCL[-c(49:58), ]
#clinical_DLBCL$samplename.x <- str_replace_all(clinical_DLBCL$samplename.x, '-', '.')

# get mRNA expression data -----------

TCGA_DLBCL_TPM_Data_full <- read.csv('DESEQ2_TPM_DLBCL2_lncRNA_mRNA.csv', header = TRUE, sep = ",",row.names=1)
TCGA_DLBCL_TPM_Data <- TCGA_DLBCL_TPM_Data_full[c(5:52)]
head(TCGA_DLBCL_TPM_Data)

##########************************************************************##########

# 1. Define the genes to be analyzed in a single list
genes_to_analyze <- c("lnc-ADAMTS1-2:1", "lnc-FXYD2-1:5", "lnc-SPARC-2:3", "ADAMTS1", "FXYD2", "SPARC", "ADAMTS5", "ATOX1", "FXYD6-FXYD2", "BCL6")

# 2. Get the mRNA expression data and reshape it to a long format
exp_long <- TCGA_DLBCL_TPM_Data %>%
  as.data.frame() %>%
  rownames_to_column(var = 'GeneSymbol') %>%
  gather(key = 'samplename', value = 'TPM', -GeneSymbol)

# 3. Create a place to store the plots if you want to view them later
median_plots <- list()
quantile_plots <- list()

# 4. Start the loop for each gene
for (gene in genes_to_analyze) {
  
  # a) Filter data for the current gene and merge with clinical data
  gene_data <- exp_long %>%
    filter(GeneSymbol == gene) %>%
    merge(clinical_DLBCL, by.x = 'samplename', by.y = 'samplename')
  
  # b) Handle median-split analysis
  if (nrow(gene_data) > 0) {
    # Denote which cases have higher or lower expression than the median
    gene_data$strata_median <- ifelse(gene_data$TPM >= median(gene_data$TPM), "HIGH", "LOW")
    
    # Fit and plot survival curve for median split
    fit_median <- survfit(Surv(overall_survival, deceased) ~ strata_median, data = gene_data)
    p_median <- ggsurvplot(fit_median,
                           data = gene_data,
                           pval = TRUE,
                           risk.table = TRUE,
                           title = paste("TCGA_DLBCL", gene, "Median_value", sep = "_"),
                           surv.median.line = "hv",
                           xlab = "Follow up in Months")
    
    median_plots[[gene]] <- p_median
    message(paste("Generated median plot for:", gene))
    print(p_median)
  }
  
  # c) Handle quantile-split analysis
  if (nrow(gene_data) > 0) {
    # Denote which cases have higher or lower expression than the 75th percentile
    gene_data$strata_quant <- ifelse(gene_data$TPM > quantile(gene_data$TPM, probs = 0.25), "Medium & High", "Low")
    
    # Fit and plot survival curve for quantile split
    fit_quant <- survfit(Surv(overall_survival, deceased) ~ strata_quant, data = gene_data)
    p_quant <- ggsurvplot(fit_quant,
                          data = gene_data,
                          pval = TRUE,
                          conf.int = TRUE,
                          risk.table = TRUE,
                          title = paste("TCGA_DLBCL", gene, "quantile", sep = "_"),
                          xlab = "Follow up in Months",
                          legend.labs = c("Medium & High", "Low"),
                          ncensor.plot = TRUE,
                          surv.median.line = "hv",
                          ggtheme = theme_bw(),
                          palette = c("#e7002f", "#2E9FDF"))
    
    quantile_plots[[gene]] <- p_quant
    message(paste("Generated quantile plot for:", gene))
    print(p_quant)
  }
}

# 5. (Optional) Arrange and save all generated plots
# This requires the 'ggpubr' package (which is part of survminer)
# To save all median plots to a PDF:
# ggsave("all_median_plots.pdf", arrange_ggsurvplots(median_plots), width = 10, height = 8)
# To save all quantile plots to a PDF:
# ggsave("all_quantile_plots.pdf", arrange_ggsurvplots(quantile_plots), width = 10, height = 8)
