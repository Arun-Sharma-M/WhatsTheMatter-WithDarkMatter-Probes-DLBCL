
library(ggplot2)
#library(htmltools)
#library(apeglm)
library(tidyverse)
#library(ggplot2)
#library(ggrepel)
#library(DEGreport)
#library(RColorBrewer)
#install.packages("ggridges")
library(ggridges)
#library(geom_flat_violin)
#load TPM and metadata

#load mRNA

DLBCL2_mRNA_NT_TPM <- read.csv('D:/PAPERS_ON_CHOL_MCL_DLBCL_NONAMER/DLBCL_2/SRA_TCGA_merge/Workflow_mRNA/RESULTS_workflow/5.All_factors_merged_TPM_DLBCL2_mRNA.csv', header = TRUE, sep = ",")
DLBCL2_mRNA_NT_TPM[1:5,1:73]
DLBCL2_mRNA_NT_TPM <- DLBCL2_mRNA_NT_TPM[!duplicated(DLBCL2_mRNA_NT_TPM$geneSymbol), ]

rownames(DLBCL2_mRNA_NT_TPM) <- DLBCL2_mRNA_NT_TPM$geneSymbol
DLBCL2_mRNA_NT_TPM[1:5,1:5]

DLBCL2_mRNA_NT_TPM <- DLBCL2_mRNA_NT_TPM[ -c(1) ] # Removes entire column be careful
DLBCL2_mRNA_NT_TPM[1:5,1:5]

#load lncRNA

DLBCL2_lncRNA_NT_TPM <- read.csv('D:/PAPERS_ON_CHOL_MCL_DLBCL_NONAMER/DLBCL_2/SRA_TCGA_merge/Workflow_mRNA/Venndia/DE_RNA/pvalue_pv_l2fc_upr_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA_AF.csv', header = TRUE, sep = ",")
DLBCL2_lncRNA_NT_TPM[1:5,1:75]
#DLBCL2_lncRNA_NT_TPM <- DLBCL2_lncRNA_NT_TPM[!duplicated(DLBCL2_lncRNA_NT_TPM$geneSymbol), ]


metaData <- read.csv('SRA_TCGA_DLBCL_for_PCA_meta.csv', header = TRUE, sep = ",")
metaData[1:5,1:5]

#load NHEJ, SLimRNA and SLiLNCRNA genelist
SLimRNA <- read.csv('String_exp_list.csv', header = FALSE, sep = ",")
SLimRNA<- SLimRNA[[1]]

################################# SLimRNA ##################################

DLBCL2_mRNA_NT_TPM_SLimRNA_int <- DLBCL2_mRNA_NT_TPM[rownames(DLBCL2_mRNA_NT_TPM) %in% SLimRNA,]      #

#add_column_name_for_first_row
DLBCL2_mRNA_NT_TPM_SLimRNA_int <- DLBCL2_mRNA_NT_TPM_SLimRNA_int %>% 
  data.frame() %>%
  rownames_to_column(var="geneSymbol") %>% 
  as_tibble()
DLBCL2_mRNA_NT_TPM_SLimRNA <- DLBCL2_mRNA_NT_TPM_SLimRNA_int[,1:73]
DLBCL2_mRNA_NT_TPM_SLimRNA[1:5,1:73]

colnames(DLBCL2_mRNA_NT_TPM_SLimRNA) <- gsub("TPM", "", colnames(DLBCL2_mRNA_NT_TPM_SLimRNA))

upr <- DLBCL2_mRNA_NT_TPM_SLimRNA %>% filter(log2FoldChange >= 2) 

# Gathering the columns to have normalized counts to a single column
gathered_SLimRNA_col_upr <- upr %>%
  gather(colnames(upr)[8:73], key = "samplename", value = "TPM")
gathered_SLimRNA_col_upr[1:5,1:9]

# remove "TPM" from column names to convert TPMC1 to C1 and so on...
gathered_top50_sigOE <- inner_join(metaData, gathered_SLimRNA_col_upr)
View(gathered_top50_sigOE)

dwnr <- DLBCL2_mRNA_NT_TPM_SLimRNA %>% filter(log2FoldChange <= -2) 
# Gathering the columns to have normalized counts to a single column
gathered_SLimRNA_col_dwnr <- dwnr %>%
  gather(colnames(dwnr)[8:73], key = "samplename", value = "TPM")
gathered_SLimRNA_col_dwnr[1:5,1:9]

# remove "TPM" from column names to convert TPMC1 to C1 and so on...
gathered_top50_sigDR <- inner_join(metaData, gathered_SLimRNA_col_dwnr)
View(gathered_top50_sigDR)
################################# SLiLNCRNA ##################################

#DLBCL2_lncRNA_NT_TPM_SLiLNCRNA_int <- DLBCL2_lncRNA_NT_TPM[sapply(DLBCL2_lncRNA_NT_TPM$geneSymbol, function(x) any(x %in% SLimRNA)), ]

DLBCL2_lncRNA_NT_TPM_SLiLNCRNA_int <- DLBCL2_lncRNA_NT_TPM %>%
  filter(geneSymbol %in% SLimRNA)

#add_column_name_for_first_row
#DLBCL2_lncRNA_NT_TPM_SLiLNCRNA_int <- DLBCL2_lncRNA_NT_TPM_SLiLNCRNA_int %>% 
 # data.frame() %>%
  #rownames_to_column(var="geneSymbol") %>% 
  #as_tibble()

DLBCL2_lncRNA_NT_TPM_SLiLNCRNA <- DLBCL2_lncRNA_NT_TPM_SLiLNCRNA_int[,1:75]
DLBCL2_lncRNA_NT_TPM_SLiLNCRNA[1:5,1:75]

colnames(DLBCL2_lncRNA_NT_TPM_SLiLNCRNA) <- gsub("TPM", "", colnames(DLBCL2_lncRNA_NT_TPM_SLiLNCRNA))

upr_lncRNA <- DLBCL2_lncRNA_NT_TPM_SLiLNCRNA %>% filter(log2FoldChange >= 2) 

# Gathering the columns to have normalized counts to a single column
gathered_SLiLNCRNA_col_upr <- upr_lncRNA %>%
  gather(colnames(upr_lncRNA)[10:75], key = "samplename", value = "TPM")

gathered_SLiLNCRNA_OE <- inner_join(metaData, gathered_SLiLNCRNA_col_upr)
View(gathered_SLiLNCRNA_OE)


dwnr_lncRNA <- DLBCL2_lncRNA_NT_TPM_SLiLNCRNA %>% filter(log2FoldChange <= -2) 

# Gathering the columns to have normalized counts to a single column
gathered_SLiLNCRNA_col_dwnr <- dwnr_lncRNA %>%
  gather(colnames(dwnr_lncRNA)[10:75], key = "samplename", value = "TPM")

gathered_SLiLNCRNA_DR <- inner_join(metaData, gathered_SLiLNCRNA_col_dwnr)
View(gathered_SLiLNCRNA_DR)

################################# Plot_BOxplot ##################################

## plot using ggplot2

ggplot(gathered_top50_sigOE) + 
  geom_boxplot(aes(x = geneSymbol, y = TPM, fill=SampleType), width=0.5, outlier.shape = NA) +
  scale_y_log10() +
  xlab("geneSymbol") +
  ylab("log10(TPM)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
ggtitle("TCGA_DLBCL2_mRNA_SLimRNA_Overexpressed") +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(gathered_SLiLNCRNA_OE) + 
  geom_boxplot(aes(x = lncRNASymbol, y = TPM, fill=SampleType), width=0.5, outlier.shape = NA) +
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(TPM)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  ggtitle("TCGA_DLBCL2_lncRNA_SLiLNCRNA_Overexpressed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(gathered_top50_sigDR) + 
  geom_boxplot(aes(x = geneSymbol, y = TPM, fill=SampleType), width=0.5, outlier.shape = NA) +
  scale_y_log10() +
  xlab("geneSymbol") +
  ylab("log10(TPM)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  ggtitle("TCGA_DLBCL2_mRNA_SLimRNA_Downregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(gathered_SLiLNCRNA_DR) + 
  geom_boxplot(aes(x = lncRNASymbol, y = TPM, fill=SampleType), width=0.5, outlier.shape = NA) +
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(TPM)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  ggtitle("TCGA_DLBCL2_lncRNA_SLiLNCRNA_Downregulated") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

# Ridge plot

library(scales)
ggplot(gathered_SLiLNCRNA) + 
  aes(x = TPM, y = geneSymbol, fill=SampleType) +
  geom_density_ridges2(scale = 2) +
  scale_x_log10() +
  scale_fill_viridis_d(alpha = 0.9) +
  xlab("log10(TPM)") +
  ylab("geneSymbol") +
  scale_x_continuous(breaks = breaks_width(5), limits = c(0, 120))
  ggtitle("TCGA_DLBCL2_mRNA_SLiLNCRNA") 
   
# Raincloudplot
  install.packages("PupillometryR")
  library(PupillometryR)
  
  ggplot(gathered_SLiLNCRNA) + 
    aes(x = TPM, y = geneSymbol, fill=SampleType) +
  geom_flat_violin(position = position_nudge(x = .5),
                   alpha = 0.7) +
    geom_point(aes(color = SampleType),
               position = position_jitter(w = .15)) +
    geom_boxplot(width = .5) + # New code added here
    xlab("log10(TPM)") +
    ylab("geneSymbol") +
    ggtitle("TCGA_DLBCL2_mRNA_SLiLNCRNA") 

  
  
  install.packages("theme_pubr")
  library(theme_pubr)
  library(ggpubr)
  ggplot(gathered_SLiLNCRNA) + 
    aes(x = geneSymbol, y = TPM, fill=SampleType) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
    geom_point(aes(color = result), 
               position = position_jitterdodge(jitter.width = .15, # to sepreate jitter into their group.
                                               dodge.width = .3), 
               size = .3, 
               alpha = .03,
               show.legend = F) +
    geom_boxplot(width = .3, 
                 outlier.shape = NA,
                 alpha = .5) +
    facet_wrap(~variable, ncol = 1, scale = "free_x") +  
    
    labs(x = "geneSymbol", 
         y = "log10(TPM)",
         title = "TCGA_DLBCL2_mRNA_SLiLNCRNA") +
    guides(fill = guide_legend(title="Result of the call")) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), na.value = "#5f5f5f") +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), na.value = "#5f5f5f") +
    scale_y_continuous(breaks = breaks_width(20)) +
    theme_pubr(base_size = 11, x.text.angle = 0) +
    theme(legend.position = "right") 
