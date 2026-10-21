#install.packages("htmltools")

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(version = "3.21")

#install.packages("locfit")
#install.packages("apeglm")
#BiocManager::install("pasilla")
#BiocManager::install("apeglm")
#BiocManager::install("DEGreport")
#install.packages("pheatmap")
#install.packages("tidyverse")

library(DESeq2)
library(ggplot2)
library(htmltools)
library(apeglm)
library(vsn)

countData <- read.csv('fordeseq_SRA_TCGA_merge_hisat2_N_T_DLBCL_lncRNA_RC.tsv', header = TRUE, sep = "\t")
head(countData)

metaData <- read.csv('SRA_TCGA_DLBCL_for_PCA_meta.csv', header = TRUE, sep = ",")
metaData

#countData <- countData[,-1]
#metaData <- factor(c(rep("normal",5),rep("cancer",29)))

#Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Age+SampleType, tidy = TRUE)
## converting counts to integer mode
#Design specifies how the counts from each gene depend on our variables in the metadata
#For this dataset the factor we care about is our treatment status (SampleType)
#tidy=TRUE argument, which tells DESeq2 to output the results table with rownames as a first #column called 'row.

#The object looks like:
dds
#
dds <- DESeq(dds)
## estimating size factors (This calculates the relative library depth of each sample)
## estimating dispersions (estimates the dispersion of counts for each gene )
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
#nbinomWaldTest (calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs)

resultsNames(dds) # lists the coefficients
# Take a look at the results table

res <- results(dds, name = "SampleType_Tumour_vs_Normal")
write.csv(as.data.frame(res), 
          file="Deseq2_results_lncRNA/DESeq2_TvsN_SRA_TCGA_merge_DLBCL_lncRNA.csv")
resLFC <- lfcShrink(dds, coef="SampleType_Tumour_vs_Normal", type="apeglm")
resLFC
#res_Tissue <- results(dds, name = "Tissue_male_vs_female")

head(results(dds, tidy=TRUE)) #let's look at the results table

#Summary of differential gene expression
summary(res) #summary of results
#summary(res_Tissue) #summary of results

#Sort summary list by p-value
res <- res[order(res$pvalue),]
head(res)


#MA Plot
  par(mfrow=c(1,2))
  pdf(file="Deseq2_results_lncRNA/1_MA_Plot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
  plotMA(res, ylim=c(-2,2))
  plotMA(resLFC, ylim=c(-2,2))
  dev.off()
 # idx <- identify(res$baseMean, res$log2FoldChange)
 #  rownames(res)[idx]

#Volcano Plot
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  pdf(file="Deseq2_results_lncRNA/2_volcano_Plot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
  #Sort L2FC in Ascending order
  res1 <- res[order(res$log2FoldChange),]
  head(res1)
  #Sort L2FC in Descending order
  res2 <- res[order(-res$log2FoldChange),]
  head(res2)
  
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-20,23)))
  
  # Add colored points: blue if pvalue<0.01, red if log2FC>1 and pvalue<0.05)
  with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="cyan"))
  with(subset(res, pvalue<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  dev.off()
  
  
 # with(res_Tissue, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,23)))
 # with(subset(res_Tissue, pvalue<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="cyan"))
 # with(subset(res_Tissue, pvalue<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))

#Effects of transformations on the variance
#PCA
  #First we need to transform the raw count data
  #vst function will perform variance stabilizing transformation
  
  vsdata <- vst(dds, blind=FALSE)
  #rld <- rlog(dds, blind=FALSE)
  pdf(file="Deseq2_results_lncRNA/3_PCA_Plot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
  plotPCA(vsdata, intgroup="SampleType") #using the DESEQ2 plotPCA fxn we can
  plotPCA(vsdata, intgroup="Tissue") #using the DESEQ2 plotPCA fxn we can
  plotPCA(vsdata, intgroup="dataset_SampleType") #using the DESEQ2 plotPCA fxn we can
  plotPCA(vsdata, intgroup="Age") #using the DESEQ2 plotPCA fxn we can
  
  pcaData <- plotPCA(vsdata, intgroup=c("dataset_SampleType", "Tissue"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=dataset_SampleType, shape=Tissue)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
#df <- data.frame(vsdata)
  pcaData <- plotPCA(vsdata, intgroup=c("Age", "dataset_SampleType"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=Age, shape=dataset_SampleType)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  dev.off()
  
#DESeq2 Normalization
ntd <- normTransform(dds)
pdf(file="Deseq2_results_lncRNA/4_meanSD_Plot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsdata))
dev.off()

#Heatmap
library("pheatmap")
#pdf(file="Deseq2_results_lncRNA/5_HmapPlot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("samplename","dataset_SampleType","Tissue")])
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsdata)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
#dev.off()

#pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
# "lnc-ACVRL1-1:3","lnc-IGFBP7-1:2","lnc-CALCRL-1:1","lnc-USP2-2:3","lnc-DIRC3-3:1","lnc-RNF26-2:1"
#BiocManager::install("GEOquery")

library(DESeq2)
library(GEOquery)
library(canvasXpress)
library(ggplot2)
library(clinfun)
library(GGally)
library(factoextra)
norm_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(norm_counts), 
          file="Deseq2_results_lncRNA/DESeq2_TvsN_normalized_lncRNA_SRA_TCGA_merge_DLBCL.csv")
pdf(file="Deseq2_results_lncRNA/6_geomPlot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
{
  lnc_IGFBP7_1 <- ggplot(NULL, aes(x=metaData$SampleType, y=log2(t(norm_counts["lnc-IGFBP7-1:2",]+1)))) +
    geom_jitter(aes(shape=metaData$Tissue, color=metaData$SampleType), size=3)+
    xlab(NULL) +
    ylab("lnc-IGFBP7-1:2 expression \n log2 (norm counts +1)") +
    theme(legend.position = "bottom") +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title =element_text(size = 25),
          legend.position = 'none') +
    stat_summary(fun=mean,
                 geom="point",
                 shape= '_',
                 size=14,
                 colour= c('#b53d35', '#066e70'))
  
  lnc_IGFBP7_1
  }
lnc_IGFBP7_1stat <- wilcox.test(norm_counts["lnc-IGFBP7-1:2",] ~ metaData$SampleType, paired = FALSE)
lnc_IGFBP7_1stat

#####

res.pca <- prcomp(t(log2(norm_counts[c("lnc-IGFBP7-1:2","lnc-USP2-2:3","lnc-DIRC3-3:1","lnc-ACVRL1-1:3","lnc-CALCRL-1:1","lnc-RNF26-2:1"),]+1)), 
                  scale = TRUE)

raw_counts <- as.matrix(countData) 
rownames(metaData) <- metaData$sample
all(rownames(metaData) %in% colnames(raw_counts))
all(colnames(raw_counts) %in% rownames(metaData))

#####
{
  p<- fviz_pca_biplot(res.pca, col.ind = metaData$SampleType,
                      geom = "point",
                      addEllipses = TRUE, 
                      palette = c('#F8766D', '#00BFC4'), 
                      title='Principal Component Analysis')
  p
}
dev.off()

## 3D Scatter Plot: "lnc-IGFBP7-1:2","lnc-USP2-2:3","lnc-DIRC3-3:1","lnc-ACVRL1-1:3","lnc-CALCRL-1:1","lnc-RNF26-2:1"
table_3Dplot <- cbind(t(log2(norm_counts[c("lnc-IGFBP7-1:2","lnc-USP2-2:3","lnc-DIRC3-3:1"),]+1)), as.data.frame(metaData$SampleType))
colnames(table_3Dplot)[4] <- "SampleType"
library(canvasXpress)
canvasXpress(
  data=t(log2(norm_counts[c("lnc-IGFBP7-1:2","lnc-USP2-2:3","lnc-DIRC3-3:1"),]+1)),
  varAnnot=as.data.frame(metaData$SampleType, row.names=rownames(metaData)),
  axisTickScaleFontFactor=0.6,
  axisTitleScaleFontFactor=0.6,
  ellipseBy="metaData$SampleType",
  colorBy="metaData$SampleType",
  colorKey=list("metaData$SampleType"=list("Tumour"="#F8766D", "Normal"="#00BFC4")),
  graphType="Scatter3D",
  title="3D scatter plot",
  xAxis=list("lnc-IGFBP7-1:2"),
  yAxis=list("lnc-USP2-2:3"),
  zAxis=list("lnc-DIRC3-3:1"),
  showLessFit = FALSE)
##Boxplot
par(mfrow=c(2,3))
#pdf(file="Deseq2_results_lncRNA/7_box_Plot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")


ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-IGFBP7-1:2",])) +
  geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
  xlab("lnc-IGFBP7-1:2") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+ 
  scale_fill_manual(values = c("#09E359", "#E31009")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-IGFBP7-1:2','pvalue'])
#[1] 2.569782e-83

  # ,,"lnc-ACVRL1-1:3","lnc-IGFBP7-1:2","lnc-CALCRL-1:1","lnc-USP2-2:3","lnc-DIRC3-3:1","lnc-RNF26-2:1"
  ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-USP2-2:3",])) +
    geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
    xlab("lnc-USP2-2:3") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+
    scale_fill_manual(values = c("#09E359", "#E31009")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-USP2-2:3','pvalue'])
  # [1] 1.051011e-55

  ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-DIRC3-3:1",])) +
    geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
    xlab("lnc-DIRC3-3:1") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+
    scale_fill_manual(values = c("#09E359", "#E31009")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-DIRC3-3:1','pvalue'])
# [1] 3.544215e-55

  ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-ACVRL1-1:3",])) +
    geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
    xlab("lnc-ACVRL1-1:3") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+
    scale_fill_manual(values = c("#09E359", "#E31009")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-ACVRL1-1:3','pvalue'])
# [1] 3.623477e-54

  ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-CALCRL-1:1",], stat = "boxplot")) + 
    geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
    xlab("lnc-CALCRL-1:1") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+
    scale_fill_manual(values = c("#09E359", "#E31009")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-CALCRL-1:1','pvalue'])
# [1] 4.076412e-54

ggplot(NULL, aes(x=metaData$SampleType, y=norm_counts["lnc-RNF26-2:1",], stat = "boxplot")) + 
  geom_boxplot(aes(fill = metaData$SampleType), position = position_dodge(0.9),  notch = FALSE,) +
  xlab("lnc-RNF26-2:1") + ylab("DESeq2_norm_counts") + geom_jitter(width = 0.2)+
  scale_fill_manual(values = c("#09E359", "#E31009")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(res['lnc-RNF26-2:1','pvalue'])
#dev.off()
#[1] 2.928397e-52

#geom_boxplot
#ggplot
  par(mfrow=c(1,1))

# load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

#add_column_name_for_first_row
norm_counts <- norm_counts %>% 
  data.frame() %>%
  rownames_to_column(var="lncRNASymbol") %>% 
  as_tibble()

# Merge DESeq and Counts data matrix 
DESseq_result <- read.csv('Deseq2_results_lncRNA/DESeq2_TvsN_SRA_TCGA_merge_DLBCL_lncRNA.csv', header = TRUE, sep = ",")
colnames(DESseq_result)[1]  <- "lncRNASymbol" 
res_tableOE_tb = merge(DESseq_result, countData, ID="lncRNASymbol", sort = F)

##
##
################################# UPREGULATED ##################################
## 
  DESeq2_pvalue_l2fc_upregulated_filter <- res_tableOE_tb %>% 
  filter(pvalue < 0.05) %>% 
  arrange(-log2FoldChange) %>% #Arrange rows by pvalue values
  pull(lncRNASymbol) %>% 		#Extract character vector of ordered genes
  head(n=50) 		#Extract the first 50 genes
  
## normalized counts for top 50 significant genes
top50_sigOE_norm_upregulated <- norm_counts %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_l2fc_upregulated_filter)
head(top50_sigOE_norm_upregulated)

#export
deseqcountsmerge =  merge(DESseq_result, norm_counts, ID="lncRNASymbol", sort = F)
top50_sigOE_norm_upregulated_export <- deseqcountsmerge %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_l2fc_upregulated_filter)

write.csv(as.data.frame(top50_sigOE_norm_upregulated_export), 
          file="Deseq2_results_lncRNA/8_upr_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.csv")

# Gathering the columns to have normalized counts to a single column
gathered_top50_sigOE <- top50_sigOE_norm_upregulated %>%
  gather(colnames(top50_sigOE_norm_upregulated)[2:67], key = "samplename", value = "DESeq2_norm_counts")

## check the column header in the "gathered" data frame
View(gathered_top50_sigOE)

gathered_top50_sigOE <- inner_join(metaData, gathered_top50_sigOE)
View(gathered_top50_sigOE)
## plot using ggplot2
pdf(file="Deseq2_results_lncRNA/8_upr_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
ggplot(gathered_top50_sigOE) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color = SampleType)) +
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 and l2fc sorted _ Upregulated_DLBCL_SRA_TCGA_merge") +
  theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
#
## plot using boxplot
ggplot(gathered_top50_sigOE) + 
  geom_boxplot(aes(x = lncRNASymbol, y = DESeq2_norm_counts, fill=SampleType), width=0.5,) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color=SampleType),position = position_dodge(width=0.5))+
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  scale_color_manual(values=c("#ff00cc", "#3300cc")) +   ####point color#####
ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 and l2fc sorted _ Upregulated_DLBCL_SRA_TCGA_merge") +
  theme_dark() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
dev.off()
#
#
################################# DOWNREGULATED #################################
#
#
DESeq2_pvalue_l2fcdownregulated_filter <- res_tableOE_tb %>% 
  filter(pvalue < 0.05) %>% 
  arrange(log2FoldChange) %>% #Arrange rows by pvalue values
  pull(lncRNASymbol) %>% 		#Extract character vector of ordered genes
  head(n=50) 		#Extract the first 50 genes

## normalized counts for top 50 significant genes

top50_sigDR_norm_downregulated <- norm_counts %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_l2fcdownregulated_filter)
head(top50_sigDR_norm_downregulated)

#export
top50_sigOE_norm_downregulated_export <- deseqcountsmerge %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_l2fcdownregulated_filter)

write.csv(as.data.frame(top50_sigOE_norm_downregulated_export), 
          file="Deseq2_results_lncRNA/9_dwnr_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.csv")

# Gathering the columns to have normalized counts to a single column
gathered_top50_sigDR <- top50_sigDR_norm_downregulated %>%
  gather(colnames(top50_sigDR_norm_downregulated)[2:67], key = "samplename", value = "DESeq2_norm_counts")

## check the column header in the "gathered" data frame
View(gathered_top50_sigDR)

gathered_top50_sigDR <- inner_join(metaData, gathered_top50_sigDR)
View(gathered_top50_sigDR)
## plot using ggplot2
pdf(file="Deseq2_results_lncRNA/9_dwnr_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
ggplot(gathered_top50_sigDR) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color = SampleType)) +
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 and l2fc sorted _ Downregulated_DLBCL_SRA_TCGA_merge") +
  theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
#
## plot using boxplot
ggplot(gathered_top50_sigDR) + 
  geom_boxplot(aes(x = lncRNASymbol, y = DESeq2_norm_counts, fill=SampleType), width=0.5,) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color=SampleType),position = position_dodge(width=0.5))+
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  scale_color_manual(values=c("#ff00cc", "#3300cc")) +   ####point color#####
ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 and l2fc sorted _ Downregulated_DLBCL_SRA_TCGA_merge") +
  theme_dark() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
dev.off()
#
#
################################## Top 50 pvalue #################################
#
#
DESeq2_pvalue_filter <- res_tableOE_tb %>% 
  filter(pvalue < 0.05) %>% #filter rows by pvalue values
  arrange(pvalue) %>% #Arrange rows by pvalue values
  pull(lncRNASymbol) %>% 		#Extract character vector of ordered genes
  head(n=50) 		#Extract the first 50 genes
#
## normalized counts for top 50 significant genes
#
top50_sigDR_norm <- norm_counts %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_filter)
head(top50_sigDR_norm)
#

#export
top50_sigOE_norm_significant_export <- deseqcountsmerge %>%
  filter(lncRNASymbol %in% DESeq2_pvalue_filter)

write.csv(as.data.frame(top50_sigOE_norm_significant_export), 
          file="Deseq2_results_lncRNA/10_pv_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.csv")

# Gathering the columns to have normalized counts to a single column
gathered_top50_sig <- top50_sigDR_norm %>%
  gather(colnames(top50_sigDR_norm)[2:67], key = "samplename", value = "DESeq2_norm_counts")
#
## check the column header in the "gathered" data frame
View(gathered_top50_sig)
#
gathered_top50_sig <- inner_join(metaData, gathered_top50_sig)
View(gathered_top50_sig)
## plot using ggplot2
pdf(file="Deseq2_results_lncRNA/10_dwnr_boxPlot_lncRNA_DLBCL_SRA_TCGA_merge.pdf")
ggplot(gathered_top50_sig) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color = SampleType)) +
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 _DLBCL_SRA_TCGA_merge") +
  theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
#
## plot using boxplot
ggplot(gathered_top50_sig) +
  geom_boxplot(aes(x = lncRNASymbol, y = DESeq2_norm_counts, fill=SampleType), width=0.5,) +
  geom_point(aes(x = lncRNASymbol, y = DESeq2_norm_counts, color=SampleType),position = position_dodge(width=0.5))+
  scale_y_log10() +
  xlab("lncRNASymbol") +
  ylab("log10(DESeq2_Normalized_Counts)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  scale_color_manual(values=c("#ff00cc", "#3300cc")) +   ####point color#####
ggtitle("Top 50 Significant DE lncRNA with pvalue < 0.05 _DLBCL_SRA_TCGA_merge") +
  theme_dark() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.1))
dev.off()
