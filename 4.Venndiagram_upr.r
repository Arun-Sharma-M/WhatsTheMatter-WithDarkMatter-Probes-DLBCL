#install.packages("ggVennDiagram")

# install.packages("devtools")
#devtools::install_github("gaospecial\\ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)

###############Padj<0.05 and L2FC >2 #############################################
highexp1057_DLBCL_SRA_lncRNA<-read.table("DE_RNA\\pvalue_pv_l2fc_upr_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA_AF.csv", header = TRUE, sep = ",")
colnames(highexp1057_DLBCL_SRA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_DLBCL_SRA_lncRNA_list <- highexp1057_DLBCL_SRA_lncRNA[,-(4:193)]
highexp1057_DLBCL_SRA_lncRNA_list <- highexp1057_DLBCL_SRA_lncRNA_list[,-(1:2)]
#rownames(highexp1057_DLBCL_SRA_lncRNA_list) <- highexp1057_DLBCL_SRA_lncRNA[,1]

highexp75_DLBCL_TCGA_mRNA<-read.table("DE_RNA\\pv_L2FC2_upr_dwnr_DESeq2_TPM_DLBCL2_T_mRNA.csv", header = TRUE, sep = ",")
colnames(highexp75_DLBCL_TCGA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp75_DLBCL_TCGA_mRNA_list <- highexp75_DLBCL_TCGA_mRNA[,-(2:191)]
#rownames(highexp75_DLBCL_TCGA_mRNA_list) <- highexp75_DLBCL_TCGA_mRNA[,1]

x_lncRNA_mRNA <- list(lncRNA=sample(highexp1057_DLBCL_SRA_lncRNA_list),
          mRNA=sample(highexp75_DLBCL_TCGA_mRNA_list))
ggVennDiagram(x_lncRNA_mRNA, label_alpha=0) + scale_fill_gradient(low="green",high = "yellow", labs(title = "Four sets Venn Diagram"))      # VennDiagram 1

Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC <- intersect(highexp1057_DLBCL_SRA_lncRNA_list, highexp75_DLBCL_TCGA_mRNA_list)
Merge_SRA_TCGA_lncRNA_Pvalue_L2FC<-data.frame(Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC)
write.csv(Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC,"Venn_Dia\\1.Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC.csv")


###############Padj<0.05 and L2FC >2_MOD #############################################
MOD_DLBCL_SRA_lncRNA<-read.table("DE_RNA\\7.pvalue_pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_AF_MOD.csv", header = TRUE, sep = ",")
colnames(MOD_DLBCL_SRA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
MOD_DLBCL_SRA_lncRNA_list <- MOD_DLBCL_SRA_lncRNA[,-(4:193)]
MOD_DLBCL_SRA_lncRNA_list <- MOD_DLBCL_SRA_lncRNA_list[,-(1:2)]
#rownames(highexp1057_DLBCL_SRA_lncRNA_list) <- highexp1057_DLBCL_SRA_lncRNA[,1]

MOD_DLBCL_TCGA_mRNA<-read.table("DE_RNA\\pv_L2FC2_upr_dwnr_DESeq2_TPM_DLBCL2_T_mRNA_MOD.csv", header = TRUE, sep = ",")
colnames(MOD_DLBCL_TCGA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
MOD_DLBCL_TCGA_mRNA_list <- MOD_DLBCL_TCGA_mRNA[,-(2:191)]
#rownames(highexp75_DLBCL_TCGA_mRNA_list) <- highexp75_DLBCL_TCGA_mRNA[,1]

x_lncRNA_mRNA_MOD <- list(lncRNA_MOD=sample(MOD_DLBCL_SRA_lncRNA_list),
                      mRNA_MOD=sample(MOD_DLBCL_TCGA_mRNA_list))
ggVennDiagram(x_lncRNA_mRNA_MOD, label_alpha=0) + scale_fill_gradient(low="blue",high = "red", labs(title = "Four sets Venn Diagram"))      # VennDiagram 1

Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC_MOD <- intersect(MOD_DLBCL_SRA_lncRNA_list, MOD_DLBCL_TCGA_mRNA_list)
Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC_MOD<-data.frame(Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC_MOD)
write.csv(Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC_MOD,"Venn_Dia\\2.Merge_SRA_TCGA_lncRNA_mRNA_Pvalue_L2FC_MOD.csv")
