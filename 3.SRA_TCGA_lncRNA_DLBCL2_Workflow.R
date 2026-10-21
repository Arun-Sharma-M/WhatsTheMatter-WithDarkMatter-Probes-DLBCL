#Input_Entire gene list
##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
countData<-read.table("fordeseq_SRA_TCGA_merge_hisat2_N_T_DLBCL_mRNA_RC.tsv",sep="\t",header=T,row.names=1)
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(countData)
#[1] 38721        34
countData[1:5,1:3]

#Input_gene length # Make sure of mRNA and gene lengths file
##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
length<-read.table("mRNA_read_length.tsv",sep="\t",header=T,row.names=1) 
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(length)
# [1] 38721         2

RC2TPM <- data.frame(countData) # Creating data
######################################################################################################################################
############################################# [READCOUNTS TO TPM CONVERSION] #########################################################
######### Introducing gene_length_in_Kb Column by performing (END - START values by 1000) to get the gene length in Kb ###############
######################################################################################################################################
######################################################################################################################################
RC2TPM[,67] = (length[,2]-length[,1])/1000
colnames(RC2TPM)[67] <- "gene_length_in_Kb"
RC2TPM[RC2TPM == 0] <- 0.1

#cHANGE: 
# 1. RC2TPM$C TO RC2TPM$CC
# 2. RC2TPM$T TO RC2TPM$T
########################### Introducing RPKT Columns ###########################
RC2TPM$RPKC1	<-	RC2TPM$C1	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC1
RC2TPM$RPKC2	<-	RC2TPM$C2	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC2
RC2TPM$RPKC3	<-	RC2TPM$C3	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC3
RC2TPM$RPKC4	<-	RC2TPM$C4	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC4
RC2TPM$RPKC5	<-	RC2TPM$C5	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC5
RC2TPM$RPKC6	<-	RC2TPM$C6	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC6
RC2TPM$RPKC7	<-	RC2TPM$C7	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC7
RC2TPM$RPKC8	<-	RC2TPM$C8	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKC8
RC2TPM$RPKT1	<-	RC2TPM$T1	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT1
RC2TPM$RPKT2	<-	RC2TPM$T2	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT2
RC2TPM$RPKT3	<-	RC2TPM$T3	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT3
RC2TPM$RPKT4	<-	RC2TPM$T4	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT4
RC2TPM$RPKT5	<-	RC2TPM$T5	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT5
RC2TPM$RPKT6	<-	RC2TPM$T6	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT6
RC2TPM$RPKT7	<-	RC2TPM$T7	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT7
RC2TPM$RPKT8	<-	RC2TPM$T8	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT8
RC2TPM$RPKT9	<-	RC2TPM$T9	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT9
RC2TPM$RPKT10	<-	RC2TPM$T10	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT10
RC2TPM$RPKT11	<-	RC2TPM$T11	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT11
RC2TPM$RPKT12	<-	RC2TPM$T12	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT12
RC2TPM$RPKT13	<-	RC2TPM$T13	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT13
RC2TPM$RPKT14	<-	RC2TPM$T14	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT14
RC2TPM$RPKT15	<-	RC2TPM$T15	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT15
RC2TPM$RPKT16	<-	RC2TPM$T16	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT16
RC2TPM$RPKT17	<-	RC2TPM$T17	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT17
RC2TPM$RPKT18	<-	RC2TPM$T18	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT18
RC2TPM$RPKT19	<-	RC2TPM$T19	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT19
RC2TPM$RPKT20	<-	RC2TPM$T20	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT20
RC2TPM$RPKT21	<-	RC2TPM$T21	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT21
RC2TPM$RPKT22	<-	RC2TPM$T22	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT22
RC2TPM$RPKT23	<-	RC2TPM$T23	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT23
RC2TPM$RPKT24	<-	RC2TPM$T24	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT24
RC2TPM$RPKT25	<-	RC2TPM$T25	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT25
RC2TPM$RPKT26	<-	RC2TPM$T26	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT26
RC2TPM$RPKT27	<-	RC2TPM$T27	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT27
RC2TPM$RPKT28	<-	RC2TPM$T28	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT28
RC2TPM$RPKT29	<-	RC2TPM$T29	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT29
RC2TPM$RPKT30	<-	RC2TPM$T30	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT30
RC2TPM$RPKT31	<-	RC2TPM$T31	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT31
RC2TPM$RPKT32	<-	RC2TPM$T32	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT32
RC2TPM$RPKT33	<-	RC2TPM$T33	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT33
RC2TPM$RPKT34	<-	RC2TPM$T34	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT34
RC2TPM$RPKT35	<-	RC2TPM$T35	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT35
RC2TPM$RPKT36	<-	RC2TPM$T36	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT36
RC2TPM$RPKT37	<-	RC2TPM$T37	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT37
RC2TPM$RPKT38	<-	RC2TPM$T38	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT38
RC2TPM$RPKT39	<-	RC2TPM$T39	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT39
RC2TPM$RPKT40	<-	RC2TPM$T40	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT40
RC2TPM$RPKT41	<-	RC2TPM$T41	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT41
RC2TPM$RPKT42	<-	RC2TPM$T42	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT42
RC2TPM$RPKT43	<-	RC2TPM$T43	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT43
RC2TPM$RPKT44	<-	RC2TPM$T44	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT44
RC2TPM$RPKT45	<-	RC2TPM$T45	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT45
RC2TPM$RPKT46	<-	RC2TPM$T46	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT46
RC2TPM$RPKT47	<-	RC2TPM$T47	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT47
RC2TPM$RPKT48	<-	RC2TPM$T48	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT48
RC2TPM$RPKT49	<-	RC2TPM$T49	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT49
RC2TPM$RPKT50	<-	RC2TPM$T50	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT50
RC2TPM$RPKT51	<-	RC2TPM$T51	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT51
RC2TPM$RPKT52	<-	RC2TPM$T52	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT52
RC2TPM$RPKT53	<-	RC2TPM$T53	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT53
RC2TPM$RPKT54	<-	RC2TPM$T54	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT54
RC2TPM$RPKT55	<-	RC2TPM$T55	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT55
RC2TPM$RPKT56	<-	RC2TPM$T56	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT56
RC2TPM$RPKT57	<-	RC2TPM$T57	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT57
RC2TPM$RPKT58	<-	RC2TPM$T58	/	RC2TPM$gene_length_in_Kb	#	Creating	new	variable	RPKT58

######################### Introducing PMRFT Columns ############################
RC2TPMC.1	<-	(sum(RC2TPM[,	'RPKC1'],	na.rm	=	TRUE))/10^6
RC2TPMC.2	<-	(sum(RC2TPM[,	'RPKC2'],	na.rm	=	TRUE))/10^6
RC2TPMC.3	<-	(sum(RC2TPM[,	'RPKC3'],	na.rm	=	TRUE))/10^6
RC2TPMC.4	<-	(sum(RC2TPM[,	'RPKC4'],	na.rm	=	TRUE))/10^6
RC2TPMC.5	<-	(sum(RC2TPM[,	'RPKC5'],	na.rm	=	TRUE))/10^6
RC2TPMC.6	<-	(sum(RC2TPM[,	'RPKC6'],	na.rm	=	TRUE))/10^6
RC2TPMC.7	<-	(sum(RC2TPM[,	'RPKC7'],	na.rm	=	TRUE))/10^6
RC2TPMC.8	<-	(sum(RC2TPM[,	'RPKC8'],	na.rm	=	TRUE))/10^6
RC2TPM.1	<-	(sum(RC2TPM[,	'RPKT1'],	na.rm	=	TRUE))/10^6
RC2TPM.2	<-	(sum(RC2TPM[,	'RPKT2'],	na.rm	=	TRUE))/10^6
RC2TPM.3	<-	(sum(RC2TPM[,	'RPKT3'],	na.rm	=	TRUE))/10^6
RC2TPM.4	<-	(sum(RC2TPM[,	'RPKT4'],	na.rm	=	TRUE))/10^6
RC2TPM.5	<-	(sum(RC2TPM[,	'RPKT5'],	na.rm	=	TRUE))/10^6
RC2TPM.6	<-	(sum(RC2TPM[,	'RPKT6'],	na.rm	=	TRUE))/10^6
RC2TPM.7	<-	(sum(RC2TPM[,	'RPKT7'],	na.rm	=	TRUE))/10^6
RC2TPM.8	<-	(sum(RC2TPM[,	'RPKT8'],	na.rm	=	TRUE))/10^6
RC2TPM.9	<-	(sum(RC2TPM[,	'RPKT9'],	na.rm	=	TRUE))/10^6
RC2TPM.10	<-	(sum(RC2TPM[,	'RPKT10'],	na.rm	=	TRUE))/10^6
RC2TPM.11	<-	(sum(RC2TPM[,	'RPKT11'],	na.rm	=	TRUE))/10^6
RC2TPM.12	<-	(sum(RC2TPM[,	'RPKT12'],	na.rm	=	TRUE))/10^6
RC2TPM.13	<-	(sum(RC2TPM[,	'RPKT13'],	na.rm	=	TRUE))/10^6
RC2TPM.14	<-	(sum(RC2TPM[,	'RPKT14'],	na.rm	=	TRUE))/10^6
RC2TPM.15	<-	(sum(RC2TPM[,	'RPKT15'],	na.rm	=	TRUE))/10^6
RC2TPM.16	<-	(sum(RC2TPM[,	'RPKT16'],	na.rm	=	TRUE))/10^6
RC2TPM.17	<-	(sum(RC2TPM[,	'RPKT17'],	na.rm	=	TRUE))/10^6
RC2TPM.18	<-	(sum(RC2TPM[,	'RPKT18'],	na.rm	=	TRUE))/10^6
RC2TPM.19	<-	(sum(RC2TPM[,	'RPKT19'],	na.rm	=	TRUE))/10^6
RC2TPM.20	<-	(sum(RC2TPM[,	'RPKT20'],	na.rm	=	TRUE))/10^6
RC2TPM.21	<-	(sum(RC2TPM[,	'RPKT21'],	na.rm	=	TRUE))/10^6
RC2TPM.22	<-	(sum(RC2TPM[,	'RPKT22'],	na.rm	=	TRUE))/10^6
RC2TPM.23	<-	(sum(RC2TPM[,	'RPKT23'],	na.rm	=	TRUE))/10^6
RC2TPM.24	<-	(sum(RC2TPM[,	'RPKT24'],	na.rm	=	TRUE))/10^6
RC2TPM.25	<-	(sum(RC2TPM[,	'RPKT25'],	na.rm	=	TRUE))/10^6
RC2TPM.26	<-	(sum(RC2TPM[,	'RPKT26'],	na.rm	=	TRUE))/10^6
RC2TPM.27	<-	(sum(RC2TPM[,	'RPKT27'],	na.rm	=	TRUE))/10^6
RC2TPM.28	<-	(sum(RC2TPM[,	'RPKT28'],	na.rm	=	TRUE))/10^6
RC2TPM.29	<-	(sum(RC2TPM[,	'RPKT29'],	na.rm	=	TRUE))/10^6
RC2TPM.30	<-	(sum(RC2TPM[,	'RPKT30'],	na.rm	=	TRUE))/10^6
RC2TPM.31	<-	(sum(RC2TPM[,	'RPKT31'],	na.rm	=	TRUE))/10^6
RC2TPM.32	<-	(sum(RC2TPM[,	'RPKT32'],	na.rm	=	TRUE))/10^6
RC2TPM.33	<-	(sum(RC2TPM[,	'RPKT33'],	na.rm	=	TRUE))/10^6
RC2TPM.34	<-	(sum(RC2TPM[,	'RPKT34'],	na.rm	=	TRUE))/10^6
RC2TPM.35	<-	(sum(RC2TPM[,	'RPKT35'],	na.rm	=	TRUE))/10^6
RC2TPM.36	<-	(sum(RC2TPM[,	'RPKT36'],	na.rm	=	TRUE))/10^6
RC2TPM.37	<-	(sum(RC2TPM[,	'RPKT37'],	na.rm	=	TRUE))/10^6
RC2TPM.38	<-	(sum(RC2TPM[,	'RPKT38'],	na.rm	=	TRUE))/10^6
RC2TPM.39	<-	(sum(RC2TPM[,	'RPKT39'],	na.rm	=	TRUE))/10^6
RC2TPM.40	<-	(sum(RC2TPM[,	'RPKT40'],	na.rm	=	TRUE))/10^6
RC2TPM.41	<-	(sum(RC2TPM[,	'RPKT41'],	na.rm	=	TRUE))/10^6
RC2TPM.42	<-	(sum(RC2TPM[,	'RPKT42'],	na.rm	=	TRUE))/10^6
RC2TPM.43	<-	(sum(RC2TPM[,	'RPKT43'],	na.rm	=	TRUE))/10^6
RC2TPM.44	<-	(sum(RC2TPM[,	'RPKT44'],	na.rm	=	TRUE))/10^6
RC2TPM.45	<-	(sum(RC2TPM[,	'RPKT45'],	na.rm	=	TRUE))/10^6
RC2TPM.46	<-	(sum(RC2TPM[,	'RPKT46'],	na.rm	=	TRUE))/10^6
RC2TPM.47	<-	(sum(RC2TPM[,	'RPKT47'],	na.rm	=	TRUE))/10^6
RC2TPM.48	<-	(sum(RC2TPM[,	'RPKT48'],	na.rm	=	TRUE))/10^6
RC2TPM.49	<-	(sum(RC2TPM[,	'RPKT49'],	na.rm	=	TRUE))/10^6
RC2TPM.50	<-	(sum(RC2TPM[,	'RPKT50'],	na.rm	=	TRUE))/10^6
RC2TPM.51	<-	(sum(RC2TPM[,	'RPKT51'],	na.rm	=	TRUE))/10^6
RC2TPM.52	<-	(sum(RC2TPM[,	'RPKT52'],	na.rm	=	TRUE))/10^6
RC2TPM.53	<-	(sum(RC2TPM[,	'RPKT53'],	na.rm	=	TRUE))/10^6
RC2TPM.54	<-	(sum(RC2TPM[,	'RPKT54'],	na.rm	=	TRUE))/10^6
RC2TPM.55	<-	(sum(RC2TPM[,	'RPKT55'],	na.rm	=	TRUE))/10^6
RC2TPM.56	<-	(sum(RC2TPM[,	'RPKT56'],	na.rm	=	TRUE))/10^6
RC2TPM.57	<-	(sum(RC2TPM[,	'RPKT57'],	na.rm	=	TRUE))/10^6
RC2TPM.58	<-	(sum(RC2TPM[,	'RPKT58'],	na.rm	=	TRUE))/10^6

####################### Introducing TPMT Columns ###############################
RC2TPM$TPMC1	<-	RC2TPM$RPKC1	/	RC2TPMC.1
RC2TPM$TPMC2	<-	RC2TPM$RPKC2	/	RC2TPMC.2
RC2TPM$TPMC3	<-	RC2TPM$RPKC3	/	RC2TPMC.3
RC2TPM$TPMC4	<-	RC2TPM$RPKC4	/	RC2TPMC.4
RC2TPM$TPMC5	<-	RC2TPM$RPKC5	/	RC2TPMC.5
RC2TPM$TPMC6	<-	RC2TPM$RPKC6	/	RC2TPMC.6
RC2TPM$TPMC7	<-	RC2TPM$RPKC7	/	RC2TPMC.7
RC2TPM$TPMC8	<-	RC2TPM$RPKC8	/	RC2TPMC.8
RC2TPM$TPMT1	<-	RC2TPM$RPKT1	/	RC2TPM.1
RC2TPM$TPMT2	<-	RC2TPM$RPKT2	/	RC2TPM.2
RC2TPM$TPMT3	<-	RC2TPM$RPKT3	/	RC2TPM.3
RC2TPM$TPMT4	<-	RC2TPM$RPKT4	/	RC2TPM.4
RC2TPM$TPMT5	<-	RC2TPM$RPKT5	/	RC2TPM.5
RC2TPM$TPMT6	<-	RC2TPM$RPKT6	/	RC2TPM.6
RC2TPM$TPMT7	<-	RC2TPM$RPKT7	/	RC2TPM.7
RC2TPM$TPMT8	<-	RC2TPM$RPKT8	/	RC2TPM.8
RC2TPM$TPMT9	<-	RC2TPM$RPKT9	/	RC2TPM.9
RC2TPM$TPMT10	<-	RC2TPM$RPKT10	/	RC2TPM.10
RC2TPM$TPMT11	<-	RC2TPM$RPKT11	/	RC2TPM.11
RC2TPM$TPMT12	<-	RC2TPM$RPKT12	/	RC2TPM.12
RC2TPM$TPMT13	<-	RC2TPM$RPKT13	/	RC2TPM.13
RC2TPM$TPMT14	<-	RC2TPM$RPKT14	/	RC2TPM.14
RC2TPM$TPMT15	<-	RC2TPM$RPKT15	/	RC2TPM.15
RC2TPM$TPMT16	<-	RC2TPM$RPKT16	/	RC2TPM.16
RC2TPM$TPMT17	<-	RC2TPM$RPKT17	/	RC2TPM.17
RC2TPM$TPMT18	<-	RC2TPM$RPKT18	/	RC2TPM.18
RC2TPM$TPMT19	<-	RC2TPM$RPKT19	/	RC2TPM.19
RC2TPM$TPMT20	<-	RC2TPM$RPKT20	/	RC2TPM.20
RC2TPM$TPMT21	<-	RC2TPM$RPKT21	/	RC2TPM.21
RC2TPM$TPMT22	<-	RC2TPM$RPKT22	/	RC2TPM.22
RC2TPM$TPMT23	<-	RC2TPM$RPKT23	/	RC2TPM.23
RC2TPM$TPMT24	<-	RC2TPM$RPKT24	/	RC2TPM.24
RC2TPM$TPMT25	<-	RC2TPM$RPKT25	/	RC2TPM.25
RC2TPM$TPMT26	<-	RC2TPM$RPKT26	/	RC2TPM.26
RC2TPM$TPMT27	<-	RC2TPM$RPKT27	/	RC2TPM.27
RC2TPM$TPMT28	<-	RC2TPM$RPKT28	/	RC2TPM.28
RC2TPM$TPMT29	<-	RC2TPM$RPKT29	/	RC2TPM.29
RC2TPM$TPMT30	<-	RC2TPM$RPKT30	/	RC2TPM.30
RC2TPM$TPMT31	<-	RC2TPM$RPKT31	/	RC2TPM.31
RC2TPM$TPMT32	<-	RC2TPM$RPKT32	/	RC2TPM.32
RC2TPM$TPMT33	<-	RC2TPM$RPKT33	/	RC2TPM.33
RC2TPM$TPMT34	<-	RC2TPM$RPKT34	/	RC2TPM.34
RC2TPM$TPMT35	<-	RC2TPM$RPKT35	/	RC2TPM.35
RC2TPM$TPMT36	<-	RC2TPM$RPKT36	/	RC2TPM.36
RC2TPM$TPMT37	<-	RC2TPM$RPKT37	/	RC2TPM.37
RC2TPM$TPMT38	<-	RC2TPM$RPKT38	/	RC2TPM.38
RC2TPM$TPMT39	<-	RC2TPM$RPKT39	/	RC2TPM.39
RC2TPM$TPMT40	<-	RC2TPM$RPKT40	/	RC2TPM.40
RC2TPM$TPMT41	<-	RC2TPM$RPKT41	/	RC2TPM.41
RC2TPM$TPMT42	<-	RC2TPM$RPKT42	/	RC2TPM.42
RC2TPM$TPMT43	<-	RC2TPM$RPKT43	/	RC2TPM.43
RC2TPM$TPMT44	<-	RC2TPM$RPKT44	/	RC2TPM.44
RC2TPM$TPMT45	<-	RC2TPM$RPKT45	/	RC2TPM.45
RC2TPM$TPMT46	<-	RC2TPM$RPKT46	/	RC2TPM.46
RC2TPM$TPMT47	<-	RC2TPM$RPKT47	/	RC2TPM.47
RC2TPM$TPMT48	<-	RC2TPM$RPKT48	/	RC2TPM.48
RC2TPM$TPMT49	<-	RC2TPM$RPKT49	/	RC2TPM.49
RC2TPM$TPMT50	<-	RC2TPM$RPKT50	/	RC2TPM.50
RC2TPM$TPMT51	<-	RC2TPM$RPKT51	/	RC2TPM.51
RC2TPM$TPMT52	<-	RC2TPM$RPKT52	/	RC2TPM.52
RC2TPM$TPMT53	<-	RC2TPM$RPKT53	/	RC2TPM.53
RC2TPM$TPMT54	<-	RC2TPM$RPKT54	/	RC2TPM.54
RC2TPM$TPMT55	<-	RC2TPM$RPKT55	/	RC2TPM.55
RC2TPM$TPMT56	<-	RC2TPM$RPKT56	/	RC2TPM.56
RC2TPM$TPMT57	<-	RC2TPM$RPKT57	/	RC2TPM.57
RC2TPM$TPMT58	<-	RC2TPM$RPKT58	/	RC2TPM.58

########################## CONFIRM TPM NORMALIZATION ###########################
if(sum(RC2TPM$TPMC1)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC1)}
if(sum(RC2TPM$TPMC2)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC2)}
if(sum(RC2TPM$TPMC3)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC3)}
if(sum(RC2TPM$TPMC4)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC4)}
if(sum(RC2TPM$TPMC5)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC5)}
if(sum(RC2TPM$TPMC6)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC6)}
if(sum(RC2TPM$TPMC7)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC7)}
if(sum(RC2TPM$TPMC8)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMC8)}
if(sum(RC2TPM$TPMT1)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT1)}
if(sum(RC2TPM$TPMT2)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT2)}
if(sum(RC2TPM$TPMT3)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT3)}
if(sum(RC2TPM$TPMT4)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT4)}
if(sum(RC2TPM$TPMT5)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT5)}
if(sum(RC2TPM$TPMT6)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT6)}
if(sum(RC2TPM$TPMT7)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT7)}
if(sum(RC2TPM$TPMT8)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT8)}
if(sum(RC2TPM$TPMT9)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT9)}
if(sum(RC2TPM$TPMT10)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT10)}
if(sum(RC2TPM$TPMT11)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT11)}
if(sum(RC2TPM$TPMT12)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT12)}
if(sum(RC2TPM$TPMT13)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT13)}
if(sum(RC2TPM$TPMT14)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT14)}
if(sum(RC2TPM$TPMT15)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT15)}
if(sum(RC2TPM$TPMT16)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT16)}
if(sum(RC2TPM$TPMT17)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT17)}
if(sum(RC2TPM$TPMT18)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT18)}
if(sum(RC2TPM$TPMT19)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT19)}
if(sum(RC2TPM$TPMT20)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT20)}
if(sum(RC2TPM$TPMT21)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT21)}
if(sum(RC2TPM$TPMT22)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT22)}
if(sum(RC2TPM$TPMT23)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT23)}
if(sum(RC2TPM$TPMT24)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT24)}
if(sum(RC2TPM$TPMT25)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT25)}
if(sum(RC2TPM$TPMT26)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT26)}
if(sum(RC2TPM$TPMT27)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT27)}
if(sum(RC2TPM$TPMT28)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT28)}
if(sum(RC2TPM$TPMT29)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT29)}
if(sum(RC2TPM$TPMT30)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT30)}
if(sum(RC2TPM$TPMT31)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT31)}
if(sum(RC2TPM$TPMT32)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT32)}
if(sum(RC2TPM$TPMT33)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT33)}
if(sum(RC2TPM$TPMT34)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT34)}
if(sum(RC2TPM$TPMT35)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT35)}
if(sum(RC2TPM$TPMT36)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT36)}
if(sum(RC2TPM$TPMT37)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT37)}
if(sum(RC2TPM$TPMT38)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT38)}
if(sum(RC2TPM$TPMT39)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT39)}
if(sum(RC2TPM$TPMT40)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT40)}
if(sum(RC2TPM$TPMT41)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT41)}
if(sum(RC2TPM$TPMT42)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT42)}
if(sum(RC2TPM$TPMT43)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT43)}
if(sum(RC2TPM$TPMT44)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT44)}
if(sum(RC2TPM$TPMT45)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT45)}
if(sum(RC2TPM$TPMT46)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT46)}
if(sum(RC2TPM$TPMT47)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT47)}
if(sum(RC2TPM$TPMT48)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT48)}
if(sum(RC2TPM$TPMT49)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT49)}
if(sum(RC2TPM$TPMT50)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT50)}
if(sum(RC2TPM$TPMT51)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT51)}
if(sum(RC2TPM$TPMT52)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT52)}
if(sum(RC2TPM$TPMT53)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT53)}
if(sum(RC2TPM$TPMT54)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT54)}
if(sum(RC2TPM$TPMT55)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT55)}
if(sum(RC2TPM$TPMT56)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT56)}
if(sum(RC2TPM$TPMT57)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT57)}
if(sum(RC2TPM$TPMT58)	==	1e+06){print("TRUE")}else{sum(RC2TPM$TPMT58)}
"red"
write.csv(RC2TPM,"RESULTS_workflow\\1.Read_counts_to_TPM_SRA_TCGA_merge_DLBCL2_mRNA.csv", row.names = TRUE)

#Making_Entire gene list without filter_DESeq result + TPM Values Control and Tumor
x <- RC2TPM
dim(x)
#[1] 38482   199
x[1:5,1:3]

x.heat <- x[134:199]############################################################
dim(x.heat)
#[1] 38482    66
"red"
write.csv(x.heat,"RESULTS_workflow\\2.Only_TPM_SRA_TCGA_merge_DLBCL2_mRNA.csv", row.names = TRUE)
"brown"
x.heat<-read.table("RESULTS_workflow\\2.Only_TPM_SRA_TCGA_merge_DLBCL2_mRNA.csv", header = TRUE, sep = ",")
colnames(x.heat)[1]  <- "geneSymbol"    # change column name for x column

##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
DESeq2_value<-read.table("DESeq2_TvsN_SRA_TCGA_merge_DLBCL_mRNA.csv", header = TRUE, sep = ",")
##              ##              ##              ##              ##              ##              ##              ##              ##
colnames(DESeq2_value)[1]  <- "geneSymbol"    # change column name for x column
dim(DESeq2_value)
#[1] 38721         7
x.1 = merge(DESeq2_value, x.heat, ID="geneSymbol", sort = F)

# Average Control value as a new Variable (mean of Columns 8 to 12 is print in last column = 52)
x.1 <- transform(x.1, Mean_Contol_TPM = rowMeans(x.1[8:15], na.rm = TRUE))  #########################################################
# Average Control value as a new Variable (mean of Columns 13 to 41 is print in last column = 53)
x.1 <- transform(x.1, Mean_Tumour_TPM = rowMeans(x.1[16:73], na.rm = TRUE)) #########################################################

# Fold change (Tumour/average control)  (Fold change of Tumor from 34 to 72)
x.1 <- data.frame(x.1)                       # Creating example data
x.1$FC1  <- x.1$TPMT1  / x.1$Mean_Contol_TPM     # Creating new variable TPM1
x.1$FC2  <- x.1$TPMT2  / x.1$Mean_Contol_TPM     # Creating new variable TPM2
x.1$FC3  <- x.1$TPMT3  / x.1$Mean_Contol_TPM     # Creating new variable TPM3
x.1$FC4  <- x.1$TPMT4  / x.1$Mean_Contol_TPM     # Creating new variable TPM4
x.1$FC5  <- x.1$TPMT5  / x.1$Mean_Contol_TPM     # Creating new variable TPM5
x.1$FC6  <- x.1$TPMT6  / x.1$Mean_Contol_TPM     # Creating new variable TPM6
x.1$FC7  <- x.1$TPMT7  / x.1$Mean_Contol_TPM     # Creating new variable TPM7
x.1$FC8  <- x.1$TPMT8  / x.1$Mean_Contol_TPM     # Creating new variable TPM8
x.1$FC9  <- x.1$TPMT9  / x.1$Mean_Contol_TPM     # Creating new variable TPM9
x.1$FC10 <- x.1$TPMT10 / x.1$Mean_Contol_TPM  # Creating new variable TPM10
x.1$FC11 <- x.1$TPMT11 / x.1$Mean_Contol_TPM   # Creating new variable TPM11
x.1$FC12 <- x.1$TPMT12 / x.1$Mean_Contol_TPM   # Creating new variable TPM12
x.1$FC13 <- x.1$TPMT13 / x.1$Mean_Contol_TPM   # Creating new variable TPM13
x.1$FC14 <- x.1$TPMT14 / x.1$Mean_Contol_TPM   # Creating new variable TPM14
x.1$FC15 <- x.1$TPMT15 / x.1$Mean_Contol_TPM   # Creating new variable TPM15
x.1$FC16 <- x.1$TPMT16 / x.1$Mean_Contol_TPM   # Creating new variable TPM16
x.1$FC17 <- x.1$TPMT17 / x.1$Mean_Contol_TPM   # Creating new variable TPM17
x.1$FC18 <- x.1$TPMT18 / x.1$Mean_Contol_TPM   # Creating new variable TPM18
x.1$FC19 <- x.1$TPMT19 / x.1$Mean_Contol_TPM   # Creating new variable TPM19
x.1$FC20 <- x.1$TPMT20 / x.1$Mean_Contol_TPM  # Creating new variable TPM20
x.1$FC21  <- x.1$TPMT21  / x.1$Mean_Contol_TPM     # Creating new variable TPM21
x.1$FC22  <- x.1$TPMT22  / x.1$Mean_Contol_TPM     # Creating new variable TPM22
x.1$FC23  <- x.1$TPMT23  / x.1$Mean_Contol_TPM     # Creating new variable TPM23
x.1$FC24  <- x.1$TPMT24  / x.1$Mean_Contol_TPM     # Creating new variable TPM24
x.1$FC25  <- x.1$TPMT25  / x.1$Mean_Contol_TPM     # Creating new variable TPM25
x.1$FC26  <- x.1$TPMT26  / x.1$Mean_Contol_TPM     # Creating new variable TPM26
x.1$FC27  <- x.1$TPMT27  / x.1$Mean_Contol_TPM     # Creating new variable TPM27
x.1$FC28  <- x.1$TPMT28  / x.1$Mean_Contol_TPM     # Creating new variable TPM28
x.1$FC29  <- x.1$TPMT29  / x.1$Mean_Contol_TPM     # Creating new variable TPM29
x.1$FC30 <- x.1$TPMT30 / x.1$Mean_Contol_TPM  # Creating new variable TPM30
x.1$FC31  <- x.1$TPMT31  / x.1$Mean_Contol_TPM     # Creating new variable TPM31
x.1$FC32  <- x.1$TPMT32  / x.1$Mean_Contol_TPM     # Creating new variable TPM32
x.1$FC33  <- x.1$TPMT33  / x.1$Mean_Contol_TPM     # Creating new variable TPM33
x.1$FC34  <- x.1$TPMT34  / x.1$Mean_Contol_TPM     # Creating new variable TPM34
x.1$FC35  <- x.1$TPMT35  / x.1$Mean_Contol_TPM     # Creating new variable TPM35
x.1$FC36  <- x.1$TPMT36  / x.1$Mean_Contol_TPM     # Creating new variable TPM36
x.1$FC37  <- x.1$TPMT37  / x.1$Mean_Contol_TPM     # Creating new variable TPM37
x.1$FC38  <- x.1$TPMT38  / x.1$Mean_Contol_TPM     # Creating new variable TPM38
x.1$FC39  <- x.1$TPMT39  / x.1$Mean_Contol_TPM     # Creating new variable TPM39
x.1$FC40 <- x.1$TPMT40 / x.1$Mean_Contol_TPM  # Creating new variable TPM40
x.1$FC41  <- x.1$TPMT41  / x.1$Mean_Contol_TPM     # Creating new variable TPM41
x.1$FC42  <- x.1$TPMT42  / x.1$Mean_Contol_TPM     # Creating new variable TPM42
x.1$FC43  <- x.1$TPMT43  / x.1$Mean_Contol_TPM     # Creating new variable TPM43
x.1$FC44  <- x.1$TPMT44  / x.1$Mean_Contol_TPM     # Creating new variable TPM44
x.1$FC45  <- x.1$TPMT45  / x.1$Mean_Contol_TPM     # Creating new variable TPM45
x.1$FC46  <- x.1$TPMT46  / x.1$Mean_Contol_TPM     # Creating new variable TPM46
x.1$FC47  <- x.1$TPMT47  / x.1$Mean_Contol_TPM     # Creating new variable TPM47
x.1$FC48  <- x.1$TPMT48  / x.1$Mean_Contol_TPM     # Creating new variable TPM48
x.1$FC49  <- x.1$TPMT49  / x.1$Mean_Contol_TPM     # Creating new variable TPM49
x.1$FC50 <- x.1$TPMT50 / x.1$Mean_Contol_TPM  # Creating new variable TPM50
x.1$FC51  <- x.1$TPMT51  / x.1$Mean_Contol_TPM     # Creating new variable TPM51
x.1$FC52  <- x.1$TPMT52  / x.1$Mean_Contol_TPM     # Creating new variable TPM52
x.1$FC53  <- x.1$TPMT53  / x.1$Mean_Contol_TPM     # Creating new variable TPM53
x.1$FC54  <- x.1$TPMT54  / x.1$Mean_Contol_TPM     # Creating new variable TPM54
x.1$FC55  <- x.1$TPMT55  / x.1$Mean_Contol_TPM     # Creating new variable TPM55
x.1$FC56  <- x.1$TPMT56  / x.1$Mean_Contol_TPM     # Creating new variable TPM56
x.1$FC57  <- x.1$TPMT57  / x.1$Mean_Contol_TPM     # Creating new variable TPM57
x.1$FC58  <- x.1$TPMT58  / x.1$Mean_Contol_TPM     # Creating new variable TPM58

"red"
write.csv(x.1,"RESULTS_workflow\\3.L2FC_SRA_TCGA_merge_DLBCL2_mRNA.csv", row.names = TRUE)
"brown"
x.1<-read.table("RESULTS_workflow\\3.L2FC_SRA_TCGA_merge_DLBCL2_mRNA.csv", header = TRUE, sep = ",", row.names = 2)
x.1[1:5,1:6]
x.1 <- x.1[ -c(1) ] # Removes entire column be careful
x.1[1:5,1:6]
#log2 transform data # pull only FOldchange data
x.log <- log(x.1[75:132], 2)#################"x.log must have 16 Variables"####################################
x.log[x.log=="-Inf"]<-NA
colnames(x.log) <- c('L2FC1','L2FC2','L2FC3','L2FC4','L2FC5','L2FC6','L2FC7','L2FC8','L2FC9','L2FC10','L2FC11','L2FC12','L2FC13','L2FC14','L2FC15','L2FC16','L2FC17','L2FC18','L2FC19','L2FC20','L2FC21','L2FC22','L2FC23','L2FC24','L2FC25','L2FC26','L2FC27','L2FC28','L2FC29','L2FC30','L2FC31','L2FC32','L2FC33','L2FC34','L2FC35','L2FC36','L2FC37','L2FC38','L2FC39','L2FC40','L2FC41','L2FC42','L2FC43','L2FC44','L2FC45','L2FC46','L2FC47','L2FC48','L2FC49','L2FC50','L2FC51','L2FC52','L2FC53','L2FC54','L2FC55','L2FC56','L2FC57','L2FC58') #Rename the columns to L2FC
head(x.log)
"red"
write.csv(x.log,"RESULTS_workflow\\4.L2FC_only_SRA_TCGA_merge_DLBCL2_mRNA.csv", row.names = TRUE)
"brown"
z<-read.table("RESULTS_workflow\\4.L2FC_only_SRA_TCGA_merge_DLBCL2_mRNA.csv", header = TRUE, sep = ",")
z[1:5,1:6]
colnames(z)[1]  <- "geneSymbol"    # change column name for x column
z[1:5,1:6]
"brown"
y<-read.table("RESULTS_workflow\\3.L2FC_SRA_TCGA_merge_DLBCL2_mRNA.csv", header = TRUE, sep = ",")
y[1:5,1:6]
y <- y[-c(1)]
y[1:5,1:6]
#******************************************************************************#
All_factors_merged_norow = merge(y, z, ID="geneSymbol", sort = F)            #
All_factors_merged <- All_factors_merged_norow[,-1]                            #
rownames(All_factors_merged) <- All_factors_merged_norow[,1]                   #
dim(All_factors_merged)                                                        #
#[1] 38482   190                                                             #
#******************************************************************************#
"red"
write.csv(All_factors_merged,"RESULTS_workflow\\5.All_factors_merged_TPM_DLBCL2_mRNA.csv", row.names = TRUE)
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
############################################## SELECTING ONLY GENES WITH pvalue <= 0.05 ####################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
library(tidyverse)
pvalue_filter <- All_factors_merged %>% filter(pvalue <= 0.05) 
dim(pvalue_filter)
#[1] 14752    190
"red"
write.csv(pvalue_filter,"RESULTS_workflow\\6.pv_DESeq2_TPM_DLBCL2_T_mRNA.csv", row.names = TRUE)
#**********************************************************************************#
#column 1 to rownames
rownames(x.heat) <- x.heat[,1]
x.heat <- x.heat[,-1]
pvalue_0.05_filtered_TPM <- x.heat[rownames(x.heat) %in% rownames(pvalue_filter),]   #
dim(pvalue_0.05_filtered_TPM)                                                        #
#[1] 14752     66  *****************************************************************#
#**********************************************************************************#
"red"
write.csv(pvalue_0.05_filtered_TPM,"RESULTS_workflow\\7.pv_DESeq2_TPM_DLBCL2_mRNA_TPM1ly.csv", row.names = TRUE)
#*#****************************************************************************************************************###########################################################################################################################################
##########################################################################################################################################
########################### SELECTING ONLY GENES WITH pvalue <= 0.05 and L2FC > 2 & -2 EXPRESSION IN ALL SAMPLES ################################
##########################################################################################################################################
##########################################################################################################################################
pvalue_0.05_L2FC_2_filter_intermediate <- pvalue_filter[c(133:190)]
dim(pvalue_0.05_L2FC_2_filter_intermediate)
#[1] 14752     58
pvalue_0.05_L2FC_2_filtered_upr <- pvalue_0.05_L2FC_2_filter_intermediate[rowSums(pvalue_0.05_L2FC_2_filter_intermediate >= 2) == ncol(pvalue_0.05_L2FC_2_filter_intermediate), ]################ L2FC values >= 2 in all samples.#################
pvalue_0.05_L2FC_2_filtered_upr <- na.omit(pvalue_0.05_L2FC_2_filtered_upr)
pvalue_0.05_L2FC_2_filtered_upr <- data.frame(pvalue_0.05_L2FC_2_filtered_upr)
dim(pvalue_0.05_L2FC_2_filtered_upr)
#[1] 205 58
pvalue_0.05_L2FC_2_filtered_dwnr <- pvalue_0.05_L2FC_2_filter_intermediate[rowSums(pvalue_0.05_L2FC_2_filter_intermediate <= -2) == ncol(pvalue_0.05_L2FC_2_filter_intermediate), ]################ L2FC values >= 2 in all samples.#################
pvalue_0.05_L2FC_2_filtered_dwnr <- na.omit(pvalue_0.05_L2FC_2_filtered_dwnr)
pvalue_0.05_L2FC_2_filtered_dwnr <- data.frame(pvalue_0.05_L2FC_2_filtered_dwnr)
dim(pvalue_0.05_L2FC_2_filtered_dwnr)
#[1]  63 58

pvalue_0.05_L2FC_2_filtered <- rbind(pvalue_0.05_L2FC_2_filtered_upr, pvalue_0.05_L2FC_2_filtered_dwnr)
dim(pvalue_0.05_L2FC_2_filtered)
#[1] 268 58
"red"
write.csv(pvalue_0.05_L2FC_2_filtered,"RESULTS_workflow\\8.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_T_mRNA.csv", row.names = TRUE)
#**********************************************************************************************************#
upr_pvalue_0.05_L2FC_2_filtered_TPM<- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered_upr),]      #
dim(upr_pvalue_0.05_L2FC_2_filtered_TPM)                                                                         #
#[1] 205 66

dwnr_pvalue_0.05_L2FC_2_filtered_TPM<- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered_dwnr),]      #
dim(dwnr_pvalue_0.05_L2FC_2_filtered_TPM)                                                                         #
#[1]  63 66

pvalue_0.05_L2FC_2_filtered_TPM <- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered),]      #
dim(pvalue_0.05_L2FC_2_filtered_TPM)                                                                         #
#[1] 268 66  *********************************************************************************************#
#**********************************************************************************************************#
"red"
write.csv(pvalue_0.05_L2FC_2_filtered_TPM,"RESULTS_workflow\\8.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_T_mRNA_TPM1ly.csv", row.names = TRUE)
#**********************************************************************************************************#
pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged <- All_factors_merged[rownames(All_factors_merged) %in% rownames(pvalue_0.05_L2FC_2_filtered),]
dim(pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged)
#[1]  268 190
"red"
write.csv(pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged,"RESULTS_workflow\\8.pv_L2FC2_DESeq2_TPM_DLBCL2_T_mRNA.csv", row.names = TRUE)
#****************************************************************************************************************#
#*##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC > 2,  Mean_Tumour_TPM > 1 and EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
only_TTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered <- upr_pvalue_0.05_L2FC_2_filtered_TPM[c(9:66)]
dim(only_TTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered)
#[1] 205    8
#cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered <- cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate[rownames(cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate) %in% rownames(TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM),]
#dim(cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered)
#[1]18  18

only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered <- only_TTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered[rowSums(only_TTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered > 1) == ncol(only_TTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered), ] 
dim(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered)
# 55    58
only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM <- x.heat[rownames(x.heat) %in% rownames(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered),]      #
dim(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM)                                                                         #
#[1] 55  66 

##########################################################################################################################################
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC > 2, and Mean_Normal_TPM < 1 EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered <- upr_pvalue_0.05_L2FC_2_filtered_TPM[c(1:8)]
dim(cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered)
#[1] 205    8
cTPM_pvalue_0.05_l2fc2_TPM_1_filtered <- cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered[rowSums(cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered < 1) == ncol(cTPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate_filtered), ] 
dim(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered)
# 192    8
cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM <- x.heat[rownames(x.heat) %in% rownames(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered),]      #
dim(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM)                                                                         #
#[1] 192  66 

##########################################################################################################################################
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC < -2 and Mean_Normal_TPM > 1 EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate <- dwnr_pvalue_0.05_L2FC_2_filtered_TPM[c(1:8)]
dim(TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate)
#[1] 63    8
#TPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered <- TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate[rownames(TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate) %in% rownames(dwnr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged),]
#dim(TPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered)
#[1]34  16

TPM_pvalue_0.05_l2fc2_cTPM_1_filtered <- TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate[rowSums(TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate > 1) == ncol(TPM_pvalue_mRNA0.055_L2FC_2_filter_intermediate), ] 
dim(TPM_pvalue_0.05_l2fc2_cTPM_1_filtered)
# 49 8
TPM_pvalue_0.05_l2fc2_cTPM_1_filtered_TPM <- x.heat[rownames(x.heat) %in% rownames(TPM_pvalue_0.05_l2fc2_cTPM_1_filtered),]      #
dim(TPM_pvalue_0.05_l2fc2_cTPM_1_filtered_TPM)                                                                         #
#[1] 49 66  

Merge_pv_l2fc2_cTPM_upr_dwnr <- rbind(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM, TPM_pvalue_0.05_l2fc2_cTPM_1_filtered_TPM)
dim(Merge_pv_l2fc2_cTPM_upr_dwnr)

#impute missing data
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("impute")
library(impute) #bioconductor package  - https://www.bioconductor.org/packages/release/bioc/html/impute.html

#1_mRNA_pvalue
mRNA_highexp_pv.heatmap.imp<-impute.knn(as.matrix(pvalue_0.05_filtered_TPM))$data 
write.csv(mRNA_highexp_pv.heatmap.imp,"RESULTS_workflow\\10.forheatmap_pv_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)
#2_mRNA_pvalue_l2fc_upr
upr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(upr_pvalue_0.05_L2FC_2_filtered_TPM))$data
write.csv(upr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow\\10.forheatmap_pv_upr_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)
#3_mRNA_pvalue_l2fc_dwnr
dwnr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(dwnr_pvalue_0.05_L2FC_2_filtered_TPM))$data 
write.csv(dwnr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow\\10.forheatmap_pv_dwnr_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)
#4_mRNA_pvalue_l2fc_upr_dwnr
upr_dwnr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(pvalue_0.05_L2FC_2_filtered_TPM))$data 
write.csv(upr_dwnr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow\\10.forheatmap_pv_upr_dwnr_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)
#5  mRNA_pvalue_l2fc_upr_cTPM<1
cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM.imp <- impute.knn(as.matrix(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM))$data
#6 mRNA_pvalue_l2fc_dwnr_cTPM>1
pvalue_0.05_L2FC_2_l2_cTPM_filtered.imp <- impute.knn(as.matrix(TPM_pvalue_0.05_l2fc2_cTPM_1_filtered_TPM))$data
#7 mRNA_pvalue_l2fc_upr_dwnr_cTPM
pvalue_0.05_L2FC2_cTPM_upr_dwnr.imp <- impute.knn(as.matrix(Merge_pv_l2fc2_cTPM_upr_dwnr))$data

##################################################
#Plot Up-regulated heatmap with column clustering#
source("arrayTools-noClanc.R")                   #
##################################################
#1_mRNA_pvalue
gc()
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\1_pvalue_0.05_mRNA_DLBCL2_heat.pdf")
myHeatmap((mRNA_highexp_pv.heatmap.imp)) ############ plot Heatmap ##################
dev.off()
##################################################
#2_mRNA_pvalue_l2fc_upr
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\2_pvalue_upr_L2FC2_mRNA_DLBCL2_heat.pdf")
myHeatmap((upr_highexp_pv_L2FC2.heatmap.imp))   ############ plot Heatmap ##################
dev.off()
##################################################
#3_mRNA_pvalue_l2fc_dwnr
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\3_pvalue_dwnr_L2FC2_mRNA_DLBCL2_heat.pdf")
myHeatmap((dwnr_highexp_pv_L2FC2.heatmap.imp)) ############ plot Heatmap ##################
dev.off()
##################################################
#4_mRNA_pvalue_l2fc_upr_dwnr
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\4_pvalue_upr_dwnr_L2FC2_mRNA_DLBCL2_heat.pdf")
myHeatmap((upr_dwnr_highexp_pv_L2FC2.heatmap.imp))   ############ plot Heatmap ##################
dev.off()
##################################################
#5_mRNA_pvalue_l2fc_upr_cTPM
pdf(file="5_pv_upr_L2FC2_cTPM_mRNA_DL.pdf")
myHeatmap((cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM.imp))   ############ plot Heatmap ##################
dev.off()
##################################################
#6_mRNA_pvalue_l2fc_dwnr_cTPM
pdf(file="5_pv_dwnr_L2FC2_cTPM_mRNA_DL.pdf")
myHeatmap((pvalue_0.05_L2FC_2_l2_cTPM_filtered.imp))   ############ plot Heatmap ##################
dev.off()
##################################################
#7 mRNA_pvalue_l2fc_upr_dwnr_cTPM
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\5_pv_upr_dwnr_L2FC2_cTPM_mRNA_DLBCL2_heat.pdf")
myHeatmap((pvalue_0.05_L2FC2_cTPM_upr_dwnr.imp))   ############ plot Heatmap ##################
dev.off()
#####################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
############################################## Gene Set Analysis_Pathway ####################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
### LBL-10725-01_Metabolic_Pathways_Panel_Gene_List was downloaded and  the pathway Annotation sheet was saved as a new .csv file
# Load required libraries
library(dplyr)
library(readr)

#oncogenes <- scan("oncokb_oncogenes.txt", what="", sep="\n")

# Step 1: Load the CSV file
#c("green","darkgreen")
#file_path <- "pathway.csv"  # Update with the correct path to your CSV file
#MetOncoTSRepair <- read_csv(file_path)
#MetOncoTSRepair <- read.table("pathway.csv", header = TRUE, sep = ",")

# Step 2: Convert "+" to respective column names
#MetOncoTSRepair <- MetOncoTSRepair %>%
#  rowwise() %>%
#  mutate(Pathways = paste(names(.)[which(c_across(-Gene) == "+")], collapse = ", ")) %>%
#  ungroup()

# Step 3: Keep only 'Gene' and 'Pathways' columns
#MetOncoTSRepair <- MetOncoTSRepair %>%
#  select(Gene, Pathways)

# Step 4: Save the cleaned DataFrame to a CSV file
#output_path <- "pathway_cleaned.csv"  # Update with the desired output path
#write_csv(MetOncoTSRepair, output_path)
#print("File saved successfully!")

#Step 4.1: The first column of "df" had values other than "+/-" Hence that has to be inserted manually in "pathway_cleaned_blanks_filled.csv" 

# Step 5: convert df name
"brown"
#MetOncoTSRepair<-read.table("RESULTS_workflow\\pathway_cleaned_blanks_filled.csv",sep=",",header=T,row.names=1)
library(tibble)

##########################################################################################################################################
"green"
MetabolicGenes<-read.table("1_pathway_Metabolicgenes.csv",sep=",",header=T,row.names=1)
MetabolicGenes[MetabolicGenes == ""] <- NA
MetabolicGenes <- tibble::rownames_to_column(MetabolicGenes, var = "geneSymbol")

"green"
ONCO_TS_Genes<-read.table("2_pathway_Onco_N_TS.csv",sep=",",header=T,row.names=1)
ONCO_TS_Genes <- tibble::rownames_to_column(ONCO_TS_Genes, var = "geneSymbol")
"green"
DNArepairGenes<-read.table("3_pathway_DNArepair.csv",sep=",",header=T,row.names=1)
DNArepairGenes <- tibble::rownames_to_column(DNArepairGenes, var = "geneSymbol")

# Step 1: Perform full joins step-by-step
metoncoTS <- full_join(MetabolicGenes, ONCO_TS_Genes, by = "geneSymbol")
MetOncoTSRepair <- full_join(metoncoTS, DNArepairGenes, by = "geneSymbol")

#rename columns
colnames(MetOncoTSRepair)[2:12] <- paste0("Metabolicpathway", 1:11)
colnames(MetOncoTSRepair)[13] <- paste0("Oncogene or Tumour Supressor")
colnames(MetOncoTSRepair)[14] <- paste0("DNARepair")
"red"
write.csv(MetOncoTSRepair,"RESULTS_workflow\\6_MetOncoTSRepair_pathway_genelist.csv", row.names = TRUE)

# column to rownames
MetOncoTSRepair <- tibble::column_to_rownames(MetOncoTSRepair, var = "geneSymbol")

##########################################################################################################################################
#pick only the genes in  MetOncoTSRepair from all filters included file/df (All_factors_merged)
filtered_df <- All_factors_merged[rownames(All_factors_merged) %in% rownames(MetOncoTSRepair),]   #
dim(filtered_df)                                                        #
# [1] 1933  190
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
############################################## SELECTING ONLY GENES WITH pvalue <= 0.05 ####################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
pvalue_filter_MOD <- filtered_df %>% filter(pvalue <= 0.05) 
dim(pvalue_filter_MOD)
#[1] 980 190
"red"
write.csv(pvalue_filter_MOD,"RESULTS_workflow\\7.pv_DESeq2_TPM_DLBCL2_T_mRNA_MOD.csv", row.names = TRUE)
#**********************************************************************************#
#column 1 to rownames

pvalue_0.05_filtered_TPM_MOD <- x.heat[rownames(x.heat) %in% rownames(pvalue_filter_MOD),]   #
dim(pvalue_0.05_filtered_TPM_MOD)                                                        #
#[1] 980  66  *****************************************************************#
#**********************************************************************************#
"red"
write.csv(pvalue_0.05_filtered_TPM_MOD,"RESULTS_workflow\\7.pv_DESeq2_TPM_DLBCL2_mRNA_MOD_TPM1ly.csv", row.names = TRUE)
#*#****************************************************************************************************************###########################################################################################################################################
##########################################################################################################################################
########################### SELECTING ONLY GENES WITH pvalue <= 0.05 and L2FC > 2 & -2 EXPRESSION IN ALL SAMPLES ################################
##########################################################################################################################################
##########################################################################################################################################
pvalue_0.05_L2FC_2_filter_intermediate_MOD <- pvalue_filter_MOD[c(133:190)]
dim(pvalue_0.05_L2FC_2_filter_intermediate_MOD)
#[1] 980  58
pvalue_0.05_L2FC_2_filtered_upr_MOD <- pvalue_0.05_L2FC_2_filter_intermediate_MOD[rowSums(pvalue_0.05_L2FC_2_filter_intermediate_MOD >= 2) == ncol(pvalue_0.05_L2FC_2_filter_intermediate_MOD), ]################ L2FC values >= 2 in all samples.#################
pvalue_0.05_L2FC_2_filtered_upr_MOD <- na.omit(pvalue_0.05_L2FC_2_filtered_upr_MOD)
pvalue_0.05_L2FC_2_filtered_upr_MOD <- data.frame(pvalue_0.05_L2FC_2_filtered_upr_MOD)
dim(pvalue_0.05_L2FC_2_filtered_upr_MOD)
#[1] 23 58
pvalue_0.05_L2FC_2_filtered_dwnr_MOD <- pvalue_0.05_L2FC_2_filter_intermediate_MOD[rowSums(pvalue_0.05_L2FC_2_filter_intermediate_MOD <= -2) == ncol(pvalue_0.05_L2FC_2_filter_intermediate_MOD), ]################ L2FC values >= 2 in all samples.#################
pvalue_0.05_L2FC_2_filtered_dwnr_MOD <- na.omit(pvalue_0.05_L2FC_2_filtered_dwnr_MOD)
pvalue_0.05_L2FC_2_filtered_dwnr_MOD <- data.frame(pvalue_0.05_L2FC_2_filtered_dwnr_MOD)
dim(pvalue_0.05_L2FC_2_filtered_dwnr_MOD)
#[1]  2 58

pvalue_0.05_L2FC_2_filtered_MOD <- rbind(pvalue_0.05_L2FC_2_filtered_upr_MOD, pvalue_0.05_L2FC_2_filtered_dwnr_MOD)
dim(pvalue_0.05_L2FC_2_filtered_MOD)
#[1] 25 58
"red"
write.csv(pvalue_0.05_L2FC_2_filtered_MOD,"RESULTS_workflow\\8.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_T_mRNA_MOD.csv", row.names = TRUE)
#**********************************************************************************************************#
upr_pvalue_0.05_L2FC_2_filtered_TPM_MOD<- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered_upr_MOD),]      #
dim(upr_pvalue_0.05_L2FC_2_filtered_TPM_MOD)                                                                         #
#[1] 23 66
dwnr_pvalue_0.05_L2FC_2_filtered_TPM_MOD<- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered_dwnr_MOD),]      #
dim(dwnr_pvalue_0.05_L2FC_2_filtered_TPM_MOD)                                                                         #
#[1]  2 66

pvalue_0.05_L2FC_2_filtered_TPM_MOD <- x.heat[rownames(x.heat) %in% rownames(pvalue_0.05_L2FC_2_filtered_MOD),]      #
dim(pvalue_0.05_L2FC_2_filtered_TPM_MOD)                                                                         #
#[1] 25 66  *********************************************************************************************#
#**********************************************************************************************************#
"red"
write.csv(pvalue_0.05_L2FC_2_filtered_TPM_MOD,"RESULTS_workflow\\8.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_T_mRNA_MOD_TPM1ly.csv", row.names = TRUE)
#**********************************************************************************************************#
pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged_MOD <- All_factors_merged[rownames(All_factors_merged) %in% rownames(pvalue_0.05_L2FC_2_filtered_MOD),]
dim(pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged_MOD)
#[1]  25 190
"red"
write.csv(pvalue_0.05_L2FC_2_filtered_TPM_all_factors_merged_MOD,"RESULTS_workflow\\8.pv_L2FC2_DESeq2_TPM_DLBCL2_T_mRNA_MOD.csv", row.names = TRUE)
#****************************************************************************************************************#
# pick gene annotations from filtered_df

pvalue_0.05_l2fc_filtered_list_MOD <- MetOncoTSRepair[rownames(MetOncoTSRepair) %in% rownames(pvalue_0.05_L2FC_2_filtered_MOD),]   #
pvalue_0.05_l2fc_filtered_list_MOD <- data.frame(pvalue_0.05_l2fc_filtered_list_MOD)
dim(pvalue_0.05_l2fc_filtered_list_MOD)
#[1] 25 13
"red"
write.csv(pvalue_0.05_l2fc_filtered_list_MOD,"RESULTS_workflow\\8.filteredgene__TPM_DLBCL2_T_mRNA_MOD.csv", row.names = TRUE)
# use venny online tool for "8.filteredgene__TPM_DLBCL2_T_mRNA.csv" and list out the pathways without duplicates
# and make string_list and paste it as below.
#****************************************************************************************************************#
string_list <- c("AMPK", "Glutamine Metabolism", "Tryptophan/Kynurenine Metabolism", 
                 "DNA Damage Repair", "Glucose Transport", "Autophagy", "Glycolysis", 
                 "Pentose Phosphate Pathway", "Amino Acid Transporters", 
                 "Transcriptional Regulation", "Epigenetic Regulation", "Cell Type", 
                 "mTOR", "Cytokine & Chemokine Signaling", "Internal Reference", 
                 "Antigen Presentation", "PI3K", "MAPK", "Fatty Acid Synthesis", 
                 "Amino Acid Synthesis", "NF-KB", "Cell Cycle", "TCR & Costimulatory Signaling", 
                 "Mitochondrial Respiration", "p53 Pathway", "TLR Signaling", 
                 "Nucleotide Synthesis", "Nucleotide Salvage", "Oncogene", "Tumour_Suppressor",
                 "Base excision repair (BER)", "Other BER and strand break joining factors",
                 "Poly(ADP-ribose) polymerase (PARP) enzymes that bind to DNA", "Direct reversal of damage",
                 "Repair of DNA-protein crosslinks", "Mismatch excision repair (MMR)", "Nucleotide excision repair (NER)",
                 "NER-related", "Homologous recombination", "Fanconi anemia", "Non-homologous end-joining (NHEJ)",
                 "Modulation of nucleotide pools", "DNA polymerases", "Editing and processing nucleases",
                 "Ubiquitination and modification", "Chromatin Structure and Modification",
                 "Genes defective in diseases associated with sensitivity to DNA damaging agents",
                 "Other identified genes with known or suspected DNA repair function", "Other conserved DNA damage response genes")

##########################################################################################################################################
##########################################################################################################################################
########################### lncRNA in CHOL_SRA_TCGA_merge ################################
##########################################################################################################################################
##########################################################################################################################################

#Input_Entire gene list
##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
countData_lncRNA <- read.table("fordeseq_SRA_TCGA_merge_hisat2_N_T_DLBCL_lncRNA_RC.tsv",sep="\t",header=T,row.names=1)
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(countData_lncRNA)
#[1] 127802     66
countData_lncRNA[1:5,1:3]

#Input_gene length # Make sure of lncRNA and gene lengths file
##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
length_lncRNA<-read.table("lncRNA_read_length.tsv",sep="\t",header=T,row.names=1) 
##              ##              ##              ##              ##              ##              ##              ##              ##
dim(length_lncRNA)
# [1] 127802      2

RC2TPM_lncRNA <- data.frame(countData_lncRNA) # Creating data
######################################################################################################################################
############################################# [READCOUNTS TO TPM CONVERSION] #########################################################
######### Introducing gene_length_in_Kb Column by performing (END - START values by 1000) to get the gene length in Kb ###############
######################################################################################################################################
######################################################################################################################################
RC2TPM_lncRNA[,67] = (length_lncRNA[,2]-length_lncRNA[,1])/1000
colnames(RC2TPM_lncRNA)[67] <- "gene_length_in_Kb"
RC2TPM_lncRNA[RC2TPM_lncRNA == 0] <- 0.1
#cHANGE: 
# 1. RC2TPM_lncRNA$C TO RC2TPM_lncRNA$CC
# 2. RC2TPM_lncRNA$T TO RC2TPM_lncRNA$TC
########################### Introducing RPKT Columns ###########################
RC2TPM_lncRNA$RPKC1	<-	RC2TPM_lncRNA$C1	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC1
RC2TPM_lncRNA$RPKC2	<-	RC2TPM_lncRNA$C2	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC2
RC2TPM_lncRNA$RPKC3	<-	RC2TPM_lncRNA$C3	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC3
RC2TPM_lncRNA$RPKC4	<-	RC2TPM_lncRNA$C4	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC4
RC2TPM_lncRNA$RPKC5	<-	RC2TPM_lncRNA$C5	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC5
RC2TPM_lncRNA$RPKC6	<-	RC2TPM_lncRNA$C6	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC6
RC2TPM_lncRNA$RPKC7	<-	RC2TPM_lncRNA$C7	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC7
RC2TPM_lncRNA$RPKC8	<-	RC2TPM_lncRNA$C8	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKC8
RC2TPM_lncRNA$RPKT1	<-	RC2TPM_lncRNA$T1	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT1
RC2TPM_lncRNA$RPKT2	<-	RC2TPM_lncRNA$T2	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT2
RC2TPM_lncRNA$RPKT3	<-	RC2TPM_lncRNA$T3	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT3
RC2TPM_lncRNA$RPKT4	<-	RC2TPM_lncRNA$T4	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT4
RC2TPM_lncRNA$RPKT5	<-	RC2TPM_lncRNA$T5	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT5
RC2TPM_lncRNA$RPKT6	<-	RC2TPM_lncRNA$T6	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT6
RC2TPM_lncRNA$RPKT7	<-	RC2TPM_lncRNA$T7	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT7
RC2TPM_lncRNA$RPKT8	<-	RC2TPM_lncRNA$T8	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT8
RC2TPM_lncRNA$RPKT9	<-	RC2TPM_lncRNA$T9	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT9
RC2TPM_lncRNA$RPKT10	<-	RC2TPM_lncRNA$T10	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT10
RC2TPM_lncRNA$RPKT11	<-	RC2TPM_lncRNA$T11	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT11
RC2TPM_lncRNA$RPKT12	<-	RC2TPM_lncRNA$T12	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT12
RC2TPM_lncRNA$RPKT13	<-	RC2TPM_lncRNA$T13	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT13
RC2TPM_lncRNA$RPKT14	<-	RC2TPM_lncRNA$T14	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT14
RC2TPM_lncRNA$RPKT15	<-	RC2TPM_lncRNA$T15	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT15
RC2TPM_lncRNA$RPKT16	<-	RC2TPM_lncRNA$T16	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT16
RC2TPM_lncRNA$RPKT17	<-	RC2TPM_lncRNA$T17	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT17
RC2TPM_lncRNA$RPKT18	<-	RC2TPM_lncRNA$T18	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT18
RC2TPM_lncRNA$RPKT19	<-	RC2TPM_lncRNA$T19	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT19
RC2TPM_lncRNA$RPKT20	<-	RC2TPM_lncRNA$T20	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT20
RC2TPM_lncRNA$RPKT21	<-	RC2TPM_lncRNA$T21	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT21
RC2TPM_lncRNA$RPKT22	<-	RC2TPM_lncRNA$T22	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT22
RC2TPM_lncRNA$RPKT23	<-	RC2TPM_lncRNA$T23	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT23
RC2TPM_lncRNA$RPKT24	<-	RC2TPM_lncRNA$T24	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT24
RC2TPM_lncRNA$RPKT25	<-	RC2TPM_lncRNA$T25	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT25
RC2TPM_lncRNA$RPKT26	<-	RC2TPM_lncRNA$T26	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT26
RC2TPM_lncRNA$RPKT27	<-	RC2TPM_lncRNA$T27	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT27
RC2TPM_lncRNA$RPKT28	<-	RC2TPM_lncRNA$T28	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT28
RC2TPM_lncRNA$RPKT29	<-	RC2TPM_lncRNA$T29	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT29
RC2TPM_lncRNA$RPKT30	<-	RC2TPM_lncRNA$T30	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT30
RC2TPM_lncRNA$RPKT31	<-	RC2TPM_lncRNA$T31	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT31
RC2TPM_lncRNA$RPKT32	<-	RC2TPM_lncRNA$T32	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT32
RC2TPM_lncRNA$RPKT33	<-	RC2TPM_lncRNA$T33	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT33
RC2TPM_lncRNA$RPKT34	<-	RC2TPM_lncRNA$T34	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT34
RC2TPM_lncRNA$RPKT35	<-	RC2TPM_lncRNA$T35	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT35
RC2TPM_lncRNA$RPKT36	<-	RC2TPM_lncRNA$T36	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT36
RC2TPM_lncRNA$RPKT37	<-	RC2TPM_lncRNA$T37	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT37
RC2TPM_lncRNA$RPKT38	<-	RC2TPM_lncRNA$T38	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT38
RC2TPM_lncRNA$RPKT39	<-	RC2TPM_lncRNA$T39	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT39
RC2TPM_lncRNA$RPKT40	<-	RC2TPM_lncRNA$T40	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT40
RC2TPM_lncRNA$RPKT41	<-	RC2TPM_lncRNA$T41	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT41
RC2TPM_lncRNA$RPKT42	<-	RC2TPM_lncRNA$T42	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT42
RC2TPM_lncRNA$RPKT43	<-	RC2TPM_lncRNA$T43	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT43
RC2TPM_lncRNA$RPKT44	<-	RC2TPM_lncRNA$T44	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT44
RC2TPM_lncRNA$RPKT45	<-	RC2TPM_lncRNA$T45	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT45
RC2TPM_lncRNA$RPKT46	<-	RC2TPM_lncRNA$T46	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT46
RC2TPM_lncRNA$RPKT47	<-	RC2TPM_lncRNA$T47	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT47
RC2TPM_lncRNA$RPKT48	<-	RC2TPM_lncRNA$T48	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT48
RC2TPM_lncRNA$RPKT49	<-	RC2TPM_lncRNA$T49	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT49
RC2TPM_lncRNA$RPKT50	<-	RC2TPM_lncRNA$T50	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT50
RC2TPM_lncRNA$RPKT51	<-	RC2TPM_lncRNA$T51	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT51
RC2TPM_lncRNA$RPKT52	<-	RC2TPM_lncRNA$T52	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT52
RC2TPM_lncRNA$RPKT53	<-	RC2TPM_lncRNA$T53	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT53
RC2TPM_lncRNA$RPKT54	<-	RC2TPM_lncRNA$T54	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT54
RC2TPM_lncRNA$RPKT55	<-	RC2TPM_lncRNA$T55	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT55
RC2TPM_lncRNA$RPKT56	<-	RC2TPM_lncRNA$T56	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT56
RC2TPM_lncRNA$RPKT57	<-	RC2TPM_lncRNA$T57	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT57
RC2TPM_lncRNA$RPKT58	<-	RC2TPM_lncRNA$T58	/	RC2TPM_lncRNA$gene_length_in_Kb	#	Creating	new	variable	RPKT58

######################### Introducing PMRFT Columns ############################
RC2TPM_lncRNAC.1	<-	(sum(RC2TPM_lncRNA[,	'RPKC1'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.2	<-	(sum(RC2TPM_lncRNA[,	'RPKC2'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.3	<-	(sum(RC2TPM_lncRNA[,	'RPKC3'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.4	<-	(sum(RC2TPM_lncRNA[,	'RPKC4'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.5	<-	(sum(RC2TPM_lncRNA[,	'RPKC5'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.6	<-	(sum(RC2TPM_lncRNA[,	'RPKC6'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.7	<-	(sum(RC2TPM_lncRNA[,	'RPKC7'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNAC.8	<-	(sum(RC2TPM_lncRNA[,	'RPKC8'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.1	<-	(sum(RC2TPM_lncRNA[,	'RPKT1'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.2	<-	(sum(RC2TPM_lncRNA[,	'RPKT2'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.3	<-	(sum(RC2TPM_lncRNA[,	'RPKT3'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.4	<-	(sum(RC2TPM_lncRNA[,	'RPKT4'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.5	<-	(sum(RC2TPM_lncRNA[,	'RPKT5'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.6	<-	(sum(RC2TPM_lncRNA[,	'RPKT6'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.7	<-	(sum(RC2TPM_lncRNA[,	'RPKT7'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.8	<-	(sum(RC2TPM_lncRNA[,	'RPKT8'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.9	<-	(sum(RC2TPM_lncRNA[,	'RPKT9'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.10	<-	(sum(RC2TPM_lncRNA[,	'RPKT10'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.11	<-	(sum(RC2TPM_lncRNA[,	'RPKT11'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.12	<-	(sum(RC2TPM_lncRNA[,	'RPKT12'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.13	<-	(sum(RC2TPM_lncRNA[,	'RPKT13'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.14	<-	(sum(RC2TPM_lncRNA[,	'RPKT14'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.15	<-	(sum(RC2TPM_lncRNA[,	'RPKT15'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.16	<-	(sum(RC2TPM_lncRNA[,	'RPKT16'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.17	<-	(sum(RC2TPM_lncRNA[,	'RPKT17'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.18	<-	(sum(RC2TPM_lncRNA[,	'RPKT18'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.19	<-	(sum(RC2TPM_lncRNA[,	'RPKT19'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.20	<-	(sum(RC2TPM_lncRNA[,	'RPKT20'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.21	<-	(sum(RC2TPM_lncRNA[,	'RPKT21'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.22	<-	(sum(RC2TPM_lncRNA[,	'RPKT22'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.23	<-	(sum(RC2TPM_lncRNA[,	'RPKT23'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.24	<-	(sum(RC2TPM_lncRNA[,	'RPKT24'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.25	<-	(sum(RC2TPM_lncRNA[,	'RPKT25'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.26	<-	(sum(RC2TPM_lncRNA[,	'RPKT26'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.27	<-	(sum(RC2TPM_lncRNA[,	'RPKT27'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.28	<-	(sum(RC2TPM_lncRNA[,	'RPKT28'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.29	<-	(sum(RC2TPM_lncRNA[,	'RPKT29'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.30	<-	(sum(RC2TPM_lncRNA[,	'RPKT30'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.31	<-	(sum(RC2TPM_lncRNA[,	'RPKT31'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.32	<-	(sum(RC2TPM_lncRNA[,	'RPKT32'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.33	<-	(sum(RC2TPM_lncRNA[,	'RPKT33'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.34	<-	(sum(RC2TPM_lncRNA[,	'RPKT34'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.35	<-	(sum(RC2TPM_lncRNA[,	'RPKT35'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.36	<-	(sum(RC2TPM_lncRNA[,	'RPKT36'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.37	<-	(sum(RC2TPM_lncRNA[,	'RPKT37'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.38	<-	(sum(RC2TPM_lncRNA[,	'RPKT38'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.39	<-	(sum(RC2TPM_lncRNA[,	'RPKT39'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.40	<-	(sum(RC2TPM_lncRNA[,	'RPKT40'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.41	<-	(sum(RC2TPM_lncRNA[,	'RPKT41'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.42	<-	(sum(RC2TPM_lncRNA[,	'RPKT42'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.43	<-	(sum(RC2TPM_lncRNA[,	'RPKT43'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.44	<-	(sum(RC2TPM_lncRNA[,	'RPKT44'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.45	<-	(sum(RC2TPM_lncRNA[,	'RPKT45'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.46	<-	(sum(RC2TPM_lncRNA[,	'RPKT46'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.47	<-	(sum(RC2TPM_lncRNA[,	'RPKT47'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.48	<-	(sum(RC2TPM_lncRNA[,	'RPKT48'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.49	<-	(sum(RC2TPM_lncRNA[,	'RPKT49'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.50	<-	(sum(RC2TPM_lncRNA[,	'RPKT50'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.51	<-	(sum(RC2TPM_lncRNA[,	'RPKT51'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.52	<-	(sum(RC2TPM_lncRNA[,	'RPKT52'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.53	<-	(sum(RC2TPM_lncRNA[,	'RPKT53'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.54	<-	(sum(RC2TPM_lncRNA[,	'RPKT54'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.55	<-	(sum(RC2TPM_lncRNA[,	'RPKT55'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.56	<-	(sum(RC2TPM_lncRNA[,	'RPKT56'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.57	<-	(sum(RC2TPM_lncRNA[,	'RPKT57'],	na.rm	=	TRUE))/10^6
RC2TPM_lncRNA.58	<-	(sum(RC2TPM_lncRNA[,	'RPKT58'],	na.rm	=	TRUE))/10^6

####################### Introducing TPMT Columns ###############################
RC2TPM_lncRNA$TPMC1	<-	RC2TPM_lncRNA$RPKC1	/	RC2TPM_lncRNAC.1
RC2TPM_lncRNA$TPMC2	<-	RC2TPM_lncRNA$RPKC2	/	RC2TPM_lncRNAC.2
RC2TPM_lncRNA$TPMC3	<-	RC2TPM_lncRNA$RPKC3	/	RC2TPM_lncRNAC.3
RC2TPM_lncRNA$TPMC4	<-	RC2TPM_lncRNA$RPKC4	/	RC2TPM_lncRNAC.4
RC2TPM_lncRNA$TPMC5	<-	RC2TPM_lncRNA$RPKC5	/	RC2TPM_lncRNAC.5
RC2TPM_lncRNA$TPMC6	<-	RC2TPM_lncRNA$RPKC6	/	RC2TPM_lncRNAC.6
RC2TPM_lncRNA$TPMC7	<-	RC2TPM_lncRNA$RPKC7	/	RC2TPM_lncRNAC.7
RC2TPM_lncRNA$TPMC8	<-	RC2TPM_lncRNA$RPKC8	/	RC2TPM_lncRNAC.8
RC2TPM_lncRNA$TPMT1	<-	RC2TPM_lncRNA$RPKT1	/	RC2TPM_lncRNA.1
RC2TPM_lncRNA$TPMT2	<-	RC2TPM_lncRNA$RPKT2	/	RC2TPM_lncRNA.2
RC2TPM_lncRNA$TPMT3	<-	RC2TPM_lncRNA$RPKT3	/	RC2TPM_lncRNA.3
RC2TPM_lncRNA$TPMT4	<-	RC2TPM_lncRNA$RPKT4	/	RC2TPM_lncRNA.4
RC2TPM_lncRNA$TPMT5	<-	RC2TPM_lncRNA$RPKT5	/	RC2TPM_lncRNA.5
RC2TPM_lncRNA$TPMT6	<-	RC2TPM_lncRNA$RPKT6	/	RC2TPM_lncRNA.6
RC2TPM_lncRNA$TPMT7	<-	RC2TPM_lncRNA$RPKT7	/	RC2TPM_lncRNA.7
RC2TPM_lncRNA$TPMT8	<-	RC2TPM_lncRNA$RPKT8	/	RC2TPM_lncRNA.8
RC2TPM_lncRNA$TPMT9	<-	RC2TPM_lncRNA$RPKT9	/	RC2TPM_lncRNA.9
RC2TPM_lncRNA$TPMT10	<-	RC2TPM_lncRNA$RPKT10	/	RC2TPM_lncRNA.10
RC2TPM_lncRNA$TPMT11	<-	RC2TPM_lncRNA$RPKT11	/	RC2TPM_lncRNA.11
RC2TPM_lncRNA$TPMT12	<-	RC2TPM_lncRNA$RPKT12	/	RC2TPM_lncRNA.12
RC2TPM_lncRNA$TPMT13	<-	RC2TPM_lncRNA$RPKT13	/	RC2TPM_lncRNA.13
RC2TPM_lncRNA$TPMT14	<-	RC2TPM_lncRNA$RPKT14	/	RC2TPM_lncRNA.14
RC2TPM_lncRNA$TPMT15	<-	RC2TPM_lncRNA$RPKT15	/	RC2TPM_lncRNA.15
RC2TPM_lncRNA$TPMT16	<-	RC2TPM_lncRNA$RPKT16	/	RC2TPM_lncRNA.16
RC2TPM_lncRNA$TPMT17	<-	RC2TPM_lncRNA$RPKT17	/	RC2TPM_lncRNA.17
RC2TPM_lncRNA$TPMT18	<-	RC2TPM_lncRNA$RPKT18	/	RC2TPM_lncRNA.18
RC2TPM_lncRNA$TPMT19	<-	RC2TPM_lncRNA$RPKT19	/	RC2TPM_lncRNA.19
RC2TPM_lncRNA$TPMT20	<-	RC2TPM_lncRNA$RPKT20	/	RC2TPM_lncRNA.20
RC2TPM_lncRNA$TPMT21	<-	RC2TPM_lncRNA$RPKT21	/	RC2TPM_lncRNA.21
RC2TPM_lncRNA$TPMT22	<-	RC2TPM_lncRNA$RPKT22	/	RC2TPM_lncRNA.22
RC2TPM_lncRNA$TPMT23	<-	RC2TPM_lncRNA$RPKT23	/	RC2TPM_lncRNA.23
RC2TPM_lncRNA$TPMT24	<-	RC2TPM_lncRNA$RPKT24	/	RC2TPM_lncRNA.24
RC2TPM_lncRNA$TPMT25	<-	RC2TPM_lncRNA$RPKT25	/	RC2TPM_lncRNA.25
RC2TPM_lncRNA$TPMT26	<-	RC2TPM_lncRNA$RPKT26	/	RC2TPM_lncRNA.26
RC2TPM_lncRNA$TPMT27	<-	RC2TPM_lncRNA$RPKT27	/	RC2TPM_lncRNA.27
RC2TPM_lncRNA$TPMT28	<-	RC2TPM_lncRNA$RPKT28	/	RC2TPM_lncRNA.28
RC2TPM_lncRNA$TPMT29	<-	RC2TPM_lncRNA$RPKT29	/	RC2TPM_lncRNA.29
RC2TPM_lncRNA$TPMT30	<-	RC2TPM_lncRNA$RPKT30	/	RC2TPM_lncRNA.30
RC2TPM_lncRNA$TPMT31	<-	RC2TPM_lncRNA$RPKT31	/	RC2TPM_lncRNA.31
RC2TPM_lncRNA$TPMT32	<-	RC2TPM_lncRNA$RPKT32	/	RC2TPM_lncRNA.32
RC2TPM_lncRNA$TPMT33	<-	RC2TPM_lncRNA$RPKT33	/	RC2TPM_lncRNA.33
RC2TPM_lncRNA$TPMT34	<-	RC2TPM_lncRNA$RPKT34	/	RC2TPM_lncRNA.34
RC2TPM_lncRNA$TPMT35	<-	RC2TPM_lncRNA$RPKT35	/	RC2TPM_lncRNA.35
RC2TPM_lncRNA$TPMT36	<-	RC2TPM_lncRNA$RPKT36	/	RC2TPM_lncRNA.36
RC2TPM_lncRNA$TPMT37	<-	RC2TPM_lncRNA$RPKT37	/	RC2TPM_lncRNA.37
RC2TPM_lncRNA$TPMT38	<-	RC2TPM_lncRNA$RPKT38	/	RC2TPM_lncRNA.38
RC2TPM_lncRNA$TPMT39	<-	RC2TPM_lncRNA$RPKT39	/	RC2TPM_lncRNA.39
RC2TPM_lncRNA$TPMT40	<-	RC2TPM_lncRNA$RPKT40	/	RC2TPM_lncRNA.40
RC2TPM_lncRNA$TPMT41	<-	RC2TPM_lncRNA$RPKT41	/	RC2TPM_lncRNA.41
RC2TPM_lncRNA$TPMT42	<-	RC2TPM_lncRNA$RPKT42	/	RC2TPM_lncRNA.42
RC2TPM_lncRNA$TPMT43	<-	RC2TPM_lncRNA$RPKT43	/	RC2TPM_lncRNA.43
RC2TPM_lncRNA$TPMT44	<-	RC2TPM_lncRNA$RPKT44	/	RC2TPM_lncRNA.44
RC2TPM_lncRNA$TPMT45	<-	RC2TPM_lncRNA$RPKT45	/	RC2TPM_lncRNA.45
RC2TPM_lncRNA$TPMT46	<-	RC2TPM_lncRNA$RPKT46	/	RC2TPM_lncRNA.46
RC2TPM_lncRNA$TPMT47	<-	RC2TPM_lncRNA$RPKT47	/	RC2TPM_lncRNA.47
RC2TPM_lncRNA$TPMT48	<-	RC2TPM_lncRNA$RPKT48	/	RC2TPM_lncRNA.48
RC2TPM_lncRNA$TPMT49	<-	RC2TPM_lncRNA$RPKT49	/	RC2TPM_lncRNA.49
RC2TPM_lncRNA$TPMT50	<-	RC2TPM_lncRNA$RPKT50	/	RC2TPM_lncRNA.50
RC2TPM_lncRNA$TPMT51	<-	RC2TPM_lncRNA$RPKT51	/	RC2TPM_lncRNA.51
RC2TPM_lncRNA$TPMT52	<-	RC2TPM_lncRNA$RPKT52	/	RC2TPM_lncRNA.52
RC2TPM_lncRNA$TPMT53	<-	RC2TPM_lncRNA$RPKT53	/	RC2TPM_lncRNA.53
RC2TPM_lncRNA$TPMT54	<-	RC2TPM_lncRNA$RPKT54	/	RC2TPM_lncRNA.54
RC2TPM_lncRNA$TPMT55	<-	RC2TPM_lncRNA$RPKT55	/	RC2TPM_lncRNA.55
RC2TPM_lncRNA$TPMT56	<-	RC2TPM_lncRNA$RPKT56	/	RC2TPM_lncRNA.56
RC2TPM_lncRNA$TPMT57	<-	RC2TPM_lncRNA$RPKT57	/	RC2TPM_lncRNA.57
RC2TPM_lncRNA$TPMT58	<-	RC2TPM_lncRNA$RPKT58	/	RC2TPM_lncRNA.58

########################## CONFIRM TPM NORMALIZATION ###########################
if(sum(RC2TPM_lncRNA$TPMC1)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC1)}
if(sum(RC2TPM_lncRNA$TPMC2)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC2)}
if(sum(RC2TPM_lncRNA$TPMC3)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC3)}
if(sum(RC2TPM_lncRNA$TPMC4)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC4)}
if(sum(RC2TPM_lncRNA$TPMC5)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC5)}
if(sum(RC2TPM_lncRNA$TPMC6)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC6)}
if(sum(RC2TPM_lncRNA$TPMC7)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC7)}
if(sum(RC2TPM_lncRNA$TPMC8)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMC8)}
if(sum(RC2TPM_lncRNA$TPMT1)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT1)}
if(sum(RC2TPM_lncRNA$TPMT2)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT2)}
if(sum(RC2TPM_lncRNA$TPMT3)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT3)}
if(sum(RC2TPM_lncRNA$TPMT4)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT4)}
if(sum(RC2TPM_lncRNA$TPMT5)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT5)}
if(sum(RC2TPM_lncRNA$TPMT6)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT6)}
if(sum(RC2TPM_lncRNA$TPMT7)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT7)}
if(sum(RC2TPM_lncRNA$TPMT8)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT8)}
if(sum(RC2TPM_lncRNA$TPMT9)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT9)}
if(sum(RC2TPM_lncRNA$TPMT10)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT10)}
if(sum(RC2TPM_lncRNA$TPMT11)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT11)}
if(sum(RC2TPM_lncRNA$TPMT12)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT12)}
if(sum(RC2TPM_lncRNA$TPMT13)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT13)}
if(sum(RC2TPM_lncRNA$TPMT14)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT14)}
if(sum(RC2TPM_lncRNA$TPMT15)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT15)}
if(sum(RC2TPM_lncRNA$TPMT16)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT16)}
if(sum(RC2TPM_lncRNA$TPMT17)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT17)}
if(sum(RC2TPM_lncRNA$TPMT18)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT18)}
if(sum(RC2TPM_lncRNA$TPMT19)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT19)}
if(sum(RC2TPM_lncRNA$TPMT20)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT20)}
if(sum(RC2TPM_lncRNA$TPMT21)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT21)}
if(sum(RC2TPM_lncRNA$TPMT22)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT22)}
if(sum(RC2TPM_lncRNA$TPMT23)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT23)}
if(sum(RC2TPM_lncRNA$TPMT24)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT24)}
if(sum(RC2TPM_lncRNA$TPMT25)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT25)}
if(sum(RC2TPM_lncRNA$TPMT26)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT26)}
if(sum(RC2TPM_lncRNA$TPMT27)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT27)}
if(sum(RC2TPM_lncRNA$TPMT28)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT28)}
if(sum(RC2TPM_lncRNA$TPMT29)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT29)}
if(sum(RC2TPM_lncRNA$TPMT30)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT30)}
if(sum(RC2TPM_lncRNA$TPMT31)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT31)}
if(sum(RC2TPM_lncRNA$TPMT32)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT32)}
if(sum(RC2TPM_lncRNA$TPMT33)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT33)}
if(sum(RC2TPM_lncRNA$TPMT34)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT34)}
if(sum(RC2TPM_lncRNA$TPMT35)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT35)}
if(sum(RC2TPM_lncRNA$TPMT36)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT36)}
if(sum(RC2TPM_lncRNA$TPMT37)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT37)}
if(sum(RC2TPM_lncRNA$TPMT38)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT38)}
if(sum(RC2TPM_lncRNA$TPMT39)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT39)}
if(sum(RC2TPM_lncRNA$TPMT40)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT40)}
if(sum(RC2TPM_lncRNA$TPMT41)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT41)}
if(sum(RC2TPM_lncRNA$TPMT42)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT42)}
if(sum(RC2TPM_lncRNA$TPMT43)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT43)}
if(sum(RC2TPM_lncRNA$TPMT44)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT44)}
if(sum(RC2TPM_lncRNA$TPMT45)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT45)}
if(sum(RC2TPM_lncRNA$TPMT46)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT46)}
if(sum(RC2TPM_lncRNA$TPMT47)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT47)}
if(sum(RC2TPM_lncRNA$TPMT48)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT48)}
if(sum(RC2TPM_lncRNA$TPMT49)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT49)}
if(sum(RC2TPM_lncRNA$TPMT50)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT50)}
if(sum(RC2TPM_lncRNA$TPMT51)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT51)}
if(sum(RC2TPM_lncRNA$TPMT52)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT52)}
if(sum(RC2TPM_lncRNA$TPMT53)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT53)}
if(sum(RC2TPM_lncRNA$TPMT54)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT54)}
if(sum(RC2TPM_lncRNA$TPMT55)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT55)}
if(sum(RC2TPM_lncRNA$TPMT56)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT56)}
if(sum(RC2TPM_lncRNA$TPMT57)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT57)}
if(sum(RC2TPM_lncRNA$TPMT58)	==	1e+06){print("TRUE")}else{sum(RC2TPM_lncRNA$TPMT58)}

"red"
write.csv(RC2TPM_lncRNA,"RESULTS_workflow_lncRNA\\1.Read_counts_to_TPM_RESULTS_workflow_lncRNA.csv", row.names = TRUE)

#Making_Entire gene list without filter_DESeq result + TPM Values Control and Tumor
XL <- RC2TPM_lncRNA
dim(XL)
#[1] 127802    103
XL[1:5,1:3]

XL.heat <- XL[134:199]############################################################
dim(XL.heat)
#[1] 127802     66
"red"
write.csv(XL.heat,"RESULTS_workflow_lncRNA\\2.Only_TPM_RESULTS_workflow_lncRNA.csv", row.names = TRUE)
"brown"
XL.heat<-read.table("RESULTS_workflow_lncRNA\\2.Only_TPM_RESULTS_workflow_lncRNA.csv", header = TRUE, sep = ",")
colnames(XL.heat)[1]  <- "lncRNASymbol"    # change column name for XL column

##              ##              ##              ##              ##              ##              ##              ##              ##
"green"
DESeq2_value_lncRNA<-read.table("DESeq2_TvsN_SRA_TCGA_merge_DLBCL_lncRNA.csv", header = TRUE, sep = ",")
##              ##              ##              ##              ##              ##              ##              ##              ##
colnames(DESeq2_value_lncRNA)[1]  <- "lncRNASymbol"    # change column name for XL column
dim(DESeq2_value_lncRNA)
#[1] 127802      7
XL.1 = merge(DESeq2_value_lncRNA, XL.heat, ID="lncRNASymbol", sort = F)

# Average Control value as a new Variable (mean of Columns 8 to 12 is print in last column = 52)
XL.1 <- transform(XL.1, Mean_Contol_TPM = rowMeans(XL.1[8:15], na.rm = TRUE))  #########################################################
# Average Control value as a new Variable (mean of Columns 13 to 41 is print in last column = 53)
XL.1 <- transform(XL.1, Mean_Tumour_TPM = rowMeans(XL.1[16:73], na.rm = TRUE)) #########################################################

# Fold change (Tumour/average control)  (Fold change of Tumor from 43 to 59)
XL.1 <- data.frame(XL.1)                       # Creating example data
XL.1$FC1  <- XL.1$TPMT1  / XL.1$Mean_Contol_TPM     # Creating new variable TPM1
XL.1$FC2  <- XL.1$TPMT2  / XL.1$Mean_Contol_TPM     # Creating new variable TPM2
XL.1$FC3  <- XL.1$TPMT3  / XL.1$Mean_Contol_TPM     # Creating new variable TPM3
XL.1$FC4  <- XL.1$TPMT4  / XL.1$Mean_Contol_TPM     # Creating new variable TPM4
XL.1$FC5  <- XL.1$TPMT5  / XL.1$Mean_Contol_TPM     # Creating new variable TPM5
XL.1$FC6  <- XL.1$TPMT6  / XL.1$Mean_Contol_TPM     # Creating new variable TPM6
XL.1$FC7  <- XL.1$TPMT7  / XL.1$Mean_Contol_TPM     # Creating new variable TPM7
XL.1$FC8  <- XL.1$TPMT8  / XL.1$Mean_Contol_TPM     # Creating new variable TPM8
XL.1$FC9  <- XL.1$TPMT9  / XL.1$Mean_Contol_TPM     # Creating new variable TPM9
XL.1$FC10 <- XL.1$TPMT10 / XL.1$Mean_Contol_TPM  # Creating new variable TPM10
XL.1$FC11 <- XL.1$TPMT11 / XL.1$Mean_Contol_TPM   # Creating new variable TPM11
XL.1$FC12 <- XL.1$TPMT12 / XL.1$Mean_Contol_TPM   # Creating new variable TPM12
XL.1$FC13 <- XL.1$TPMT13 / XL.1$Mean_Contol_TPM   # Creating new variable TPM13
XL.1$FC14 <- XL.1$TPMT14 / XL.1$Mean_Contol_TPM   # Creating new variable TPM14
XL.1$FC15 <- XL.1$TPMT15 / XL.1$Mean_Contol_TPM   # Creating new variable TPM15
XL.1$FC16 <- XL.1$TPMT16 / XL.1$Mean_Contol_TPM   # Creating new variable TPM16
XL.1$FC17 <- XL.1$TPMT17 / XL.1$Mean_Contol_TPM   # Creating new variable TPM17
XL.1$FC18 <- XL.1$TPMT18 / XL.1$Mean_Contol_TPM   # Creating new variable TPM18
XL.1$FC19 <- XL.1$TPMT19 / XL.1$Mean_Contol_TPM   # Creating new variable TPM19
XL.1$FC20 <- XL.1$TPMT20 / XL.1$Mean_Contol_TPM  # Creating new variable TPM20
XL.1$FC21  <- XL.1$TPMT21  / XL.1$Mean_Contol_TPM     # Creating new variable TPM21
XL.1$FC22  <- XL.1$TPMT22  / XL.1$Mean_Contol_TPM     # Creating new variable TPM22
XL.1$FC23  <- XL.1$TPMT23  / XL.1$Mean_Contol_TPM     # Creating new variable TPM23
XL.1$FC24  <- XL.1$TPMT24  / XL.1$Mean_Contol_TPM     # Creating new variable TPM24
XL.1$FC25  <- XL.1$TPMT25  / XL.1$Mean_Contol_TPM     # Creating new variable TPM25
XL.1$FC26  <- XL.1$TPMT26  / XL.1$Mean_Contol_TPM     # Creating new variable TPM26
XL.1$FC27  <- XL.1$TPMT27  / XL.1$Mean_Contol_TPM     # Creating new variable TPM27
XL.1$FC28  <- XL.1$TPMT28  / XL.1$Mean_Contol_TPM     # Creating new variable TPM28
XL.1$FC29  <- XL.1$TPMT29  / XL.1$Mean_Contol_TPM     # Creating new variable TPM29
XL.1$FC30 <- XL.1$TPMT30 / XL.1$Mean_Contol_TPM  # Creating new variable TPM30
XL.1$FC31  <- XL.1$TPMT31  / XL.1$Mean_Contol_TPM     # Creating new variable TPM31
XL.1$FC32  <- XL.1$TPMT32  / XL.1$Mean_Contol_TPM     # Creating new variable TPM32
XL.1$FC33  <- XL.1$TPMT33  / XL.1$Mean_Contol_TPM     # Creating new variable TPM33
XL.1$FC34  <- XL.1$TPMT34  / XL.1$Mean_Contol_TPM     # Creating new variable TPM34
XL.1$FC35  <- XL.1$TPMT35  / XL.1$Mean_Contol_TPM     # Creating new variable TPM35
XL.1$FC36  <- XL.1$TPMT36  / XL.1$Mean_Contol_TPM     # Creating new variable TPM36
XL.1$FC37  <- XL.1$TPMT37  / XL.1$Mean_Contol_TPM     # Creating new variable TPM37
XL.1$FC38  <- XL.1$TPMT38  / XL.1$Mean_Contol_TPM     # Creating new variable TPM38
XL.1$FC39  <- XL.1$TPMT39  / XL.1$Mean_Contol_TPM     # Creating new variable TPM39
XL.1$FC40 <- XL.1$TPMT40 / XL.1$Mean_Contol_TPM  # Creating new variable TPM40
XL.1$FC41  <- XL.1$TPMT41  / XL.1$Mean_Contol_TPM     # Creating new variable TPM41
XL.1$FC42  <- XL.1$TPMT42  / XL.1$Mean_Contol_TPM     # Creating new variable TPM42
XL.1$FC43  <- XL.1$TPMT43  / XL.1$Mean_Contol_TPM     # Creating new variable TPM43
XL.1$FC44  <- XL.1$TPMT44  / XL.1$Mean_Contol_TPM     # Creating new variable TPM44
XL.1$FC45  <- XL.1$TPMT45  / XL.1$Mean_Contol_TPM     # Creating new variable TPM45
XL.1$FC46  <- XL.1$TPMT46  / XL.1$Mean_Contol_TPM     # Creating new variable TPM46
XL.1$FC47  <- XL.1$TPMT47  / XL.1$Mean_Contol_TPM     # Creating new variable TPM47
XL.1$FC48  <- XL.1$TPMT48  / XL.1$Mean_Contol_TPM     # Creating new variable TPM48
XL.1$FC49  <- XL.1$TPMT49  / XL.1$Mean_Contol_TPM     # Creating new variable TPM49
XL.1$FC50 <- XL.1$TPMT50 / XL.1$Mean_Contol_TPM  # Creating new variable TPM50
XL.1$FC51  <- XL.1$TPMT51  / XL.1$Mean_Contol_TPM     # Creating new variable TPM51
XL.1$FC52  <- XL.1$TPMT52  / XL.1$Mean_Contol_TPM     # Creating new variable TPM52
XL.1$FC53  <- XL.1$TPMT53  / XL.1$Mean_Contol_TPM     # Creating new variable TPM53
XL.1$FC54  <- XL.1$TPMT54  / XL.1$Mean_Contol_TPM     # Creating new variable TPM54
XL.1$FC55  <- XL.1$TPMT55  / XL.1$Mean_Contol_TPM     # Creating new variable TPM55
XL.1$FC56  <- XL.1$TPMT56  / XL.1$Mean_Contol_TPM     # Creating new variable TPM56
XL.1$FC57  <- XL.1$TPMT57  / XL.1$Mean_Contol_TPM     # Creating new variable TPM57
XL.1$FC58  <- XL.1$TPMT58  / XL.1$Mean_Contol_TPM     # Creating new variable TPM58
"red"
write.csv(XL.1,"RESULTS_workflow_lncRNA\\3.L2FC_RESULTS_workflow_lncRNA.csv", row.names = TRUE)
"brown"
XL.1<-read.table("RESULTS_workflow_lncRNA\\3.L2FC_RESULTS_workflow_lncRNA.csv", header = TRUE, sep = ",", row.names = 2)
XL.1[1:5,1:6]
XL.1 <- XL.1[ -c(1) ] # Removes entire column be careful
XL.1[1:5,1:6]
#log2 transform data
XL.log <- log(XL.1[75:132], 2)#################"XL.log must have 16 Variables"####################################
XL.log[XL.log=="-Inf"]<-NA
colnames(XL.log) <- c('L2FC1','L2FC2','L2FC3','L2FC4','L2FC5','L2FC6','L2FC7','L2FC8','L2FC9','L2FC10','L2FC11','L2FC12','L2FC13','L2FC14','L2FC15','L2FC16','L2FC17','L2FC18','L2FC19','L2FC20','L2FC21','L2FC22','L2FC23','L2FC24','L2FC25','L2FC26','L2FC27','L2FC28','L2FC29','L2FC30','L2FC31','L2FC32','L2FC33','L2FC34','L2FC35','L2FC36','L2FC37','L2FC38','L2FC39','L2FC40','L2FC41','L2FC42','L2FC43','L2FC44','L2FC45','L2FC46','L2FC47','L2FC48','L2FC49','L2FC50','L2FC51','L2FC52','L2FC53','L2FC54','L2FC55','L2FC56','L2FC57','L2FC58') #Rename the columns to L2FC
head(XL.log)
"red"
write.csv(XL.log,"RESULTS_workflow_lncRNA\\4.L2FC_only_RESULTS_workflow_lncRNA.csv", row.names = TRUE)
"brown"
zL<-read.table("RESULTS_workflow_lncRNA\\4.L2FC_only_RESULTS_workflow_lncRNA.csv", header = TRUE, sep = ",")
zL[1:5,1:6]
colnames(zL)[1]  <- "lncRNASymbol"    # change column name for XL column
zL[1:5,1:6]
"brown"
yL<-read.table("RESULTS_workflow_lncRNA\\3.L2FC_RESULTS_workflow_lncRNA.csv", header = TRUE, sep = ",")
yL[1:5,1:6]
yL <- yL[-c(1)]
yL[1:5,1:6]
#******************************************************************************#
lncRNA_All_facotrs_merged_norow = merge(yL, zL, ID="lncRNASymbol", sort = F)            #
lncRNA_All_facotrs_merged <- lncRNA_All_facotrs_merged_norow[,-1]                            #
rownames(lncRNA_All_facotrs_merged) <- lncRNA_All_facotrs_merged_norow[,1]                   #
dim(lncRNA_All_facotrs_merged)                                                        #
#[1] 127802    190                                                             #
#******************************************************************************#
"red"
write.csv(lncRNA_All_facotrs_merged,"RESULTS_workflow_lncRNA\\5.All_factors_merged_TPM_DLBCL2_lncRNA.csv", row.names = TRUE)
#
##########
##################################################################################
#########################Pick lncRNAs corresponding to shortlisted mRNA #########################################################
##################################################################################
######### manually stip the lncRNA annotaion and keep the HUGO gene Symbol ###############
search4_lncRNA_using_mRNA <- c(rownames(pvalue_0.05_L2FC_2_filtered_MOD)) # creating values for the 9 mRNAs to filter its corresponding lncRNAs

#lncRNA_All_facotrs_merged<-read.table("RESULTS_workflow_lncRNA\\5.All_factors_merged_TPM_DLBCL2_lncRNA_strip.csv", header = TRUE, sep = ",")
"green"
#file_path_1 <- read_csv("RESULTS_workflow_lncRNA\\5.All_factors_merged_TPM_DLBCL2_lncRNA_striplist.csv")  # Update with the correct path to your CSV file
#lncRNA_All_facotrs_merged <- read_csv(file_path_1)

lncRNA_All_facotrs_merged <- read.table("RESULTS_workflow_lncRNA\\5.All_factors_merged_TPM_DLBCL2_lncRNA_strip.csv", header = TRUE, sep = ",")

library(data.table)
lncRNA_All_facotrs_merged_all <- data.table(lncRNA_All_facotrs_merged)
df_lncRNA_pathway_Annotaion_genesenlisted_filtered <- lncRNA_All_facotrs_merged_all[lncRNA_All_facotrs_merged_all$geneSymbol %in% search4_lncRNA_using_mRNA,]   #

dim(df_lncRNA_pathway_Annotaion_genesenlisted_filtered)  
#[1]  224 193

df_lncRNA_pathway_Annotaion_genesenlisted_filtered <- data.frame(df_lncRNA_pathway_Annotaion_genesenlisted_filtered)
rownames(df_lncRNA_pathway_Annotaion_genesenlisted_filtered)<- df_lncRNA_pathway_Annotaion_genesenlisted_filtered[,1]
df_lncRNA_pathway_Annotaion_genesenlisted_filtered <- df_lncRNA_pathway_Annotaion_genesenlisted_filtered[,-1]
dim(df_lncRNA_pathway_Annotaion_genesenlisted_filtered)
#[1] 224 192
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
############################################## SELECTING ONLY GENES WITH pvalue <= 0.055 ####################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
## 
library(dplyr)
pvalue_filter_lncRNA_MOD <- df_lncRNA_pathway_Annotaion_genesenlisted_filtered %>% filter(pvalue <= 0.055) 
dim(pvalue_filter_lncRNA_MOD)
#[1] 140 192
"red"
write.csv(pvalue_filter_lncRNA_MOD,"RESULTS_workflow_lncRNA\\6.pvalue_sig_DLBCL2_SRA_TCGA_merge_lncRNA_MOD.csv", row.names = TRUE)
#**********************************************************************************#
rownames(XL.heat)<- XL.heat[,1]
XL.heat <- XL.heat[,-1]
#rownames(pvalue_filter_lncRNA)<- pvalue_filter_lncRNA[,1]
pvalue_0.055_filtered_TPM_lncRNA_MOD <- XL.heat[rownames(XL.heat) %in% rownames(pvalue_filter_lncRNA_MOD),]   #
dim(pvalue_0.055_filtered_TPM_lncRNA_MOD)                                                        #
#[1] 140 66  *****************************************************************#
#**********************************************************************************#
"red"
write.csv(pvalue_0.055_filtered_TPM_lncRNA_MOD,"RESULTS_workflow_lncRNA\\6.pvalue_sig_upr_DLBCL2_SRA_TCGA_merge_lncRNA_TPM1ly_MOD.csv", row.names = TRUE)
#*#****************************************************************************************************************###########################################################################################################################################
##########################################################################################################################################
########################### SELECTING ONLY GENES WITH pvalue <= 0.055 and L2FC >2 EXPRESSION IN ALL SAMPLES ################################
##########################################################################################################################################

##########################################################################################################################################

#**********************************************************************************************************#
pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_MOD <- pvalue_filter_lncRNA_MOD[c(135:192)]
dim(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_MOD)
#[1] 140 58
upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD <- pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_MOD[rowSums(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_MOD >= 2) == ncol(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_MOD), ]################ L2FC values >= 2 in all samples.#################
upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD <- na.omit(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD)
dim(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD)
#[1] 2 58
"red"
write.csv(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD,"RESULTS_workflow_lncRNA\\7.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_MOD.csv", row.names = TRUE)

upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_MOD <- XL.heat[rownames(XL.heat) %in% rownames(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD),]      #
dim(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_MOD)                                                                         #
#[1] 2 66  *********************************************************************************************#
#**********************************************************************************************************#
"red"
write.csv(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_MOD,"RESULTS_workflow_lncRNA\\7.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_TPM1ly_MOD.csv", row.names = TRUE)
#**********************************************************************************************************#
upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged_MOD <- df_lncRNA_pathway_Annotaion_genesenlisted_filtered[rownames(df_lncRNA_pathway_Annotaion_genesenlisted_filtered) %in% rownames(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_MOD),]
dim(upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged_MOD)
#[1] 2 192
"red"
write.csv(upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged_MOD,"RESULTS_workflow_lncRNA\\7.pvalue_pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_AF_MOD.csv", row.names = TRUE)
#****************************************************************************************************************#
##########################################################################################################################################
########################### SELECTING ONLY GENES WITH pvalue <= 0.055 and L2FC < -2 EXPRESSION IN ALL SAMPLES ################################
##########################################################################################################################################
pvalue_filter_lncRNA <- lncRNA_All_facotrs_merged %>% filter(pvalue <= 0.055)
dim(pvalue_filter_lncRNA)
# [1] 44267   193
"red"
write.csv(pvalue_filter_lncRNA,"RESULTS_workflow_lncRNA\\6.pv_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE)

search4_lncRNA_all <- c(rownames(pvalue_0.05_L2FC_2_filtered_TPM))

pvalue_filter_lncRNA_merged <- data.table(pvalue_filter_lncRNA)
df_lncRNA_all_genesenlisted_filtered <- pvalue_filter_lncRNA_merged[pvalue_filter_lncRNA_merged$geneSymbol %in% search4_lncRNA_all,]   #

dim(df_lncRNA_all_genesenlisted_filtered)  
# [1]   677 192

rownames(pvalue_filter_lncRNA)<- pvalue_filter_lncRNA[,1]
pvalue_filter_lncRNA <- pvalue_filter_lncRNA[,-1]
dim(pvalue_filter_lncRNA)
# [1] 44267   192

df_lncRNA_all_genesenlisted_filtered <- data.frame(df_lncRNA_all_genesenlisted_filtered)
rownames(df_lncRNA_all_genesenlisted_filtered)<- df_lncRNA_all_genesenlisted_filtered[,1]
df_lncRNA_all_genesenlisted_filtered <- df_lncRNA_all_genesenlisted_filtered[,-1]
dim(df_lncRNA_all_genesenlisted_filtered)
# [1] 677   192

pvalue_0.055_filtered_TPM_lncRNA <- XL.heat[rownames(XL.heat) %in% rownames(df_lncRNA_all_genesenlisted_filtered),]   #
dim(pvalue_0.055_filtered_TPM_lncRNA)                                                        #
#[1] 677 66  *****************************************************************#
#**********************************************************************************#
"red"
write.csv(pvalue_0.055_filtered_TPM_lncRNA,"RESULTS_workflow_lncRNA\\6.pvalue_DLBCL2_SRA_TCGA_merge_lncRNA_TPM1ly.csv", row.names = TRUE)
#*#****************************************************************************************************************###########################################################################################################################################

#**********************************************************************************************************#

pvalue_lncRNA_0.055_L2FC_2_filter_intermediate <- pvalue_filter_lncRNA[c(135:192)]
dim(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate)
#[1] 44267    58
upr_pvalue_lncRNA_0.055_L2FC_2_filtered <- pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rowSums(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate >= 2) == ncol(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate), ]################ L2FC values >= 2 in all samples.#################
upr_pvalue_lncRNA_0.055_L2FC_2_filtered <- na.omit(upr_pvalue_lncRNA_0.055_L2FC_2_filtered)
dim(upr_pvalue_lncRNA_0.055_L2FC_2_filtered)
#[1] 486  58
"red"
write.csv(upr_pvalue_lncRNA_0.055_L2FC_2_filtered,"RESULTS_workflow_lncRNA\\7.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE)


upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM <- XL.heat[rownames(XL.heat) %in% rownames(upr_pvalue_lncRNA_0.055_L2FC_2_filtered),]      #
dim(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM)                                                                         #
#[1] 486  66  *********************************************************************************************#
#**********************************************************************************************************#
"red"
write.csv(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM,"RESULTS_workflow_lncRNA\\7.pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_TPM1ly.csv", row.names = TRUE)
#**********************************************************************************************************#
upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged <- pvalue_filter_lncRNA[rownames(pvalue_filter_lncRNA) %in% rownames(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM),]
dim(upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged)
#[1] 486 192
"red"
write.csv(upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged,"RESULTS_workflow_lncRNA\\7.pvalue_pv_l2fc_upr_DLBCL2_SRA_TCGA_merge_lncRNA_AF.csv", row.names = TRUE)
#****************************************************************************************************************#


dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered <- pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rowSums(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate <= -2) == ncol(pvalue_lncRNA_0.055_L2FC_2_filter_intermediate), ]################ L2FC values >= 2 in all samples.#################
dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered <- na.omit(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered)
dim(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered)
# [1] 366 58
"red"
write.csv(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered,"RESULTS_workflow_lncRNA\\9.pv_l2fc_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE)


dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM <- XL.heat[rownames(XL.heat) %in% rownames(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered),]      #
dim(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM)                                                                         #
#[1] 366 66  *********************************************************************************************#
#**********************************************************************************************************#
"red"
write.csv(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM,"RESULTS_workflow_lncRNA\\9.pv_l2fc_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA_TPM1ly.csv", row.names = TRUE)
#**********************************************************************************************************#
dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_all_factors_merged <- pvalue_filter_lncRNA[rownames(pvalue_filter_lncRNA) %in% rownames(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM),]
dim(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_all_factors_merged)
#[1] 366 192
"red"
write.csv(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_all_factors_merged,"RESULTS_workflow_lncRNA\\9.pvalue_pv_l2fc_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA_AF.csv", row.names = TRUE)


pvalue_0.05_L2FC_2_filtered_lncRNA <- rbind(upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged, dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_all_factors_merged)
dim(pvalue_0.05_L2FC_2_filtered_lncRNA)
#[1] 852 192

"red"
write.csv(pvalue_0.05_L2FC_2_filtered_lncRNA,"RESULTS_workflow_lncRNA\\10.pv_l2fc_upr_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE)
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC > 2,  Mean_Tumour_TPM > 1 and EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
only_TTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered <- upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged[c(17:74)]
dim(only_TTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered)
#[1] 486    8
#cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered <- cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rownames(cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate) %in% rownames(TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM),]
#dim(cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered)
#[1]18  18

only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered <- only_TTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered[rowSums(only_TTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered > 1) == ncol(only_TTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered), ] 
dim(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered)
# 486    8
only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM <- XL.heat[rownames(XL.heat) %in% rownames(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered),]      #
dim(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM)                                                                         #
#[1] 436  66 



##########################################################################################################################################
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC > 2,  Mean_Tumour_TPM < 1 and Mean_Normal_TPM > 1 EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered <- upr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged[c(9:16)]
dim(cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered)
#[1] 486    8
#cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered <- cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rownames(cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate) %in% rownames(TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM),]
#dim(cTPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered)
#[1]18  18

cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_lnc <- cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered[rowSums(cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered < 1) == ncol(cTPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate_filtered), ] 
dim(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_lnc)
# 486    8
cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM_lnc <- XL.heat[rownames(XL.heat) %in% rownames(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_lnc),]      #
dim(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM)                                                                         #
#[1] 436  66 

##########################################################################################################################################
##########################################################################################################################################
############### SELECTING  GENES WITH pvalue <= 0.05 and L2FC>2, Normal_TPM < 1 and Mean_Tumour_TPM < 1 EXPRESSION IN ALL SAMPLES #######################
##########################################################################################################################################
##########################################################################################################################################
TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate <- cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM[c(9:66)]
dim(TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate)
#[1] 436    58
#TPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered <- TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rownames(TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate) %in% rownames(dwnr_pvalue_0.055_L2FC_2_filtered_TPM_all_factors_merged),]
#dim(TPM_pvalue_0.05_l2fc_2_TPM_intermidiate_filtered)
#[1]34  16

TPM_pvalue_0.05_l2fc2_TPM_1_filtered <- TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate[rowSums(TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate > 1) == ncol(TPM_pvalue_lncRNA_0.055_L2FC_2_filter_intermediate), ] 
dim(TPM_pvalue_0.05_l2fc2_TPM_1_filtered)
# 78 58
TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM <- XL.heat[rownames(XL.heat) %in% rownames(TPM_pvalue_0.05_l2fc2_TPM_1_filtered),]      #
dim(TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM)                                                                         #
#[1] 78 66  
##############***************&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$#########################@@@@@@@@@@@@@@@@@@@@@@@

##################################################
#impute missing data                             #
##################################################
#1_lncRNA_pvalue
pvalue_0.055_filtered_TPM_lncRNA.heatmap.imp<-impute.knn(as.matrix(pvalue_0.055_filtered_TPM_lncRNA))$data 
write.csv(pvalue_0.055_filtered_TPM_lncRNA.heatmap.imp,"RESULTS_workflow_lncRNA\\10.forheatmap_pv_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE, col.names = TRUE)
#2_lncRNA_pvalue_l2fc_upr
upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp<-impute.knn(as.matrix(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM))$data            #### $data    add if see Error in data.imputed[1:features,]: incorrect number of dimensions
write.csv(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp,"RESULTS_workflow_lncRNA\\10.forheatmap__upr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE, col.names = TRUE)
#3_lncRNA_pvalue_l2fc_dwnr
dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp<-impute.knn(as.matrix(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM))$data            #### $data    add if see Error in data.imputed[1:features,]: incorrect number of dimensions
write.csv(dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp,"RESULTS_workflow_lncRNA\\10.forheatmap_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE, col.names = TRUE)
#5_lncRNA_pvalue_l2fc_upr_dwnr_tumourTPM
TPM_upr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(TPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM))$data 
write.csv(TPM_upr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow_lncRNA\\10.forheatmap_pv_dupr_TPM_DLBCL2_SRA_TCGA_merge_lncRNA.csv", row.names = TRUE, col.names = TRUE)
#6_lncRNA_pvalue_l2fc_upr_dwnr_tumourTPM_cTPM
cTPM_dwnr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(cTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM))$data 
write.csv(cTPM_dwnr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow_lncRNA\\10.forheatmap_pv_dwnr_TPM_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)
#7_lncRNA_pvalue_l2fc_upr_dwnr_tumour_tumour_TPM_only
only_TTPM_upr_highexp_pv_L2FC2.heatmap.imp<-impute.knn(as.matrix(only_TTPM_pvalue_0.05_l2fc2_TPM_1_filtered_TPM))$data 
write.csv(only_TTPM_upr_highexp_pv_L2FC2.heatmap.imp,"RESULTS_workflow_lncRNA\\11.forheatmap_pv_upr_only_TTPM_DLBCL2_SRA_TCGA_merge_mRNA.csv", row.names = TRUE, col.names = TRUE)


##################################################
#Plot Up-regulated heatmap with column clustering#
##################################################
#1_lncRNA_pvalue
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\5_pvalue_lncRNA_DLBCL2_heat.pdf")
myHeatmap((pvalue_0.055_filtered_TPM_lncRNA.heatmap.imp))     ############ plot Heatmap ##################
dev.off()
##################################################
#2_lncRNA_pvalue_l2fc_upr
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\6_pvalue_upr_l2fc_lncRNA_DLBCL2_heat.pdf")
myHeatmap((upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##################################################
#3_lncRNA_pvalue_l2fc_dwnr
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\7_pvalue_dwnr_l2fc_lncRNA_DLBCL2_heat.pdf")
myHeatmap((dwnr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##################################################
#4_lncRNA_pvalue_l2fc_dwnr_TPM
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\8_pvalue_upr_l2fc_TPM_lncRNA_DLBCL2_heat.pdf")
myHeatmap((TPM_upr_highexp_pv_L2FC2.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##################################################
#5_lncRNA_pvalue_l2fc_dwnr_TPM
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\9_pvalue_dwnr_l2fc_c_T_TPM_lncRNA_DLBCL2_heat.pdf")
myHeatmap((cTPM_dwnr_highexp_pv_L2FC2.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##################################################
#6_lncRNA_pvalue_l2fc_upr_only_TTPM
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\10_pvalue_upr_l2fc_only_T_TPM_lncRNA_DLBCL2_heat.pdf")
myHeatmap((only_TTPM_upr_highexp_pv_L2FC2.heatmap.imp))  ############ plot Heatmap ##################
dev.off() 
##############***************&&&&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$#########################@@@@@@@@@@@@@@@@@@@@@@@
pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM <- tibble(upr_pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap.imp)
pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM <- data.frame(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM)
#write.csv(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM,"RESULTS_workflow\\10.for_heatmap_mRNA.csv", row.names = TRUE)

#Plot Up-regulated heatmap without column clustering
"brown"
pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap<-read.table("RESULTS_workflow_lncRNA\\10.forheatmap_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.csv", header = TRUE, sep = ",")
colnames(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap)[1]  <- "lncRNASymbol"    # change column name for x column
#pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap <- data.matrix(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap) 10.forheatmap_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.csv
#####################################
# Basic Heatmap                     #
library(gplots)                     #
library(pheatmap)                   #
#####################################
################################################################################################################################################
pdf(file="RESULTS_workflow\\Figures_SRA_TCGA_merge_DLBCL2_mRNA\\7_highexp_dwnr_DLBCL2_SRA_TCGA_merge_lncRNA.pdf")
pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered <- data.frame(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC1,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC2,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC3,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC4,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC5,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC6,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC7,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC8,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC9,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC10,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC11,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC12,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC13,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC14,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC15,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC16,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC17,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMC18,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT1,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT2,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT3,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT4,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT5,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT6,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT7,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT8,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT9,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT10,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT11,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT12,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT13,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT14,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT15,pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$TPMT16)
colnames(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered) <- c("TPMC1","TPMC2","TPMC3","TPMC4","TPMC5","TPMC6","TPMC7","TPMC8","TPMC9","TPMC10","TPMC11","TPMC12","TPMC13","TPMC14","TPMC15","TPMC16","TPMC17","TPMC18","TPMT1","TPMT2","TPMT3","TPMT4","TPMT5","TPMT6","TPMT7","TPMT8","TPMT9","TPMT10","TPMT11","TPMT12","TPMT13","TPMT14","TPMT15","TPMT16")
rownames(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered) <- pvalue_lncRNA_0.055_L2FC_2_filtered_TPM.heatmap$lncRNASymbol
pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered <- data.matrix(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered)
pheatmap(pvalue_lncRNA_0.055_L2FC_2_filtered_TPM_TPM_filtered ,scale = "row", cluster_rows=FALSE, cluster_cols=FALSE, color = colorRampPalette(c("green","red"))(50),margins=c(3,25),fontsize =8)
dev.off()
