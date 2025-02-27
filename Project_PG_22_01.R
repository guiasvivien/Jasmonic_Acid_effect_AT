setwd("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data")
getwd()

#------ loading libraries ------
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(genefilter) 
library(ggplot2)
library(corrplot)
library(UpSetR)
library(tidyverse)
library(dplyr)
library(tidyr)
library(VennDiagram)
library(ggVennDiagram)
library(Glimma)
library(DataCombine)
library(viridisLite)
library(UpSetR)
library(viridis)



######### ----------------------- DESeq2 ---------------------------- ##############
###1. Calculating dds

#loading the sample table
sampleTable <- read.csv("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/sample_table.csv")
#View(sampleTable)

#transforming the Condition into a factor
sampleTable$Condition <- as.factor(sampleTable$Condition)

#reading the raw matrix
raw_counts <-read.csv("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/data_1.csv")
#View(raw_counts)

#transforming the GeneID column into rownames
raw_counts<-column_to_rownames(raw_counts, var='GeneID')
#View(raw_counts)
raw_counts <- round(raw_counts)
nrow(raw_counts)
head(raw_counts)
unique=unique(raw_counts)



nrow(unique)
#Calculating the dds
dds<-DESeqDataSetFromMatrix(countData = unique, colData = sampleTable,
                            design = ~Condition)
#View(dds)




############## --------------- Apply regularized algorithm -------------------###########

## we calculate the sums of rows 
row_sums<-rowSums(counts(dds))

## we eliminate the ones with rowcounts 0 because it cannot be logaritmised 
dds_noZero=dds[row_sums>0]
rld<-vst(dds_noZero, blind=TRUE) 
rld

#extracting the vallues as assay
rld_values<-assay(rld)
rld_values

ggplot(rld_values, aes(x=T0_1 , y=Mock_T6_1))+
  geom_point() +
  theme_bw()




############### -------------------- 3. NORMALIZE DATA WITH DEseq Normalization Factor-------------------################

#estimate the size factor = SF for normalized counts table
dds_SF<-estimateSizeFactors(dds)
#dds_SF

#we use the SF to normalize the raw counts , (normalized=T)
normalized_counts<-counts(dds_SF, normalized=T)
normalized_counts

ggplot(normalized_counts, aes(x=T0_1 , y=Mock_T6_1))+
  geom_point() +
  theme_bw()
nrow(normalized_counts) ### resulted in 26638 differentiated genes

#######---------------------- OVERVIEW PLOTS -------------------------- 

# 1------- HEATMAP

### sample distances REPLICATES 
sampleDists <-dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
#View(sampleDistMatrix)

#setting colors for the heatmaps 
colors<-colorRampPalette(rev(brewer.pal(9,"Greens")))(250)

#plotting the heatmaps 
graphics.off()
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors, border_color = 'white')

# 2 ---------- PCA 

## getting the PCA DATA

pcaData<-plotPCA(rld, intgroup=c('Condition', 'Treatment'), returnData=T)

#transform as percentage 
percentVar<-round(100*attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Condition, shape = Treatment)) +
  geom_point(size=3.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw(base_size = 12)





############## 3 ------- Top Variable Genes ---------------------#############

# 3.1 Finding most variable Genes 
topvarGenes<-head(order(rowVars(assay(rld)), decreasing= T), 50)
#print(topvarGenes) 

# to be able to visualize it we need to use the assay function 
topVarGenes_arabidopsis<-assay(rld)[topvarGenes, ]
#View(topVarGenes_arabidopsis)
## 3.2 Pheatmap 

## purple-upregulated, green is downregulated

pheatmap( topVarGenes_arabidopsis, scale = 'row',
          color = colorRampPalette( rev(brewer.pal(9, "PRGn")) )(255))


#### ------- Gene addonation ------
adno_table <- read.delim('D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mercator_1.csv', sep = ",")
#View(adno_table)
#head(adno_table)
# we extract only the needed columns 
adno_GeneName<-select(adno_table, GenID, Gene_Name, Gene_Family)
head(adno_GeneName)


# MERGING THE ANNOTATION TABLE WITH TOPVARGENES 
### first we transform the topVarGenes intro a data frame so we can
### use the tibble function to transform the rownames into a new column, we need this
#to then mergebthe two tables

topVarGenes_arabidopsis<-as.data.frame(topVarGenes_arabidopsis)
topVarGenes_arabidopsis<-tibble::rownames_to_column(topVarGenes_arabidopsis, "GenID")
#View(topVarGenes_arabidopsis)
rld_values_df <- tibble::rownames_to_column(as.data.frame(rld_values))

### now we can merge the mercator table and the topVarGenes table
topVarGenes_arabidopsis<-merge(adno_GeneName, topVarGenes_arabidopsis, by="GenID", all.y=T )

View(topVarGenes_arabidopsis)
topVarGenes_arabidopsis_anno<-topVarGenes_arabidopsis%>% tidyr:: unite("Gene_Family",
                                                                       GenID:Gene_Family,
                                                                       remove=T, sep="/")
View(topVarGenes_arabidopsis_anno)
topVarGenes_arabidopsis_anno<-tibble::column_to_rownames(topVarGenes_arabidopsis_anno,
                                                         "Gene_Family")
View(topVarGenes_arabidopsis_anno)

pheatmap(topVarGenes_arabidopsis_anno, scale='row', 
         color=colorRampPalette(rev(brewer.pal(9, "PRGn")))(255))


###########-------------------- PAIRWISE COMPARISON -----------------###########


##-----------------
## 1) DEGs analysis
##-----------------

#calculating Differential Expression table for normalized dataset (dds_SF)

dds_DEGs<-DESeq(dds_SF)

# we need to extract the results 
res<-results(dds_DEGs)
summary(res)
# what we see printed in the console are the results from a "random" pairwise 
#comparison" -> later we have to retrieve the comparisons we are interested in

res

# notice that we only have the baseMean (which is the mean between the two 
#conditions), not the baseMean for each condition -> this is why we have to 
#calculate baseMeanPerLvl
#View(dds_DEGs$Condition)
baseMeanPerLvl<-sapply(levels(dds_DEGs$Condition), function(lvl)
  rowMeans (counts(dds_DEGs, normalized = T)[,dds_DEGs$Condition == lvl]))

head(baseMeanPerLvl)

# for further usage we transform it in a data frame
baseMeanPerLvl<-as.data.frame(baseMeanPerLvl)

#for merginf we need the rownames 
baseMeanPerLvl_tab<-tibble::rownames_to_column(baseMeanPerLvl, "DmID")
baseMeanPerLvl
#View(baseMeanPerLvl)
colnames(baseMeanPerLvl_tab)
ncol(baseMeanPerLvl)
colnames(baseMeanPerLvl)[1]

#"DmID"
#"MeJA_0.5"
#"MeJA_1"  
#"MeJA_10" 
#"MeJA_16" 
#"MeJA_2"  
#"MeJA_3"
#"MeJA_6"  
#"Mock_0.5"
#"Mock_1"  
#"Mock_10" 
#"Mock_16" 
#"Mock_2"
#"Mock_3" 
#"Mock_6" 
#"time_0"

# if we want to have positive values as UPREGULATED GENES after the application of treatment, we have to divide a bigger number to a smaller number (for calculating the Fold Chnage value)
# CondA = treatment (the bigger number) -> expression goes up
# CondB = control (the smaller number) -> the ase where we start from (our comparison)
# therefore we divide CondA to CondB = treatment / control (eg: Trap_APs-1h_ERR / Trap_Control_ERR)


#---------------------------------
#	MetJA_16h vs  Mock_16h 
#----------------------------------

Mt_16_vs_Mck_16=results(dds_DEGs, contrast=c("Condition","MeJA_16", "Mock_16" 
))
summary(Mt_16_vs_Mck_16)
Mt_16_vs_Mck_16

#baseMeanA- is Mt_16
colnames(baseMeanPerLvl)[4]
Mt_16_vs_Mck_16$baseMeanA=baseMeanPerLvl[,4]

#baseMeanB - is Mck_16
colnames(baseMeanPerLvl)[11]
Mt_16_vs_Mck_16$baseMeanB=baseMeanPerLvl[,11]

Mt_16_vs_Mck_16
# to add suffix we need to change it into a data frame
Mt_16_vs_Mck_16<-data.frame(Mt_16_vs_Mck_16)
Mt_16_vs_Mck_16_colName<- Mt_16_vs_Mck_16 %>% dplyr::rename_all(paste0, "_Mt_16_vs_Mck_16") # add suffix 
#View(Mt_16_vs_Mck_16_colName)





#-----------------------------------
#	MetJA_16h vs T0 
#-----------------------------------


Mt_16_vs_T_0=results(dds_DEGs, contrast=c("Condition","MeJA_16","time_0"))
summary(Mt_16_vs_T_0)
Mt_16_vs_T_0

#baseMeanA- is control_15min
Mt_16_vs_T_0$baseMeanA=baseMeanPerLvl[,4]
#baseMeanB - is SPB00525_60min
colnames(baseMeanPerLvl)[15]
Mt_16_vs_T_0$baseMeanB=baseMeanPerLvl[,15]
Mt_16_vs_T_0

# to add suffix we need to change it into a data frame
Mt_16_vs_T_0<-data.frame(Mt_16_vs_T_0)
Mt_16_vs_T_0_colName<- Mt_16_vs_T_0 %>% dplyr::rename_all(paste0, "_Mt_16_vs_T_0") # add suffix 
#View(Mt_16_vs_T_0_colName)








#-----------------------------------
# MetJA_1h_vs_Mock_16h 
#-----------------------------------


Mt_1_vs_Mck_16=results(dds_DEGs, contrast=c("Condition","MeJA_1","Mock_16"))
summary(Mt_1_vs_Mck_16)
Mt_1_vs_Mck_16

#baseMeanA- is control_15min
colnames(baseMeanPerLvl)[2]
Mt_1_vs_Mck_16$baseMeanA=baseMeanPerLvl[,2]
#baseMeanB - is SPB00525_60min
Mt_1_vs_Mck_16$baseMeanB=baseMeanPerLvl[,11]
Mt_1_vs_Mck_16

# to add suffix we need to change it into a data frame
Mt_1_vs_Mck_16<-data.frame(Mt_1_vs_Mck_16)
Mt_1_vs_Mck_16_colName<- Mt_1_vs_Mck_16 %>% dplyr::rename_all(paste0, "_Mt_1_vs_Mck_16") # add suffix 
View(Mt_1_vs_Mck_16_colName)








#-----------------------------------
#MetJA_1h vs T0  
#-----------------------------------


Mt_1_vs_T_0=results(dds_DEGs, contrast=c("Condition","MeJA_1","time_0"))
summary(Mt_1_vs_T_0)
Mt_1_vs_T_0

#baseMeanA- is control_15min
Mt_1_vs_T_0$baseMeanA=baseMeanPerLvl[,2]
#baseMeanB - is SPB00525_60min
Mt_1_vs_T_0$baseMeanB=baseMeanPerLvl[,15]
Mt_1_vs_T_0

# to add suffix we need to change it into a data frame
Mt_1_vs_T_0<-data.frame(Mt_1_vs_T_0)
Mt_1_vs_T_0_colName<- Mt_1_vs_T_0 %>% dplyr::rename_all(paste0, "_Mt_1_vs_T_0") # add suffix 
#View(Mt_1_vs_T_0_colName)



################################################################################

##########--------------------MERGING THE TABLES -------------------------######

Mt_16_vs_Mck_16_colName<-tibble::rownames_to_column(Mt_16_vs_Mck_16_colName, "DmID")
Mt_16_vs_T_0_colName<-tibble::rownames_to_column(Mt_16_vs_T_0_colName, "DmID")
Mt_1_vs_Mck_16_colName<-tibble::rownames_to_column(Mt_1_vs_Mck_16_colName, "DmID")
Mt_1_vs_T_0_colName<-tibble::rownames_to_column(Mt_1_vs_T_0_colName, "DmID")
##########   NOW WE MERGE THEM

merged_DEGs_tabs<-merge(Mt_16_vs_Mck_16_colName, Mt_16_vs_T_0_colName, by="DmID",
                        all.y=TRUE)

merged_DEGs_tabs_1<-merge(Mt_1_vs_Mck_16_colName, Mt_1_vs_T_0_colName, by="DmID",
                        all.y=TRUE)

merged_DEGs_tabs<-merge(merged_DEGs_tabs, merged_DEGs_tabs_1, by="DmID",
                          all.y=TRUE)
View(merged_DEGs_tabs)

#we extract the log2FC and padjusted
merged_DEGs_tabs_log2FC_padj<-select(merged_DEGs_tabs, DmID, log2FoldChange_Mt_16_vs_Mck_16, 
                                     padj_Mt_16_vs_Mck_16, log2FoldChange_Mt_16_vs_T_0,
                                     padj_Mt_16_vs_T_0,log2FoldChange_Mt_1_vs_Mck_16,
                                     padj_Mt_1_vs_Mck_16,log2FoldChange_Mt_1_vs_T_0,
                                     padj_Mt_1_vs_T_0)
head(merged_DEGs_tabs_log2FC_padj)



#############---------------ADDING ADDNOTATIONS ---------------######################


View(adno_GeneName)
View(merged_DEGs_tabs_log2FC_padj)

#renaming 
adno_GeneName <- adno_GeneName %>%
  rename(DmID = GenID)
View(adno_GeneName)

#merge the big table with the gene Name mercator
merged_DEGs_tabs_anno<-merge(adno_GeneName,merged_DEGs_tabs, by='DmID', all.y = T)
view(merged_DEGs_tabs_anno)

#merge the small table with the mercator 
merged_DEGs_tabs_sel_anno<-merge(adno_GeneName, merged_DEGs_tabs_log2FC_padj,
                                 by='DmID', all.y=T)
view(merged_DEGs_tabs_sel_anno)
# we can export this and inspect it in excel




#######################################################################
# --------------- DEGS HEATMAP--------------------------------------

# we have the tables with all genes -> now we have to extract the one which we consider STATISTICALLY SIGNIFICANT
# the thresholds could be: 
# p value < 0.05
# or stricter: p-adj < 0.05
# even more stric: p-adj < 0.01 or p-adj < 0.001
# we can add threshold at the log2FC level: log2FC > 1 or log2FC > 2 (to get upregualated genes after application of treatment in "treatment vs control" comparison)
# we can add threshold at the log2FC level: log2FC < -1 or log2FC < -2 (to get downregulated genes after application of treatment in "treatment vs control" comparison)
# we can add threshold at baseMean level to explude genes with small counts (eg: baseMean > 10)

# do we want the filteres applied at the same time or we need one or the other: try | , &
# in this case we apply filter for log2FC and padj both at the same time, so we use & sign:



########-------------------------------------------


#---------------------------------
#	 DEGS----------------  MetJA_16h vs  Mock_16h 
#----------------------------------

#Upregulated Mt_16_vs_Mck_16

Up_DEGs_Mt_16_vs_Mck_16<-filter(Mt_16_vs_Mck_16_colName, padj_Mt_16_vs_Mck_16<0.05, log2FoldChange_Mt_16_vs_Mck_16>1)
head(Up_DEGs_Mt_16_vs_Mck_16)
View(Up_DEGs_Mt_16_vs_Mck_16)

#Downregulated Mt_16_vs_Mck_16
Down_DEGs_Mt_16_vs_Mck_16<-filter(Mt_16_vs_Mck_16_colName, padj_Mt_16_vs_Mck_16<0.05, log2FoldChange_Mt_16_vs_Mck_16< -1)
head(Down_DEGs_Mt_16_vs_Mck_16)
View(Down_DEGs_Mt_16_vs_Mck_16)

## sometimes we need them together so we merge them 
Degs_Mt_16_vs_Mck_16<-rbind(Up_DEGs_Mt_16_vs_Mck_16, Down_DEGs_Mt_16_vs_Mck_16)
View(Degs_Mt_16_vs_Mck_16)


##########---------------------------------------



#------------------------------------------------------------------------
# DEGS-----------------------	MetJA_16h vs  T_0 
#------------------------------------------------------------------------

#Upregulated Mt_16_vs_T0

Up_DEGs_Mt_16_vs_T_0<-filter(Mt_16_vs_T_0_colName, padj_Mt_16_vs_T_0<0.05, log2FoldChange_Mt_16_vs_T_0>1)
head(Up_DEGs_Mt_16_vs_T_0)
View(Up_DEGs_Mt_16_vs_T_0)

#Downregulated Mt_16_vs_T0
Down_DEGs_Mt_16_vs_T_0<-filter(Mt_16_vs_T_0_colName, padj_Mt_16_vs_T_0<0.05, log2FoldChange_Mt_16_vs_T_0< -1)
head(Down_DEGs_Mt_16_vs_T_0)
View(Down_DEGs_Mt_16_vs_T_0)

## sometimes we need them together so we merge them 
Degs_Mt_16_vs_T_0<-rbind(Up_DEGs_Mt_16_vs_T_0, Down_DEGs_Mt_16_vs_T_0)
View(Degs_Mt_16_vs_T_0)

##########---------------------------------------



#------------------------------------------------------------------------
# DEGS-----------------------	o	MetJA_1h vs  Mock_16h  
#------------------------------------------------------------------------

#Upregulated Mt_1_vs_M_16

Up_DEGs_Mt_1_vs_Mck_16<-filter(Mt_1_vs_Mck_16_colName, padj_Mt_1_vs_Mck_16<0.05, log2FoldChange_Mt_1_vs_Mck_16>1)
head(Up_DEGs_Mt_1_vs_Mck_16)
View(Up_DEGs_Mt_1_vs_Mck_16)

#Downregulated Mt_16_vs_T0
Down_DEGs_Mt_1_vs_Mck_16<-filter(Mt_1_vs_Mck_16_colName, padj_Mt_1_vs_Mck_16<0.05, log2FoldChange_Mt_1_vs_Mck_16< -1)
head(Down_DEGs_Mt_1_vs_Mck_16)
View(Down_DEGs_Mt_1_vs_Mck_16)

## sometimes we need them together so we merge them 
Degs_Mt_1_vs_Mck_16<-rbind(Up_DEGs_Mt_1_vs_Mck_16, Down_DEGs_Mt_1_vs_Mck_16)
View(Degs_Mt_1_vs_Mck_16)


#------------------------------------------------------------------------
# DEGS-----------------------	o	MetJA_1h vs T0   
#------------------------------------------------------------------------

#Upregulated Mt_1_vs_T_0

Up_DEGs_Mt_1_vs_T_0<-filter(Mt_1_vs_T_0_colName, padj_Mt_1_vs_T_0<0.05, log2FoldChange_Mt_1_vs_T_0>1)
head(Up_DEGs_Mt_1_vs_T_0)
View(Up_DEGs_Mt_1_vs_T_0)

#Downregulated Mt_16_vs_T0
Down_DEGs_Mt_1_vs_T_0<-filter(Mt_1_vs_T_0_colName, padj_Mt_1_vs_T_0<0.05, log2FoldChange_Mt_1_vs_T_0< -1)
head(Down_DEGs_Mt_1_vs_T_0)
View(Down_DEGs_Mt_1_vs_T_0)

## sometimes we need them together so we merge them 
Degs_Mt_1_vs_T_0<-rbind(Up_DEGs_Mt_1_vs_T_0, Down_DEGs_Mt_1_vs_T_0)
View(Degs_Mt_1_vs_T_0)




#####--------------- MERGING TO GET ONLY LOGFOLD FOR HEATMAP -----------------############
#Degs_Mt_16_vs_Mck_16
#Degs_Mt_16_vs_T_0
#Degs_Mt_1_vs_Mck_16
#Degs_Mt_1_vs_T_0

Degs_Mt_16_vs_T_0_Degs_Mt_16_vs_T_0 <-merge(Degs_Mt_16_vs_Mck_16, Degs_Mt_16_vs_T_0, by='DmID')
Degs_Mt_1_vs_Mck_16_Degs_Mt_1_vs_T_0 <-merge(Degs_Mt_1_vs_Mck_16, Degs_Mt_1_vs_T_0, by='DmID')
Degs<-merge(Degs_Mt_16_vs_T_0_Degs_Mt_16_vs_T_0, Degs_Mt_1_vs_Mck_16_Degs_Mt_1_vs_T_0, by='DmID')

colnames(Degs)

#### now only logfold

# we only select the log2FC 
Degs_matrix<-select(Degs, DmID, log2FoldChange_Mt_16_vs_Mck_16,
                    log2FoldChange_Mt_16_vs_T_0,log2FoldChange_Mt_1_vs_Mck_16,log2FoldChange_Mt_1_vs_T_0)
head(Degs_matrix)

## to use it for a heatmap we put column as rowname 
Degs_matrix<-column_to_rownames(Degs_matrix, var='DmID')

#plot the heatmaps 
View(Degs_matrix)
plot1 <- pheatmap( Degs_matrix,
                   scale = "column",
                   color = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
                   cluster_cols = F,
                   cluster_rows = T,
                   legend = T,
                   border_color = "NA",
                   labels_col = c("Mt_16_vs_Mck_16", "Mt_16_vs_T_0", 'Mt_1_vs_Mck_16', 'Mt_1_vs_T_0'), # the column names are very long, we rename them here with labels_col -> make sure the order is correct
                   show_rownames = F)
plot1






########################## ------------- Venn Diagram________________________________
########_________________________upreg_MT_16_vs_Mck_16_T_0 --------------------------

upreg_MT_16_vs_Mck_16_T_0 <-list(Mt_16_vs_Mck_16= Up_DEGs_Mt_16_vs_Mck_16$DmID,
             Mt_16_vs_T_0=Up_DEGs_Mt_16_vs_T_0$DmID)
ggVennDiagram(upreg_MT_16_vs_Mck_16_T_0)+scale_fill_gradient(low="grey90", high='red')+
  scale_x_continuous(expand = expansion( mult =0.2))




###########---------- downreg_MT_16_VS_MCK_16_T_0__________________________

downreg_MT_16_vs_Mck_16_T_0 <-list(Mt_16_vs_Mck_16= Down_DEGs_Mt_16_vs_Mck_16$DmID,
                                 Mt_16_vs_T_0=Down_DEGs_Mt_16_vs_T_0$DmID)
ggVennDiagram(downreg_MT_16_vs_Mck_16_T_0)+scale_fill_gradient(low="grey90", high='blue')+
  scale_x_continuous(expand = expansion( mult =0.2))


########### ------- Advanced Venn Diagram UP-REGULATED  -------

# we use a matrix where 1 suggests pressence of DEG s and 0 non-DEGs

##we merge all the upregulated genes
upreg_MT_16_vs_Mck_16_T_0_2<-merge(Up_DEGs_Mt_16_vs_Mck_16,Up_DEGs_Mt_16_vs_T_0, by='DmID', all=T)

#we choose the logfold2
upreg_MT_16_vs_Mck_16_T_0_2_MATRIX<-select(upreg_MT_16_vs_Mck_16_T_0_2, DmID, log2FoldChange_Mt_16_vs_Mck_16,
                           log2FoldChange_Mt_16_vs_T_0)

head(upreg_MT_16_vs_Mck_16_T_0_2_MATRIX)
# for venn we should not have any characters so we transform dmid into rownames
upreg_MT_16_vs_Mck_16_T_0_2_MATRIX<-column_to_rownames(upreg_MT_16_vs_Mck_16_T_0_2_MATRIX, var='DmID')
#renaming because of awfully long name I gave
venn_matrix<-upreg_MT_16_vs_Mck_16_T_0_2_MATRIX
#because we cannot use NA in a matrix we transform it into 0
venn_matrix[is.na(venn_matrix)]<-0
#all the other values will be 1 
venn_matrix[venn_matrix!=0]<-1
head(venn_matrix)

#draw the pairwise venn
graphics.off()
draw.pairwise.venn(area1=(nrow(subset(venn_matrix, log2FoldChange_Mt_16_vs_Mck_16==1))),
                   area2 = (nrow(subset(venn_matrix, log2FoldChange_Mt_16_vs_T_0==1))),
                   cross.area = nrow(subset(venn_matrix,log2FoldChange_Mt_16_vs_Mck_16==1 & log2FoldChange_Mt_16_vs_T_0==1)),
                   category= c("Mt_16_vs_T0 ", "Mt_16_vs_Mck_16"),
                   lty="blank",
                   fill=c("dodgerblue", "red"))


########### ------- Advanced Venn Diagram DOWN-REGULATED  -------

# we use a matrix where 1 suggests pressence of DEG s and 0 non-DEGs

##we merge all the upregulated genes
downreg_MT_16_vs_Mck_16_T_0_2<-merge(Down_DEGs_Mt_16_vs_Mck_16,Down_DEGs_Mt_16_vs_T_0, by='DmID', all=T)

#we choose the logfold2
downreg_MT_16_vs_Mck_16_T_0_2_MATRIX<-select(downreg_MT_16_vs_Mck_16_T_0_2, DmID, log2FoldChange_Mt_16_vs_Mck_16,
                                           log2FoldChange_Mt_16_vs_T_0)

head(downreg_MT_16_vs_Mck_16_T_0_2_MATRIX)
# for venn we should not have any characters so we transform dmid into rownames
downreg_MT_16_vs_Mck_16_T_0_2_MATRIX<-column_to_rownames(downreg_MT_16_vs_Mck_16_T_0_2_MATRIX, var='DmID')
#renaming because of awfully long name I gave
venn_matrix_down<-downreg_MT_16_vs_Mck_16_T_0_2_MATRIX
#because we cannot use NA in a matrix we transform it into 0
venn_matrix_down[is.na(venn_matrix_down)]<-0
#all the other values will be 1 
venn_matrix_down[venn_matrix_down!=0]<-1
head(venn_matrix_down)

#draw the pairwise venn
graphics.off()
draw.pairwise.venn(area1=(nrow(subset(venn_matrix_down, log2FoldChange_Mt_16_vs_Mck_16==1))),
                   area2 = (nrow(subset(venn_matrix_down, log2FoldChange_Mt_16_vs_T_0==1))),
                   cross.area = nrow(subset(venn_matrix_down,log2FoldChange_Mt_16_vs_Mck_16==1 & log2FoldChange_Mt_16_vs_T_0==1)),
                   category= c("Mt_16_vs_Mck_16 ", "Mt_16_vs_T_0"),
                   lty="blank",
                   fill=c("dodgerblue", "red"))






















########################## ------------- Venn Diagram________________________________
########_________________________Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0 --------------------------




upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0 <-list(Mt_1_vs_Mck_16= Up_DEGs_Mt_1_vs_Mck_16$DmID,
                                           Mt_1_vs_T_0=Up_DEGs_Mt_1_vs_T_0$DmID)
ggVennDiagram(upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0)+scale_fill_gradient(low="grey90", high='red')+
  scale_x_continuous(expand = expansion( mult =0.2))




downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0 <-list(Mt_1_vs_Mck_16= Down_DEGs_Mt_1_vs_Mck_16$DmID,
                                             Mt_1_vs_T_0=Down_DEGs_Mt_1_vs_T_0$DmID)
ggVennDiagram(downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0)+scale_fill_gradient(low="grey90", high='blue')+
  scale_x_continuous(expand = expansion( mult =0.2))





########### ------------------------------- Advanced Venn Diagram UP-REGULATED  -------

# we use a matrix where 1 suggests pressence of DEG s and 0 non-DEGs
##we merge all the upregulated genes

upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0<-merge(Up_DEGs_Mt_1_vs_Mck_16,Up_DEGs_Mt_1_vs_T_0, by='DmID', all=T)

#we choose the logfold2
upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX<-select(upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0, DmID, log2FoldChange_Mt_1_vs_Mck_16,
                                           log2FoldChange_Mt_1_vs_T_0)
head(upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX)

# for venn we should not have any characters so we transform dmid into rownames
upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX<-column_to_rownames(upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX, var='DmID')

#renaming because of awfully long name I gave
venn_matrix<-upreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX

#because we cannot use NA in a matrix we transform it into 0
venn_matrix[is.na(venn_matrix)]<-0
#all the other values will be 1 
venn_matrix[venn_matrix!=0]<-1
head(venn_matrix)

#draw the pairwise venn
graphics.off()
draw.pairwise.venn(area1=(nrow(subset(venn_matrix, log2FoldChange_Mt_1_vs_Mck_16==1))),
                   area2 = (nrow(subset(venn_matrix, log2FoldChange_Mt_1_vs_T_0==1))),
                   cross.area = nrow(subset(venn_matrix,log2FoldChange_Mt_1_vs_Mck_16==1 & log2FoldChange_Mt_1_vs_T_0==1)),
                   category= c("Mt_1_vs_Mck_16  ", "Mt_1_vs_T_0"),
                   lty="blank",
                   fill=c("dodgerblue", "red"))


########### ------- Advanced Venn Diagram DOWN-REGULATED  -------

# we use a matrix where 1 suggests pressence of DEG s and 0 non-DEGs

##we merge all the upregulated genes
downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2<-merge(Down_DEGs_Mt_1_vs_Mck_16,Down_DEGs_Mt_1_vs_T_0, by='DmID', all=T)

#we choose the logfold2
downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX<-select(downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2, DmID, log2FoldChange_Mt_1_vs_Mck_16,
                                             log2FoldChange_Mt_1_vs_T_0)

head(downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX)
# for venn we should not have any characters so we transform dmid into rownames
downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX<-column_to_rownames(downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX, var='DmID')
#renaming because of awfully long name I gave
venn_matrix_down<-downreg_Mt_1_vs_Mck_16_vs_Mt_1_vs_T_0_2_MATRIX
#because we cannot use NA in a matrix we transform it into 0
venn_matrix_down[is.na(venn_matrix_down)]<-0
#all the other values will be 1 
venn_matrix_down[venn_matrix_down!=0]<-1
head(venn_matrix_down)

#draw the pairwise venn
graphics.off()
draw.pairwise.venn(area1=(nrow(subset(venn_matrix_down, log2FoldChange_Mt_1_vs_Mck_16==1))),
                   area2 = (nrow(subset(venn_matrix_down, log2FoldChange_Mt_1_vs_T_0==1))),
                   cross.area = nrow(subset(venn_matrix_down,log2FoldChange_Mt_1_vs_Mck_16==1 & log2FoldChange_Mt_1_vs_T_0==1)),
                   category= c("Mt_1_vs_Mck_16 ", "Mt_1_vs_T_0"),
                   lty="blank",
                   fill=c("dodgerblue", "red"))








########################------------------------------ COMBINED VENN- UPREGULATED -------------------------##############################



# Assuming you have four sets of genes
set1 <- Up_DEGs_Mt_16_vs_Mck_16$DmID
set2 <- Up_DEGs_Mt_16_vs_T_0$DmID
set3 <- Up_DEGs_Mt_1_vs_Mck_16$DmID
set4 <- Up_DEGs_Mt_1_vs_T_0$DmID

# Combine them into a list
venn_data <- list(
  Mt_16_vs_Mck_16 = set1,
  Mt_16_vs_T_0 = set2,
  Mt_1_vs_Mck_16 = set3,
  Mt_1_vs_T_0 = set4
)

# Generate a Quatro Venn Diagram
library(ggVennDiagram)

ggVennDiagram(venn_data, label_alpha = 0.5) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Quatro Venn Diagram for DEGs") +
  theme_minimal()





########################------------------------------ COMBINED VENN- DOWNREGULATED -------------------------##############################
#Mt_16_vs_Mck_16
#Mt_16_vs_T_0
#Mt_1_vs_Mck_16
#Mt_1_vs_T_0


# Assuming you have four sets of genes
set1 <- Down_DEGs_Mt_16_vs_Mck_16$DmID
set2 <- Down_DEGs_Mt_16_vs_T_0$DmID
set3 <- Down_DEGs_Mt_1_vs_Mck_16$DmID
set4 <- Down_DEGs_Mt_1_vs_T_0$DmID
set5<- Down_DEGs

# Combine them into a list
venn_data <- list(
  Mt_16_vs_Mck_16 = set1,
  Mt_16_vs_T_0 = set2,
  Mt_1_vs_Mck_16 = set3,
  Mt_1_vs_T_0 = set4
)

# Generate a Quatro Venn Diagram
library(ggVennDiagram)

graphics.off()
ggVennDiagram(venn_data, label_alpha = 0.5) +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Quatro Venn Diagram for DEGs") +
  theme_minimal()







############----------------------- UPSET PLOT - UPREGULATED ----------------######################
upreg_all <- full_join(Up_DEGs_Mt_16_vs_Mck_16, Up_DEGs_Mt_16_vs_T_0, by = "DmID") %>%
  full_join(Up_DEGs_Mt_1_vs_Mck_16, by = "DmID") %>%
  full_join(Up_DEGs_Mt_1_vs_T_0, by = "DmID")

upreg_matrix <- select(upreg_all, DmID, 
                       log2FoldChange_Mt_16_vs_Mck_16, 
                       log2FoldChange_Mt_16_vs_T_0, 
                       log2FoldChange_Mt_1_vs_Mck_16, 
                       log2FoldChange_Mt_1_vs_T_0)

upreg_matrix <- column_to_rownames(upreg_matrix, var = "DmID")
upreg_matrix[is.na(upreg_matrix)] <- 0
upreg_matrix[upreg_matrix != 0] <- 1
colnames(upreg_matrix) <- gsub("log2FoldChange_", "upreg_", colnames(upreg_matrix))
view(upreg_matrix)


# Create the UpSet 
graphics.off()
upset(
  upreg_matrix,
  sets = c("upreg_Mt_16_vs_Mck_16", "upreg_Mt_16_vs_T_0", "upreg_Mt_1_vs_Mck_16", "upreg_Mt_1_vs_T_0"),
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "darkred",
  text.scale = 1.5
)








############----------------------- UPSET PLOT - DOWNREGULATED ----------------######################


downreg_all <- full_join(Down_DEGs_Mt_16_vs_Mck_16, Down_DEGs_Mt_16_vs_T_0, by = "DmID") %>%
  full_join(Down_DEGs_Mt_1_vs_Mck_16, by = "DmID") %>%
  full_join(Down_DEGs_Mt_1_vs_T_0, by = "DmID")
downreg_matrix <- select(downreg_all, DmID, 
                         log2FoldChange_Mt_16_vs_Mck_16, 
                         log2FoldChange_Mt_16_vs_T_0, 
                         log2FoldChange_Mt_1_vs_Mck_16, 
                         log2FoldChange_Mt_1_vs_T_0)
downreg_matrix <- column_to_rownames(downreg_matrix, var = "DmID")
downreg_matrix[is.na(downreg_matrix)] <- 0
downreg_matrix[downreg_matrix != 0] <- 1
colnames(downreg_matrix) <- gsub("log2FoldChange_", "downreg_", colnames(downreg_matrix))
view(downreg_matrix)


# Create the UpSet 
graphics.off()
upset(
  downreg_matrix,
  sets = c("downreg_Mt_16_vs_Mck_16", "downreg_Mt_16_vs_T_0", "downreg_Mt_1_vs_Mck_16", "downreg_Mt_1_vs_T_0"),
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "darkred",
  text.scale = 1.5
)





##### --------------- GIGA_TABLE UPREGULATED---------------------############
#### Creating gigantic table FOR UPREG


View(as_data_frame(Mt_1_vs_T_0_colName))
head(Mt_16_vs_T_0_colName)
head(Mt_1_vs_Mck_16_colName)
head(Mt_1_vs_T_0_colName)

# Merge all tables by DmID
mega_table<- full_join(Mt_16_vs_Mck_16_colName, Mt_16_vs_T_0_colName, by = "DmID") %>%
  full_join(Mt_1_vs_Mck_16_colName, by = "DmID") %>%
  full_join(Mt_1_vs_T_0_colName, by = "DmID")

View(mega_table)
# Add matrix 
upreg_matrix <- rownames_to_column(upreg_matrix, var = "DmID")
#view(upreg_matrix)
mega_table_upreg <- merge(mega_table, upreg_matrix , by='DmID')
#view(mega_table_upreg)

# Check 
head(mega_table_upreg)

# Save the mega-table to a CSV file
write.csv(mega_table_upreg, "mega_table_upreg.csv", row.names = FALSE)


##### --------------- GIGA_TABLE DOWNREGULATED---------------------############
#### Creating gigantic table FOR UPREG


View(as_data_frame(Mt_1_vs_T_0_colName))
head(Mt_16_vs_T_0_colName)
head(Mt_1_vs_Mck_16_colName)
head(Mt_1_vs_T_0_colName)

downreg_matrix <- rownames_to_column(downreg_matrix, var = "DmID")
#view(upreg_matrix)
mega_table_downreg<- merge(mega_table, downreg_matrix , by='DmID')
#view(mega_table_upreg)

# Check 
head(mega_table_downreg)

# Save the mega-table to a CSV file
write.csv(mega_table_downreg, "mega_table_downreg.csv", row.names = FALSE)








############--------------------- ADD MERCATOR TO MEGA TABLE ---------------#################

mercator_results<-read.delim('D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/Arabidopsis.results.txt',
                             header = T, sep ="\t")

head(mercator_results)
str(mercator_results)
View(mercator_results)


colnames(mercator_results)[colnames(mercator_results)== "IDENTIFIER"] <- "DmID"
mercator_results$DmID <- toupper(mercator_results$DmID)
mercator_results$DmID <- gsub("\\..*", "", mercator_results$DmID)
mercator_results$DmID <- gsub("^'", "", mercator_results$DmID) 
mercator_results$DmID <- gsub("^'|'$", "", mercator_results$DmID)
mercator_results$DmID <- gsub("'", "", mercator_results$DmID)
mercator_results$DESCRIPTION <- gsub(".*description:", "", mercator_results$DESCRIPTION)  # Remove everything before "description:"
#mercator_results$DESCRIPTION <- gsub("\\[.*", "", mercator_results$DESCRIPTION_cleaned)  # Remove everything after "["

head(mercator_results$DESCRIPTION )

head(mercator_results)

write.csv(mercator_results, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mercator_results_modified.csv", row.names = FALSE)
mercator_results=read.csv("mercator_results_modified.csv")
mega_table_upreg=read.csv("mega_table_upreg.csv")

View(mega_table_upreg)
View(mercator_results)

mega_table_upreg<-merge(mega_table_upreg, mercator_results, by= "DmID", all.x= T)
head(mega_table_upreg)
View(mega_table_upreg)


# Remove the DESCRIPTION_cleaned and NAME columns
mega_table_upreg <- mega_table_upreg[, !(colnames(mega_table_upreg) %in% c("TYPE"))]
# Verify the columns have been removed
colnames(mega_table_upreg)
View(mega_table_upreg)

write.csv(mega_table_upreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_upreg.csv", row.names = FALSE)




#########---------- NOW LETS MAKE THE DOWNREG GIGA TABLE ---------#########

mega_table_downreg=read.csv("mega_table_downreg.csv")

View(mega_table_downreg)
View(mercator_results)

mega_table_downreg<-merge(mega_table_downreg, mercator_results, by= "DmID", all.x= T)
head(mega_table_downreg)
View(mega_table_downreg)


# Remove the DESCRIPTION_cleaned and NAME columns
mega_table_downreg <- mega_table_downreg[, !(colnames(mega_table_downreg) %in% c("TYPE"))]
# Verify the columns have been removed
colnames(mega_table_downreg)


write.csv(mega_table_downreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_downreg.csv", row.names = FALSE)










####################-------------- ADDING TF  TO THE MEGA TABLE ------------

tfdb_data<-read.csv("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/Ath_TF_list.txt", sep='\t')
head(tfdb_data)

colnames(tfdb_data)[colnames(tfdb_data)== "Gene_ID"] <- "DmID"

tfdb_data <- tfdb_data[, !(colnames(tfdb_data) %in% c("TF_ID"))]
head(tfdb_data)

mega_table_upreg<-read.csv("mega_table_mercator_upreg.csv")
mega_table_upreg<-merge(mega_table_upreg, tfdb_data, by= "DmID", all.x= T)
head(mega_table_upreg)
View(mega_table_upreg)

write.csv(mega_table_upreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_upreg.csv", row.names = FALSE)
mega_table_upreg=read.csv("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_upreg.csv")
View(mega_table_upreg)

########--------------------------Same for downreg-----------------------#########


mega_table_downreg<-merge(mega_table_downreg, tfdb_data, by= "DmID", all.x= T)
head(mega_table_downreg)
View(mega_table_downreg)

write.csv(mega_table_downreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_downreg.csv", row.names = FALSE)
mega_table_downreg=read.csv("D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_downreg.csv")
View(mega_table_downreg)




# DELETING DUPLICATES APPEARING BECAUSE MERCATOR HAS .X AFTER ITS DMID!!!!!!!!

mega_table_upreg <- read.csv("mega_table_mercator_upreg.csv")
# Retain the first occurrence of each unique DmID
mega_table_upreg <- mega_table_upreg %>%
  arrange(DmID) %>%  # Optional: Sort by DmID for consistency
  distinct(DmID, .keep_all = TRUE)
View(mega_table_upreg)
write.csv(mega_table_upreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_upreg.csv", row.names = FALSE)
mega_table_upreg <- read.csv("mega_table_mercator_upreg.csv")


mega_table_downreg <- read.csv("mega_table_mercator_downreg.csv")
# Retain the first occurrence of each unique DmID
mega_table_downreg <- mega_table_downreg %>%
  arrange(DmID) %>%  # Optional: Sort by DmID for consistency
  distinct(DmID, .keep_all = TRUE)
View(mega_table_downreg)
write.csv(mega_table_downreg, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/mega_table_mercator_downreg.csv", row.names = FALSE)
mega_table_downreg <- read.csv("mega_table_mercator_downreg.csv")





######################---------------- EXTRACT FROM UPSET MT_16 T1 ----------------------#################
mega_table_upreg <- read.csv("mega_table_mercator_upreg.csv")
View(mega_table_upreg)
# Genes different from T0 in MetJa
diff_metja <- mega_table_upreg %>%
  filter(upreg_Mt_16_vs_T_0 == 1) %>%  # Significantly upregulated in MetJa vs T0
  select(DmID)

cat("Total genes different in MetJa (diff_metja):", nrow(diff_metja), "\n")

unique_genes <- mega_table_upreg %>%
  filter(upreg_Mt_16_vs_T_0 == 1 & upreg_Mt_16_vs_Mck_16 == 1 &
           upreg_Mt_1_vs_Mck_16 == 0 & upreg_Mt_1_vs_T_0 == 0) %>%
  select(DmID)


# Debugging: Check the count of genes different in MetJa
cat("Total genes different in MetJa (unique_genes):", nrow(unique_genes), "\n")
write.csv(unique_genes, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/very_unique.csv", row.names = FALSE)




target_genes_metja_specific <- mega_table_upreg %>%
  filter(upreg_Mt_16_vs_T_0 == 1 & upreg_Mt_16_vs_Mck_16 == 0) %>%
  select(DmID)

cat("Total MetJa16-specific upregulated genes compared to T0:", nrow(target_genes_metja_specific), "\n")



# Filter strictly for MetJa16 vs T0 and MetJa16 vs Mock
metja16_vs_t0 <- mega_table_upreg %>%
  filter(upreg_Mt_16_vs_T_0 == 1) %>%
  pull(DmID)

metja16_vs_mock <- mega_table_upreg %>%
  filter(upreg_Mt_16_vs_Mck_16 == 1) %>%
  pull(DmID)

# Create a list with filtered sets
upset_data <- list(
  "MetJa16 vs T0" = metja16_vs_t0,
  "MetJa16 vs Mock" = metja16_vs_mock
)

# Generate the UpSet plot
upset(fromList(upset_data),
      sets = c("MetJa16 vs T0", "MetJa16 vs Mock"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkred",
      text.scale = 1.5)

View(target_genes_metja_specific)

# Save the 2137 MetJa16-specific gene IDs to a file
write.csv(target_genes_metja_specific, "D:/Faculty/Master_2/Plant_Genomics/Project_Pl_Genom/Plant_Genomics_Project/data/target_genes_metja_specific.csv", row.names = FALSE)








##########------------TOP TRANSCRIPTION FACTORS  ----------------------#################



# Load the mega table and 2137 genes
mega_table <- read.csv("mega_table_mercator_upreg.csv")
metja16_genes <- read.csv("target_genes_metja_specific.csv", header = FALSE)
colnames(metja16_genes) <- c("DmID")  # Ensure the column name matches

# Subset the mega table for the 2137 specific genes
mega_table_filtered <- mega_table %>%
  filter(DmID %in% metja16_genes$DmID)

View(mega_table_filtered)

# Save the filtered table for further analysis
write.csv(mega_table_filtered, "mega_table_filtered_TF.csv", row.names = FALSE)

count_data<-table(mega_table_filtered$Family)

ggplot() +
  geom_bar(
    aes(x = names(count_data), y = as.vector(count_data), fill = names(count_data)), 
    stat = 'identity'
  ) +
  xlab('Transcription Factors') + 
  ylab('Count') +
  ggtitle("Barplot of Transcription Factors") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, option = "D") +  # Use discrete viridis palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





