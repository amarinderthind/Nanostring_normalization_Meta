### This scripts is Merge/Integrate various Packages and functions for NANostring Normalization using RUVg (Where no replicates are available b/w different batches) and RUVIII (Coming Soon)

setwd('/Data_drectory_of_RCC_files')
datadir <- "/Users/Data_drectory_of_RCC_files"

##################################
source('functions.r')

## Quality control of .rcc files
library(NanoStringNCTools) ## Bioconductor packages for Nanostring QC
library(ggthemes)
library(ggiraph)
library(PCAtools)

files.RCC <- dir(datadir, pattern = "*\\.RCC$", full.names = TRUE)
files.RCC
sample_annotation <- file.path(datadir, "meta_info.csv")  ##Conditions/treatment/normal are defined in the coulmn "Group"
sample_annotation
demoData <- readNanoStringRccSet(files.RCC, rlfFile = NULL, phenoDataFile = sample_annotation, phenoDataRccColName = "RCC")

design( demoData ) <- ~ `Group`
design( demoData )

dimLabels( demoData )
protocolData(demoData)[["Sample ID"]] <- sampleNames(demoData)
dimLabels( demoData )[2] <- "Sample ID"
dimLabels( demoData )

## how many condition in total 
unique( sData( demoData )$"Group" )
with( controlSubset( demoData ) , table( CodeClass ) )
with( nonControlSubset( demoData ) , table( CodeClass ) )


## extracting the expression data
exprs_df <- transform( assayData( demoData )[["exprs"]] )

## QC plots 
girafe( ggobj = autoplot( demoData , "lane-bindingDensity" ) )
girafe( ggobj = autoplot( demoData , "lane-fov" ) )
girafe( ggobj = autoplot( demoData, "housekeep-geom" ) ) 
girafe( ggobj = autoplot( demoData , "ercc-linearity" ) )
girafe( ggobj = autoplot( demoData , "ercc-lod" ) )

girafe( ggobj = autoplot( demoData , "boxplot-feature" , index = featureNames(demoData)[3] , elt = "exprs" ) )

##normalize the expression values using nSOLVER default and save to exprs_norm

demoData <- normalize( demoData , type="nSolver", fromELT = "exprs" , toELT = "exprs_norm" )
autoplot( demoData , "heatmap-genes" , elt = "exprs_norm" )

#######################################################
require(EnvStats)

raw_expression <- as.data.frame(demoData@assayData[["exprs"]])
pData <- as.data.frame(demoData@phenoData@data)
fData <- as.data.frame(demoData@featureData@data)
protocolData <- as.data.frame(demoData@protocolData@data)


##### RUVg BASED analysis
## https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa163/5891144
## removed this step because RUVg doesnot consider the technical replicates and after removal of calibrator, not found any gene associated with biology of interest.

cIdx <- fData$GeneName[fData$CodeClass == "Housekeeping"]
pData$HK_Gene_Miss = colSums(raw_expression[cIdx,] == 0)
rownames(fData) = fData$Gene
rownames(raw_expression) = fData$Gene
rownames(pData) = colnames(raw_expression)


#### CHECK IF HK ARE ASSOCIATED WITH PRIMARY PHENOTYPE
## exclude hk genes which are associated with phenotype

hk_raw = raw_expression[cIdx,]
pval = vector(length = nrow(hk_raw))

require(MASS)

for (i in 1:nrow(hk_raw)){
  
  reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(pData$Group)) # negative binomial 
  pval[i] = coef(summary(reg))[2,4]
}

sum(pval <= .05)
summary(reg)
idx <- pval <= .05
#View(hk_raw[idx,])
exc <- row.names(hk_raw[idx, ]) ## save these genes into exc and pass to RUV_total 

################

######### RLE plot based on raw_expression

data <- as.matrix(raw_expression)
data <- log2(data+1)

boxplot(data,las=2)

library(plotrix)
RLEplot_mod(data,pData$Run) ### function is defined below


### grid of datasets
library(RUVSeq)
library(DESeq2)
library(limma)
library(matrixStats)


k = 3  ## Please check the original paper for more details ## iteration of plots using different K values are recommended

vsd = RUV_total(raw_expression,pData,fData,k = k, exclude = exc)$vsd ## put exclude if there any hk associated with biology of interest 
set = RUV_total(raw_expression,pData,fData,k = k, exclude = exc )$set # ,exclude = exc


normalizedcount <- set@assayData[["normalizedCounts"]]

write.csv(normalizedcount,"normalized_count_Ruvg_based_.csv")

nordata <- log2(normalizedcount)
RLEplot_mod(nordata,pData$Run) ### function is defined below

p <-pca(nordata, metadata = pData, removeVar = 0.01) ## -- removing the lower 1% of variables based on variance
#biplot(p)

biplot(p,
       #lab =NULL,
       colby= 'Group',#'BATCH.ID',
       #shape = 'Alt_Batch_ID',
       hline = 0, vline = 0,
       legendPosition = 'right', encircle = T, 
       encircleLineSize=4,encircleFill = F,
       title='PCA plot')


####

dds <- DESeqDataSetFromMatrix(countData = counts(set),
                              colData = pData(set),
                              design = ~ Group) #counts(set)[1:652,], ##use design with W_1 or W_2 based on k value

dds <- DESeq(dds)


#temp_count <- as.data.frame(counts(dds,normalized=TRUE))

firstC <- 'CL'
SecondC <- 'MET'
contrast<- c("Group",firstC,SecondC)


#### Valcano Plot #####
plotMA(results(dds,contrast =contrast))
res_deseq2 <- as.data.frame(results(dds,contrast =contrast))

res_valcanoplot <- results(dds,contrast =contrast)

library(EnhancedVolcano)
EnhancedVolcano(res_valcanoplot,
                lab = rownames(res_valcanoplot),
                x = 'log2FoldChange',
                y = 'pvalue')


EnhancedVolcano(res_valcanoplot,
                lab = rownames(res_valcanoplot),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Condition 1 vs Condition 2',
                pCutoff = 0.01,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 6.0)


nam <- paste('up_in',firstC, sep = '_')
res_deseq2[, nam] <- as.logical(res_deseq2$log2FoldChange > 0)

res_deseq2$threshold <- as.logical(res_deseq2$padj < 0.05)  #Threshold defined earlier
row.names(res_deseq2)[which(res_deseq2$threshold)]

filename <- paste0('Deseq_ruvg_k_',k,'_DEA_.csv')

norm_mean <- sapply( levels(dds$Group), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$Group == lvl, drop=F] ) )
colnames(norm_mean) <- paste('Rowmean_exp_',levels(dds$Group),sep='')

res3 <- cbind(norm_mean,res_deseq2)
res3 <- res3[,-3]

write.csv(res3,filename)
