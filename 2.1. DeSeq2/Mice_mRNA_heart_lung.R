#################################################################################################
#################################################################################################
#### Analysis of RNA-Seq data for the Wildfire Project
####
#### Using DESeq2 for sample normalization and statistical analysis
#### For this example: to identify mRNAs associated with various biomass exposures in mice tissue (heart or lung)
####
#### Because this is just a toxicity analysis, code does not include covariate imputation, and just an optional SVA step
#### 
#### 
#### Code drafted by Alexis Payton and Julia Rager
#### Lasted updated: January 7, 2021
#################################################################################################
#################################################################################################


sessionInfo()
rm(list=ls())

#################################################################################################
#################################################################################################
#### Installing appropriate packages in R (if already have these installed, SKIP THIS STEP)
#################################################################################################
#################################################################################################
# install.packages("data.table")

# Note that BiocManager allows you to access the Bioconductor repo in which DeSeq2 package is. DeSeq2 pacakge has a "tibble" "RcppArmadillo" and "rlang" dependency.
# to install RcppArmadillo, you must install a Fortan compiler (I used gfortran-6.1.pkg, https://cran.r-project.org/bin/macosx/tools/).

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# install.packages("tibble")
# install.packages("rlang")
# install.packages("RcppArmadillo")

# Install DESeq and SVA:
# BiocManager::install("DESeq2", version = "3.11")
# BiocManager::install("sva", version ="3.11")
# BiocManager::install("BiocManager")

# Install plot packages
# install.packages("ggbeeswarm") # note that I had to manually install this package using tools --> install packages
# install.packages("gridExtra")
# install.packages("tidyverse") # Install tidyverse packages 
# install.packages("gplots")
# install.packages("RColorBrewer")
# install.packages("scales") # Install scales (needed for tidyverse activation):



#################################################################################################
#################################################################################################
#### Activating appropriate packages in R
#################################################################################################
#################################################################################################

# Activating appropriate libraries
library(data.table)
library(rlang)
library(DESeq2)
library(limma)
library(sva)
library(ggplot2)
library(stats)
library(scales)
library(ggbeeswarm)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(tidyverse)


#################################################################################################
#################################################################################################
#### Set working directory, output folder path, etc
#################################################################################################
#################################################################################################

# Set working directory
setwd('/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/Experiments/2. Wildfire Analysis/2.1. DeSeq2/2.1.1. Mice mRNA/1. Input')
getwd()

# Create an output folder (make sure to make the folder first, and then point to it here)
Output <- "/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/Experiments/2. Wildfire Analysis/2.1. DeSeq2/2.1.1. Mice mRNA/2. Output"
Output_StatResults = "/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/Experiments/2. Wildfire Analysis/2.1. DeSeq2/2.1.1. Mice mRNA/2. Output/1. StatResults"
Output_StatResults_SVA = "/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/Experiments/2. Wildfire Analysis/2.1. DeSeq2/2.1.1. Mice mRNA/2. Output/2. StatResults_w_SVA"
cur_date = "010721"

#################################################################################################
#################################################################################################
#### Loading count data and metadata (sample information data)
#### Organizing and filtering to keep samples that are present all files
#################################################################################################
#################################################################################################


# Read in count data
countdata <- read.csv(file = 'SP0222_gene_counts.csv', check.names = FALSE)
dim(countdata)
# 171 samples, 30146 mRNAs

# visualize this data quickly by viewing top left corner:
countdata[1:3,1:6]

# Check for duplicate mRNA IDs in the countdata 
Dups <- duplicated(countdata[,1])
summary(Dups)
# Not an issue for this wildfire dataset


# Read in metadata (sample information file)
subjectinfo <- read.csv(file = "Sample_Info_112520.csv", check.names = FALSE)
dim(subjectinfo) #170 rows, 9 col

# Visualize this data quickly by viewing top left corner:
subjectinfo[1:3,1:7]

#reading in table of contrasts 
contrasts = read_csv("Table_of_Contrasts_121120.csv")
head(contrasts)

#################################################################################################
#################################################################################################
#### Create dataframes that are formatted for proceeding code, as well as DESeq2 functions
#### countdata and coldata
#################################################################################################
#################################################################################################

# Make the gene variable the rowname 
countdata <- countdata %>% 
  column_to_rownames(var="Gene") 


# Creating the coldata object, based on information in the subjectinfo file
coldata <- subjectinfo


# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$SampleID_BioSpyderCountFile))
#replacing the biospyder ids in the countdata file with the ids
colnames(countdata) <- coldata$ID

# Double checking that the same variables appear between the two dataframes
setequal(as.character(coldata$ID), colnames(countdata))

# Additionally checking that not only the sets of variables are the same, but that they are in the same order
identical(as.character(coldata$ID), colnames(countdata))


#################################################################################################
#################################################################################################
#### RNASeq QA/QC on raw count data to identify potential outlier samples
#### This may or may not result in another filter step
#################################################################################################
#################################################################################################
# One way to evaluate sequencing data quality / identify potential outliers is through Principal Component Analysis (PCA)
# PCA helps in identifying outlying samples for quality control, and gives a feeling for the principal causes of variation in a dataset


# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))
screeplot(pca)
# this scree plot indicates that nearly all variation is explained in PC1,2 and 3, 
# thus PCA is a useful exercise for this dataset 


# Make dataframe for PCA plot generation using first three components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], Sample=colnames(countdata))
# Add attributes (covariates from coldata) to pca df
pca_df <- merge(pca_df, coldata, by.x="Sample", by.y="ID")


# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)


# Generating PCA plot annotated by plate batch 
ggplot(pca_df, aes(PC1,PC2, color = PlateBatch))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output,"/", cur_date,"PCA_PlateBatch.png"))
dev.off()


# Generating PCA plot annotated by treatment 
ggplot(pca_df, aes(PC1,PC2, color = Treatment))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output,"/", cur_date,"PCA_Treatment.png"))
dev.off()


# Generating PCA plot annotated by tissue
ggplot(pca_df, aes(PC1,PC2, color = Tissue))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output,"/", cur_date,"PCA_Tissue.png"))
dev.off()
# all the heart samples seem to be grouped together, which makes sense
# lung samples display more variance than heart - so we may want to capture unwanted variance through SVA or RUV techniques


# Lets identify which samples are which
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output,"/", cur_date,"PCA_IDsLabeled.png"))
dev.off()

# This would suggest that M112_RedOakFlame may be a potential outlier, but we can investigate further
# We will conduct heirachical clustering to investigate further and to determine if we should remove these samples


# create a dataframe with samples as rows and genes as columns to input into hclust
countdata_forclustering <- t(countdata)
countdata_forclustering[1:5,1:10]

# run hierarchical clustering
hc <-hclust(dist(countdata_forclustering))
jpeg(filename= paste0(Output,"/",cur_date, "Hierachical_clustering_outliers.jpeg"),width = 4500, height = 350)
plot(hc)
dev.off()


# Here, the one aforementioned sample (M112_RedOakFlame) was separate from the rest, so should be removed



#################################################################################################
#################################################################################################
#### Removing sample(s) with potential QA/QC issues
#################################################################################################
#################################################################################################

# Since this is the first time we've filtered out a sample (unlike the ELGAN code) we need to set up vectors of IDs to keep, and then run a filter
# This step wasn't needed for this data set

# First, pull all the IDs from the count data
keep_samples <- colnames(countdata)
# Note that in this case, this represents all 170 IDs

# Then, remove the one sample outlier from the "keep_samples" vectors
keep_samples <- keep_samples[! keep_samples %in% c("M112_RedOakFlame")]
# Note that this represents all remaining 169 IDs

# Filter the countdata dataframe for these IDs
countdata <- countdata[, colnames(countdata) %in% keep_samples]
countdata[1:5, 1:10]

# Filter the coldata dataframe for these IDs
coldata <- coldata[coldata$ID %in% keep_samples, ]


# Export the raw data for just the included samples to potentially generate plots outside of R
write.csv(countdata, paste0(Output,"/", cur_date, "RawCounts_SamplesIncluded.csv"), row.names= TRUE)
write.csv(coldata, paste0(Output,"/", cur_date, "SampleInfo_SamplesIncluded.csv"), row.names= FALSE)





#################################################################################################
#################################################################################################
#### Setting up the actual DESeq2 function and associated algorithm ("experiment")
#################################################################################################
#################################################################################################
# More information can be found about the DeSeq2 package at https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8, as well as many other resources, as this package is well documented


# Ensuring that the appropriate variable types are recognized within the DESeq2 algorithm (e.g., factor vs numeric)
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Tissue <- as.factor(coldata$Tissue)
coldata$Group <- as.factor(coldata$Group)


# At this point, we need to make sure that the coldata is in the following format:
# rownames = the subject ID that matches the countdata
# all other columns indicate the covariate data that we may/may not include in the final experiment
coldata[1:5, 1:8]   # viewing data


# We also need to make sure that the countdata is in the following format:
# rownames = the gene identifiers
# all other columns indicate individual samples that match (and are in the same order as) the coldata
countdata[1:5, 1:10]   # viewing data



# All values should be integers (DESeq2 requires integer count values), so we'll need to convert values here if they aren't
# Let's check the class of the values:
sapply(countdata, class)
# for this dataset, all values are recognized as integers, so no need to convert




#################################################################################################
#### Creating a subset of the data to contrast
#################################################################################################


  #coldata_sub = coldata[coldata$Group %in% c(contrasts$Group[i], contrasts$Control[i]),]
  #countdata_sub = countdata[, colnames(countdata) %in% coldata_sub$ID]


# Create the final DESeq2 experiment, with appropriate experimental design:
# Note in the design: last factor indicates the factor we want to test the effect of; the first factor(s) = factors we want to control for
# Note that in this particular experiment, it gives us a warning that we are using integer data as numeric variables - this is what we want for some variables

dds = DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~Group)

# View what the experiment contains
dds
    


# Let's also make sure that we have the main contrast is in the order we want to calculate appropriate fold change values
# was no longer able to do this
#dds$Group <- relevel (dds$Group, contrasts$Control[i]) 


# Need to next, estimate the size factors, since size factors are used to normalize the counts (next step)
# The "iterate" estimator iterates between estimating the dispersion with a design of ~1, and finding a size factor vector by numerically optimizing the likelihood of the ~1 model.
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #check size factors

    
#################################################################################################
#################################################################################################
#### Pulling and exporting normalized count, for internal records and for potential future plots/figures
#################################################################################################
#################################################################################################
# Note that code is drafted to output several forms of normalized data - so users can choose which applies to their external purposes the best
# Also note that variance stabilized (vsd) values are commonly used for figures, which requires the experiment to run below, so will export those in the code below


# normalized counts:
normcounts<- counts(dds, normalized=TRUE)
    
write.csv(normcounts, paste0(Output,"/",cur_date, "_NormCounts_.csv"), row.names=TRUE)
    
    
# log2 pseudocounts (y=log2(n+1))
log2normcounts <- log2(normcounts+1)
write.csv(log2normcounts, paste0(Output,"/",cur_date, "_NormCounts_pslog2_.csv"), row.names=TRUE)
    

  


# Background filter: note that this is a different approach than what we usually apply to human epi-based analyses
# This is because of our variable contrasting conditions, cross-tissue analyses, etc
# I also needed this here, because genes were being retained that were expressed at values of zero throughout this subset, which didn't allow for SVA
# But this is actually the more commonly applied approach within DESeq2 examples
# Here, remove rows with only zeros, or only a single count across all samples:
idx <- rowSums(normcounts) > 1     # note that I've also seen a rowMeans > 1 filter applied here
CountsAboveBack <- normcounts[idx,]
nrow(CountsAboveBack)
    
write.csv(CountsAboveBack, paste0(Output,"/", cur_date, "NormCounts_AboveBack_.csv"), row.names=TRUE)
    
    
  # Also need to filter in the entire DESeq2 experiment:
dds <- dds[ rowSums(counts(dds, normalized=TRUE)) > 1, ]




#################################################################################################
#################################################################################################
#### Statistical analysis to detect significantly differentially expressed genes, first without SVA
#################################################################################################
#################################################################################################

  # Running the differential expression statistical pipeline
dds <- DESeq(dds, betaPrior=FALSE)      # because we used a user-defined model matrix, need to set betaPrior=FALSE
  # This ran a Wald test p-value
    
resultsNames(dds)
  
#start of running loop to compare all samples
significant_genes = c()
for (i in 1:length(contrasts$Group)){


  # Pulling statistical results
  res <- results(dds, pAdjustMethod = "BH", contrast = c("Group", contrasts$Group[i], contrasts$Control[i]))  #Statistical output with multiple test correction by the default, BH (aka FDR)
  #head(res)
    
  # Exporting statistical results:
  res_df = as.data.frame(res)[order(res$padj),] 
  filtered_res_df = res_df %>% filter(padj < 0.1) #filtering for only sig genes
  write.csv(filtered_res_df, paste0(Output_StatResults,"/", cur_date, "Stat_Results_", contrasts$Group[i] ,".csv"))
  
  #count number of significantly associated genes
  significant_genes = c(significant_genes, contrasts$Group[i], length(filtered_res_df$padj))
  
}

#putting significant genes into a table and exporting
dim(significant_genes) = c(2, length(significant_genes)/2)
sig_genes = data.frame(t(significant_genes))
colnames(sig_genes) = c('Treatment', 'Gene Count')
Model = rep('No', times = length(sig_genes$Treatment)) #adding a col denoting whether SVA was used or not
sig_genes = cbind(sig_genes, Model)
#write.csv(sig_genes, paste0(Output,"/", cur_date, "Total_Sig_Genes.csv"))

#################################################################################################
#################################################################################################
#### Surrogate variable analysis (SVA) to capture potential variances between samples (including sample / cell type heterogeneity)
#### Ending up not including in the final wildfire analysis
#################################################################################################
#################################################################################################
## For more information: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030161
## "Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis", Jeffrey T Leek,  John D Storey

# First creating an additional experiment to run SVA through
dds_SVA <- dds

# Set the model matrix for what's being used to fit the data (minus the SVA variables)
mod <- model.matrix(~Group, colData(dds_SVA))


# Set the null model matrix being compared against when fitting the data for the SVA analysis
mod0 <- model.matrix( ~1, colData(dds_SVA))


# Calculated the number of significant surrogate variables, using the following code:
svseq <- svaseq(CountsAboveBack, mod, mod0, n.sv=NULL)  
# Here, we derive 4 significant SVs

# Running SVA, here, using number of SVs = 4
svseq <- svaseq(CountsAboveBack, mod, mod0, n.sv=4)


#################
#### Adding surrogate variables into the DESeq2 analysis
#################


dds_SVA$SV1 <- svseq$sv[,1]
dds_SVA$SV2 <- svseq$sv[,2]
dds_SVA$SV3 <- svseq$sv[,3]
dds_SVA$SV4 <- svseq$sv[,4]
design(dds_SVA) <- ~ SV1 + SV2 + SV3 + SV4 + Group


# Export variance stabilized counts, which will be influenced by the new design above
# vsd will be useful for plots and future figure generation
vsd <- varianceStabilizingTransformation(dds_SVA, blind=FALSE)
vsd_matrix <-as.matrix(assay(vsd))
write.csv(vsd_matrix, paste0(Output,"/", cur_date, "VSDCounts_w_SVAs.csv"), row.names=TRUE)



#################
#### Evaluating variance stabilized values to guage whether SVs accounted for unwanted sources of variation between samples
#################

plotPCA(vsd, intgroup="Treatment", returnData=FALSE)
plotPCA(vsd, intgroup="Tissue", returnData=FALSE)
plotPCA(vsd, intgroup="Group", returnData=FALSE)
  
plotPCA(vsd, intgroup="SV1", returnData=FALSE)
plotPCA(vsd, intgroup="SV2", returnData=FALSE)
plotPCA(vsd, intgroup="SV3", returnData=FALSE)
plotPCA(vsd, intgroup="SV4", returnData=FALSE)


#In sum, it looks like the SVs are potentially accounting for a lot of the unwanted variation



#################################################################################################
#################################################################################################
#### Statistical analysis to detect significantly differentially expressed genes, now including surrogate variables
#################################################################################################
#################################################################################################

# Running the differential expression statistical pipeline
dds_SVA <- DESeq(dds_SVA, betaPrior=FALSE)
  
resultsNames(dds_SVA) #check the available comparisons

#start of running loop to compare all samples
significant_genes_w_SVA = c()
for (i in 1:length(contrasts$Group)){
  
  
  # Pulling statistical results
  res <- results(dds_SVA, pAdjustMethod = "BH", contrast = c("Group", contrasts$Group[i], contrasts$Control[i]))  #Statistical output with multiple test correction by the default, BH (aka FDR)
  #head(res)
  
  # Exporting statistical results:
  res_df = as.data.frame(res)[order(res$padj),] 
  filtered_res_df = res_df %>% filter(padj < 0.1) #filtering for only sig genes_w_SVA
  write.csv(filtered_res_df, paste0(Output_StatResults_SVA,"/", cur_date, "Stat_Results_w_SVA", contrasts$Group[i] ,".csv"))
  
  #count number of significantly associated genes_w_SVA
  significant_genes_w_SVA = c(significant_genes_w_SVA, contrasts$Group[i], length(filtered_res_df$padj))
  
}

#putting significant genes_w_SVA into a table and exporting
dim(significant_genes_w_SVA) = c(2, length(significant_genes_w_SVA)/2)
sig_genes_w_SVA = data.frame(t(significant_genes_w_SVA))
colnames(sig_genes_w_SVA) = c('Treatment', 'Gene Count')
Model = rep('Yes', times = length(sig_genes_w_SVA$Treatment)) #adding a col denoting whether SVA was used or not
sig_genes_w_SVA = cbind(sig_genes_w_SVA, Model)
final_sig_gene_total = rbind(sig_genes, sig_genes_w_SVA)
write.csv(final_sig_gene_total, paste0(Output,"/", cur_date, "Total_Sig_Genes.csv"), row.names = FALSE)
