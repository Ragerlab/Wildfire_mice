# Installing DESeq2:
# Note that BiocManager allows you to access the Bioconductor repo in which DeSeq2 package is. DeSeq2 pacakge has a "tibble" "RcppArmadillo" and "rlang" dependency.
# to install RcppArmadillo, you must install a Fortan compiler (I used gfortran-6.1.pkg, https://cran.r-project.org/bin/macosx/tools/).

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("tibble")
install.packages("rlang")
install.packages("RcppArmadillo")

# Install DESeq and SVA:
BiocManager::install("DESeq2", version = "3.11")
BiocManager::install("sva", version ="3.11")
BiocManager::install("BiocManager")

# Install plot packages
install.packages("ggbeeswarm") # note that I had to manually install this package using tools --> install packages
#install.packages("gridExtra")

# Install tidyverse packages 
#install.packages("tidyverse")

# Install gplots packages
install.packages("gplots")

# Install RColorBrewer package
install.packages("RColorBrewer")

# Install scales (needed for tidyverse activation):
install.packages("scales")

# Install data.tables
install.packages("data.table")

# Activating appropriate libraries (skip step above if everything is already installed)
# library(data.table)
library(rlang)
# library(DESeq2)
# library(limma)
# library(sva)
library(ggplot2)
library(stats)
library(scales)
library(ggbeeswarm)
library(gridExtra)
library(tidyverse)
library(dplyr)
# library(gplots)
library(RColorBrewer)
# library(data.table)
library(readxl)

# Set working directory
setwd('/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/2_Wildfire_Analysis/Input')

# Create an output folder (make sure to make the folder first, and then point to it here)
Output <- ("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/2_Wildfire_Analysis/Output")

# Read in count data
# countdata <- read.table(file = 'Wildfire_GeneCounts_102320.csv', header = TRUE)
countdata <- read_csv(file = 'Wildfire_GeneCounts_102320.csv') 

#only interested in comparing mouse heart peat flame (MH112-MH118) vs. mouse heart saline samples (MH101-MH106)
subsetted_colnames = substr(colnames(countdata),1,5)
samples_of_interest = c('Gene','MH112', 'MH113', 'MH114', 'MH115', 'MH116', 'MH117', 'MH118','MH101',
                                           'MH102', 'MH103', 'MH104', 'MH105','MH106') 
countdata_subset = countdata[,(substr(colnames(countdata),1,5) %in% samples_of_interest)]

dim(countdata_subset) #14 samples (plus Gene ID column), 30146 genes
# visualize this data quickly by viewing top left corner:
countdata_subset[1:3,1:6]

# Check for duplicate gene IDs in the countdata_subset 
Dups <- duplicated(countdata_subset[,1])
summary(Dups)
# Not an issue for this wildfire dataset

# Read in metadata (subject info file)
subjectinfo <- read_excel("SampleID_Information_102320.xlsx")
dim(subjectinfo) #170 rows, 7 col

# Visualize this data quickly by viewing top left corner:
subjectinfo[1:3,1:6]

##########
##########
# (1) Filtering to first keep subjects present in both files
##########
##########

# Pull all mRNAseq IDs that overlap between countdata_subset & subjectinfo
keep_samples <- colnames(countdata_subset[, colnames(countdata_subset) %in% subjectinfo$FullID]) 
# Note that in this case, this represents 13 IDs
# Filter the countdata dataframe for these IDs
# Because the first column is the gene identifier information, this will get lost if we don't also grab that 
# separately (hence, a few steps here)
countdata_temp <- countdata[, colnames(countdata) %in% keep_samples] # pulling countdata, without gene identifier col
countdata <- cbind(countdata[,1], countdata_temp)    # adding the gene identifier column back in
colnames(countdata)[1] <- "Gene"   # renaming the gene identifer column
countdata_temp <- NULL    # removing temporary dataframe

# Viewing merged, filtered count data
countdata[1:3,1:10]
# Filter the subjectinfo dataframe for these IDs
subjectinfo <- subjectinfo[subjectinfo$FullID %in% keep_samples, ]

##########
##########
# (2) Background filter, to filter out genes that were universally lowly expressed
##########
##########


# Create variable for number of samples, assuming that one variable is Gene: 
numberofsamples = (ncol(countdata_subset)-1)

# Create a variable that is the median of all expression levels across all samples and genes
totalmedian = median(as.matrix(countdata_subset[,2:ncol(countdata_subset)]))

# Create a list of genes that pass the background filter:
Genes_backgroundfiltered <- countdata_subset %>% 
  #gather so that there is one row for each gene-sample pairing
  gather(key = "sampleID", value = "expression", 2:ncol(countdata_subset)) %>% 
  #add a variable that takes the value 1 is the expression value is above the median 
  mutate(abovemedian = ifelse(expression > totalmedian,1,0)) %>% 
  #sum the abovemedian variable created above by gene so that totalabovemedian calcualtes the number of samples for each gene
  #that have expression levels above the median
  group_by(Gene) %>% 
  summarize(totalabovemedian = sum(abovemedian)) %>% 
  #remove all genes that have a totalabovemedian equal to less than a specified percentage of the total number of samples 
  filter(totalabovemedian > (0.20 * numberofsamples)) %>% 
  #select only gene so that you essentially have a list of the genes to include 
  select(Gene) 
# Select genes that passed background filter only to proceed with analysis
countdata_subset <- left_join(Genes_backgroundfiltered, countdata_subset, by = "Gene")
# This resulted in filtering the original list of 30,146 genes down to 16,163 genes

##########
##########
#### (3) Check subject filter, to check that we have filtered out subjects that had samples present all values of zero
#### (and thus not detected)
##########
##########

# Note that this code was originally drafted using a transposed version of the countdata_subset dataframe
# So here, rather than modifying, we will transpose and then transpose back at the end

#transpose count data 
countdata_subset_t <- countdata_subset %>% 
  gather(key = "sampleID", value = "expression", 2:ncol(countdata_subset)) %>% 
  spread(Gene, expression)  
# view transposed count data
# countdata_subset_t[1:4,1:10]

countdata_subset_t <- countdata_subset_t %>% 
  #create variable that is the sum of count values across all genes, for each sample
  mutate(rowsum =  rowSums(countdata_subset_t[,2:nrow(countdata_subset_t)],na.rm  =  FALSE, dims = 1L)) 

#identify the samples that will be removed, for records 
samples_zerocounts <- countdata_subset_t %>% 
  filter(rowsum  == 0) %>% 
  select("sampleID") 

countdata_subset_t <- countdata_subset_t %>% 
  #remove samples that have a total sum of all counts of zero 
  filter(rowsum != 0) %>% 
  select(-rowsum)

# SAMPLES REMOVED: No samples were removed 


# Transposing the countdata_subset back such that the genes are the rownames and the samples are the columnnames
countdata_subset <- countdata_subset_t %>% 
  gather(Gene, value, 2:ncol(countdata_subset_t)) %>% 
  spread(sampleID, value) 

##########
##########
# (4) Re-filter to include subjects meeting the above filters in (2) and (3)
##########
##########

# Pull all mRNAseq IDs that overlap between countdata_subset & subjectinfo
keep_samples <- colnames(countdata_subset[, colnames(countdata_subset) %in% subjectinfo$FullID]) 
# Note that in this case, this represents 372 IDs

# Filter the countdata_subset dataframe for these IDs
# Because the first column is the gene identifier information, this will get lost if we don't also grab that 
# separately (hence, a few steps here)
# pulling all countdata_subset, without gene identifier column
countdata_subset_temp <- countdata_subset[, colnames(countdata_subset) %in% keep_samples]   
# adding the gene identifier column back in
countdata_subset <- cbind(countdata_subset$Gene, countdata_subset_temp)    
colnames(countdata_subset)[1] <- "Gene"   # renaming the gene identifer column
# Viewing merged, filtered count data
# countdata_subset[1:3,1:10]

# Filter the subjectinfo dataframe for these IDs
subjectinfo <- subjectinfo[subjectinfo$FullID %in% keep_samples, ]

# Make the gene variable the rowname 
countdata_subset <- countdata_subset %>% 
  column_to_rownames(var = "Gene") 

# Creating the coldata object, based on information in the subjectinfo file
coldata <- subjectinfo
# Let's change the name of the primary ID column, to use from here on out
names(coldata)[names(coldata) == 'mRNAseq_ID'] <- 'ID'


# Set the rownames of coldata and column names of countdata_subset to be in the same order 
countdata_subset <- setcolorder(countdata_subset, as.character(coldata$ID))

# Double checking that the same variables appear between the two dataframes
setequal(as.character(coldata$ID), colnames(countdata_subset))

# Additionally checking that not only the sets of variables are the same, but that they are in the same order
identical(as.character(coldata$ID), colnames(countdata_subset))

# One way to evaluate sequencing data quality / identify potential outliers is through Principal Component Analysis (PCA):
# PCA helps in identifying outlying samples for quality control, and gives a feeling for the principal causes of variation in a dataset

# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))
screeplot(pca)
#this scree plot indicates that nearly all variation is explained in PC1,2 and 3, 
#thus PCA is a useful exercise 

# Make dataframe for PCA plot generation using first three components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], Sample=colnames(countdata))
# Add attributes (covariates from coldata) to pca df
pca_df <- merge(pca_df, coldata, by.x="Sample", by.y="ID")



# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)

# Generating PCA plot annotated by sex 
ggplot(pca_df, aes(PC1,PC2, color = sex))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by maternal age 
ggplot(pca_df, aes(PC1,PC2, color = mage))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by gestational age 
ggplot(pca_df, aes(PC1,PC2, color = gadays))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by race 
ggplot(pca_df, aes(PC1,PC2, color = race1))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Lets identify the outliers
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))


# This would suggest that LS18017 and LS18097 are potential outliers, as well as potentially LS18183 and LS18076
# We will conduct heirachical clustering to investigate further and to determine if we should also remove these samples


# create a dataframe with samples as rows and genes as columns to input into hclust
countdata_forclustering <- t(countdata)
countdata_forclustering[1:5,1:10]


# run heirarchical clustering
hc <-hclust(dist(countdata_forclustering))
jpeg(filename= paste0(Output_Folder,"/",cur_date, "_heirachicalclusteringforoutliers.jpeg"),width = 4500, height = 350)
plot(hc)
dev.off()

# Here, the four aforementioned samples were clearly separate from the rest, so should be removed

# Remove these subject IDs from our "keep_samples" vector
keep_samples <- keep_samples[! keep_samples %in% c("LS18017", "LS18097", "LS18183", "LS18076")]
# Note that in this case, this represents 368 IDs

# Filter the countdata dataframe for these IDs
countdata <- countdata[, colnames(countdata) %in% keep_samples]
countdata[1:5, 1:10]

# Filter the coldata dataframe for these IDs
coldata <- coldata[coldata$ID %in% keep_samples, ]



# Export the raw data for just the included samples to potentially generate plots outside of R
write.csv(countdata, paste0(Output_Folder,"/", cur_date, "_RawCounts_SubjectsIncludedinModel.csv"), row.names= TRUE)
write.csv(coldata, paste0(Output_Folder,"/", cur_date, "_SubjectInfo_SubjectsIncludedinModel.csv"), row.names= FALSE)
# Note that if, for future studies, we exclude subjects with missing covariate data, that step will have to incorporated prior to this export
# But for this analysis, we're imputing missing covariate data, so went ahead and exported raw data here



# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))
screeplot(pca)
#this scree plot indicates that nearly all variation is explained in PC1,2 and 3, 
#thus PCA is a useful exercise 

# Make dataframe for PCA plot generation using first three components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], Sample=colnames(countdata))
# Add attributes (covariates from coldata) to pca df
pca_df <- merge(pca_df, coldata, by.x="Sample", by.y="ID")

# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)

# Generating PCA plot annotated by sex 
ggplot(pca_df, aes(PC1,PC2, color = sex))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by maternal age 
ggplot(pca_df, aes(PC1,PC2, color = mage))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by gestational age 
ggplot(pca_df, aes(PC1,PC2, color = gadays))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Generating PCA plot annotated by race 
ggplot(pca_df, aes(PC1,PC2, color = race1))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# Lets identify potential outliers (although here, since this is the 2nd run, we would be much stricter on this; meaning that there would have to be a super obvious outlier, or else the cycle would just continue)
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))

# At this point, we need to make sure that the coldata is in the following format:
# rownames = the subject ID that matches the countdata
# all other columns indicate the covariate data that we may/may not include in the final experiment
coldata[1:5, 1:7]   # viewing data


# At this point, we need to make sure that the countdata is in the following format:
# rownames = the gene identifiers
# all other columns indicate individual samples that match (and are in the same order as) the coldata
countdata[1:5, 1:10]   # viewing data


# all values should be integers (DESeq2 requires integer count values), need to convert values here if haven't already
sapply(countdata, class)   # checking the class of the values
countdata_2 <- as.data.frame(sapply(countdata, as.integer))   # converting into integers if showing numeric
sapply(countdata_2, class)   # re-checking the class of the values
row.names(countdata_2) <- row.names(countdata)  # this removed the row names (gene identifiers), so had to merge these back
countdata <- countdata_2



# Create the final DESeq2 experiment, with appropriate experimental design:
# Note in the design: last factor indicates the factor we want to test the effect of; the first factor(s) = factors we want to control for
# Note that in this particular experiment, it gives us a warning that we are using integer data as numeric variables - this is what we want for some variables

dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~race1 + mage + sex + gadays + zcat4 + medu + ASD)

# View what the experiment contains
dds


# Let's also make sure that we have the main contrast (last term) in the order we want to calculate appropriate fold change values
# In this case, we want ASD=0 as denominators (first level, aka, the "control" level), and ASD=1, cases as numerators
dds$ASD <- relevel (dds$ASD, "0")



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
# Also note that variance stabilized (vsd) values are commonly used for figures, which requires the experiment to run and adjust for covariates below, so will export these in the code below


# normalized counts:
normcounts<- counts(dds, normalized=TRUE)
write.csv(normcounts, paste0(Output_Folder,"/", cur_date, "_NormCounts_IncludedSubjects.csv"), row.names=TRUE)

# pseudocounts to make plots easier:
ps_normcounts <- normcounts + 1
write.csv(ps_normcounts, paste0(Output_Folder,"/",cur_date, "_NormCounts_ps_IncludedSubjects.csv"),row.names=TRUE)

# log2 pseudocounts (y=log2(n+1))
log2normcounts <- log2(normcounts+1)
write.csv(log2normcounts, paste0(Output_Folder,"/", cur_date, "_NormCounts_pslog2_IncludedSubjects.csv"), row.names=TRUE)




#################################################################################################
#################################################################################################
#### Surrogate variable analysis (SVA) to capture potential variances between samples (including batch, but also sample / cell type heterogeneity)
#################################################################################################
#################################################################################################
## For more information: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030161
## "Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis", Jeffrey T Leek,  John D Storey


# Set the model matrix for what's being used to fit the data (minus the SVA variables)
mod <- model.matrix(~race1 + mage + sex + gadays + zcat4 + medu + ASD, colData(dds))


# Set the null model matrix being compared against when fitting the data for the SVA analysis
mod0 <- model.matrix( ~1, colData(dds))


# Running SVA, here, using number of SVs = 3
# Note that this number of significant SVs selection was based on several iterations, where different numbers of SVs were tested for the ELGAN dataset
svseq <- svaseq(normcounts, mod, mod0, n.sv=3)

# I also tested an iteration where svaseq calculated the number of significant SVs, using the following code:
# svseq <- svaseq(normcounts, mod, mod0, n.sv=NULL)  
# Number of significant surrogate variables is much higher, but this muted global signal responses


#################
#################
#### Adding surrogate variables into the DESeq2 analysis
#################
#################

dds$SV1 <- svseq$sv[,1]
dds$SV2 <- svseq$sv[,2]
dds$SV3 <- svseq$sv[,3]
design(dds) <- ~ SV1 + SV2 + SV3 + race1 + mage + sex + gadays + zcat4 + medu + ASD


# Export variance stabilized counts, which will be influenced by the new design above
# vsd will be useful for plots and future figure generation
# CAUTION: THIS STEP TAKES SOME TIME! RUN TIME IS ABOUT 20 MINUTES FOR MRNASEQ DATA ON JULIA'S COMPUTER
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_matrix <-as.matrix(assay(vsd))
write.csv(vsd_matrix, paste0(Output_Folder,"/", cur_date, "_VSDCounts_IncludedSubjects.csv"), row.names=TRUE)



#################
#################
#### Evaluating variance stabilized values to guage whether SVs accounted for unwanted sources of variation between samples
#################
#################

plotPCA(vsd, intgroup="sex", returnData=FALSE)
plotPCA(vsd, intgroup="mage", returnData=FALSE)
plotPCA(vsd, intgroup="SV1", returnData=FALSE)
#SV1 accounts for a lot of variation along PC1 as shown by the steady gradient in color across the PC1 axis
plotPCA(vsd, intgroup="SV2", returnData=FALSE)
plotPCA(vsd, intgroup="SV3", returnData=FALSE)
#SV2 and 3 also account for a lot of variation across PCs

#In sum, it looks like the SVs are accounting for a lot of the unwanted variation which is great!



#################################################################################################
#################################################################################################
#### Statistical analysis to detect significantly differentially expressed genes
#################################################################################################
#################################################################################################

# Running the differential expression statistical pipeline
# CAUTION: THIS STEP TAKES A LOT OF TIME! RUN TIME IS ABOUT 1 HOUR FOR MRNASEQ DATA ON JULIA'S COMPUTER
dds <- DESeq(dds, betaPrior=FALSE)      # because we used a user-defined model matrix, need to set betaPrior=FALSE
# This ran a Wald test p-value (ASD 1 vs. 0)

resultsNames(dds) #check the available comparisons

# Pulling statistical results
res <- results(dds, pAdjustMethod = "BH")  #The  with multiple test correction by the default, BH (aka FDR)
head(res)


# Exporting statistical results:
write.csv(as.data.frame(res)[order(res$padj),], paste0(Output_Folder,"/", cur_date,"_AllStatResults.csv"))


#################################################################################################
#################################################################################################
#### Overview of statistical results
#################################################################################################
#################################################################################################
summary(res)

# Basic filters
table(res$padj<0.10)

table(res$padj<0.10 & res$log2FoldChange > log2(1.5))
table(res$padj<0.10 & res$log2FoldChange < -log2(1.5))

table(res$padj<0.10 & res$log2FoldChange > log2(1.3))
table(res$padj<0.10 & res$log2FoldChange < -log2(1.3))

table(res$padj<0.10 & (res$log2FoldChange > log2(1.5) | res$log2FoldChange < -log2(1.5)))


# Sort results and create a dataframe (df) for exporting & certain plots in the code below
resdf <- as.data.frame(res)[order(res$padj),]



# Also create a df for filtered data for certain plots in the code below
resdf_filter <- resdf[which(resdf$padj<0.10 & (resdf$log2FoldChange > log2(1.5) | resdf$log2FoldChange < -log2(1.5))),]
resdf_pfilter <- resdf[which(resdf$padj<0.10),]

###Heatmap###

#make dataset for heatmap in partek
countsforhm <- vsd_matrix

resdf_pfilter <- resdf_pfilter %>% 
  rownames_to_column("gene") %>% 
  filter(abs(log2FoldChange)>log2(1.3)) %>% 
  arrange(log2FoldChange)

genenames <- as.data.frame(resdf_pfilter$gene) %>% 
  rename("gene"="resdf_pfilter$gene")

countsforhm <- as.data.frame(countsforhm) %>% 
  rownames_to_column(var="gene") 

countsforhm <- left_join(genenames, countsforhm, by="gene") %>% 
  column_to_rownames(var="gene")

countsforhm <- as.data.frame(t(countsforhm)) %>% 
  rownames_to_column(var="ID") 

countsforhm <- full_join(countsforhm, coldata, by="ID") %>% 
  select(-"race",-"gadays",- "mage")

countsforhm <- countsforhm %>% 
  column_to_rownames(var="ID")

write.csv(countsforhm, paste0(Output_Folder,"/", cur_date,"_countdata_forheatmap_mrnas.csv"))

### MA plot ###

###################################
####  MA plot in ggplot2
###################################
MA <- resdf
MA[is.na(MA)] <- 1
MA1 <- MA[ which(MA$padj>=0.1),]  # all the rest
MA2 <- MA[ which(MA$padj<0.1 & MA$log2FoldChange > log2(1.3)),]   # sig up-regulated w/ FC filter-FIREBRICK
MA3 <- MA[ which(MA$padj<0.1 & MA$log2FoldChange < -log2(1.3)),]  # sig down-regulated w/ FC filter-dodgerblue4
MA4 <- MA[ which(MA$padj<0.1 & MA$log2FoldChange > 0 & MA$log2FoldChange <= log2(1.3)),]  # sig up-regulated w/o FC filter- darkorchid3
MA5 <- MA[ which(MA$padj<0.1 & MA$log2FoldChange < 0 & MA$log2FoldChange >= -log2(1.3)), ] # sig down-regulated w/o FC filter- chartreuse4

# get names of probes if annotating the graphic with them, according to specific criteria (not used here):
namesmmu=rownames(MA)[which((MA$log2FoldChange > 2 | MA$log2FoldChange < -2) & MA$padj < 0.1 & MA$baseMean > 0)]
namesmmu2=rownames(MA2)[which((MA2$log2FoldChange > 2 | MA2$log2FoldChange < -2) & MA2$padj < 0.1 & MA2$padj >= 0.05 & MA2$baseMean > 0)]
namesmmu3=rownames(MA3)[which((MA3$log2FoldChange > 2 | MA3$log2FoldChange < -2) & MA3$padj < 0.05 & MA3$baseMean > 0)]                                               

plot <-
  ggplot(subset(MA1, baseMean>=0), aes(x = baseMean, y = log2FoldChange)) + #can subset to exclude the very low counts (<1 average count) if we want
  geom_point(color="gray56", size = 2) + 
  geom_point(data = MA2, color="firebrick", size=3, show.legend = TRUE) + # colors used: firebrick, dodgerblue4, chocolate1, darkorchid3, chartreuse4
  geom_point(data = MA3, color="dodgerblue4", size=3, show.legend = TRUE) +
  geom_point(data = MA4, color="darkorchid3", size=3, show.legend = TRUE) +
  geom_point(data = MA5, color="chartreuse4", size=3, show.legend = TRUE) +
  theme_bw() +
  scale_x_continuous(trans = "log10", breaks=c(1,10,100, 1000, 10000, 100000, 1000000), labels=c("1","10","100", "1000", "10000", "100000", "1000000")) + 
  scale_y_continuous(limits=c(-1.5, 1.5), breaks=c(1.5,1,0.5,0,-0.5,-1,-1.5), labels=c(1.5,1,0.5,0,-0.5,-1,-1.5)) +
  xlab("Expression (Normalized Count)") + ylab("Fold Change (log2) Cases/Control") +
  labs(title="ASD vs. Control") +   
  #coord_fixed(ratio=0.5) +
  geom_hline(yintercept=0) 
#annotate("text",label="# DEPs higher in Cr(VI): 78\n(representing 77 DEGs)", x=1.1,y=10, size=4, hjust=0, color="darkorchid3") +
#annotate("text",label="# DEPs higher in Captan: 86\n(representing 85 DEGs)", x=1.1,y=-10, size=4, hjust=0, color="chartreuse4")+
#geom_text(data=MA2[namesmmu2,], size=3, hjust=-.05, vjust=.7, col="deepskyblue", aes(label=namesmmu2,size=5, angle=0)) + #annotate specific probes with text
#geom_text(data=MA3[namesmmu3,], size=3, hjust=-.05, vjust=.7, col="firebrick3", aes(label=namesmmu3,size=5, angle=0))

plot

#export single plots
png(paste0(Output_Folder, "/", cur_date, "_MAplot_v1.png"),
    width = 1600, height = 1200, units = "px", res = 250)
plot
dev.off()
