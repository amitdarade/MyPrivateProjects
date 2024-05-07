# MyPrivateProjects
Hi, This is Private projects. Dont Waste your Time.

#Clear Memory and Residues
rm(list = ls(all.names = TRUE))
gc()
#Packages 

#Libraries 
library(affy)
library(GEOquery)
library(tidyverse)
library(oligo)
library(WGCNA)
library(NOISeq)
library(impute)
library(arrayQualityMetrics)
library(dplyr)
library(tidyr)
library(genefilter)
library(Biobase)
library(hgu133plus2.db)
library(pheatmap)
library(limma)
library(ggplot2)
library(ggrepel)
library(readr)
library(limma)
library(org.Hs.eg.db)
library(EnhancedVolcano)

#Once all packages have been loaded, run sessionInfo().
sessionInfo()

#Importing and Prepare Data.

# Import Raw Data

getGEOSuppFiles("GSE18090")

untar("your path", exdir = 'data/')

raw.data <- ReadAffy(celfile.path = "data/")

gse <- getGEO("GSE18090", GSEMatrix = TRUE)

#Import GSE matrix for phenodata.

gse <- getGEO("GSE18090", GSEMatrix = TRUE)

feature.data <- gse$GSE18090_series_matrix.txt.gz@featureData@data

pheno.data <- gse$GSE18090_series_matrix.txt.gz@phenoData@data

# Extract expression data from raw data

class(raw.data)

exprs_r.d1 <- exprs(raw.data)

# Extract clinical variables from pheno.data.
sampleInfo3 <- pheno.data
sampleInfo3 <- dplyr::select(sampleInfo3, title)
sampleInfo3_split <- sampleInfo3 %>%
  separate(title, into = c("Group", "Patients"), sep = " ", extra = "merge")
print(sampleInfo3_split)

# Visualize data & Check Quality Before Normalization.

# Histogram
hist(raw.data)

#MA Plot
ma.plot(rowMeans(log2(exprs_r.d)), log2(exprs_r.d[,2])-log2(exprs_r.d[,3]), cex=1)

# Dendrogram
htree <- hclust(dist(t(exprs_r.d)), method = "average")
plot(htree)

#Boxplot
# Obtain group information and corresponding colors
group <- sampleInfo3_split$Group
group_colors <- c("pink", "red", "darkgreen")
sample_colors <- group_colors[match(group, unique(group))]

# Set up plot parameters
par(mar = c(8, 4, 4, 6)) # Adjust margins

# Create boxplot with colored boxes and legend
{boxplot(exprs_r.d1, outline = FALSE, col = sample_colors, 
        names.arg = colnames(exprs_r.d1), las = 2, cex.axis = 0.6, boxwex = 0.6,
        main = "Boxplot Before Normalization", # Add title
        ylab = "Expression Levels") # Add axis labels
  legend("topright", legend = unique(group), fill = group_colors, title = "Group", cex = 0.6)
}







# NORMALIZATION

normalized.data <- rma(raw.data)

normalized.expr <- as.data.frame(exprs(normalized.data))

#Check Missing Values.

missing_values <- sum(is.na(normalized.expr) )

head(normalized.expr)

# Check for Outliers with goodSampleGenes() function.

gsg <- goodSamplesGenes(normalized.expr)

summary(gsg)

gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# Boxplot After Normalization.
boxplot(normalized.expr,outline=FALSE, col = "skyblue")
{boxplot(normalized.expr, outline = FALSE, col = sample_colors, 
         names.arg = colnames(exprs_r.d1), las = 2, cex.axis = 0.6, boxwex = 0.6,
         main = "Boxplot After Normalization", 
         ylab = "Expression Levels") 
  legend("topright", legend = unique(group), fill = group_colors, title = "Group", cex = 0.6)
}

# PCA before Batch correction 
pca <- prcomp(t(normalized.expr))
pca.dat <- pca$x
pca_df <- cbind(sampleInfo3_split, as.data.frame(pca$x))

{ggplot(pca_df, aes(x = PC1, y = PC2, col = Group, label = paste("Patient", Patients))) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") + 
    geom_point(size = 3, alpha = 0.8) +  
    geom_text_repel(size = 3, segment.size = 0.2, segment.color = "grey60") +  
    stat_ellipse(aes(geom = "polygon", alpha = .5))+
    theme_grey() +  
    labs(x = "PC1", y = "PC2") +  # Customize axis labels
    scale_color_manual(values = c("purple", "red", "darkgreen")) + 
    ggtitle("PCA Plot of Gene Expression Data Showing Batch Effects") + 
    theme(plot.title = element_text(hjust = 0.5))}

# BATCH CORRECTION 

BC_PCA <- readData(normalized.data , factor = sampleInfo3_split)
myPCA = ARSyNseq(BC_PCA, factor="Group", batch = FALSE, norm = "n", logtransf = TRUE)

# Rename dataset 

Normalize.batchdata <- myPCA

# Extract expression data
Batch.expr <- as.data.frame(exprs(Normalize.batchdata ))

# PCA After Batch correction 

pca_BC <- prcomp(t(Batch.expr)
pca.dat_BC <- pca_BC$x
pca_df_BC <- cbind(sampleInfo3_split, as.data.frame(pca_BC$x))

{ggplot(pca_df_BC, aes(x = PC1, y = PC2, col = Group, label = paste("Patient", Patients))) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") + 
    geom_point(size = 3, alpha = 0.8) +  
    geom_text_repel(size = 3, segment.size = 0.2, segment.color = "grey60") +  
    stat_ellipse(aes(geom = "polygon", alpha = .5))+
    theme_grey() +  
    labs(x = "PC1", y = "PC2") +  # Customize axis labels
    scale_color_manual(values = c("purple", "red", "darkgreen")) + 
    ggtitle("PCA Plot of Gene Expression Data") + 
    theme(plot.title = element_text(hjust = 0.5))}

# other plots

hist(Normalize.batchdata)

ma.plot(rowMeans(log2(Batch.expr)), log2(Batch.expr[,2])-log2(Batch.expr[,3]), cex=1)

htree <- hclust(dist(t(Batch.expr)), method = "average")
plot(htree)







# GENES FILTERING

# Detecting Most Variable genes.

sds <- apply ((Batch.expr), 1, sd)
sdsO<- sort(sds)

# Visualization 
plot(1:length(sdsO), sdsO, 
     main = "Distribution of variability for all genes",
     sub = "Vertical lines represent 90% and 95% percentiles",
     xlab = "Gene index (from least to most variable)", 
     ylab = "Standard deviation")
abline(v=length(sds)*c(0.9,0.95))


#Filtering least variable genes.
annotation(Normalize.batchdata) <- "hgu133plus2.db"

filtered_1 <- nsFilter(Normalize.batchdata, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.50, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")

print(filtered_1$filter.log)
print(filtered_1$eset)

# Save filtered gene for DGE analysis.

eset_filtered <-filtered$eset

# Visualization 

corMatrix <- cor(eset_filtered, use="c")
pheatmap(corMatrix) 

rownames(sampleInfo3_split)
colnames(corMatrix)

rownames(sampleInfo3_split) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo3_split)






# DGE Analysis with Limma.

sampleInfo3_split$Group

# Model matrix.
design <- model.matrix(~0+sampleInfo3_split$Group)

dim(design)
head(design)
design

# Rename the column names 
colnames(designbc) <- c("DF","DHF","ND")
designbc

# Model fitting
fit <- lmFit(exprs_filtered, designbc)
head(fit$coefficients)

# Contrast matrix
contrasts <- makeContrasts(ND - DF, ND - DHF, DF - DHF, levels=design)

# Model fitting 2
fit2 <- contrasts.fit(fit, contrasts)

# Ebayes 
fit2 <- eBayes(fit2)

#Results

topTable(fit2)

decideTests(fit2)

table(decideTests(fit2))

# Obtaining lists of differentially expressed genes

topTab_NDvsDF <- topTable (fit2, number=nrow(fit2), coef=1, adjust="fdr") 
 head(topTab_NDvsDF)

topTab_NDvsDHF <- topTable (fit2, number=nrow(fit2), coef=2, adjust="fdr") 
head(topTab_NDvsDHF)

topTab_DFvsDHF <- topTable (fit2, number=nrow(fit2), coef=3, adjust="fdr") 
head(topTab_DFvsDHF)

# GENE  ANNOTATION

annotatedTopTable <- function(topTab, anotPackage){
       topTab <- cbind(PROBEID=rownames(topTab), topTab)
       myProbes <- rownames(topTab)
       thePackage <- eval(parse(text = anotPackage))
       geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
       annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
     return(annotatedTopTab)}

topAnnotated_NDvsDF <- annotatedTopTable(topTab_NDvsDF, anotPackage="hgu133plus2.db")

topAnnotated_NDvsDHF<- annotatedTopTable(topTab_NDvsDHF, anotPackage="hgu133plus2.db")

topAnnotated_DFvsDHF <- annotatedTopTable(topTab_DFvsDHF, anotPackage="hgu133plus2.db")

# Save Results

write.csv(topAnnotated_NDvsDF, file = "NDvsDF.csv", row.names = TRUE)

write.csv(topAnnotated_NDvsDHF, file = "NDvsDHF.csv", row.names = TRUE)

write.csv(topAnnotated_DFvsDHF, file = "DFvsDHF.csv", row.names = TRUE)


# VIZUALIZATION OF DEG OF ND-DF,ND-DHF,DF-DHF


# Theme
theme_set(theme_bw(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1),  
              color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#The Volcano Plot for NDvsDF coef=1

{
NDvsDF <- topTable(fit2, coef=1, number=Inf)
NDvsDF <- tibble::rownames_to_column(NDvsDF,"ID")

NDvsDF$diffexpress <- 'NO'
NDvsDF$diffexpress[NDvsDF$logFC > 0.5 & NDvsDF$P.Value < 0.05] <- 'UP'
NDvsDF$diffexpress[NDvsDF$logFC < -0.5 & NDvsDF$P.Value < 0.05] <- 'DOWN'

top30deg <- head(NDvsDF[order(NDvsDF$P.Value),"Gene.Symbol"], 30)

NDvsDF$delabel <- ifelse(NDvsDF$Gene.Symbol %in% top30deg, NDvsDF$Gene.Symbol, NA)

head(NDvsDF)

write.csv(NDvsDF, file = "NDvsDF.csv", row.names = TRUE)

ggplot(NDvsDF, aes(x = logFC, y = -log10(P.Value))) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.5), col = "gray", linetype = 'dashed') + 
  geom_point() 

ggplot(NDvsDF, aes(x = logFC, y = -log10(P.Value), col = diffexpress, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = c(0.1), col = "black", linetype = 'dashed') + 
  geom_point() +
  scale_colour_manual(values = c("#00AFBB", "grey", "red"),
                      labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 5), xlim = c(-5, 5)) +
  scale_x_continuous(breaks = seq(-4, 4, 1))+
  labs(color = 'NDvsDF', 
       x = expression("LogFC"), y = expression("P.Value")) +
  ggtitle("DGE for NDvsDF")+
  geom_text_repel(max.overlaps = Inf)
}


######## REPEAT FOR FULL_RESULT: ND-DHF AND DF-DHF.




