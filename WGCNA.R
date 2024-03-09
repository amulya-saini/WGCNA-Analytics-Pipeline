# loading the required packages
# install.packages("WGCNA")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")

library("WGCNA")
library("tidyverse")
library("dplyr")
library("magrittr")
library("tidyr")

# loading the expression data
exp_data <- read.csv("C:\\Users\\saini\\Downloads\\HW03_expression.csv", header = TRUE)
View(exp_data)
dim(exp_data)

# loading the Traits data
traits_data <- read.csv("C:\\Users\\saini\\Downloads\\HW03_traits.csv", header = TRUE)
View(traits_data)
dim(traits_data)

######## Removing NaN values and Pre-processing #########

# dropping the NA values from traits data
traits_data <- na.omit(traits_data)
dim(traits_data)

# setting the "ID" column as row names
traits_data <- transform(traits_data, row.names = ID)

# removing the original column "ID"
traits_data <- traits_data[, -1]
View(traits_data)

# extracting the row names from the dataset
traits_data_rownames <- rownames(traits_data)

# getting the column names of exp_data besides gene_symbol
exp_data_colnames <- colnames(exp_data)[-1]

# extracting the column names that are not present in traits data
columns_to_drop <- setdiff(exp_data_colnames, traits_data_rownames)

# dropping the columns that we donot have traits data for
exp_data <- exp_data[, !(colnames(exp_data) %in% columns_to_drop)]
dim(exp_data)
View(exp_data)

# Drop rows with NA values
exp_data <- na.omit(exp_data)

################## Checking if the data is normalized #####################

# Get the column names 
variables <- names(exp_data)

# Loop through each variable and create a histogram
for (variable in variables) {
  
  # Exclude the "gene_symbol" column
  if (variable != "gene_symbol") {
    
    # Check if the variable is numeric
    if (is.numeric(exp_data[[variable]])) {
      
      # Create a new plot for each numeric variable
      hist(exp_data[[variable]], main = paste("Histogram of", variable), xlab = variable)
    } else {
      
      # Print a message if the variable is not numeric
      cat("Variable", variable, "is not numeric and cannot be plotted as a histogram.\n")
    }
  }
}

# From the histogram plots, the data looks normalized

##################### Tranposing the data ###############

# making the gene_symbol values unique using make.names
exp_data$gene_symbol <- make.names(exp_data$gene_symbol, unique = TRUE)

# setting the first column as row names
row.names(exp_data) <- exp_data[, 1]

# Remove the first column from the dataframe
exp_data <- exp_data[, -1]
View(exp_data)

# Transposing the data
input_mat = t(exp_data)
View(input_mat)

##################### End of Pre-processing ##############

# Choosing a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20 , by = 2))
powers

# Call the network topology analysis function
sft <- pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)
sft_data <- sft$fitIndices
View(sft_data)

# we need a soft threshold that has maximum rsquare and minimum mean connectivity
###### Plotting  the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))

# Setting some parameters
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# the selected soft threshold power
picked_power = 10

#################### Gene clustering based on TOM based dissimilarity ############### 

# Turn adjacency into topological overlap matrix (TOM)
adjacency <- adjacency(input_mat, power = picked_power)
TOMadj <- TOMsimilarity(adjacency)


dissTOMadj <- 1- TOMadj

# Clustering using TOM
# Call the hierarchical clustering function 
hclustGeneTree <- hclust(as.dist(dissTOMadj), method = "average")

# Plot the resulting clustering tree (dendogram)
sizeGrWindow(12, 9)
plot(hclustGeneTree, xlab = "", sub = "", 
     main = "Gene Clustering on TOM-based disssimilarity", 
     labels = FALSE, hang = 0.04)

######################## Creating Modules ###########################

# determining the minimum size a module should have. Setting it higher means smaller modules are not considered
minModuleSize <- 30

# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, 
                             distM = dissTOMadj,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Converting numeric lables obtained from dynamic tree cut into colors
# makes it easier to visulaize and interpret
dynamicColors <- labels2colors(dynamicMods)

# printing the dynamic colors and number of gene in each module
table(dynamicColors)

# Plotting the dendrogram with the identified modules highlighted by different colors
# the modules are not merged yet
sizeGrWindow(8,6)
plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
# epigengenes represent the gene expression of their module across all samples
dynamic_MEList <- moduleEigengenes(input_mat, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes

# Calculating dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs)

# creating a new hierarchial clustering tree based on the Dissimilarity between module eigengenes
dynamic_METree <- hclust(as.dist(dynamic_MEDiss))

# Plot the hclust
sizeGrWindow(7,6)
plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
     xlab = "", sub = "")

######################## MERGE SIMILAR MODULES ############################

# defining a module dissimilarty threshold below which the modulesare merged
dynamic_MEDissThres <- 0.25

# Plot the cut line
#abline(h = dynamic_MEDissThres, col = "red")

# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(input_mat, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)

# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors

# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs
mergedMEs

table(dynamic_mergedColors)

# plotting a dendrogram with merged module colors
plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


############################ Relating to traits ##########################

# defining the number of genes and samples
n_genes <- ncol(input_mat)
n_samples <- nrow(input_mat)

# recalculate MEs with color labels
MEs0 <- moduleEigengenes(input_mat, dynamic_mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
names(MEs) <- substring(names(MEs), 3)

temp_cor <- cor       
cor <- WGCNA::cor
cor <- temp_cor

moduleTraitCor <- cor(MEs, traits_data, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, n_samples)

############################ HEATMAP ###################################
#sizeGrWindow(10,6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(traits_data),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

