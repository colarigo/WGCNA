### Langfelder and Horvath WGCNA tutorial adapted ###
# Load the WGCNA package
library(WGCNA);
library(dplyr)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in  data set
expData = read.csv("C_ComBat_seq_vst.csv")
setName = "C_Pv_Trinity"
shortLabels = "C"


#Transpose to gene colname and sample rowname
datExpr0 = as.data.frame(t(expData[,-1]));
names(datExpr0) = expData[,1];
rownames(datExpr0) = names(expData)[-1];

#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes and samples from the data
if (!gsg$allOK)
{
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gsg$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
	if (sum(!gsg$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
	# Remove the offending genes and samples from the data:
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#Cluster sample to check for outlier
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

###To modify if outlier
# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

### Link trait to dataset
traitData = read.csv(list.files(getwd(), "*info.csv"));
###If necessary remove columns
# remove columns that hold information we do not need.
allTraits = traitData[,];
# Form a data frame analogous to expression data that will hold the clinical traits.
#allTraits = traitData
Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
#Trait of interest 
TraitI = as.data.frame(datTraits$Plasmopara)
names(TraitI) = "Plasmopara"
rownames(TraitI) = rownames(datTraits)
collectGarbage();

# Re-cluster samples
#sampleTree2 = hclust(dist(datExpr), method = "average")

save(datExpr, datTraits, file = paste0(setName, "-dataInput.RData"))

enableWGCNAThreads()
lnames = load(file = paste0(setName, "-dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
#pdf(paste0("Plots_scalefree_topology_", setName,".pdf"), wi = 12, he =12)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 8
deepsplit = 0
start = Sys.time()
net = blockwiseModules(datExpr, power = softPower,deepSplit = deepsplit,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       loadTOM = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3,
                       maxBlockSize = ncol(datExpr))
end = Sys.time()
as.numeric(end- start, units= "mins")

table(net$colors)
# If too many modules
# merge = mergeCloseModules(datExpr, net$colors, cutHeight = 0.35, verbose = 3)


# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]];

# Plot the dendrogram and the module colors underneath
pdf(paste0(setName, "_dendro_module.pdf"), width=12, height =9)
plotDendroAndColors(geneTree, moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()


#Alternative dendrogramn with GS

GS = cor(datExpr, TraitI, use = 'p')
colnames(GS) = colnames(TraitI)
# sizeGrWindow(12,9)
pdf(paste0(setName, "_dendro_module_GS.pdf"), width=12, height =9)
plotDendroAndColors(geneTree, cbind(moduleColors, 
numbers2colors(GS,signed=TRUE)), 
c("Modules", spaste("GS_", shortLabels )), 
main = paste("Consensus dendrogram module colors", 
"and gene significance for response to Plasmopara"), 
dendroLabels = FALSE, hang = 0.02, guideHang = 0.05, addGuide = TRUE,
marAll = c(1,9,2,1), cex.lab = 1.2, cex.colorLabels = 1.1)
dev.off()


save(datExpr, datTraits, TraitI, setName, MEs, moduleLabels, moduleColors, geneTree, GS, shortLabels, file = paste0(setName, "-networkConstruction-auto.RData"))


