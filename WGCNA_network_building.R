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

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, TraitI, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf(paste(setName, "module_signif.pdf", sep="_"), width = 6, height = 12)
# sizeGrWindow(24,12)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.8, 3, 2.2));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(TraitI),
yLabels = colnames(MEs),
ySymbols = colnames(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 1,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

dev.off()

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, TraitI))
# Plot the relationships among the eigengenes and the trait
# sizeGrWindow(5,7.5);
pdf(paste0("Plot_", setName, "_eigengenesSignificance.pdf"), wi = 5, he = 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)

# Calculating Module Membership, Gene significance and network Screening
modNames = substring(names(MEs), 3)
datkME=signedKME(datExpr, MEs, outputColumnName="", corFnc = "bicor", corOptions = "use = 'p'")
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(datkME), nSamples));
# names(datkME) = paste("MM.", modNames, sep="");
names(MMPvalue) = paste("p.", modNames, sep="");
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(GS), nSamples));
names(GSPvalue) = paste("p.", names(TraitI), sep="");
NS=networkScreening(y=TraitI[,1], datME=MEs, datExpr=datExpr, oddPower=3, blockSize=42000, minimumSampleSize=4, addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)


# ### define modules of interest by list
# colorlevels = c("grey60","purple", "midnightblue","salmon" )
# nmodules = length(colorlevels)

# Get modules of interest p < 0.05 and > 0.3
nSamples = nrow(datExpr);
moduleTraitCor = cor(MEs, TraitI, use = "p")
modTrait = as.data.frame(moduleTraitCor)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
colnames(moduleTraitPvalue) = "Traitpvalue"
modTrait["value"] = moduleTraitPvalue
passmodules = modTrait %>% slice_max(abs(modTrait[,names(TraitI)]), n = 3) %>% filter (value < 0.05)
colorlevels = gsub("ME", "",rownames(passmodules))
nmodules = length(colorlevels)

if (12 > nmodules & nmodules > 6)
{
# sizeGrWindow(6*nmodules, 12);
pdf(paste(setName,"modules_signif.pdf",sep="_"), width= 24, height = 12)
par(mfrow = c(2,6));
} else {
# sizeGrWindow(6*nmodules, 6);
pdf(paste(setName,"modules_signif.pdf",sep="_"), width= 24, height = 6)
par(mfrow = c(1,nmodules));
}
for (nb in 1:nmodules)
{
column = match(colorlevels[nb], modNames);
moduleGenes = moduleColors==colorlevels[nb];
verboseScatterplot(abs(datkME[moduleGenes, column]),
abs(GS[moduleGenes, 1]),
xlab = paste("Module Membership in", colorlevels[nb], "module"),
ylab = "Gene significance for Plasmopara",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = colorlevels[nb])
abline(v = 0.8, h = 0.6)
}
dev.off()


# Read in the probe annotation
#annot = read.csv(file = list.files( pattern= ".*[Aa]nnot.*"));
probes = names(datExpr)
#probes2annot = match(probes, annot$Gene)

#geneInfo0 = data.frame(Gene = probes, geneSymbol = annot$Annot[probes2annot], moduleColor = moduleColors, GS, GSPvalue, NS$cor.Weighted, NS$p.Weighted)
geneInfo0 = data.frame(Gene = probes, moduleColor = moduleColors, GS, GSPvalue, NS$cor.Weighted, NS$p.Weighted)
old = names(geneInfo0)
modOrder = order(abs(moduleTraitCor));
for (mod in 1:ncol(datkME))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, datkME[, modOrder[mod]], MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$Plasmopara));
geneInfo = geneInfo0[geneOrder, ]
# Full dataframe
write.csv(geneInfo, file = paste(setName, "geneInfo.csv", sep="_"))

# Limited dataframe with interest modules
info0 = geneInfo %>% select(starts_with(old) | ends_with(colorlevels))
info = info0 %>% filter(info0$moduleColor %in% colorlevels)
write.csv(info, file = paste(setName, "modint_info.csv", sep="_"))


save(datkME, MMPvalue, GSPvalue, NS, colorlevels, file= paste0(setName, "_hub_selection.RData"))


