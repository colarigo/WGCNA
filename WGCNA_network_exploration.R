perlib = "/home/genome/R/x86_64-pc-linux-gnu-library/3.6"
.libPaths(c(perlib, .libPaths()))
# Load the WGCNA package
library(WGCNA)
library(dplyr)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
lnames = load(file = list.files(getwd(), "-networkConstruction-auto.RData$"));

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