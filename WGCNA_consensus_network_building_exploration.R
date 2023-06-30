library(WGCNA)
library(dplyr)
library(purrr)
library(expss, lib.loc= perlib)
library(mgsub, lib.loc= perlib)
source("/path/to/function.R")
options(stringsAsFactors = FALSE);

setName = "Full_Pv_all_ComBat_seq_vst"
dd = list.files(pattern= "*.txt")
Dataname = lapply(dd, read.csv)
# TofData = read.csv("./Tof_ComBat_seq_vst1.txt")
#Rpv12 = read.csv("Rpv12_ComBat_seq_vst.txt")
# PNThData = read.csv("./PNTh_vst.txt")
#TB = read.csv("./TB_ComBat_seq_vst.txt")
#TM = read.csv("./TM_ComBat_seq_vst.txt")
#TP = read.csv("./TP_ComBat_seq_vst.txt")
#PN = read.csv("./PNTh_sens_vst.txt")
#TM = TM[,-13]

#Check if dataframe of different size
if (length(Reduce(intersect, lapply(Dataname, nrow))) < 1)
{
template = Reduce(intersect, lapply(Dataname, function(x) x[,1]))
Dataname = map(Dataname, function(x) filter(x, x[,1] %in% template))
}

nSets = length(dd)
setLabels = c("V.Amurensis","TBianca", "TMgaloblishvilii", "TPinot")
setNames = c("Rpv_Vv_VAG","TB_Rpv3_Trinity", "TM_Trinity", "TP")
# shortLabels = c("Tof", "Rpv12", "PNTh")
shortLabels = c("TB", "TM", "TP", "Rpv", "PN")
#Dataname = list(TB, TM,TP, Rpv12, PN)
multiExpr = list()
for (set in 1:nSets)
{
multiExpr[[set]] = list(data = as.data.frame(t(Dataname[[set]][,-1])))
names(multiExpr[[set]]$data) = Dataname[[set]][,1]
rownames(multiExpr[[set]]$data)= names(Dataname[[set]][,-1])
}
exprSize = checkSets(multiExpr)

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
# Print information about the removed genes:
if (sum(!gsg$goodGenes) > 0)
printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
collapse = ", ")))
for (set in 1:exprSize$nSets)
{
if (sum(!gsg$goodSamples[[set]]))
printFlush(paste("In set", setLabels[set], "removing samples",
paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
# Remove the offending genes and samples
multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
}
# Update exprSize
exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets)
{
sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

#pdf(file = "Plots_SampleClustering.pdf", width = 12, height = 12);
sizeGrWindow(12,12)
par(mfrow=c(3,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
#dev.off();

# Remove outliers, might need some fine-tuning for specific data
cutHeights = Reduce(c,lapply(sampleTrees, function(x) unname(quantile(x$height, probs=0.8)+ sd(x$height)*2)))
# Choose the "base" cut height depending on max height of dendrogram
#baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
#cutHeights = rep(baseHeight,nSets);
# Re-plot the dendrograms including the cut lines
sizeGrWindow(12,12)
par(mfrow=c(3,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();

for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels = cutree(sampleTrees[[set]], h = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

traitData = read.csv(list.files(getwd(), "*_batch.csv"))
allTraits = traitData

Traits = list();
for (set in 1:nSets)
{
setSamples = rownames(multiExpr[[set]]$data);
traitRows = match(setSamples, allTraits$Sample);
Traits[[set]] = list(data = allTraits[traitRows, -1]);
rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
### One trait of interest
TraitI = list();
for (set in 1:nSets)
{
setSamples = rownames(multiExpr[[set]]$data);
traitRows = match(setSamples, rownames(Traits[[set]]$data));
TraitI[[set]] = list(data = as.data.frame(Traits[[set]]$data[traitRows, "Plasmopara"]));
rownames(TraitI[[set]]$data) = rownames(Traits[[set]]$data)
names(TraitI[[set]]$data) = "Plasmopara"
}
collectGarbage();

nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

#save(multiExpr, Traits, TraitI, nGenes, nSamples, setLabels, shortLabels, exprSize, nSets, file = "Consensus_vst-dataInput.RData");

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
{
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, RsquaredCut = 0.8,
verbose = 5));
collectGarbage();
}
# Retrieve the softpower with > 0.8 and common in  all data set
softP = Reduce(intersect,Reduce(c,lapply(powerTables, function(x) x$data$fitIndices %>% filter(SFT.R.sq > 0.8, slope < 0 , truncated.R.sq >0.8)%>% select(Power))))
softP = softP[1:6]
softP
# Plot the results:
#Adapted number of colors depending on nSets
colors = rainbow(nSets)
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
for (col in 1:length(plotCols))

{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data$fitIndices[, plotCols[col]], na.rm = TRUE);
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data$fitIndices[, plotCols[col]], na.rm = TRUE);
}
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
#sizeGrWindow(8, 6)
pdf(file = paste0("Plots_scaleFreeAnalysis_", setName, ".pdf"), wi = 12, he = 9)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data$fitIndices[,1], -sign(powerTables[[set]]$data$fitIndices[,3])*powerTables[[set]]$data$fitIndices[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col]);
addGrid();
}
if (col==1)
{
text(powerTables[[set]]$data$fitIndices[,1], -sign(powerTables[[set]]$data$fitIndices[,3])*powerTables[[set]]$data$fitIndices[,2],
labels=powers,cex=cex1,col=colors[set]);
} else
text(powerTables[[set]]$data$fitIndices[,1], powerTables[[set]]$data$fitIndices[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set]);
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
} else
legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();

softPower = softP[1]
deepsplit = 0
net = blockwiseConsensusModules(multiExpr, power = softPower, minModuleSize = 30, deepSplit = deepsplit,
                                pamRespectsDendro = FALSE,
                                mergeCutHeight = 0.25, numericLabels = TRUE,
                                minKMEtoStay = 0,
                                verbose = 5,maxBlockSize = nGenes,
                                TOMType = "signed",
                                saveConsensusTOMs = FALSE)

consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];

pdf(file = paste0("Plots_", setName,"_Dendrogram_auto.pdf"), wi = 12, he = 9)
plotDendroAndColors(consTree, moduleColors,
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Consensus gene dendrogram and module colors")
dev.off()

#Gene significance
GS.TC = multiData.mapply(bicor, multiExpr, TraitI, MoreArgs = list( use = "p"), mdmaSimplify = TRUE)
GS.TC = GS.TC[, , , drop = TRUE]
GScolors = numbers2colors(GS.TC, signed = TRUE, colors = blueWhiteRed(100))

pdf(paste0("Plot_", setName,"_module_GS.pdf"), width = 12, height = 9)
# sizeGrWindow(12,9)
plotDendroAndColors(net$dendrograms[[1]], cbind(labels2colors(moduleLabels), GScolors), 
c("Consensus\nModules", setLabels), rowText = spaste(moduleLabels, ": ", labels2colors(moduleLabels)),
textPositions = nSets + 1, main = "Consensus dendrogram ",
dendroLabels = FALSE, hang = 0.01, guideHang = 0.03, abHeight = 0.997, adddGuide = TRUE, marAll = c(1,5,1,1), cex.lab= 1.2, cex.colorLabels = 1.1)

dev.off()

save(nSets, setLabels, shortLabels, multiExpr, softPower, consMEs, moduleLabels, moduleColors, consTree, Traits, TraitI, setName, deepsplit,
 file = paste0(setName, "-NetworkConstruction-auto.RData"))

loading = load(file = list.files(pattern = "*NetworkConstruction-auto.RData$"))

### Second Part
exprSize = checkSets(multiExpr);
# nSets = exprSize$nSets;
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
nTraits = checkSets(TraitI)$nGenes

moduleTraitCor = list();
moduleTraitPvalue = list();

MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");

# Calculate the correlations
for (set in 1:nSets)
{
moduleTraitCor[[set]] = cor(consMEs[[set]]$data, TraitI[[set]]$data, use = "p");
moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = eval(parse(text = paste0("moduleTraitCor[[", 1:nSets,"]] < 0 ", collapse ="& ")))
consensusCor[negative] = eval(parse(text = paste0("pmax(",paste0("moduleTraitCor[[", 1:nSets, "]][negative]", collapse= ", "), ")")))
consensusPvalue[negative] = eval(parse(text = paste0("pmax(",paste0("moduleTraitPvalue[[", 1:nSets, "]][negative]", collapse= ", "), ")")))
# # Find consensus positive correlations
positive = eval(parse(text = paste0("moduleTraitCor[[", 1:nSets,"]] > 0 ", collapse ="& ")))
consensusCor[positive] = eval(parse(text = paste0("pmax(",paste0("moduleTraitCor[[", 1:nSets, "]][positive]", collapse= ", "), ")")))
consensusPvalue[positive] = eval(parse(text = paste0("pmax(",paste0("moduleTraitPvalue[[", 1:nSets, "]][positive]", collapse= ", "), ")")))

modsign = Reduce(intersect,lapply(lapply(moduleTraitPvalue, function(x) as.data.frame(x) %>% filter(Plasmopara < 0.1)), function(x) rownames(x)))
if (length(modsign) == 0)
{}
### Plot each set module Trait relationships 
for (set in 1:nSets)
{
pdf(file = paste0("Plots_ModuleTraitRelationships_ds", deepsplit, shortLabels[[set]], ".pdf", collapse = "_"), wi = 6, he = 12);
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(TraitI[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships in", setLabels[set]))
dev.off();
}

###Plot consensus module Trait relationship
textMatrix = paste(signif(consensusCor, 2), "\n(",
signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
# sizeGrWindow(10,7)
pdf(file = paste0("Plots_ModuleTraitRelationships_", setName,".pdf"), wi = 6, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
xLabels = names(TraitI[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Consensus module--trait relationships across\n",
paste(setLabels, collapse = " and ")))
dev.off()


# Recalculate MEs
MEs = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
nMEs = ncol(MEs[[1]]$data)
# Module Eigengene significance as list
### Think to load function.R and networkFunctions-extra-05.R
MEsignif0 = multiData.mapply(bicorAndPvalue, TraitI, MEs, mdmaSimplify = FALSE)
# Re-format
MEsignif = pMEsignif = array(0, dim = c(nSets, 1 , nMEs))
for (set in 1:nSets)
{
MEsignif[set,,] = MEsignif0[[set]]$data$bicor;
pMEsignif[set,,] = MEsignif0[[set]]$data$p;
}
# For one trait
MEMAs = metaAnalysis(MEs, TraitI, useRankPvalue = FALSE, corFnc = bicor)

# for several traits
# MEMAs = list()
# colName = "Z.RootDoFWeights"
# for (t in 1:nTraits)
# {
# printFlush(paste("Working on ", traitNames[t]))
# MEMAs[[t]]= metaAnalysis(MEs,multiData.subset(Traits, t), useRankPvalue = FALSE, corFnc = bicor)
# }

colName = "Z.RootDoFWeights"
labLevels = sort(unique(moduleLabels[moduleLabels!=0]))
modSizes = table(moduleLabels[moduleLabels!=0])
xLabels = spaste("ME", labels2colors(labLevels))
xSymbols = labels2colors(labLevels)
# text matrix
mat = MEsignif
dim(mat)= c(nSets*nTraits, nMEs)
textMat = spaste(round(MEsignif, 2), "\n", signif(pMEsignif, 1))
dim(textMat) = dim(mat)
#  Maximum value for the color legends
max = max(abs(MEsignif), na.rm = TRUE)

### Meta-analysis of module eigengene significance
# For one trait
index = ( 1) : (1*nSets)
Z = MEMAs [ , match(colName, colnames(MEMAs))]
Z1 = Z/max(abs(Z))*max
mat1 = rbind(Z1, mat[index,])
pColName = sub("Z", "p", colName)
p = MEMAs[, match(pColName, colnames(MEMAs))]
text1 = rbind(spaste(signif(Z,2), "\n", signif(p,1)),textMat[index,])
yLabels1 = c("meta-analysis Z and p ", setLabels)
pdf(file=paste0("Plots_", setName, "_eigengenesSignificance.pdf"),  width = 6, height = 18)
par(mar = c( 7,7,2.2,2))
labeledHeatmap(t(mat1),
	xLabels = yLabels1,
	yLabels = xLabels,
	ySymbols = xSymbols,
	textMat = t(text1),
	cex.tex  = 1,
	main = paste0("Module eigegene significance for", colnames(TraitI[[1]])),
	colors = blueWhiteRed(100),
	zlim = c(-max,max),  setStdMargins = FALSE)
dev.off()

###################################################################################################################################################################################

																		############################
																		### Network exploration ####
																		############################

### Find hub genes
annot = read.csv(file = list.files(getwd(), pattern = "*annotation.csv"));
# Match probes in the data set to the probe IDs in the annotation file
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot$Gene)

# Find the top 3 Z.meta significative modules
# Get vector of top 3 modules
col = match(colName, colnames(MEMAs))
bms = as.vector((MEMAs %>% slice_max(abs(MEMAs[,col]), n= 3) %>% dplyr::select(ID))[,1])
modCol = labels2colors(match(gsub("ME", "", bms), labLevels))

save(colName, MEMAs, labLevels, col, bms, modCol, file = paste0(setName, "_module_cytoscape.RData"))

																		###########################################	
																		### Simple meta-analysis for GS and kME ###
																		###########################################

#Calculating Gene significance and module membership
GS = list();
kME = list();
NS = list();
for (set in 1:nSets)
{
GS[[set]] = bicorAndPvalue(multiExpr[[set]]$data, TraitI[[set]]$data);
kME[[set]] = bicorAndPvalue(multiExpr[[set]]$data, MEs[[set]]$data);
NS[[set]] = networkScreening(y=TraitI[[set]]$data[,1], datME = MEs[[set]]$data, datExpr=multiExpr[[set]]$data,
 oddPower = 3, blockSize = 24000, minimumSampleSize = 4, addMEy = TRUE, removeDiag = FALSE, weightESy=0.5, corFnc = "bicor")
}

### Create a string like (GS[[1]]$Z + GS[[2]]$Z + GS[[3]]$Z)/sqrt(nSets)
### eval(parse) interpret the string as variables
GS.metaZ = eval(parse(text = paste0("(", paste0("GS[[", 1:nSets, "]]$Z", collapse=" + "), ")", "/sqrt(nSets)")))
kME.metaZ = eval(parse(text = paste0("(", paste0("kME[[", 1:nSets, "]]$Z", collapse=" + "), ")", "/sqrt(nSets)")))
NS.metaZ = eval(parse(text = paste0("(", paste0("NS[[", 1:nSets, "]]$Z.Weighted", collapse=" + "), ")", "/sqrt(nSets)")))
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);
NS.metaP = 2*pnorm(abs(NS.metaZ), lower.tail = FALSE);

#Create table for output GS
GSmat = eval(parse(text = paste0("rbind(", paste0(paste0("GS[[", 1:nSets, "]]$bicor", collapse=", "), ", ", paste0("GS[[", 1:nSets, "]]$p" , collapse=', ')), ", GS.metaZ, ", "GS.metaP)")))
nTraits = length(colnames(TraitI[[1]]$data))
traitNames = colnames(TraitI[[1]]$data)
dim(GSmat) = c(nGenes, (nSets*2+2)*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(c(paste0("GS.set", 1:nSets), paste0("p.GS.set", 1:nSets),  "Z.GS.meta.", "p.GS.meta"),
rep(traitNames, rep((nSets*2+2), nTraits)))

# Same code for kME:
kMEmat = eval(parse(text = paste0("rbind(", paste0(paste0("kME[[", 1:nSets, "]]$bicor", collapse=", "), ", ", paste0("kME[[", 1:nSets, "]]$p" , collapse=', ')), ", kME.metaZ, ", "kME.metaP)")))
MEnames = colnames(MEs[[1]]$data);
nMEs = checkSets(MEs)$nGenes
dim(kMEmat) = c(nGenes, (nSets*2+2)*nMEs)
# rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(c(paste0("kME.set",1:nSets), paste0("p.kME.set", 1:nSets),"Z.kME.meta.", "p.kME.meta"),
rep(MEnames, rep((nSets*2+2), nMEs)))

# Same code for NS
NSmat = eval(parse(text = paste0("rbind(", paste0(paste0("NS[[", 1:nSets, "]]$bicor.Weighted", collapse=", "), ", ", paste0("NS[[", 1:nSets, "]]$p.Weighted" , collapse=', ')), ", NS.metaZ, ", "NS.metaP)")))
traitNames = colnames(TraitI[[1]]$data)
dim(NSmat) = c(nGenes, (nSets*2+2)*nTraits)
rownames(NSmat) = probes;
colnames(NSmat) = spaste(c(paste0("NS.set",1:nSets), paste0("p.NS.set", 1:nSets),"Z.NS.meta.", "p.NS.meta"),
rep(traitNames, rep((nSets*2+2), nTraits)))

																		###############################
																		#### Complete network info ####
																		

#Write output with GS and kME for consensus module
info = data.frame(Probe = probes, ModuleColor = labels2colors(moduleLabels), GSmat, kMEmat, NSmat);
#Filter non significative genes 
#Limit to top 3 modules' columns
pmod = c("p.GS.metaPlasmopara", spaste("p.kME.meta", bms))
# Alt select all p.value meta columns
# pmod = c("p.GS.metaPlasmopara", spaste("p.kME.metaME", seq(from=1, to= length(labLevels))))
#info = info %>% filter(across(pmod,lt(0.05)))
colnames(info) <- mgsub::mgsub(colnames(info),paste0("ME", 1:nMEs), paste0("ME", labels2colors(match(1:nMEs, labLevels))))
write.csv(info, file = paste0(setName,"-simple_meta_long.csv"),row.names = FALSE, quote = FALSE);

																		##############################
																		#### Limited network info ####


# Short version
# Filter columns corresponding to top 3 modules 
subkMEmat = as.data.frame(kMEmat) %>% dplyr::select(matches(paste0(bms, "$", collapse="|")))
info0 = data.frame(Probe = probes, ModuleColor = labels2colors(moduleLabels), GSmat, subkMEmat, NSmat)
#pmod0 = c("p.GS.metaPlasmopara", spaste("p.kME.meta", bms))
#info0 = info0 %>% filter(across(pmod0,lt(0.05)))
colnames(info0) <-  mgsub::mgsub(colnames(info0), bms, paste0("ME", labels2colors(match(gsub("ME", "", bms), labLevels))))
write.csv(info0, file = paste0(setName,"-simple_meta_short.csv"),row.names = FALSE, quote = FALSE);

																		###################################
																		#### Alternative meta-analysis ####
																		###################################

#Marginal meta-analysis
metaGS = metaAnalysis(multiExpr, TraitI, useRankPvalue = FALSE, corFnc = bicor)

#Module membership analysis
bestMElabel = bestMEindex = bestMEcolor = rep(0, nSets)
bestME = list()
for (a in 1:nSets)
{
bestMEindex[a] = which.max(abs(MEMAs[,col]))
bestMEcolor[a] = labels2colors(match(bestMEindex[a], labLevels))
bestMElabel[a]= MEMAs$ID[bestMEindex[a]]
bestME[[a]] = MEs[[a]]$data[bestMElabel[a]]
}
allsame = function(x) length(unique(x)) == 1
if (allsame(bestMEindex) == TRUE) {
conKME = consensusKME(multiExpr, moduleLabels, consensusQuantile = 0.25, corAndPvalue = bicorAndPvalue,	useModules = bestMEindex[1])
} 
# Limit table to best module
colkeep = c("ID",paste0("consensus.kME", bestMEindex[1]),paste0("meta.Z.RootDoFWeights.kME",bestMEindex[1]), 
	paste0("meta.p.RootDoFWeights.kME",bestMEindex[1]), paste0(".Set_", 1:nSets), paste0(".Set_", 1:nSets, ".1"))
collectGarbage();

# Code for top modules
bestmodules = gsub("ME", "", bms)
# Get module colors for better visualization
conKME = consensusKME(multiExpr, moduleLabels, consensusQuantile = 0.25,corAndPvalue = bicorAndPvalue, useModules = bestmodules)
collectGarbage();
colkeep = c("ID",paste0("consensus.kME", bestmodules),paste0("meta.Z.RootDoFWeights.kME",bestmodules), paste0("meta.p.RootDoFWeights.kME",bestmodules),
 c(with(expand.grid(a=bestmodules, i=1:nSets), paste0("kME", a, ".Set_", i)), with(expand.grid(a=bestmodules, i=1:nSets), paste0("p.kME", a, ".Set_", i))))

#For all modules
conKME = consensusKME(multiExpr, moduleLabels, consensusQuantile = 0.25,corAndPvalue = bicorAndPvalue)
colkeep = c("ID",paste0("consensus.kME", 1:nMEs),paste0("meta.Z.RootDoFWeights.kME",1:nMEs), 
	paste0("meta.p.RootDoFWeights.kME",1:nMEs), c(with(expand.grid(a = 1:nMEs, b = 1:nSets), paste0("kME", a, ".Set_", b))),
	 c(with(expand.grid(a = 1:nMEs, b = 1:nSets), paste0("p.kME", a, ".Set_", b))))

#colkeepGS = c("ID", "Z.RootDoFWeights", "p.RootDoFWeights", paste0("cor.in.Set_", 1:nSets), paste0("p.Student.in.Set_", 1:nSets))
# If corPearson used
colkeepGS = c("ID", "Z.RootDoFWeights", "p.RootDoFWeights", paste0("corPearson.Set_", 1:nSets), paste0("pvalueStudent.0.vs.1.Set_", 1:nSets))
conkME = conKME[, colkeep]
metaGS = metaGS[, colkeepGS]
# Rename columns
modColors = labels2colors(match(1:nMEs, labLevels))
# for best module
#colnames(conkME) <- mgsub::mgsub(colnames(conkME),c(paste0("kME",bestMEindex[1]),paste0(".Set_", 1:nSets), paste0(".Set_", 1:nSets, ".1")), c(paste0("kME.",bestMEcolor[1]),paste0("kME.Set", 1:nSets), paste0("pkME.Set", 1:nSets)))
# for top modules
#colnames(conkME) <- mgsub::mgsub(colnames(conkME),c(paste0("kME", bestmodules)), c(paste0("kME.", labels2colors(match(bestmodules, labLevels)))))
# for all modules 
colnames(conkME) <- mgsub::mgsub(colnames(conkME),c(paste0("kME",1:nMEs), c(with(expand.grid(a = 1:nMEs, b = 1:nSets), paste0("kME", a, ".Set_", b))), c(with(expand.grid(a = 1:nMEs, b = 1:nSets), paste0("p.kME", a, ".Set_", b)))),
 c(paste0("kME.",modColors),c(with(expand.grid(a = modColors, b = 1:nSets), paste0("kME", a, ".Set_", b))), c(with(expand.grid(a = modColors, b = 1:nSets), paste0("p.kME", a, ".Set_", b)))))
# Merge metaGS and conkME dataframe
merged = merge(metaGS, conkME, by = "ID")

# Format dataframe
meta_info = data.frame(Probe = probes, ModuleColor = labels2colors(moduleLabels), merged);
# Limit dataframe to significative genes 
#meta_info = meta_info %>% filter(p.RootDoFWeights < 0.05, eval(parse(text = paste0("meta.p.RootDoFWeights.kME.", bestMEcolor[1]))) < .05)
write.csv(meta_info, file = "metaAnalysis_consensusGS_kME_all_mod.csv",row.names = FALSE, quote = FALSE);

