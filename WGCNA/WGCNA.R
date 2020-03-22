library(WGCNA)
options(stringsAsFactors = FALSE)
setwd('D:\\PAAD\\WGCNA')


###### 
univariate_data <- read.table('diffmRNAExp.csv',header = T,row.names = 1,
                              sep = ',')

#univariate_data  <- log2(univariate_data +0.001)
univariate_data <- univariate_data[,-c(1:7)]


###### 
metadata <- data.frame(colnames(univariate_data))

for (i in 1:length(metadata[,1])) {
  num <- as.numeric(as.character(substring(metadata[i,1],14,15)))
  if (num == 1 ) {metadata[i,2] <- "T"}
  if (num != 1) {metadata[i,2] <- "N"}
}

names(metadata) <- c("id","group")
metadata$group <- as.factor(metadata$group)
metadata <- subset(metadata,metadata$group == "T")
exprSet <- univariate_data[,which(colnames(univariate_data) %in% metadata$id)]
colnames(exprSet)  <- substr(x=colnames(exprSet),start = 1,stop = 12)
colnames(exprSet)  <- chartr(old='.',new = '-',x=colnames(exprSet) )
femData = exprSet



#Take a quick look at what is in the data set:
dim(femData)
names(femData)
head(femData)
datExpr0 = as.data.frame(t(femData));




gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}




if (T) {
  sampleTree = hclust(dist(datExpr0), method = "average");
  sizeGrWindow(12,9)
  pdf("Plots_sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
      cex.axis = 1.5, cex.main = 2)
  abline(h = 80, col = "red");
  clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
  dev.off()
  table(clust)
}


keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


traitData = read.csv("survival.csv");

traitData$T_stage[grepl('T1',traitData$TNM)] <- 1
traitData$T_stage[grepl('T2',traitData$TNM)] <- 2
traitData$T_stage[grepl('T3',traitData$TNM)] <- 3
traitData$T_stage[grepl('T4',traitData$TNM)] <- 4
traitData$T_stage[grepl('TX',traitData$TNM)] <- 0
traitData$T_stage[traitData$T_stage == NA] <- 0

traitData$N_stage[grepl('N1',traitData$TNM)] <- 1
traitData$N_stage[grepl('N2',traitData$TNM)] <- 2
traitData$N_stage[grepl('N3',traitData$TNM)] <- 3
traitData$N_stage[grepl('N0',traitData$TNM)] <- 0
traitData$N_stage[grepl('NX',traitData$TNM)] <- 4

traitData$M_stage[grepl('M1',traitData$TNM)] <- 1
traitData$M_stage[grepl('M2',traitData$TNM)] <- 2
traitData$M_stage[grepl('M3',traitData$TNM)] <- 3
traitData$M_stage[grepl('M0',traitData$TNM)] <- 0
traitData$M_stage[grepl('MX',traitData$TNM)] <- 4
traitData$M_stage[traitData$M_stage == NA] <- 0
traitData$TNM <- NULL

traitData$stage[traitData$Stage == 'Stage I'] <- 1
traitData$stage[traitData$Stage == 'Stage IA'] <- 1
traitData$stage[traitData$Stage == 'Stage IB'] <- 1
traitData$stage[traitData$Stage == 'Stage II'] <- 2
traitData$stage[traitData$Stage == 'Stage IIA'] <- 2
traitData$stage[traitData$Stage == 'Stage IIB'] <- 2
traitData$stage[traitData$Stage == 'Stage III'] <- 3
traitData$stage[traitData$Stage == 'Stage IIIA'] <- 3
traitData$stage[traitData$Stage == 'Stage IIIB'] <- 3
traitData$stage[traitData$Stage == 'Stage IV'] <- 4
traitData$stage[traitData$Stage == 0] <- 0
traitData$Stage <- NULL

traitData$cancer_status[traitData$cancer_status == 'TUMOR FREE'] <- 2
traitData$cancer_status[traitData$cancer_status == 'WITH TUMOR'] <- 1
traitData$cancer_status[traitData$cancer_status == '0'] <- 0

traitData$Age[traitData$Age < 60] <- 0
traitData$Age[traitData$Age >= 60] <- 1

traitData$OS.Time <- NULL
traitData$cancer_status <- as.numeric(traitData$cancer_status)

names(traitData)





allTraits <- traitData
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Barcode);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();





pdf("Sample dendrogram and trait heatmap.pdf",height=10,width=20)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,cex.dendroLabels = 0.7,
                    cex.colorLabels = 1.5,
                    groupLabels = names(datTraits), 
                    main = "")
dev.off()




powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)# 设定窗口大小
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




softPower = 5;

adjacency = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 30;

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")





MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);

METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2

abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors

colorOrder = c("grey", standardColors(50));

moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;




nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);




sizeGrWindow(10,6)

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


pdf('plot1.pdf',10,10)
par(mar = c(10, 10, 4, 4));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               font.lab.x = 2,
               cex.lab = 1.5,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



weight = as.data.frame(datTraits$OS);
names(weight) = "OS"

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
write.csv(GSPvalue, file = "GSPvalue")





module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
names(datExpr)
names(datExpr)[moduleColors=="green"]




nSelect = 400

set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];

selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];

sizeGrWindow(9,9)

plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")





MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

weight = as.data.frame(datTraits$OS);
names(weight) = "OS"

MET = orderMEs(cbind(MEs, weight))

sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)


sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)





TOM = TOMsimilarityFromExpr(datExpr, power = 5);

modules = c('green');

probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
