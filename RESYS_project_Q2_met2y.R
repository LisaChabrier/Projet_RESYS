library(rstudioapi) # Set working directory
library(RTN) # Infer regulatory network based on the downloaded gene expression data
library(snow) # Parallel processing in RTN
library(RedeR) # Better visualization of the network
library(CoRegNet) # Get the human TFs
library(lumi) # Convert nuID > entrezID
library(lumiHumanIDMapping) # Convert nuID > entrezID
library(org.Hs.eg.db) # Convert TFs SYMBOLS into ENTREZ

data("HumanTF_entrezgene")

# Design matrix : G0 = no mestatases.

#### Getting the path of your current open file ####
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())



################################################################
####   Loading & preprocessing the data  ####
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE21257", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10295", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "X1XXX1XXX1XXXXXXX1X110000X1X10X10X000X0100X1XX0000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }



################################################################
####    Network inference     ####
################################################################

#### Preprocessing ####
exEntrez <- ex
rowAnnotation <- as.matrix(nuID2EntrezID(nuID = rownames(ex), lib.mapping = 'lumiHumanIDMapping'))
rowAnnotation <- cbind(rowAnnotation, rownames(rowAnnotation))
rowAnnotation <- rowAnnotation[which(rowAnnotation[,1] != ""),] # Filter nuID without entrezID
colnames(rowAnnotation) <- c("entrez", "nuID")
rowAnnotation <- rowAnnotation[!duplicated(rowAnnotation[,1]),]

exEntrez <- cbind(exEntrez, as.numeric(nuID2EntrezID(nuID = rownames(ex), lib.mapping = 'lumiHumanIDMapping')))
exEntrez <- subset(exEntrez, exEntrez[,31] %in% rowAnnotation[,1]) # Subset exEntrez with nuID in both exEntrez & nuIDtoEntrez
rownames(exEntrez) <- as.numeric(exEntrez[,31])
exEntrez <- exEntrez[,-31]
exEntrez <- exEntrez[!duplicated(rownames(exEntrez)),] # Remove duplicated rownames

#### TFs from immune cell types of interest ####
macrophage <- c("ATM", "JUN", "CEBPB", "CREB1", "IRF1", "IRF4", "IRF5", 
                "IRF7", "IRF8", "MKI67", "MNDA", "NFKB1", "TP53", "PCNA", 
                "PPARG", "SPI1", "SSRP1", "STAT1", "STAT3", "STAT6", "XBP1")
macrophage <- mapIds(org.Hs.eg.db, macrophage, "ENTREZID", "SYMBOL")

Bcell <- c("BACH2", "BCL6", "PRDM1", "IRF4", "IRF8", "PAX5", "SPI1", 
           "STAT5A", "XBP1")
Bcell <- mapIds(org.Hs.eg.db, Bcell, "ENTREZID", "SYMBOL")

Tcell <- c("IKZF3", "ATM", "BATF", "BCL6", "PRDM1", "FOXP3", "GATA3",
           "IKZF2", "IKZF1", "IRF4", "IRF7", "MKI67", "NFKB1", "TP53",
           "PCNA", "SSRP1", "STAT1", "STAT4", "TBX21", "TCL1A", "ZBTB7B",
           "ZAP70")
Tcell <- mapIds(org.Hs.eg.db, Tcell, "ENTREZID", "SYMBOL")

#### Transcriptional Network Inference ####
# For the macrophage TFs
rtni_macro <- tni.constructor(expData = as.matrix(exEntrez),
                              regulatoryElements = macrophage,
                              rowAnnotation = rowAnnotation)
rtni_macro <- tni.permutation(rtni_macro, nPermutations = 100)
rtni_macro <- tni.bootstrap(rtni_macro)
rtni_macro <- tni.dpi.filter(rtni_macro)
tni.regulon.summary(rtni_macro)
g_macro <- tni.graph(rtni_macro, tfs = macrophage)

# Visualization on RedeR
rdp <- RedPort()
calld(rdp, checkcalls=T)
addGraph(rdp, g_macro, layout=NULL)
addLegend.color(rdp, g_macro, type="edge")
addLegend.shape(rdp, g_macro)
relax(rdp, ps = TRUE)

# For the B cell TFs
rtni_Bcell <- tni.constructor(expData = as.matrix(exEntrez),
                              regulatoryElements = Bcell,
                              rowAnnotation = rowAnnotation)
rtni_Bcell <- tni.permutation(rtni_Bcell, nPermutations = 100)
rtni_Bcell <- tni.bootstrap(rtni_Bcell)
rtni_Bcell <- tni.dpi.filter(rtni_Bcell)
tni.regulon.summary(rtni_Bcell)
g_Bcell <- tni.graph(rtni_Bcell, tfs = Bcell)

# Visualization on RedeR
rdp <- RedPort()
calld(rdp, checkcalls=T)
addGraph(rdp, g_Bcell, layout=NULL)
addLegend.color(rdp, g_Bcell, type="edge")
addLegend.shape(rdp, g_Bcell)
relax(rdp, ps = TRUE)

# For the T cell TFs
rtni_Tcell <- tni.constructor(expData = as.matrix(exEntrez),
                              regulatoryElements = Tcell,
                              rowAnnotation = rowAnnotation)
rtni_Tcell <- tni.permutation(rtni_Tcell, nPermutations = 100)
rtni_Tcell <- tni.bootstrap(rtni_Tcell)
rtni_Tcell <- tni.dpi.filter(rtni_Tcell)
tni.regulon.summary(rtni_Tcell)
g_Tcell <- tni.graph(rtni_Tcell, tfs = Tcell)

# Visualization on RedeR
rdp <- RedPort()
calld(rdp, checkcalls=T)
addGraph(rdp, g_Tcell, layout=NULL)
addLegend.color(rdp, g_Tcell, type="edge")
addLegend.shape(rdp, g_Tcell)
relax(rdp, ps = TRUE)



################################################################
####    Differential gene expression analysis   ####
################################################################

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=48701)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SEQUENCE"))
write.table(tT, file=stdout(), row.names=F, sep="\t")



################################################################
####    Transcriptional Network Analysis    ####
################################################################

### Preprocessing ###
tTb <- subset(tT, tT$B >= 1.9)
tTb.entrez <- nuID2EntrezID(tTb$ID, lib.mapping = 'lumiHumanIDMapping')
tTb <- cbind(tTb, tTb.entrez)
tT.entrez <- nuID2EntrezID(tT$ID, lib.mapping = 'lumiHumanIDMapping')
tT <- cbind(tT, tT.entrez)

# Building different possibilities for the phenotype argument
b.stat <- tT$B
names(b.stat) <- rownames(tT)

logFC <- tT$logFC
names(logFC) <- rownames(tT)

# Building hits argument
tTb.hits <- rownames(tTb)
tTb.hits <- as.character(tTb.hits)

# Building phenoIDs argument
phenoID <- rownames(tT)
phenoID <- cbind(phenoID, tT.entrez)

#### TNA for the macrophage TFs ####

# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
rtna_macro <- tni2tna.preprocess(object = rtni_macro,
                                 phenotype = logFC,
                                 hits = tTb.hits,
                                 phenoIDs = phenoID)

# Run the MRA method
rtna_macro <- tna.mra(rtna_macro)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_macro <- tna.get(rtna_macro, what="mra", ntop = -1)
head(mra_macro)

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_macro <- tna.gsea1(rtna_macro, stepFilter=FALSE, nPermutations=100)

# Get GSEA results
gsea1_macro <- tna.get(rtna_macro, what="gsea1", ntop = -1)
head(gsea1_macro)

# Plot GSEA results
tna.plot.gsea1(rtna_macro, labPheno="abs(log2 fold changes)", ntop = -1,
               file = "tna_gsea_macro")

# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna_macro <- tna.gsea2(rtna_macro, stepFilter = FALSE, nPermutations = 100)

# Get GSEA-2T results
gsea2_macro <- tna.get(rtna_macro, what = "gsea2", ntop = -1)
head(gsea2_macro$differential)

# Plot GSEA-2T results
tna.plot.gsea2(rtna_macro, labPheno="log2 fold changes", tfs="IRF7",
               file = "tna_gsea2_macro")


#### TNA for the B-cell TFs ####

# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
rtna_Bcell <- tni2tna.preprocess(object = rtni_Bcell,
                                 phenotype = logFC,
                                 hits = tTb.hits,
                                 phenoIDs = phenoID)

# Run the MRA method
rtna_Bcell <- tna.mra(rtna_Bcell)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_Bcell <- tna.get(rtna_Bcell, what="mra", ntop = -1)
head(mra_Bcell)

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_Bcell <- tna.gsea1(rtna_Bcell, stepFilter=FALSE, nPermutations=100)

# Get GSEA results
gsea1_Bcell <- tna.get(rtna_Bcell, what="gsea1", ntop = -1)
head(gsea1_Bcell)

# Plot GSEA results
tna.plot.gsea1(rtna_Bcell, labPheno="abs(log2 fold changes)", ntop = -1,
               file = "tna_gsea_bcell")

# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna_Bcell <- tna.gsea2(rtna_Bcell, stepFilter = FALSE, nPermutations = 100)

# Get GSEA-2T results
gsea2_Bcell <- tna.get(rtna_Bcell, what = "gsea2", ntop = -1)
head(gsea2_Bcell$differential)

# Plot GSEA-2T results
tna.plot.gsea2(rtna_Bcell, labPheno="log2 fold changes", tfs="IRF8",
               file = "tna_gsea2_bcell")

#### TNA for the T-cell TFs ####

# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
rtna_Tcell <- tni2tna.preprocess(object = rtni_Tcell,
                                 phenotype = logFC,
                                 hits = tTb.hits,
                                 phenoIDs = phenoID)

# Run the MRA method
rtna_Tcell <- tna.mra(rtna_Tcell)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_Tcell <- tna.get(rtna_Tcell, what="mra", ntop = -1)
head(mra_Tcell)

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_Tcell <- tna.gsea1(rtna_Tcell, stepFilter=FALSE, nPermutations=100)

# Get GSEA results
gsea1_Tcell <- tna.get(rtna_Tcell, what="gsea1", ntop = -1)
head(gsea1_Tcell)

# Plot GSEA results
tna.plot.gsea1(rtna_Tcell, labPheno="abs(log2 fold changes)", ntop = -1,
               file = "tna_gsea_tcell")

# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna_Tcell <- tna.gsea2(rtna_Tcell, stepFilter = FALSE, nPermutations = 100)

# Get GSEA-2T results
gsea2_Tcell <- tna.get(rtna_Tcell, what = "gsea2", ntop = -1)
head(gsea2_Tcell$differential)

# Plot GSEA-2T results
tna.plot.gsea2(rtna_Tcell, labPheno="log2 fold changes", tfs="IRF7",
               file = "tna_gsea2_tcell")



################################################################
####   Boxplot for selected GEO samples   ####
################################################################
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE21257", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10295", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "X1XXX1XXX1XXXXXXX1X110000X1X10X10X000X0100X1XX0000000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  # set group names

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("no+metastases","metastases2y")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE21257", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
