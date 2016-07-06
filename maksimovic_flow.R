# ssh -Y biowulf2.nih.gov
# sinteractive --mem=20G --x11=first
# module load R/3.3.0_gcc-4.9.1
# package loading
library(data.table)
library(dplyr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(Gviz)
library(DMRcate)
library(stringr)
library(data.table)
library(dplyr)
# closely follows the workflow laid out by Maksimovic et al.
# http://f1000research.com/articles/5-1281/v1

# pull annotation
ann450k <-  getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# bring in raw data (Illumina .idat)
setwd('/data/mcgaugheyd/projects/nhgri/brody/fetal_alcohol_methylation_study/')
metadata <- fread('sentrix_id_conversion_table_demographics_v6.txt')
targets <- read.metharray.sheet("/data/mcgaugheyd/projects/nhgri/brody/fetal_alcohol_methylation_study/Tsai_Sample_083",pattern="Sheet.csv$")
rgSet <- read.metharray.exp(targets=targets)
sampleNames(rgSet) <- targets$Sample_Name

# add sample group info to targets
full_sample_names <- targets$Sample_Name
targets <- left_join(data.frame(metadata),targets,by=c("Sample"="Sample_Name"))
targets$Sample_Group <- targets$Case.Control
targets$Sample_Name <- targets$Sample

# remove samples based on previous outlier analysis / ethnicity
keep <-  full_sample_names %in% targets$Sample_Name
rgSet <- rgSet[,keep]

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

# plot detection p-values
# everything looks fine in this data-set
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")


# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalization
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)


# MDS plots to look at largest sources of variation
par(mfrow=c(2,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Gender)])
legend("top", legend=levels(factor(targets$Gender)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Ethnicity)])
legend("top", legend=levels(factor(targets$Ethnicity)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Plate)])
legend("top", legend=levels(factor(targets$Sample_Plate)), text.col=pal,
       bg="white", cex=0.7)

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in%
                                                      c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
table(keep)

# pull in my cg sites that I'm skipping due to overlap with SNPs and 
# are cross-reactive and are on chrs X Y
cg_remove_1 <- scan('450k_cg_sites_that_LITERALLY_are_on_SNPS.maf005_somepop_plus1bp.txt',what='character')
cg_remove_2 <- scan('450k_cg_sites_that_overlap_dbsnp138.maf005_somepop_plus1bp.txt',what='character')
xReactive <- scan('illumina450k_positions_to_exclude.not_including_dbsnp_overlapping.withChen.dat',what='character')
cg_remove <- c(cg_remove_1, cg_remove_2, xReactive)
cg_remove <- unique(cg_remove)

keep <- !(featureNames(mSetSqFlt) %in% cg_remove)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]           


# MDS plots to look at largest sources of variation
# now with the filtered set
par(mfrow=c(2,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Gender)])
legend("top", legend=levels(factor(targets$Gender)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Ethnicity)])
legend("top", legend=levels(factor(targets$Ethnicity)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Plate)])
legend("top", legend=levels(factor(targets$Sample_Plate)), text.col=pal,
       bg="white", cex=0.7)


# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")




case.control <- factor(targets$Case.Control)
ethnicity <- factor(targets$Ethnicity)
design <- model.matrix(~0+case.control+ethnicity, data=targets)
colnames(design) <- c("Case","Control","Hispanic","OtherEth","White")

fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Case-Control,Hispanic-White,levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)


# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)