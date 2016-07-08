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
# mostly follows the workflow laid out by Maksimovic et al.
# http://f1000research.com/articles/5-1281/v1

# sva correction
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#batch-effects-correction-with-sva

# pull annotation
ann450k <-  getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# bring in raw data (Illumina .idat)
metadata <- fread('~/git/fetal_alcohol/sentrix_id_conversion_table_demographics_v6.txt')
targets <- read.metharray.sheet("/Volumes/ThunderBay/PROJECTS/brody/fetal_alcohol/Tsai_Sample_083/",pattern="Sheet.csv$")
rgSet <- read.metharray.exp(targets=targets)
sampleNames(rgSet) <- rgSet$Sample_Name

# add sample group info to targets, fewer rows in metadata, already parsed down to remove outliers
targets <- left_join(data.frame(metadata),targets,by=c("Sample"="Sample_Name"))
targets$Sample_Group <- targets$Case.Control
targets$Sample_Name <- targets$Sample

# remove samples based on previous outlier analysis / ethnicity
keep <-  rgSet$Sample_Name %in% targets$Sample_Name
rgSet <- rgSet[,keep]

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

# plot detection p-values
# everything looks fine in this data-set
# first match up targets to detP
targets <- targets[match(colnames(detP), targets$Sample),]

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
set.seed(23345230974)
mSetSq <- preprocessSWAN(rgSet)

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalization
targets <- targets[match(mSetSq$Sample_Name, targets$Sample),]
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)


# MDS plots to look at largest sources of variation
# first, line up targets and mSetSq
targets <- targets[match(mSetSq$Sample_Name, targets$Sample),]

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
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in%  c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
table(keep)

# pull in my cg sites that I'm skipping due to overlap with SNPs and 
# are cross-reactive and are on chrs X Y
#mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps=c("SBE","CpG"), maf=0)
cg_remove_1 <- scan('~/git/fetal_alcohol/450k_cg_sites_that_LITERALLY_are_on_SNPS.maf005_somepop_plus1bp.txt',what='character')
cg_remove_2 <- scan('~/git/fetal_alcohol/450k_cg_sites_that_overlap_dbsnp138.maf005_somepop_plus1bp.txt',what='character')
xReactive <- scan('~/git/fetal_alcohol/illumina450k_positions_to_exclude.not_including_dbsnp_overlapping.withChen.dat',what='character')
cg_remove <- unique(c(cg_remove_1,cg_remove_2,xReactive))

keep <- !(featureNames(mSetSqFlt) %in% cg_remove)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]           


# MDS plots to look at largest sources of variation
# now with the filtered set
# again, ensure targets and mSetSqFlt line up
targets <- targets[match(mSetSqFlt$Sample_Name, targets$Sample),]
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
#head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
#head(bVals[,1:5])

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")

# explicitly match targets sample order (rows) with mVars
targets <- targets[match(colnames(mVals), targets$Sample),]

# setup for limma
case_control <- factor(targets$Case.Control)
ethnicity <- factor(targets$Ethnicity)
gender <- factor(targets$Gender)
#education <- as.numeric(targets$Education)
cd8t <- as.numeric(targets$CD8T)
cd4t <- as.numeric(targets$CD4T)
nk <- as.numeric(targets$NK)
mono <- as.numeric(targets$Mono)
gran <- as.numeric(targets$Gran)
#smoking <- factor(targets$Smoking)
 
# limma run
design <- model.matrix(~0+case_control+gender+cd8t+cd4t+nk+mono+gran+ethnicity, data=targets)
design_ethnicity <- model.matrix(~0+ethnicity+case_control+gender+cd8t+cd4t+nk+mono+gran, data=targets)

colnames(design)<-c("Case","Control","Gender","CD8T","CD4","NK","MONO","GRAN","Hispanic","Other","White")
colnames(design_ethnicity)<-c("Black","Hispanic","Other","White","Case.Control","Gender","CD8T","CD4","NK","MONO","GRAN")

cmtx <- makeContrasts( "Case-Control", levels=design)
cmtx_ethnicity <- makeContrasts( "Black-White", levels=design_ethnicity)

fit <- lmFit(mVals, design)
fit2 <- contrasts.fit(fit, cmtx)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

fit_eth <- lmFit(mVals, design_ethnicity)
fit2_eth <- contrasts.fit(fit_eth, cmtx_ethnicity)
fit2_eth <- eBayes(fit2_eth)
summary(decideTests(fit2_eth))



# get the table of results for the first contrast 
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)

DMPs_eth <- topTable(fit2_eth, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs_eth)

myAnnotation <- cpg.annotate(mVals, datatype = "array",
                             analysis.type="differential", design=design,
                             contrasts = TRUE, cont.matrix = cmtx,
                             coef="Case-Control")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# grab ranges for plotting
results.ranges <- extractRanges(DMRs, genome = "hg19")
groups <- c(Case='magenta',Control='forestgreen')
cols <- groups[as.character(targets$Case.Control)]
DMR.plot(ranges=results.ranges,dmr=1,CpGs=mVals,phen.col=as.factor(targets$Case.Control),genome='hg19')
DMR.plot(ranges=results.ranges,dmr=1,CpGs=mVals,phen.col=cols,genome='hg19',samps=c(1:10))




###########################
# don't use, results in a confusing mess
#
# RUV to improve signal by removing spurious patterns
# https://www.bioconductor.org/packages/devel/bioc/vignettes/missMethyl/inst/doc/missMethyl.pdf
# first, ID positions that do not contribute to case-control 
# this part is kind of fuzzy and I selected 0.9 because that resulted in ~ 1/3 false and 2/3 true
ctl <- rownames(mVals) %in% rownames(DMPs[DMPs$adj.P.Val > 0.9,])
table(ctl)
grp <- factor(targets$Case.Control,levels=c("Case","Control"))
des <- model.matrix(~grp)
rfit1 <- RUVfit(data=mVals, design=des, coef=2, ctl=ctl) # Stage 2 analysis
rfit2 <- RUVadj(rfit1)
nrow(topRUV(rfit2,number = Inf) %>% filter(p.ebayes.BH<0.05))
##############################



# variance analysis
design <- model.matrix(~factor(targets$Sample_Group))
fitvar <- varFit(mVals, design = design, coef = c(1,2))
topDV <- topVar(fitvar, coef=2)