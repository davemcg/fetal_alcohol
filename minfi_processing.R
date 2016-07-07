require(minfi) #1.8.7
require(IlluminaHumanMethylation450kmanifest)
targets=read.450k.sheet("/data/mcgaugheyd/projects/nhgri/brody/fetal_alcohol_methylation_study/Tsai_Sample_083/",pattern="Sheet.csv$")
RGset=read.450k.exp(targets=targets)
sampleNames(RGset) <- targets$Sample_Name
Mset.swan=preprocessSWAN(RGset)
Mset.swan.fixOutliers=fixMethOutliers(Mset.swan,verbose=TRUE)
Beta.swan.fixOutliers=getBeta(Mset.swan.fixOutliers,type="Illumina")
write.table(Beta.swan.fixOutliers,file="/cluster/ifs/projects/brody/fetal_alcohol_methylation_study/Illumina_beta.SWAN.fixOutliers.dat",quote=FALSE,sep="\t")
# Illumina_beta.SWAN.fixOutliers.t.IDs.dat created by hand massaging information from the sheets 
# (above) to add the Sentrix and Sample IDs to the transposed data above
 