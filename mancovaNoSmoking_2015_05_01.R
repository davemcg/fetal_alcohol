#!/home/mcgaugheyd/R/Rscript


args <- commandArgs(TRUE)

require(data.table)
all_data<-data.frame(fread("Illumina_beta.SWAN.fixOutliers.t.IDs.dat",header=T))
row.names(all_data)<-all_data$Sample
cgs_to_remove <- scan("illumina450k_positions_to_exclude.not_including_dbsnp_overlapping.withChen.dat",what="character")
samples_to_remove <- c("B50856","B55188","A41404","A32697")
more_cgs_to_remove <- scan("450k_cg_sites_that_overlap_dbsnp138.maf005_somepop_plus1bp.txt",what="character")
even_more <- scan("450k_cg_sites_that_LITERALLY_are_on_SNPS.maf005_somepop_plus1bp.txt",what="character")

all_cgs_to_remove <- unique(c(cgs_to_remove,more_cgs_to_remove,even_more))

spot_trimmed <- all_data[,!(colnames(all_data) %in% all_cgs_to_remove)]
almost_ready <- spot_trimmed[!(rownames(spot_trimmed) %in% samples_to_remove),]
almost_ready$Type <- substr(almost_ready$Sample,1,1)
almost_ready <- subset(almost_ready,Type=="A" | Type=="B")
rownames(almost_ready) <- almost_ready$Sample
ready <- almost_ready[!(rownames(almost_ready) %in% samples_to_remove),]
ready <-ready[order(row.names(ready)),]

id <- read.table("sentrix_id_conversion_table_demographics_v6.txt",header=T,sep="\t")
row.names(id) <- id$Sample
readyInfo <- merge(id,ready,by="row.names")

row.names(readyInfo) <- readyInfo$Row.names

#readyInfo$BMI <- (readyInfo$Prepregnancy.Weight..pounds./(readyInfo$Height..inches.)^2)*703


readyInfo.1<-readyInfo[complete.cases(readyInfo[ncol(readyInfo)]),]
readyInfo.1<-readyInfo.1[complete.cases(readyInfo.1[ncol(readyInfo.1)-1]),]
readyInfo.2<-subset(readyInfo.1,Ethnicity=="Black, NonHispanic"|Ethnicity=="White, NonHispanic")

readyInfo.2<-readyInfo.2[,order(names(readyInfo.2))]

#readyInfo.2<-subset(readyInfo.2,Case.Control=="Control")



# identifies the position of the first array spot in the data frame
positions <- grep("^cg",names(readyInfo.2))

start <- positions[1]
end <- positions[length(positions)]


# function to split numeric vector x into n sized parts
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
# splitting positions of all of the cg columns into 20 parts
splits <- chunk2(seq(start,end,1),20)
# selecting which of the 20 parts to use, based on argument from bash
indices <- unlist(splits[args[1]])

# logit
logit <- function(x) log(x/(1-x))

# car package for Mancova function ("Anova")
require(car)
mancova  <- sapply(logit(readyInfo.2[,indices]), function(x) Anova(lm(x~

# factor(readyInfo.2$Smoking) +
 factor(readyInfo.2$Case.Control) +
# factor(readyInfo.2$Cocaine) +
 factor(readyInfo.2$Ethnicity) +
 as.numeric(readyInfo.2$Maternal.Age) +
 factor(readyInfo.2$Gender) +
# factor(readyInfo.2$nyc) +
 as.numeric(readyInfo.2$Birth.Weight..grams.) +
 as.numeric(readyInfo.2$Number.of.Live.Births) +
 factor(readyInfo.2$Education) + 
 as.numeric(readyInfo.2$CD8T) +
 as.numeric(readyInfo.2$CD4T) +
 as.numeric(readyInfo.2$NK) +
 as.numeric(readyInfo.2$Bcell) +
 as.numeric(readyInfo.2$Mono) +
 as.numeric(readyInfo.2$Gran)))) 


output <- apply(mancova,2,unlist)
filename <- paste("mancova_noSmoking_2015-05-01_logit_",args[1],".dat",sep="")


 write.table(output, filename,quote=FALSE,sep="\t",row.names=FALSE)
