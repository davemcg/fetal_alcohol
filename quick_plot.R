library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(minfi)

# my old data
tBeta_1 <- fread('Illumina_beta.SWAN.fixOutliers.t.IDs.dat')
# my new data
load('mSetSqFlt.Rdata')
bVals <- getBeta(mSetSqFlt)
name <- colnames(bVals)
bVals <- data.table(t(bVals))
bVals$Sample <- name

cg_positions <- DMPs %>% filter(adj.P.Val < 0.05) %>% .[['Name']]
plot_data <- bVals %>% select(one_of(c("Sample",cg_positions))) %>% left_join(.,metadata) %>% select(one_of(c("Sample","Case.Control",cg_positions))) %>% gather(key=Sample,Case.Control)
colnames(plot_data)<-c("Sample","Case.Control","cg","Value")
ggplot(data=plot_data,aes(x=Case.Control,y=Value)) + geom_jitter() + geom_boxplot(alpha=0.5) + ylim(c(0,1)) + facet_wrap(~cg)


# ethnicity
cg_positions <- DMPs_eth %>% filter(adj.P.Val < 0.05) %>% head(20) %>% .[['Name']]
plot_data <- bVals %>% select(one_of(c("Sample",cg_positions))) %>% left_join(.,metadata) %>% select(one_of(c("Sample","Ethnicity",cg_positions))) %>% gather(key=Sample,Ethnicity)
colnames(plot_data)<-c("Sample","Ethnicity","cg","Value")
ggplot(data=plot_data,aes(x=Ethnicity,y=Value)) + geom_jitter() + geom_boxplot(alpha=0.5) + ylim(c(0,1)) + facet_wrap(~cg)


# individual plot
ggplot(data=subset(plot_data,cg=='cg10831285'),aes(x=Case.Control,y=Value)) + geom_jitter() + geom_boxplot(alpha=0.5) + ylim(c(0,1))
       

#tsne

set.seed(935489)
# before filtering probes
tsne_out <- Rtsne(as.matrix(t(getM(mSetSq))),perplexity = 5, check_duplicates = FALSE)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$Sample <- row.names(t(mVals))
tsne_plot <- left_join(tsne_plot,metadata)
ggplot(tsne_plot,aes(x=X1,y=X2,label=Sample,shape=Ethnicity,colour=Gender)) +
geom_text_repel() + geom_point() + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")
# after probe filtering
tsne_out <- Rtsne(as.matrix(t(getM(mSetSqFlt))),perplexity = 5, check_duplicates = FALSE)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$Sample <- row.names(t(mVals))
tsne_plot <- left_join(tsne_plot,metadata)
ggplot(tsne_plot,aes(x=X1,y=X2,label=Sample,shape=Gender,colour=Ethnicity)) +
  geom_text_repel() + geom_point() + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")
