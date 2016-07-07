library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

tBeta <- fread('Illumina_beta.SWAN.fixOutliers.t.IDs.dat')

cg_positions <- DMPs %>% filter(adj.P.Val < 0.05) %>% .[['Name']]
plot_data <- tBeta %>% select(one_of(c("Sample",cg_positions))) %>% left_join(.,metadata) %>% select(one_of(c("Sample","Case.Control",cg_positions))) %>% gather(key=Sample,Case.Control)
colnames(plot_data)<-c("Sample","Case.Control","cg","Value")
ggplot(data=plot_data,aes(x=Case.Control,y=Value)) + geom_jitter() + geom_boxplot(alpha=0.5) + ylim(c(0,1)) + facet_wrap(~cg)


# individual plot
ggplot(data=subset(plot_data,cg=='cg13220109'),aes(x=Case.Control,y=Value)) + geom_jitter() + geom_boxplot(alpha=0.5) + ylim(c(0,1))
       
