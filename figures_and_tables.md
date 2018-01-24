# Figures and Tables
January, 2018  


```r
library(ggbio)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: ggplot2
```

```
## Need specific help about ggbio? try mailing 
##  the maintainer or visit http://tengfei.github.com/ggbio/
```

```
## 
## Attaching package: 'ggbio'
```

```
## The following objects are masked from 'package:ggplot2':
## 
##     geom_bar, geom_rect, geom_segment, ggsave, stat_bin,
##     stat_identity, xlim
```

```r
library(EnsDb.Hsapiens.v75)
```

```
## Loading required package: ensembldb
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicFeatures
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: AnnotationFilter
```

```r
ensdb <- EnsDb.Hsapiens.v75
library(tidyverse)
```

```
## Warning: package 'tidyverse' was built under R version 3.4.2
```

```
## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ tibble  1.4.2     ✔ purrr   0.2.4
## ✔ tidyr   0.7.2     ✔ dplyr   0.7.4
## ✔ readr   1.1.1     ✔ stringr 1.2.0
## ✔ tibble  1.4.2     ✔ forcats 0.2.0
```

```
## Warning: package 'tibble' was built under R version 3.4.3
```

```
## Warning: package 'tidyr' was built under R version 3.4.2
```

```
## Warning: package 'purrr' was built under R version 3.4.2
```

```
## Warning: package 'dplyr' was built under R version 3.4.2
```

```
## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ ggbio::autoplot()      masks ggplot2::autoplot()
## ✖ dplyr::collapse()      masks IRanges::collapse()
## ✖ dplyr::combine()       masks Biobase::combine(), BiocGenerics::combine()
## ✖ dplyr::desc()          masks IRanges::desc()
## ✖ tidyr::expand()        masks S4Vectors::expand()
## ✖ dplyr::filter()        masks stats::filter()
## ✖ dplyr::first()         masks S4Vectors::first()
## ✖ stringr::fixed()       masks ggbio::fixed()
## ✖ ggbio::geom_bar()      masks ggplot2::geom_bar()
## ✖ ggbio::geom_rect()     masks ggplot2::geom_rect()
## ✖ ggbio::geom_segment()  masks ggplot2::geom_segment()
## ✖ ggbio::ggsave()        masks ggplot2::ggsave()
## ✖ dplyr::lag()           masks stats::lag()
## ✖ ggplot2::Position()    masks BiocGenerics::Position(), base::Position()
## ✖ purrr::reduce()        masks GenomicRanges::reduce(), IRanges::reduce()
## ✖ dplyr::rename()        masks S4Vectors::rename()
## ✖ dplyr::select()        masks ensembldb::select(), AnnotationDbi::select()
## ✖ dplyr::slice()         masks IRanges::slice()
## ✖ ggbio::stat_bin()      masks ggplot2::stat_bin()
## ✖ ggbio::stat_identity() masks ggplot2::stat_identity()
## ✖ ggbio::xlim()          masks ggplot2::xlim()
```

```r
library(ggbeeswarm)
library(ggthemes)
library(ggrepel)
```

```
## Warning: package 'ggrepel' was built under R version 3.4.2
```

```r
library(ggsci)
```

```
## Warning: package 'ggsci' was built under R version 3.4.2
```

```r
library(cowplot)
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:ggbio':
## 
##     ggsave
```

```
## The following object is masked from 'package:ggplot2':
## 
##     ggsave
```

```r
load('reanalysis_processed.Rdata')
```
# Figure 1

```r
# custom fuction to make data avail for ggplot
grabber <- function(cg_id){
  data <- bVals[cg_id,]
  data <- t(data) %>% data.frame() %>% tibble::rownames_to_column('Sample')
  data$Group <- targets$Sample_Group 
  data <- data %>% gather(Probe, Beta, -Group, -Sample)
  data <- data %>% mutate(Group = case_when(Group == 'Case' ~ 'FAS',
                                            TRUE ~ 'Control'))
  return(data)
}

# volcano
volcano <- DMP_e %>% 
  mutate(Significant = case_when(adj.P.Val < 0.05 ~ 'FDR < 0.05',
                                 TRUE ~ 'Not Significant')) %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value))) +
  scale_color_manual(values = c("#BC3C29FF", "grey")) +
  geom_point(aes(colour=Significant)) +
  geom_label_repel(
    data = DMP_e %>% 
      mutate(Label= case_when(Symbol=='' ~ Probe,
                              TRUE ~ paste(Probe, Symbol, sep=' | '))) %>% 
      filter(adj.P.Val < 0.1 & grepl('MAP3K13', Symbol) | adj.P.Val < 0.1 & grepl('GFI1', Symbol)),
    aes(label = Label),
    alpha=0.5,
    force=10,
    max.iter = 100000) +
  coord_cartesian(xlim=c(-0.75,0.75)) +
  ylab('-log10(P value)')

# plot most significantly differentially methylated CpGs
top_CG_up <- DMP_e  %>% filter(adj.P.Val < 0.1 & grepl('MAP3K13', Symbol)) %>% arrange(P.Value) %>%  pull(Probe)
top_CG_down <- DMP_e  %>% filter(adj.P.Val < 0.1 & grepl('GFI1', Symbol)) %>% arrange(P.Value) %>%  pull(Probe)
up <- grabber(top_CG_up)
down <- grabber(top_CG_down)

up_plot <- ggplot(up, aes(x=Group, y=Beta, colour=Group)) + 
  facet_wrap(~Probe, ncol=10, strip.position = 'bottom') + 
  geom_boxplot(width=0.4, outlier.shape = NA) +
  geom_quasirandom(size=0.5, alpha=0.4) +
  theme_minimal() +
  ggsci::scale_colour_nejm() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(angle=45)) + 
  coord_cartesian(ylim=c(0,1))

down_plot <- ggplot(down, aes(x=Group, y=Beta, colour=Group)) + 
  facet_wrap(~Probe, ncol=10, strip.position = 'bottom') + 
  geom_boxplot(width=0.4, outlier.shape = NA) +
  geom_quasirandom(size=0.5, alpha=0.4) +
  theme_minimal() +
  ggsci::scale_colour_nejm() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(angle=45)) + 
  coord_cartesian(ylim=c(0,1))

legend <- get_legend(up_plot + theme(legend.position="left"))
box_plots <- plot_grid(up_plot + theme(legend.position="none"), 
                       down_plot + theme(legend.position="none"), 
                       nrow=2, align = 'h', labels = c('B','C'))
plot_grid(volcano, box_plots, legend, nrow=1, labels=c('A',NULL,NULL),rel_widths = c(2,1,0.2))
```

![](figures_and_tables_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

# Figure 2

```r
gfi1 <- DMP_e %>% filter(chr=='chr1', pos > 92940000, pos < 92955000) %>% arrange(chr, pos) %>% pull(Probe)
gfi1_b <- grabber(gfi1)
gfi1_b <- left_join(gfi1_b, DMP_e %>% filter(grepl('GFI1$', Symbol)), by=c('Probe')) 
gfi1_meaned <- gfi1_b %>% group_by(Group, Probe) %>% summarise(Beta = mean(Beta), pos=mean(pos))
gfi1_plot <- gfi1_b %>% 
  ggplot(aes(x=pos, y=Beta, group=Group, colour=Group))  + 
  geom_jitter(width=100, alpha=0.01) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_line(data = gfi1_meaned %>% data.frame(), aes(x=pos,y=Beta,colour=Group), size=0.5) + 
  scale_color_nejm() +
  coord_cartesian(xlim=c(92940000,92955000)) + xlab('chr1')

map13_start <- 185000000
map13_stop <- 185010000
map13 <- DMP_e %>% filter(chr=='chr3', pos > map13_start, pos < (map13_stop+1000000)) %>% arrange(chr, pos) %>% pull(Probe)
map13 <- grabber(map13)
map13<- left_join(map13, DMP_e %>% filter(chr=='chr3', pos > map13_start, pos < (map13_stop+1000000)), by=c('Probe')) 
map13_meaned <- map13 %>% group_by(Group, Probe) %>% summarise(Beta = mean(Beta), pos=mean(pos))
map13_plot <- map13 %>% 
  ggplot(aes(x=pos, y=Beta, group=Group, colour=Group))  + 
  geom_jitter(width=10, alpha=0.01) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_line(data = map13_meaned %>% data.frame(), aes(x=pos,y=Beta,colour=Group), size=0.5) + 
  scale_color_nejm() +
  coord_cartesian(xlim=c(map13_start,map13_stop)) + xlab('chr3')


gfi_gene <- ggplot() + geom_alignment(ensdb, which = TxNameFilter("ENST00000427103"), names.expr=NULL, stat='identity') +
  coord_cartesian(xlim=c(92940000,92955000)) + xlab('GFI1')
```

```
## Fetching data...OK
## Parsing exons...OK
## Defining introns...OK
## Defining UTRs...OK
## Defining CDS...OK
## aggregating...
## Done
## Constructing graphics...
```

```r
map3_gene <- ggplot() + geom_alignment(ensdb, which = TxNameFilter("ENST00000424227"), names.expr=NULL, stat='identity') +
  coord_cartesian(xlim=c(map13_start,map13_stop)) + xlab('MAP3K13')
```

```
## Fetching data...OK
## Parsing exons...OK
## Defining introns...OK
## Defining UTRs...OK
## Defining CDS...OK
## aggregating...
## Done
## Constructing graphics...
```

```r
plot_grid(gfi1_plot + theme(legend.position="none"), map13_plot, gfi_gene, map3_gene, labels = c('A','B',NULL,NULL), ncol=2, align='v', axis = 'lr', rel_heights = c(1,0.3))
```

```
## Warning: Removed 1 rows containing missing values (geom_text).

## Warning: Removed 1 rows containing missing values (geom_text).
```

![](figures_and_tables_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

# GO enrichment
Not much for GO enrichment when using Enrichr (only uses direct-ish overlaps between probes and genes)

LOADS of stuff (but anything useful?) when using GREAT (which takes all probes and assigns closest gene)

```r
# paste these lists into enrichr (http://amp.pharm.mssm.edu/Enrichr/)
DMP_e %>% filter(adj.P.Val < 0.05, logFC > 0) %>% pull(Symbol) %>% unique() %>% paste(., collapse = ',')
```

```
## [1] ",GPD1,SPSB1,SLFN12L,PURA,LINC01024,TMEM200C,RBM20,HAPLN4,DOCK3,EMX1,MAST1,CHMP1B,GNAL,TPGS1,IGFBP2,NGFR,MCOLN3,PPT2,PPT2-EGFL8,PRRT1,LOC100507547,CACNA1A,SV2B,WT1,IGFBP5,RIIAD1,ATP1A2,HOXC6,HOXC5,HOXC4,TMEM184B,HSPB2-C11orf52,HSPB2,RGS20,CLEC4F,NOTCH4,PLCD4,MIR9-3HG,THBS4,ALKAL1,NTN1,SEZ6,TNXB,ATF6B,CCNF,LINC00028,REM1,LHPP,MYO1G,DFNA5,AQP1,INMT-MINDY4,MINDY4,ZDHHC24,ACTN3,USP44,FGFBP1,ADORA2B,BTBD11,MFSD1,SLC22A23,COL13A1,C1orf204,VSIG8,SLC17A8,PRICKLE1,SH2D4A,UNC5B,LBX2,LBX2-AS1"
```

```r
DMP_e %>% filter(adj.P.Val < 0.05, logFC < 0) %>% pull(Symbol) %>% unique() %>% paste(., collapse = ',')
```

```
## [1] "MAP3K13,GFI1,RPH3A,RPRD1B,TTI1,SEC61G,PTGR2,ZNF410,TSACC,CCT3,SURF4,STKLD1,DMRT2,SNORD52,C6orf48,NCAPH,SNX19,"
```

```r
# for GREAT (http://great.stanford.edu)
sig <- DMP_e %>% dplyr::select(Probe, chr, pos, Gene=Symbol, Type, logFC,FDR = adj.P.Val) %>% filter(FDR < 0.05) %>% mutate(Type = gsub('_',' ', Type)) %>% mutate(Type=gsub('hg19 genes ','',Type)) %>% mutate(Type=gsub('hg19 cpg ','',Type)) %>% mutate(pos2=pos+1) %>% dplyr::select(chr, pos, pos2) %>% data.frame()
write_tsv(sig, path='sig.bed', col_names = F)
background <- DMP_e %>% mutate(pos2=pos+1) %>% dplyr::select(chr, pos, pos2)
write_tsv(background %>% arrange(chr, pos) %>% mutate(pos2=as.integer(pos2)),path = 'background.bed', col_names = FALSE)
```


# Table 2
Currently filtering at FDR < 0.1.
In paper with FDR < 0.01 (11 probes)?

```r
table2_1 <- DMP_e %>% dplyr::select(Probe, chr, pos, Gene=Symbol, Type, logFC,FDR = adj.P.Val) %>% filter(FDR < 0.1) %>% mutate(Type = gsub('_',' ', Type)) %>% mutate(Type=gsub('hg19 genes ','',Type)) %>% mutate(Type=gsub('hg19 cpg ','',Type))
table2_1
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Probe"],"name":[1],"type":["chr"],"align":["left"]},{"label":["chr"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pos"],"name":[3],"type":["int"],"align":["right"]},{"label":["Gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Type"],"name":[5],"type":["chr"],"align":["left"]},{"label":["logFC"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"cg09249404","2":"chr3","3":"185001026","4":"MAP3K13","5":"introns, shores","6":"-0.31696083","7":"0.0009452857"},{"1":"cg26250086","2":"chr1","3":"6830972","4":"","5":"intergenic, inter","6":"0.28952853","7":"0.0009452857"},{"1":"cg16092017","2":"chr5","3":"134827643","4":"","5":"intergenic, islands","6":"0.23282585","7":"0.0009452857"},{"1":"cg18901140","2":"chr17","3":"8908225","4":"","5":"intergenic, shores","6":"0.28398190","7":"0.0009452857"},{"1":"cg09935388","2":"chr1","3":"92947588","4":"GFI1","5":"introns, islands","6":"-0.67990449","7":"0.0009452857"},{"1":"cg00636368","2":"chr3","3":"185000586","4":"MAP3K13","5":"promoters, islands","6":"-0.60189683","7":"0.0014548343"},{"1":"cg06168149","2":"chr12","3":"50497778","4":"GPD1","5":"5UTRs, exons, inter","6":"0.19971845","7":"0.0022728845"},{"1":"cg09256832","2":"chr1","3":"9382878","4":"SPSB1","5":"introns, inter","6":"0.36945903","7":"0.0027097446"},{"1":"cg21390082","2":"chr17","3":"33842255","4":"SLFN12L","5":"introns, inter","6":"0.25842455","7":"0.0059649991"},{"1":"cg11067714","2":"chr5","3":"139492230","4":"PURA,LINC01024","5":"1to5kb, shores","6":"0.22151331","7":"0.0083603415"},{"1":"cg00058163","2":"chr18","3":"5891948","4":"TMEM200C","5":"exons, shores","6":"0.27230991","7":"0.0083603415"},{"1":"cg05567511","2":"chr7","3":"151113760","4":"","5":"intergenic, inter","6":"0.29379268","7":"0.0109362749"},{"1":"cg01257889","2":"chr12","3":"113230065","4":"RPH3A","5":"5UTRs, exons, inter","6":"-0.28211237","7":"0.0115376387"},{"1":"cg22623319","2":"chr10","3":"112427911","4":"RBM20","5":"introns, shelves","6":"0.31270774","7":"0.0115376387"},{"1":"cg18146737","2":"chr1","3":"92946700","4":"GFI1","5":"introns, islands","6":"-0.45644750","7":"0.0143669043"},{"1":"cg03775416","2":"chr19","3":"19368775","4":"HAPLN4","5":"exons, islands","6":"0.54350349","7":"0.0143669043"},{"1":"cg25415695","2":"chr5","3":"134828036","4":"","5":"intergenic, shores","6":"0.23262314","7":"0.0143669043"},{"1":"cg25959131","2":"chr3","3":"50855436","4":"DOCK3","5":"introns, inter","6":"0.29507065","7":"0.0143669043"},{"1":"cg07140194","2":"chr2","3":"73150236","4":"EMX1","5":"introns, shores","6":"0.24081803","7":"0.0143669043"},{"1":"cg16800939","2":"chr14","3":"100631798","4":"","5":"intergenic, inter","6":"0.17495037","7":"0.0143669043"},{"1":"cg07633618","2":"chr3","3":"185000208","4":"MAP3K13","5":"promoters, shores","6":"-0.35563716","7":"0.0143669043"},{"1":"cg06307939","2":"chr19","3":"12984645","4":"MAST1","5":"exons, islands","6":"0.57912627","7":"0.0143669043"},{"1":"cg09662411","2":"chr1","3":"92946132","4":"GFI1","5":"introns, islands","6":"-0.44348323","7":"0.0144212215"},{"1":"cg27093273","2":"chr18","3":"5892213","4":"TMEM200C","5":"promoters, shores","6":"0.33219057","7":"0.0144212215"},{"1":"cg18848419","2":"chr18","3":"11850153","4":"CHMP1B,GNAL","5":"1to5kb, introns, shores","6":"0.29696901","7":"0.0144212215"},{"1":"cg21177456","2":"chr19","3":"519611","4":"TPGS1","5":"exons, 3UTRs, islands","6":"0.38505605","7":"0.0153467225"},{"1":"cg03792653","2":"chr20","3":"36661448","4":"RPRD1B,TTI1","5":"promoters, introns, shores","6":"-0.31763669","7":"0.0158173961"},{"1":"cg22979368","2":"chr2","3":"217501581","4":"IGFBP2","5":"introns, shelves","6":"0.22691361","7":"0.0158173961"},{"1":"cg09163442","2":"chr7","3":"54827431","4":"SEC61G","5":"promoters, shores","6":"-0.28540276","7":"0.0172949234"},{"1":"cg04613258","2":"chr17","3":"47589594","4":"NGFR","5":"introns, shores","6":"0.21572373","7":"0.0172949234"},{"1":"cg22500876","2":"chr7","3":"45187902","4":"","5":"intergenic, inter","6":"0.26718511","7":"0.0174336495"},{"1":"cg16299091","2":"chr18","3":"11850139","4":"CHMP1B,GNAL","5":"1to5kb, introns, shores","6":"0.25319342","7":"0.0174336495"},{"1":"cg12037450","2":"chr1","3":"85513124","4":"MCOLN3","5":"5UTRs, exons, introns, shores","6":"0.22391695","7":"0.0174336495"},{"1":"cg06211243","2":"chr18","3":"12907446","4":"","5":"intergenic, shelves","6":"0.29764566","7":"0.0190867231"},{"1":"cg23261715","2":"chr12","3":"54144380","4":"","5":"intergenic, shores","6":"0.16609154","7":"0.0190867231"},{"1":"cg27629948","2":"chr15","3":"89285170","4":"","5":"intergenic, islands","6":"0.15085767","7":"0.0190867231"},{"1":"cg18316974","2":"chr1","3":"92947035","4":"GFI1","5":"introns, islands","6":"-0.36952309","7":"0.0190867231"},{"1":"cg25426302","2":"chr6","3":"32120826","4":"PPT2,PPT2-EGFL8,PRRT1,LOC100507547","5":"promoters, 1to5kb, 5UTRs, exons, shores","6":"0.38053470","7":"0.0190867231"},{"1":"cg26579653","2":"chr14","3":"74319300","4":"PTGR2,ZNF410","5":"introns, shores","6":"-0.53526531","7":"0.0190867231"},{"1":"cg00536080","2":"chr19","3":"13409931","4":"CACNA1A","5":"exons, islands","6":"0.21643817","7":"0.0190867231"},{"1":"cg24516362","2":"chr15","3":"91838184","4":"SV2B","5":"exons, 3UTRs, shores","6":"0.20638024","7":"0.0190867231"},{"1":"cg10244666","2":"chr11","3":"32421808","4":"WT1","5":"introns, inter","6":"0.42862800","7":"0.0190867231"},{"1":"cg05324273","2":"chr10","3":"60086738","4":"","5":"intergenic, inter","6":"0.25708554","7":"0.0190867231"},{"1":"cg19008649","2":"chr2","3":"217560292","4":"IGFBP5","5":"promoters, shores","6":"0.18170663","7":"0.0190867231"},{"1":"cg08893853","2":"chr1","3":"151701251","4":"RIIAD1","5":"exons, inter","6":"0.28131044","7":"0.0205996388"},{"1":"cg22932677","2":"chr19","3":"12983876","4":"MAST1","5":"introns, shores","6":"0.28151224","7":"0.0207234886"},{"1":"cg22171758","2":"chr12","3":"54146133","4":"","5":"intergenic, shores","6":"0.26290478","7":"0.0213916176"},{"1":"cg19741660","2":"chr14","3":"100631897","4":"","5":"intergenic, inter","6":"0.31752176","7":"0.0213916176"},{"1":"cg18191092","2":"chr8","3":"64331042","4":"","5":"promoters, inter","6":"0.25033544","7":"0.0213916176"},{"1":"cg20632887","2":"chr1","3":"160085581","4":"ATP1A2","5":"5UTRs, exons, inter","6":"0.26217008","7":"0.0214285062"},{"1":"cg25122233","2":"chr12","3":"54412506","4":"HOXC6,HOXC5,HOXC4","5":"introns, shores","6":"0.19773618","7":"0.0215399369"},{"1":"cg08575722","2":"chr11","3":"32421845","4":"WT1","5":"introns, inter","6":"0.26651303","7":"0.0216583153"},{"1":"cg14736365","2":"chr5","3":"173947646","4":"","5":"intergenic, inter","6":"0.31811621","7":"0.0219198001"},{"1":"cg02755938","2":"chr18","3":"11850146","4":"CHMP1B,GNAL","5":"1to5kb, introns, shores","6":"0.23578405","7":"0.0228112076"},{"1":"cg09851989","2":"chr17","3":"55235","4":"","5":"intergenic, shores","6":"0.26237602","7":"0.0230101527"},{"1":"cg05565052","2":"chr22","3":"38669495","4":"TMEM184B","5":"promoters, shores","6":"0.21182639","7":"0.0257834612"},{"1":"cg00079023","2":"chr5","3":"139492535","4":"PURA,LINC01024","5":"1to5kb, shores","6":"0.15007000","7":"0.0259579293"},{"1":"cg09065113","2":"chr16","3":"23798532","4":"","5":"intergenic, inter","6":"0.32680840","7":"0.0259579293"},{"1":"cg11845417","2":"chr11","3":"111789613","4":"HSPB2-C11orf52,HSPB2","5":"5UTRs, exons, inter","6":"0.15753960","7":"0.0281874316"},{"1":"cg14285533","2":"chr7","3":"63386328","4":"","5":"intergenic, shores","6":"0.79409541","7":"0.0294015587"},{"1":"cg06932776","2":"chr8","3":"54852232","4":"RGS20","5":"exons, introns, inter","6":"0.24574344","7":"0.0294015587"},{"1":"cg12748639","2":"chr2","3":"71035898","4":"CLEC4F","5":"exons, 3UTRs, inter","6":"0.19755615","7":"0.0294015587"},{"1":"cg15654136","2":"chr2","3":"71047710","4":"CLEC4F","5":"5UTRs, exons, inter","6":"0.23760551","7":"0.0294015587"},{"1":"cg24210717","2":"chr12","3":"50497827","4":"GPD1","5":"5UTRs, exons, inter","6":"0.16187876","7":"0.0294015587"},{"1":"cg17918558","2":"chr1","3":"46945595","4":"","5":"intergenic, inter","6":"0.26429141","7":"0.0298173793"},{"1":"cg00285620","2":"chr11","3":"102147694","4":"","5":"intergenic, inter","6":"0.24113937","7":"0.0298173793"},{"1":"cg24378227","2":"chr6","3":"1396771","4":"","5":"intergenic, shelves","6":"0.19237925","7":"0.0298173793"},{"1":"cg10585648","2":"chr1","3":"156307318","4":"TSACC,CCT3","5":"1to5kb, 5UTRs, exons, introns, shores","6":"-0.29869219","7":"0.0298173793"},{"1":"cg20816361","2":"chr1","3":"2467100","4":"","5":"intergenic, islands","6":"0.28494122","7":"0.0298173793"},{"1":"cg04671145","2":"chr6","3":"157041096","4":"","5":"intergenic, islands","6":"0.56850891","7":"0.0303898583"},{"1":"cg27532171","2":"chr1","3":"46952054","4":"","5":"intergenic, shores","6":"0.20263663","7":"0.0309150766"},{"1":"cg03507824","2":"chr11","3":"64176993","4":"","5":"intergenic, inter","6":"0.21857230","7":"0.0309150766"},{"1":"cg15090036","2":"chr17","3":"48706075","4":"","5":"intergenic, inter","6":"0.19091583","7":"0.0309150766"},{"1":"cg03257743","2":"chr2","3":"11123616","4":"","5":"intergenic, inter","6":"0.26106589","7":"0.0309150766"},{"1":"cg15318627","2":"chr6","3":"32181764","4":"NOTCH4","5":"introns, inter","6":"0.16481723","7":"0.0309150766"},{"1":"cg18887096","2":"chr2","3":"219472410","4":"PLCD4","5":"promoters, inter","6":"0.16492253","7":"0.0309150766"},{"1":"cg24360241","2":"chr2","3":"233370823","4":"","5":"intergenic, shelves","6":"0.15105338","7":"0.0313213107"},{"1":"cg07038187","2":"chr15","3":"89922989","4":"MIR9-3HG","5":"introns, shores","6":"0.24476287","7":"0.0316590879"},{"1":"cg22244122","2":"chr9","3":"136243698","4":"SURF4,STKLD1","5":"promoters, 1to5kb, introns, shores","6":"-0.15025068","7":"0.0335117719"},{"1":"cg20714328","2":"chr5","3":"79331850","4":"THBS4","5":"introns, shores","6":"0.22737554","7":"0.0336611354"},{"1":"cg16868095","2":"chr3","3":"197090774","4":"","5":"intergenic, inter","6":"0.15783484","7":"0.0341757214"},{"1":"cg12046602","2":"chr6","3":"157041162","4":"","5":"intergenic, islands","6":"0.66157222","7":"0.0341757214"},{"1":"cg16198103","2":"chr8","3":"53471884","4":"ALKAL1","5":"introns, inter","6":"0.27020855","7":"0.0348173888"},{"1":"cg03740162","2":"chr9","3":"1050680","4":"DMRT2","5":"promoters, 5UTRs, exons, introns, shores","6":"-0.22299226","7":"0.0348173888"},{"1":"cg02481697","2":"chr17","3":"8928012","4":"NTN1","5":"introns, shores","6":"0.37023935","7":"0.0366332614"},{"1":"cg13541527","2":"chr6","3":"31804078","4":"SNORD52,C6orf48","5":"promoters, introns, shores","6":"-0.30999096","7":"0.0369018262"},{"1":"cg05697249","2":"chr11","3":"111789693","4":"HSPB2-C11orf52,HSPB2","5":"5UTRs, exons, inter","6":"0.18372231","7":"0.0379053835"},{"1":"cg20480274","2":"chr17","3":"27313125","4":"SEZ6","5":"1to5kb, introns, islands","6":"0.30446238","7":"0.0379053835"},{"1":"cg02657865","2":"chr6","3":"32077744","4":"TNXB,ATF6B","5":"promoters, introns, inter","6":"0.17996402","7":"0.0379053835"},{"1":"cg25337631","2":"chr10","3":"103500093","4":"","5":"intergenic, inter","6":"0.22181830","7":"0.0380471870"},{"1":"cg05870645","2":"chr16","3":"2503590","4":"CCNF","5":"introns, inter","6":"0.24202060","7":"0.0380471870"},{"1":"cg02725313","2":"chr20","3":"30071913","4":"LINC00028,REM1","5":"1to5kb, introns, islands","6":"0.17535648","7":"0.0382667021"},{"1":"cg13220109","2":"chr10","3":"126150491","4":"LHPP","5":"exons, islands","6":"0.36465814","7":"0.0382667021"},{"1":"cg04180046","2":"chr7","3":"45002736","4":"MYO1G","5":"introns, islands","6":"0.31189560","7":"0.0396225455"},{"1":"cg08356445","2":"chr12","3":"67461558","4":"","5":"intergenic, shores","6":"0.17288785","7":"0.0396440263"},{"1":"cg20764575","2":"chr7","3":"24797884","4":"DFNA5","5":"promoters, shores","6":"0.22880172","7":"0.0396440263"},{"1":"cg14370520","2":"chr7","3":"30891837","4":"AQP1,INMT-MINDY4,MINDY4","5":"1to5kb, exons, inter","6":"0.21383594","7":"0.0396440263"},{"1":"cg08158233","2":"chr8","3":"37351056","4":"","5":"intergenic, inter","6":"0.23511581","7":"0.0396440263"},{"1":"cg25497796","2":"chr13","3":"20876280","4":"","5":"intergenic, shores","6":"0.36039959","7":"0.0396440263"},{"1":"cg04663194","2":"chr11","3":"66315239","4":"ZDHHC24,ACTN3","5":"1to5kb, introns, shores","6":"0.20578565","7":"0.0396440263"},{"1":"cg14565151","2":"chr12","3":"95945384","4":"USP44","5":"promoters, 1to5kb, shelves","6":"0.18922403","7":"0.0396440263"},{"1":"cg03913456","2":"chr2","3":"97000924","4":"NCAPH","5":"promoters, shores","6":"-0.24989119","7":"0.0396440263"},{"1":"cg00800679","2":"chr4","3":"15940058","4":"FGFBP1","5":"1to5kb, 5UTRs, exons, inter","6":"0.17520673","7":"0.0396440263"},{"1":"cg07563400","2":"chr17","3":"15849556","4":"ADORA2B","5":"introns, shores","6":"0.25129935","7":"0.0396440263"},{"1":"cg14141102","2":"chr12","3":"107974557","4":"BTBD11","5":"promoters, introns, shores","6":"0.17191115","7":"0.0396440263"},{"1":"cg18165065","2":"chr1","3":"8130665","4":"","5":"intergenic, inter","6":"0.29449614","7":"0.0413647687"},{"1":"cg13839439","2":"chr12","3":"54144521","4":"","5":"intergenic, shores","6":"0.17876482","7":"0.0423545342"},{"1":"cg01971181","2":"chr3","3":"158518888","4":"MFSD1","5":"promoters, shores","6":"0.23975704","7":"0.0455989051"},{"1":"cg10091752","2":"chr6","3":"3455459","4":"SLC22A23","5":"promoters, introns, shores","6":"0.17516286","7":"0.0459078651"},{"1":"cg06322432","2":"chr10","3":"112289853","4":"","5":"intergenic, shores","6":"0.20214313","7":"0.0461909963"},{"1":"cg27268279","2":"chr10","3":"71672083","4":"COL13A1","5":"introns, inter","6":"0.18661871","7":"0.0464463907"},{"1":"cg00303378","2":"chr1","3":"159825552","4":"C1orf204,VSIG8","5":"promoters, introns, shores","6":"0.19036064","7":"0.0474053499"},{"1":"cg24446178","2":"chr12","3":"100750702","4":"SLC17A8","5":"promoters, inter","6":"0.20092436","7":"0.0477936885"},{"1":"cg02264943","2":"chr11","3":"130785842","4":"SNX19","5":"5UTRs, exons, introns, shores","6":"-0.13003060","7":"0.0478579049"},{"1":"cg21903286","2":"chr12","3":"42876128","4":"PRICKLE1","5":"introns, shores","6":"0.27493755","7":"0.0482443167"},{"1":"cg18132363","2":"chr6","3":"166260572","4":"","5":"introns, inter","6":"-0.25916233","7":"0.0487868991"},{"1":"cg03138955","2":"chr8","3":"19170250","4":"SH2D4A","5":"promoters, 1to5kb, shores","6":"0.27274163","7":"0.0492576637"},{"1":"cg23622162","2":"chr10","3":"73059249","4":"UNC5B","5":"exons, 3UTRs, inter","6":"0.27034701","7":"0.0495221477"},{"1":"cg02100410","2":"chr2","3":"74731354","4":"LBX2,LBX2-AS1","5":"promoters, 1to5kb, exons, 3UTRs, shores","6":"0.18847304","7":"0.0495990577"},{"1":"cg25200182","2":"chr1","3":"158119395","4":"","5":"intergenic, shores","6":"-0.22205929","7":"0.0495990577"},{"1":"cg24690588","2":"chr1","3":"57285399","4":"C1orf168","5":"promoters, inter","6":"0.21885776","7":"0.0502513728"},{"1":"cg11367633","2":"chr20","3":"40247916","4":"CHD6","5":"promoters, shores","6":"0.17841227","7":"0.0502764637"},{"1":"cg04694209","2":"chr16","3":"89228970","4":"LINC02138,LINC00304","5":"1to5kb, introns, shelves","6":"0.20284730","7":"0.0504191566"},{"1":"cg13986670","2":"chr1","3":"11708792","4":"FBXO2","5":"exons, islands","6":"0.31236985","7":"0.0506055061"},{"1":"cg17310600","2":"chr4","3":"165305334","4":"MARCH1","5":"promoters, shores","6":"0.13079640","7":"0.0509636651"},{"1":"cg19590578","2":"chr5","3":"77266117","4":"","5":"intergenic, shelves","6":"0.31588337","7":"0.0512748800"},{"1":"cg18590502","2":"chr3","3":"49203081","4":"CCDC71","5":"1to5kb, introns, shores","6":"-0.21725068","7":"0.0525364092"},{"1":"cg13734833","2":"chr11","3":"68779817","4":"MRGPRF","5":"promoters, introns, shores","6":"0.15769606","7":"0.0525364092"},{"1":"cg01952185","2":"chr5","3":"134813213","4":"","5":"intergenic, inter","6":"0.29091021","7":"0.0549098967"},{"1":"cg00769843","2":"chr11","3":"63754453","4":"OTUB1","5":"5UTRs, exons, introns, islands","6":"-0.33110623","7":"0.0567455528"},{"1":"cg24158160","2":"chr1","3":"151693222","4":"RIIAD1,CELF3","5":"promoters, 1to5kb, introns, shores","6":"0.17244062","7":"0.0575942456"},{"1":"cg16347256","2":"chr1","3":"57285160","4":"C1orf168","5":"5UTRs, exons, inter","6":"0.21682984","7":"0.0579343162"},{"1":"cg02343335","2":"chr2","3":"113785852","4":"IL36B","5":"introns, inter","6":"0.22312588","7":"0.0580597372"},{"1":"cg08757448","2":"chr10","3":"104613767","4":"BORCS7,BORCS7-ASMT","5":"promoters, shores","6":"-0.16643203","7":"0.0580597372"},{"1":"cg15620905","2":"chr1","3":"44024150","4":"PTPRF","5":"introns, inter","6":"0.21183064","7":"0.0580597372"},{"1":"cg07506407","2":"chr8","3":"121385077","4":"","5":"intergenic, inter","6":"0.15319561","7":"0.0580597372"},{"1":"cg24437523","2":"chr1","3":"160085568","4":"ATP1A2","5":"5UTRs, exons, inter","6":"0.22317630","7":"0.0580597372"},{"1":"cg12791065","2":"chr15","3":"101419269","4":"ALDH1A3","5":"promoters, islands","6":"-0.20232177","7":"0.0580597372"},{"1":"cg17986992","2":"chr1","3":"159825761","4":"C1orf204,VSIG8","5":"promoters, exons, shores","6":"0.38101505","7":"0.0583031545"},{"1":"cg26069230","2":"chr17","3":"29250073","4":"ADAP2","5":"exons, shores","6":"0.22952564","7":"0.0587261519"},{"1":"cg02491754","2":"chr12","3":"53773040","4":"SP1","5":"promoters, 1to5kb, shores","6":"-0.24066257","7":"0.0587743562"},{"1":"cg25191241","2":"chr3","3":"124448994","4":"UMPS","5":"promoters, shores","6":"-0.23930798","7":"0.0588269788"},{"1":"cg09795979","2":"chr11","3":"45308098","4":"SYT13","5":"promoters, shores","6":"0.18106658","7":"0.0588269788"},{"1":"cg06397322","2":"chr19","3":"48102118","4":"","5":"intergenic, shores","6":"0.20340782","7":"0.0588269788"},{"1":"cg13870494","2":"chr9","3":"72658358","4":"MAMDC2","5":"promoters, shores","6":"0.30419800","7":"0.0600072098"},{"1":"cg06210630","2":"chr1","3":"160085433","4":"ATP1A2","5":"promoters, inter","6":"0.19837089","7":"0.0600072098"},{"1":"cg18932699","2":"chr8","3":"64518293","4":"","5":"intergenic, inter","6":"0.20461351","7":"0.0600072098"},{"1":"cg21787078","2":"chr1","3":"34099118","4":"CSMD2","5":"exons, islands","6":"-0.11709807","7":"0.0600072098"},{"1":"cg11632438","2":"chr3","3":"158519176","4":"MFSD1","5":"promoters, shores","6":"0.21125335","7":"0.0600072098"},{"1":"cg13352518","2":"chr11","3":"124708433","4":"","5":"intergenic, shores","6":"0.15965378","7":"0.0600072098"},{"1":"cg11855943","2":"chr17","3":"47964664","4":"","5":"intergenic, shelves","6":"0.12302991","7":"0.0600072098"},{"1":"cg10831246","2":"chr17","3":"18553663","4":"","5":"intergenic, inter","6":"0.21578578","7":"0.0600072098"},{"1":"cg06595479","2":"chr19","3":"519035","4":"TPGS1","5":"exons, islands","6":"0.23932240","7":"0.0600072098"},{"1":"cg21205031","2":"chr2","3":"32502699","4":"YIPF4","5":"promoters, shores","6":"-0.16894807","7":"0.0600191664"},{"1":"cg06108383","2":"chr6","3":"32120899","4":"PPT2,PPT2-EGFL8,PRRT1,LOC100507547","5":"promoters, 1to5kb, exons, shores","6":"0.28396724","7":"0.0600191664"},{"1":"cg15012027","2":"chr2","3":"99351550","4":"MGAT4A","5":"1to5kb, shelves","6":"0.18868337","7":"0.0608076422"},{"1":"cg20687616","2":"chr1","3":"56883736","4":"","5":"intergenic, inter","6":"0.29265686","7":"0.0608076422"},{"1":"cg11429111","2":"chr5","3":"134813329","4":"","5":"intergenic, inter","6":"0.33229123","7":"0.0608076422"},{"1":"cg21151061","2":"chr12","3":"50498016","4":"GPD1","5":"introns, inter","6":"0.18249576","7":"0.0614133298"},{"1":"cg26645709","2":"chr15","3":"89987792","4":"","5":"intergenic, inter","6":"0.20256801","7":"0.0614133298"},{"1":"cg14585415","2":"chr10","3":"102778683","4":"PDZD7","5":"exons, islands","6":"0.14798590","7":"0.0614133298"},{"1":"cg21207436","2":"chr14","3":"74815316","4":"VRTN","5":"5UTRs, exons, inter","6":"0.18205458","7":"0.0614133298"},{"1":"cg21073930","2":"chr19","3":"12984170","4":"MAST1","5":"exons, islands","6":"0.48688892","7":"0.0614133298"},{"1":"cg05881762","2":"chr15","3":"25684849","4":"UBE3A","5":"promoters, shores","6":"-0.31035179","7":"0.0614133298"},{"1":"cg17200137","2":"chr10","3":"48416896","4":"GDF2","5":"promoters, shelves","6":"0.28396956","7":"0.0621939390"},{"1":"cg00806138","2":"chr1","3":"47881404","4":"FOXE3","5":"promoters, shores","6":"0.17543907","7":"0.0627627702"},{"1":"cg10906284","2":"chr12","3":"63544430","4":"AVPR1A","5":"exons, islands","6":"0.30187570","7":"0.0627627702"},{"1":"cg05575921","2":"chr5","3":"373378","4":"AHRR","5":"introns, shores","6":"-0.29800069","7":"0.0629884196"},{"1":"cg16781992","2":"chr4","3":"20985623","4":"KCNIP4","5":"5UTRs, exons, introns, inter","6":"0.33357153","7":"0.0631661083"},{"1":"cg01176715","2":"chr2","3":"211037310","4":"KANSL1L","5":"1to5kb, shores","6":"0.19142941","7":"0.0638849183"},{"1":"cg04155830","2":"chr3","3":"120314864","4":"NDUFB4","5":"promoters, shores","6":"0.26134460","7":"0.0642663740"},{"1":"cg08814800","2":"chr4","3":"174421377","4":"","5":"intergenic, islands","6":"0.24096629","7":"0.0642663740"},{"1":"cg02017109","2":"chr5","3":"76788595","4":"WDR41","5":"promoters, shores","6":"0.15521015","7":"0.0642663740"},{"1":"cg21535366","2":"chr12","3":"107974612","4":"BTBD11","5":"promoters, introns, shores","6":"0.19981443","7":"0.0643865916"},{"1":"cg07780348","2":"chr11","3":"118842413","4":"FOXR1","5":"promoters, islands","6":"0.21343448","7":"0.0643865916"},{"1":"cg21625386","2":"chr4","3":"20985696","4":"KCNIP4","5":"promoters, introns, inter","6":"0.23057560","7":"0.0650799856"},{"1":"cg26181818","2":"chr2","3":"233323935","4":"ALPI","5":"exons, 3UTRs, islands","6":"0.28728349","7":"0.0653851470"},{"1":"cg06205614","2":"chr16","3":"31128351","4":"KAT8","5":"promoters, shores","6":"-0.28084229","7":"0.0662330715"},{"1":"cg22941086","2":"chr12","3":"40019883","4":"C12orf40","5":"promoters, inter","6":"0.22422354","7":"0.0662330715"},{"1":"cg10442355","2":"chr8","3":"61326412","4":"","5":"promoters, inter","6":"0.19636372","7":"0.0664810672"},{"1":"cg00063111","2":"chr3","3":"39448952","4":"SNORA6,SNORA62,RPSA","5":"promoters, 1to5kb, introns, shores","6":"-0.21778652","7":"0.0666681642"},{"1":"cg16273979","2":"chr4","3":"122376685","4":"","5":"intergenic, inter","6":"0.25735937","7":"0.0670850315"},{"1":"cg11600596","2":"chr2","3":"2617199","4":"","5":"intergenic, inter","6":"-0.22944694","7":"0.0670850315"},{"1":"cg20467658","2":"chr15","3":"89933347","4":"MIR9-3HG","5":"introns, inter","6":"0.31905195","7":"0.0670991604"},{"1":"cg15484406","2":"chr14","3":"94461913","4":"LINC00521","5":"1to5kb, inter","6":"0.21462810","7":"0.0674193449"},{"1":"cg01189038","2":"chr3","3":"183994516","4":"ECE2","5":"promoters, introns, shores","6":"0.13158845","7":"0.0674193449"},{"1":"cg02953960","2":"chr7","3":"55756785","4":"FKBP9P1","5":"1to5kb, 5UTRs, exons, introns, inter","6":"0.14721098","7":"0.0675008524"},{"1":"cg12221621","2":"chr1","3":"9381317","4":"SPSB1","5":"introns, inter","6":"0.21009445","7":"0.0678547948"},{"1":"cg21450784","2":"chr19","3":"519609","4":"TPGS1","5":"exons, 3UTRs, islands","6":"0.39667523","7":"0.0678547948"},{"1":"cg01293485","2":"chr3","3":"42847044","4":"ACKR2,HIGD1A","5":"1to5kb, exons, shores","6":"0.16160675","7":"0.0678547948"},{"1":"cg01440864","2":"chr17","3":"44892087","4":"WNT3","5":"introns, shelves","6":"0.18731226","7":"0.0688386487"},{"1":"cg11565911","2":"chr12","3":"72233249","4":"TBC1D15","5":"promoters, shores","6":"-0.25495032","7":"0.0688386487"},{"1":"cg06840491","2":"chr10","3":"48416966","4":"GDF2","5":"promoters, shelves","6":"0.30704911","7":"0.0693571786"},{"1":"cg13736131","2":"chr11","3":"8832907","4":"ST5","5":"promoters, introns, inter","6":"0.17667395","7":"0.0696922987"},{"1":"cg05113927","2":"chr2","3":"27531244","4":"UCN","5":"promoters, islands","6":"0.29282116","7":"0.0705301258"},{"1":"cg12803068","2":"chr7","3":"45002919","4":"MYO1G","5":"introns, shores","6":"0.47279519","7":"0.0710025217"},{"1":"cg24684168","2":"chr6","3":"132271588","4":"CTGF","5":"exons, islands","6":"-0.20650679","7":"0.0710025217"},{"1":"cg09244277","2":"chr2","3":"12538247","4":"","5":"introns, inter","6":"0.18327308","7":"0.0727138124"},{"1":"cg00020474","2":"chr1","3":"8214963","4":"","5":"intergenic, inter","6":"0.16630216","7":"0.0727138124"},{"1":"cg25598308","2":"chr11","3":"120440821","4":"GRIK4","5":"introns, inter","6":"0.13514229","7":"0.0727138124"},{"1":"cg00313914","2":"chr1","3":"201618901","4":"NAV1","5":"introns, islands","6":"-0.36160791","7":"0.0738606646"},{"1":"cg22821355","2":"chr4","3":"10117479","4":"WDR1","5":"1to5kb, exons, introns, shores","6":"-0.20019305","7":"0.0763804705"},{"1":"cg06956786","2":"chr7","3":"3713348","4":"SDK1","5":"introns, inter","6":"0.15709972","7":"0.0763804705"},{"1":"cg24699976","2":"chr11","3":"66384069","4":"RBM14,RBM14-RBM4","5":"5UTRs, exons, islands","6":"-0.29865844","7":"0.0763804705"},{"1":"cg25540816","2":"chr5","3":"148928885","4":"CSNK1A1","5":"introns, shores","6":"-0.24118663","7":"0.0763804705"},{"1":"cg14847483","2":"chr15","3":"34516640","4":"EMC4","5":"promoters, inter","6":"-0.37933229","7":"0.0763804705"},{"1":"cg03268942","2":"chr2","3":"74148310","4":"","5":"intergenic, inter","6":"0.18952739","7":"0.0763804705"},{"1":"cg01656853","2":"chr19","3":"49199172","4":"FUT2","5":"promoters, shores","6":"0.13741141","7":"0.0763886112"},{"1":"cg21940331","2":"chr2","3":"131451866","4":"","5":"introns, shores","6":"-0.18101020","7":"0.0763886112"},{"1":"cg05641656","2":"chr19","3":"58728610","4":"","5":"intergenic, islands","6":"0.16483922","7":"0.0763886112"},{"1":"cg20703122","2":"chr15","3":"28339563","4":"OCA2","5":"introns, shores","6":"0.29429003","7":"0.0763886112"},{"1":"cg09867953","2":"chr16","3":"23706666","4":"ERN2","5":"exons, inter","6":"0.16919088","7":"0.0763886112"},{"1":"cg10831285","2":"chr19","3":"17918927","4":"B3GNT3","5":"exons, islands","6":"0.39405447","7":"0.0763886112"},{"1":"cg04470557","2":"chr5","3":"131630995","4":"SLC22A4","5":"promoters, exons, introns, shores","6":"0.18751904","7":"0.0763886112"},{"1":"cg21488617","2":"chr1","3":"38022586","4":"SNIP1,DNALI1","5":"1to5kb, 5UTRs, exons, islands","6":"0.13812249","7":"0.0763886112"},{"1":"cg24035898","2":"chr1","3":"159825595","4":"C1orf204,VSIG8","5":"promoters, introns, shores","6":"0.26986434","7":"0.0763886112"},{"1":"cg14557236","2":"chr10","3":"104613771","4":"BORCS7,BORCS7-ASMT","5":"promoters, shores","6":"-0.15569293","7":"0.0763886112"},{"1":"cg05480594","2":"chr3","3":"185000760","4":"MAP3K13","5":"5UTRs, exons, islands","6":"-0.51795712","7":"0.0772277395"},{"1":"cg14576824","2":"chr1","3":"213224402","4":"RPS6KC1","5":"promoters, shores","6":"-0.23480841","7":"0.0772277395"},{"1":"cg01038207","2":"chr3","3":"50330361","4":"IFRD2,HYAL3","5":"promoters, exons, 3UTRs, shores","6":"-0.16289873","7":"0.0772277395"},{"1":"cg02652018","2":"chr7","3":"45188275","4":"","5":"intergenic, inter","6":"0.15452449","7":"0.0773183501"},{"1":"cg07008193","2":"chr11","3":"46318527","4":"CREB3L1","5":"introns, shores","6":"0.15667713","7":"0.0800818671"},{"1":"cg18470427","2":"chr17","3":"33842301","4":"SLFN12L","5":"introns, inter","6":"0.25164733","7":"0.0800818671"},{"1":"cg08985785","2":"chr6","3":"39196340","4":"KCNK5","5":"introns, shores","6":"-0.37951923","7":"0.0809304966"},{"1":"cg19227924","2":"chr17","3":"28565709","4":"SLC6A4","5":"1to5kb, shelves","6":"0.15933567","7":"0.0817770283"},{"1":"cg19961800","2":"chr14","3":"62550978","4":"SYT16","5":"exons, inter","6":"0.21010079","7":"0.0817770283"},{"1":"cg18221429","2":"chr7","3":"150020136","4":"LRRC61,ACTR3C","5":"promoters, introns, islands","6":"-0.34566308","7":"0.0818052631"},{"1":"cg22697136","2":"chr1","3":"11714248","4":"FBXO44,FBXO2","5":"promoters, exons, introns, islands","6":"0.14778985","7":"0.0820579183"},{"1":"cg19398269","2":"chr6","3":"28678460","4":"","5":"intergenic, inter","6":"0.20357719","7":"0.0826898762"},{"1":"cg20893180","2":"chr8","3":"36957379","4":"","5":"intergenic, inter","6":"0.17624957","7":"0.0826898762"},{"1":"cg07068187","2":"chr22","3":"25349743","4":"","5":"intergenic, islands","6":"0.15351291","7":"0.0826898762"},{"1":"cg08862291","2":"chr2","3":"158648039","4":"ACVR1","5":"introns, inter","6":"-0.17906857","7":"0.0831190325"},{"1":"cg14096207","2":"chr1","3":"223407010","4":"SUSD4","5":"introns, inter","6":"0.14074224","7":"0.0831952016"},{"1":"cg23435671","2":"chr8","3":"29211078","4":"DUSP4","5":"1to5kb, shores","6":"0.31313043","7":"0.0838829939"},{"1":"cg17586988","2":"chr18","3":"5892245","4":"TMEM200C","5":"promoters, shores","6":"0.27274305","7":"0.0838829939"},{"1":"cg24632139","2":"chr1","3":"159832353","4":"VSIG8","5":"5UTRs, exons, inter","6":"0.16386037","7":"0.0838829939"},{"1":"ch.18.1546242R","2":"chr18","3":"76891500","4":"ATP9B","5":"introns, inter","6":"0.25142677","7":"0.0838829939"},{"1":"cg03476948","2":"chr14","3":"78715384","4":"NRXN3","5":"introns, inter","6":"0.19257079","7":"0.0838829939"},{"1":"cg18517898","2":"chr19","3":"36435675","4":"LRFN3","5":"exons, islands","6":"0.19781261","7":"0.0838829939"},{"1":"cg07534186","2":"chr2","3":"113915259","4":"","5":"intergenic, islands","6":"0.13729328","7":"0.0838829939"},{"1":"cg18418457","2":"chr5","3":"135364328","4":"TGFBI","5":"promoters, shores","6":"0.20749590","7":"0.0844489670"},{"1":"cg25577130","2":"chr17","3":"58226188","4":"CA4","5":"1to5kb, shores","6":"0.17255323","7":"0.0853885748"},{"1":"cg13568106","2":"chr17","3":"26835115","4":"FOXN1","5":"introns, inter","6":"0.25530450","7":"0.0853885748"},{"1":"cg24059075","2":"chr12","3":"49688316","4":"PRPH","5":"promoters, shores","6":"0.20965093","7":"0.0855974484"},{"1":"cg25607643","2":"chr15","3":"40730161","4":"BAHD1","5":"1to5kb, shores","6":"0.20643073","7":"0.0855974484"},{"1":"cg01461862","2":"chr10","3":"27548591","4":"","5":"introns, shores","6":"0.14223676","7":"0.0861839521"},{"1":"cg16968985","2":"chr17","3":"27313254","4":"SEZ6","5":"1to5kb, introns, islands","6":"0.21079569","7":"0.0867589860"},{"1":"cg16489360","2":"chr7","3":"129845574","4":"TMEM209,SSMEM1","5":"promoters, 1to5kb, shores","6":"-0.15705942","7":"0.0867589860"},{"1":"cg00497735","2":"chr16","3":"86364736","4":"","5":"intergenic, inter","6":"0.24818770","7":"0.0867589860"},{"1":"cg01090161","2":"chr6","3":"34032747","4":"GRM4","5":"1to5kb, introns, inter","6":"0.15589465","7":"0.0867589860"},{"1":"cg09461185","2":"chr1","3":"186649631","4":"PTGS2","5":"promoters, islands","6":"-0.13736981","7":"0.0867589860"},{"1":"cg18505593","2":"chr14","3":"100059479","4":"CCDC85C","5":"introns, inter","6":"0.34079991","7":"0.0867589860"},{"1":"cg19831077","2":"chr7","3":"151107924","4":"WDR86,WDR86-AS1","5":"promoters, 1to5kb, introns, shores","6":"0.21322672","7":"0.0868606398"},{"1":"cg06415891","2":"chr1","3":"2467096","4":"","5":"intergenic, islands","6":"0.24356416","7":"0.0868606398"},{"1":"cg25112848","2":"chr3","3":"39470521","4":"","5":"intergenic, inter","6":"0.25104558","7":"0.0868606398"},{"1":"cg08812802","2":"chr12","3":"49686619","4":"PRPH","5":"1to5kb, shelves","6":"0.17477265","7":"0.0868606398"},{"1":"cg04404543","2":"chr12","3":"54019307","4":"ATF7","5":"introns, shores","6":"-0.13773445","7":"0.0874815781"},{"1":"cg05204104","2":"chr2","3":"235403141","4":"ARL4C","5":"exons, 3UTRs, shores","6":"0.24968230","7":"0.0874815781"},{"1":"cg08626436","2":"chr10","3":"88555339","4":"BMPR1A","5":"introns, inter","6":"0.32853004","7":"0.0875238274"},{"1":"cg03083228","2":"chr11","3":"43660436","4":"","5":"intergenic, inter","6":"0.22172444","7":"0.0875542711"},{"1":"cg25425390","2":"chr12","3":"42763958","4":"PPHLN1","5":"introns, inter","6":"-0.19336960","7":"0.0875542711"},{"1":"cg26438105","2":"chr18","3":"5892119","4":"TMEM200C","5":"promoters, shores","6":"0.28314670","7":"0.0875542711"},{"1":"cg18716679","2":"chr8","3":"39835944","4":"IDO2","5":"promoters, introns, inter","6":"0.14743217","7":"0.0875542711"},{"1":"cg04807108","2":"chr15","3":"78697933","4":"","5":"intergenic, inter","6":"0.23656916","7":"0.0875542711"},{"1":"cg11067786","2":"chr2","3":"74780467","4":"LOXL3,DOK1","5":"promoters, 1to5kb, introns, shores","6":"-0.17389990","7":"0.0881950664"},{"1":"cg03115048","2":"chr1","3":"152077402","4":"","5":"intergenic, shelves","6":"0.25216707","7":"0.0881950664"},{"1":"cg09203312","2":"chr13","3":"20805196","4":"GJB6","5":"5UTRs, exons, introns, shores","6":"0.32541925","7":"0.0882237108"},{"1":"cg00166722","2":"chr3","3":"10149974","4":"FANCD2OS","5":"promoters, 1to5kb, islands","6":"0.13125683","7":"0.0895681109"},{"1":"cg26218577","2":"chr15","3":"69744466","4":"RPLP1","5":"promoters, shores","6":"-0.31529944","7":"0.0897258158"},{"1":"cg00499707","2":"chr20","3":"34679665","4":"EPB41L1","5":"promoters, introns, shores","6":"0.17717031","7":"0.0897258158"},{"1":"cg20075683","2":"chr5","3":"132083910","4":"CCNI2","5":"introns, islands","6":"0.17759039","7":"0.0904374999"},{"1":"cg14190151","2":"chr2","3":"14778017","4":"FAM84A","5":"exons, introns, 3UTRs, shelves","6":"0.20625739","7":"0.0904374999"},{"1":"cg16135695","2":"chr10","3":"1889970","4":"","5":"intergenic, inter","6":"0.26873168","7":"0.0904374999"},{"1":"cg01519094","2":"chr11","3":"125034181","4":"PKNOX2","5":"promoters, shores","6":"0.22087307","7":"0.0904374999"},{"1":"cg02975463","2":"chr2","3":"45178953","4":"","5":"intergenic, shores","6":"0.24972867","7":"0.0904374999"},{"1":"cg22777724","2":"chr17","3":"46622516","4":"HOXB2,HOXB-AS1","5":"promoters, exons, shores","6":"-0.17990926","7":"0.0912746811"},{"1":"cg22145172","2":"chr19","3":"19174151","4":"SLC25A42","5":"promoters, shores","6":"-0.17115844","7":"0.0912746811"},{"1":"cg04017512","2":"chr10","3":"35195138","4":"","5":"intergenic, inter","6":"0.25685773","7":"0.0913700943"},{"1":"cg25879142","2":"chr7","3":"4671391","4":"","5":"intergenic, inter","6":"0.24140303","7":"0.0913700943"},{"1":"cg03276149","2":"chr6","3":"11651032","4":"","5":"intergenic, inter","6":"0.25388156","7":"0.0919340977"},{"1":"cg15482422","2":"chr16","3":"66956240","4":"CDH16,RRAD","5":"1to5kb, exons, islands","6":"0.23203856","7":"0.0919340977"},{"1":"cg20566484","2":"chr14","3":"94489832","4":"OTUB2","5":"1to5kb, shelves","6":"0.13883916","7":"0.0919340977"},{"1":"cg12711760","2":"chr11","3":"20132923","4":"NAV2","5":"introns, inter","6":"-0.14289688","7":"0.0919340977"},{"1":"cg26890706","2":"chr6","3":"13924425","4":"RNF182","5":"promoters, shores","6":"0.17535457","7":"0.0919340977"},{"1":"cg04143626","2":"chr10","3":"85924004","4":"","5":"intergenic, inter","6":"0.17079222","7":"0.0919340977"},{"1":"cg23035209","2":"chr5","3":"72795181","4":"BTF3","5":"exons, introns, shores","6":"-0.18995324","7":"0.0920373878"},{"1":"cg04873861","2":"chr3","3":"12525454","4":"TSEN2","5":"promoters, shores","6":"0.15479892","7":"0.0920373878"},{"1":"cg25383593","2":"chr5","3":"134827820","4":"","5":"intergenic, shores","6":"0.19875431","7":"0.0920488774"},{"1":"cg09066361","2":"chr7","3":"126890254","4":"GRM8","5":"introns, shores","6":"0.23869161","7":"0.0923610374"},{"1":"cg09835878","2":"chr14","3":"74815131","4":"VRTN","5":"promoters, inter","6":"0.21035946","7":"0.0928098981"},{"1":"cg19815818","2":"chr17","3":"7037365","4":"","5":"intergenic, inter","6":"0.18432953","7":"0.0932784741"},{"1":"cg23025459","2":"chr1","3":"3134420","4":"PRDM16","5":"introns, inter","6":"-0.23278276","7":"0.0937009876"},{"1":"cg07870237","2":"chr5","3":"131348431","4":"ACSL6","5":"promoters, 1to5kb, shores","6":"0.20403997","7":"0.0948566543"},{"1":"cg08246428","2":"chr1","3":"160160106","4":"CASQ1","5":"promoters, inter","6":"0.14670496","7":"0.0952932163"},{"1":"cg15556873","2":"chr2","3":"78263008","4":"","5":"introns, inter","6":"-0.20409136","7":"0.0952932163"},{"1":"cg05223459","2":"chr11","3":"45928854","4":"C11orf94","5":"promoters, inter","6":"0.20807859","7":"0.0952932163"},{"1":"cg08762819","2":"chr11","3":"11757920","4":"","5":"intergenic, inter","6":"0.16039696","7":"0.0952932163"},{"1":"cg23462052","2":"chr20","3":"2452871","4":"SNRPB","5":"1to5kb, shores","6":"-0.14112667","7":"0.0952932163"},{"1":"cg16049246","2":"chr11","3":"76464026","4":"","5":"intergenic, inter","6":"0.09622486","7":"0.0952932163"},{"1":"cg18844118","2":"chr7","3":"26191489","4":"NFE2L3","5":"promoters, shores","6":"-0.24529078","7":"0.0957327773"},{"1":"cg05789250","2":"chr6","3":"31804306","4":"SNORD52,C6orf48","5":"promoters, introns, shores","6":"-0.25815286","7":"0.0957327773"},{"1":"cg26607337","2":"chr5","3":"145758881","4":"","5":"intergenic, islands","6":"0.19086282","7":"0.0957327773"},{"1":"cg03502812","2":"chr19","3":"36435729","4":"LRFN3","5":"exons, islands","6":"0.30215588","7":"0.0957327773"},{"1":"cg06581938","2":"chr7","3":"64403228","4":"","5":"intergenic, inter","6":"0.13425564","7":"0.0957327773"},{"1":"cg20892847","2":"chr12","3":"58011875","4":"SLC26A10","5":"promoters, 1to5kb, shores","6":"0.20917242","7":"0.0975548561"},{"1":"cg23038338","2":"chr1","3":"46996546","4":"","5":"intergenic, shelves","6":"0.16873097","7":"0.0975548561"},{"1":"cg06171908","2":"chr19","3":"57630340","4":"USP29","5":"1to5kb, islands","6":"0.23663515","7":"0.0984771339"},{"1":"cg05893709","2":"chr8","3":"39836064","4":"IDO2","5":"promoters, introns, inter","6":"0.22550580","7":"0.0984771339"},{"1":"cg25718438","2":"chr20","3":"43743921","4":"WFDC5","5":"promoters, inter","6":"0.11264565","7":"0.0984771339"},{"1":"cg21918313","2":"chr2","3":"73496203","4":"FBXO41","5":"exons, islands","6":"0.33972658","7":"0.0984771339"},{"1":"cg20722210","2":"chr5","3":"134827136","4":"","5":"intergenic, shores","6":"0.19166712","7":"0.0994890442"},{"1":"cg04211115","2":"chr17","3":"13504008","4":"HS3ST3A1","5":"exons, islands","6":"0.16228496","7":"0.0998237058"},{"1":"cg26500033","2":"chr6","3":"166585954","4":"T","5":"1to5kb, shelves","6":"0.25346710","7":"0.0998262351"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
write.csv(table2_1, file='figures/table2.csv')
```

# Kobor validation 
## Table 3

```r
kobor_sig <- readxl::read_xlsx('13072_2016_74_MOESM2_ESM.xlsx')
```

```
## Warning in strptime(x, format, tz = tz): unknown timezone 'default/America/
## New_York'
```

```r
kobor_probes <- kobor_sig %>% pull(Probe)

kobor_sig2 <- readxl::read_xlsx('13148_2018_439_MOESM1_ESM.xlsx', sheet = 'Supplementary table 2')
kobor_probes_validate2 <- kobor_sig2 %>% pull(CpG)
# number of Kobor probes in our dataset
overlap_num <- DMPs %>% tibble::rownames_to_column('Probe') %>% filter(Probe %in% kobor_probes) %>% nrow()

DMP_e %>% filter(Probe %in% kobor_probes, P.Value < (0.05/overlap_num)) %>% dplyr::select(Probe, chr, pos, Gene=Symbol, Type, logFC,FDR = adj.P.Val)  %>% mutate(Type = gsub('_',' ', Type)) %>% mutate(Type=gsub('hg19 genes ','',Type)) %>% mutate(Type=gsub('hg19 cpg ','',Type))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Probe"],"name":[1],"type":["chr"],"align":["left"]},{"label":["chr"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pos"],"name":[3],"type":["int"],"align":["right"]},{"label":["Gene"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Type"],"name":[5],"type":["chr"],"align":["left"]},{"label":["logFC"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"cg04180046","2":"chr7","3":"45002736","4":"MYO1G","5":"introns, islands","6":"0.3118956","7":"0.03962255"},{"1":"cg10831285","2":"chr19","3":"17918927","4":"B3GNT3","5":"exons, islands","6":"0.3940545","7":"0.07638861"},{"1":"cg09066361","2":"chr7","3":"126890254","4":"GRM8","5":"introns, shores","6":"0.2386916","7":"0.09236104"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Supplementary Figure 1

```r
z <- DMPs %>% 
  tibble::rownames_to_column('Probe') %>% 
  mutate(Kobor_sig = case_when(Probe %in% kobor_probes_validate2 ~ 'Validated',
                               TRUE ~ 'No')) 


t.test(z$logFC ~ z$Kobor_sig)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  z$logFC by z$Kobor_sig
## t = -2.2041, df = 131.06, p-value = 0.02926
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.026188002 -0.001414446
## sample estimates:
##        mean in group No mean in group Validated 
##            9.328511e-05            1.389451e-02
```

```r
z <- DMPs %>% 
  tibble::rownames_to_column('Probe') %>% 
  mutate(Kobor_sig = case_when(Probe %in% kobor_probes ~ 'Validated',
                               TRUE ~ 'No')) 


t.test(z$logFC ~ z$Kobor_sig)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  z$logFC by z$Kobor_sig
## t = -3.3781, df = 541.97, p-value = 0.0007822
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.016811501 -0.004448867
## sample estimates:
##        mean in group No mean in group Validated 
##            8.237309e-05            1.071256e-02
```

```r
DMP_e %>% 
 mutate(`Kobor Probes` = case_when(Probe %in% kobor_probes_validate2 ~ 'Validated',
                               Probe %in% kobor_probes ~ 'Significant',
                               TRUE ~ 'Not Significant')) %>% 
  ggplot(aes(x=abs(logFC), colour=`Kobor Probes`)) + 
  geom_density() +
  theme_bw() +
  scale_color_nejm() + coord_cartesian(xlim=c(0,0.5))
```

![](figures_and_tables_files/figure-html/unnamed-chunk-7-1.png)<!-- -->
