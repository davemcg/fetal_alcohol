# compare https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-016-0074-4
# with this set

library(dplyr)
library(tidyr)
library(readxl)
library(data.table)
library(minfi)

# my data
load('mSetSqFlt.Rdata')
mVals <- getM(mSetSqFlt)
name <- colnames(mVals)
mVals <- data.table(t(mVals))
mVals$Sample <- name
# my metadata
metadata <- fread('sentrix_id_conversion_table_demographics_v6.txt')
# their data (Portales-Casamar, Lussier, et al.)
PC_L <- read_excel('validation/13072_2016_74_MOESM2_ESM.xlsx')

x<-left_join(mVals, metadata) %>% select(one_of(c("Sample","Case.Control",PC_L)))
                                         