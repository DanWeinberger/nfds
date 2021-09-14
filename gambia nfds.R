#Data from GPS, figshare: https://figshare.com/articles/dataset/Roary_gene_presence_and_absence_of_the_whole_collection/11357837
library(lme4)
library(fastDummies)
library(glmnet)
library(stringr)
library(readxl)
library(quadprog)
source('easyNF.R')
source('plotNF.R')
#install.packages('fastDummies')
#install.packages("glmnet")

#https://figshare.com/projects/Gladstone_et_al_MGEN_2019/69173
#freq of cloud genes, etc https://figshare.com/articles/dataset/Roary_count_of_core_and_accessory_genes_of_the_whole_collection/11357828
# d1<-read.csv('./Data/gene_presence_absence_minimised.csv')
# # d1[,-c(1:14)] <- apply(d1[,-c(1:14)] ,1,function(x){
# #       x[is.na(x)] <- 0
# #       return(x)
# #     }  )
# d1$Gene <- gsub('group_', 'group', d1$Gene)
# gene_name <-  str_extract(d1$Gene, "[^_]+")
# 
# #d1.spl <- split(d1[-c(1:14)], gene_name)
# library(reshape2)
# library(dplyr)
# 
# ##STEP 1: Filter to COGs that are 5-95% frequency among the whole GPS collection (as in croucher and Colijn)
# d2 <- cbind.data.frame(gene_name,d1[,-c(1:14)])
# #d2.m <- melt(d2, id.vars=c('gene_name'))
# #d2.m <- d2.m[!is.na(d2.m$value),]
# 
# cog_freq <- apply(d2,1, function(x) sum(as.numeric(as.character(x[-1])) , na.rm=T))/ncol(d2)
# 
# ###Gladstone: A further 1957 genes were classified as shell genes (>=15 to <95%) and 24219 as cloud genes (>=1 to <15%). The average number of genes defined as core (>=95%) within lineages (GPSCs) was 1276.
# d3 <- d2[cog_freq>0.05 & cog_freq<0.95,] #filter to just intermediate-frequency COGs
# d3.t <- as.data.frame(t(d3[,-1]))
# names(d3.t) <- paste0('gene',names(d3.t))
# saveRDS(d3.t, './Data/int_freq_GPS.rds')

d3.t <-readRDS('./Data/int_freq_GPS.rds')

d3.t$id2 <- gsub('.velvet','',row.names(d3.t))

## STEP 2: BRING IN GPS META-DATA, PULL OUT CARRIAGE/IPD ISOLATE PRE-PCV
meta1a<-read_excel('./Data/gladstone gps.xlsx')
meta1a$id2 <- gsub('#','.', meta1a$ID)

meta1b<-read_excel('./Data/gladstone gps.xlsx', sheet='T2-GPSC assignment dataset')
meta1b$id2 <- gsub('#','.', meta1b$Taxon)

meta2 <- merge(meta1a, meta1b, by='id2')
meta2 <- meta2[,c('id2',"Clinical Manifest",'Year','Country.x','CC',"In_Silico_Serotype","GPSC.x","Vaccine_Period","Age_group")]
meta2$id2 <- paste0('X', meta2$id2)
d4 <- merge(d3.t, meta2, by='id2')

table(meta2$Country.x[meta2$`Clinical Manifest`=='Carriage'])
table(meta2$Country.x[meta2$`Clinical Manifest`=='Disease']) #S Africa, US most represented


d4.gamb.pre <- d4[d4$Country.x=='The Gambia' & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV',]
d4.gamb.post <- d4[d4$Country.x=='The Gambia' & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Post-PCV7',]

d4.israel.pre.IPD <- d4[d4$Country.x=='Israel' & d4$`Clinical Manifest`== "Disease" & d4$Vaccine_Period=='Pre-PCV' & d4$Age_group %in% c('<=2','>2<=5'),]

##################
country="The Gambia"
G <- d4[,grep('gene',names(d4))]
G <- G[,-which(names(gene.mat)=='geneid2')]
gene.mat <- G[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV',]
gene.mat <- apply(gene.mat, 2, function(x){ x[is.na(x)] <-0 
return(as.numeric(as.character(x)))
})

sample.serotypes <- d4$In_Silico_Serotype[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV']
sample.GPSCs <- d4$GPSC[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV']

serotypes.period2 <-d4$In_Silico_Serotype[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Post-PCV7']
set.vts <- c('4','6A','6B','9V','14','18C','19F','23F')
set.rr <-0.2

nf.obj <- list('sample.GPSCs'=sample.GPSCs,'sample.serotypes'=sample.serotypes,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr, 'post.serotypes'=serotypes.period2, 'country'='The Gambia')

#Run model
preds <- easyNF(nf.obj)

#Generate plots of observed and expected and calculate correlation for NVTs
plots <- plotNF(preds)

#Print plot
plots$ObsExpPlot

#Print correlations
plots$correlations

