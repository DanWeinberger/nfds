#Data from GPS, figshare: https://figshare.com/articles/dataset/Roary_gene_presence_and_absence_of_the_whole_collection/11357837
library(lme4)
library(fastDummies)
library(glmnet)
library(stringr)
library(readxl)
library(quadprog)
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
# cog_freq1 <- apply(d2[1:20000,-1],1, function(x) sum(x , na.rm=T))/(ncol(d2)-1)
# cog_freq2 <- apply(d2[20001:nrow(d2),-1],1, function(x) sum(x , na.rm=T))/(ncol(d2)-1)
# cog_freq <- c(cog_freq1,cog_freq2)
# 
# ###Gladstone: A further 1957 genes were classified as shell genes (>=15 to <95%) and 24219 as cloud genes (>=1 to <15%). The average number of genes defined as core (>=95%) within lineages (GPSCs) was 1276.
# d3 <- d2[cog_freq>0.05 & cog_freq<0.95,] #filter to just intermediate-frequency COGs
# d3.t <- as.data.frame(t(d3[,-1]))
# names(d3.t) <- paste0('gene_',d3$gene_name)
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
for(country in c("The Gambia", 'South Africa')){
G <- d4[,grep('gene',names(d4))]
gene.mat <- as.matrix(G[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV',])
gene.mat[is.na(gene.mat)] <- 0

sample.serotypes <- d4$In_Silico_Serotype[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV']
sample.GPSCs <- d4$GPSC[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Pre-PCV']

serotypes.period2 <-d4$In_Silico_Serotype[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Post-PCV7']
serotypes.period3 <-d4$In_Silico_Serotype[d4$Country.x==country & d4$`Clinical Manifest`=='Carriage' & d4$Vaccine_Period=='Post-PCV13']
set.vts <- c('4','6A','6B','9V','14','18C','19F','23F')
set.rr <-0.3

nf.obj <- list('sample.GPSCs'=sample.GPSCs,'sample.serotypes'=sample.serotypes,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr)

#Run model
preds <- easyNF(nf.obj)

st.freq.post2 <- as.data.frame(table(serotypes.period2) /sum(table(serotypes.period2)))
names(st.freq.post2) <- c('st', 'obs.post.prev2')

st.freq.post3 <- as.data.frame(table(serotypes.period3)/sum(table(serotypes.period3)))
names(st.freq.post3) <- c('st', 'obs.post.prev3')

preds2 <- merge(preds, st.freq.post2, by='st')
preds2 <- merge(preds2, st.freq.post3, by='st')

preds2$rr.st2 <- preds2$obs.post.prev2/preds2$obs.pre.prev
preds2$rr.st3 <- preds2$obs.post.prev3/preds2$obs.pre.prev

par(mfrow=c(1,1))
yrange <-range(c( preds2$obs.post.prev2,preds2$pred.prev.post))
plot( preds2$obs.post.prev2,preds2$pred.prev.post, col='white', xlim=yrange, ylim=yrange)
text(  preds2$obs.post.prev2,preds2$pred.prev.post, preds2$st)
abline(a=0, b=1)

yrange <-range(c( preds2$obs.post.prev3,preds2$pred.prev.post))
plot( preds2$obs.post.prev3,preds2$pred.prev.post, col='white', xlim=yrange, ylim=yrange)
text(  preds2$obs.post.prev3,preds2$pred.prev.post, preds2$st)
abline(a=0, b=1)

}

