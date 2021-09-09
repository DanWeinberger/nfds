#NOTE: for Maela, a subset of 674 isolates are sampled (see Colijn); all unvaccinated population
#Can look at Cambodia carriage survey post-vax for comparison:https://pubmed.ncbi.nlm.nih.gov/31175819/#&gid=article-figures&pid=figure-2-uid-1

#Data from GPS, figshare: https://figshare.com/articles/dataset/Roary_gene_presence_and_absence_of_the_whole_collection/11357837
library(lme4)
library(fastDummies)
library(glmnet)
library(rmatio)
library(readxl)
library(quadprog)
source('easyNF.R')
#install.packages('fastDummies')
#install.packages("glmnet")

#Effect of pcv against carriage
set.rr <-0.3

set.vts <- c('1','3','4','5','6A','6B','7F','9V','14','18C','19A','19F','23F')

d1<-read.mat('./Data/redmaela.mat')
G <- as.data.frame(d1[[1]]$G[[1]]) #already filtered for intermediate frequency loci
names(G) <- paste0('gene', 1:ncol(G))
  
strain.name <- unlist(d1$redmaela$taxon, )
st <- unlist(d1$redmaela$serotype, )

g2 <- cbind.data.frame(st, strain.name, G)
g2$strain.name <-gsub('_0','',g2$strain.name)

#Meta data from Chewapreecha's original  study
x1a <- read_excel('./Data/NIHMS56749-supplement-2.xls')
x1b <- read_excel('./Data/gladstone gps.xlsx', sheet='T2-GPSC assignment dataset')
x1b <- x1b[x1b$Study=="Chewapreecha et al",]
x1b <- x1b[, c('Taxon','GPSC', 'ERR')]
x2 <- merge(x1a, g2, by.x='isolate_id', by.y='strain.name')
x2 <- merge(x2, x1b, by.x='ENA_accession_no', by.y='ERR')

mass.pre <- x2

gene.mat <- as.matrix(mass.pre[,grep('gene',names(mass.pre))])
gene.mat[is.na(gene.mat)] <- 0

sample.serotypes <- mass.pre$serotype
sample.GPSCs <- mass.pre$GPSC

nf.obj <- list('sample.GPSCs'=sample.GPSCs,'sample.serotypes'=sample.serotypes,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr)

#Run model
preds <- easyNF(nf.obj)



 yrange <-range(c( preds$obs.pre.prev,preds$pred.prev.post))
 plot( preds$obs.pre.prev,preds$pred.prev.post, col='white', xlim=yrange, ylim=yrange, main='Observed pre vs expected post')
 text(  preds$obs.pre.prev,preds$pred.prev.post, preds$st)
 abline(a=0, b=1)

# yrange <-range(c( preds$obs.post.prev,preds$pred.prev.post))
# plot( preds$obs.post.prev,preds$pred.prev.post, col='white', xlim=yrange, ylim=yrange)
# text(  preds$obs.post.prev,preds$pred.prev.post, preds$st)
# abline(a=0, b=1)

