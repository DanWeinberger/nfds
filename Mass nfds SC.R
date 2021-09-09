#Issues: if 6A is not included as a VT, the model fits very well for all sts (except 6A)
#If 6A is included, fit is worse. Is this due to small numbers issues?
##Follow-up: if we use GPSC instead of the clonal grouping in the Colijn paper,
#it is much more robust (33 groups of GPSC vs 14 for the strain clusters)

#Data from GPS, figshare: https://figshare.com/articles/dataset/Roary_gene_presence_and_absence_of_the_whole_collection/11357837
library(lme4)
library(fastDummies)
library(glmnet)
library(rmatio)
library(patchwork)
library(quadprog)
library(readxl)
source('plotNF.R')
source('easyNF_migration.R')

#install.packages('fastDummies')
#install.packages("glmnet")

#Effect of PCV7 against carriage
set.rr <-0.1
set.vts <- c('4','6A','6B','9V','14','18C','19F','23F')
vt.eff <-  c(0.0, 0.1,0.1, 0, 0, 0, 0.2, 0.1 )


d1<-read.mat('./Data/massdata.mat')
G <- as.data.frame(d1[[1]]$G[[1]]) #already filtered for intermediate frequency loci
names(G) <- paste0('gene', 1:ncol(G))
locus_weight <- d1$massdata$locusweights[[1]]  


strain.name <- unlist(d1$massdata$taxon, )
st <- unlist(d1$massdata$serotype, )

g2 <- cbind.data.frame(st, strain.name, G)
g2$strain.name <-gsub('_0','',g2$strain.name)

#Meta data from croucher's original Nat Gen study
x1a <- read_excel("C:/Users/dmw63/Desktop/My documents h/LAB/pneumo metabolic genes/NIHMS474991-supplement-2.xlsx")
x1b <-read_excel('./Data/gladstone gps.xlsx', sheet='T2-GPSC assignment dataset')
x1c <- x1b[x1b$Study=="Croucher et al",]
x1c <- x1c[, c('Taxon','GPSC', 'ERR')]
x1a$Serotype <- NULL
x2 <- merge(x1a, g2, by.x='Taxon ID', by.y='strain.name')
x2 <- merge(x2, x1c, by.x='Accession', by.y='ERR')

meta1 <- read_excel('./Data/gladstone gps.xlsx')
meta2 <-merge(meta1, x1b, by='ERR')
meta2.carr <- meta2[meta2$`Clinical Manifest`=='Carriage' & meta2$Vaccine_Period=='Pre-PCV',]

table(meta2$Country.x[meta2$`Clinical Manifest`=='Disease' & meta2$Vaccine_Period=='Post-PCV13'])
post.ipd <- meta2[meta2$`Clinical Manifest`=='Disease' & meta2$Vaccine_Period=='Post-PCV13',]

mass.pre <- x2[x2$`Year of Isolation`==2001,]
mass.post <- x2[x2$`Year of Isolation`==2007,]
gene.mat <- as.matrix(mass.pre[,grep('gene',names(mass.pre))])
gene.mat[is.na(gene.mat)] <- 0

sample.serotypes <- mass.pre$st
sample.GPSCs <- mass.pre$GPSC
sample.SC <- mass.pre$`Strain Cluster (SC)`
sample.SC.post <- mass.post$`Strain Cluster (SC)`

sample.GPSCs.post <- mass.post$GPSC
sample.serotypes.post <- mass.post$st

serotypes.post <- mass.post$st

df <-cbind.data.frame(x2$`Strain Cluster (SC)`, x2$GPSC, x2$`Consensus serotype`, 'year'=x2$`Year of Isolation`)
names(df) <-c('sc','gpsc','st','year')
test1 <-df[df$gpsc==1,]
table(test1$gpsc, test1$st)

nf.obj <- list('sample.GPSCs'=sample.SC,'sample.serotypes'=sample.serotypes,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr, 'post.serotypes'=serotypes.post, country='US_MA')

#Run model
preds <- easyNF(nf.obj)

#Generate plots of observed and expected and calculate correlation for NVTs
plots <- plotNF(preds)

#Print plots
plots$ObsExpPlot + plots$PrePost

#Print correlations
plots$correlations

#Invasiveness
inv1 <- read.csv('./Data/mcmc_invasive_single_stage.csv')
inv2 <- merge(inv1,preds$val1, by='st', all=T)
inv2$expected.IPD.nfds <- inv2$pred.prev.post * exp(inv2$log.inv.age1)
inv2$carr.nvt <- inv2$obs.pre.prev
inv2$carr.nvt[inv2$st %in% preds$nf.obj$set.vts] <- NA
inv2$carr.nvt <- inv2$carr.nvt/sum(inv2$carr.nvt, na.rm=T)
inv2$expected.IPD.nfds <- inv2$pred.prev.post * exp(inv2$log.inv.age1*1.0)
inv2$expected.IPD.simple <- inv2$carr.nvt * exp(inv2$log.inv.age1*1.0)

plot(inv2$expected.IPD.nfds, inv2$expected.IPD.simple, col='white', bty='l')
text(inv2$expected.IPD.nfds, inv2$expected.IPD.simple, inv2$st)

abline(a=0, b=1)

#Compare SC predictions from simple approach with Colijn's results
pred.orig <- cbind.data.frame('SC'=d1$massdata$finalSC[[1]], 'pred.freq.colijn'=d1$massdata$finalfreqSC[[1]])
names(pred.orig) <- c('SC','pred.freq.colijn')

pred.here <- cbind.data.frame('pred.freq.here'=preds$gpsfreq, 'SC'=preds$gpsc.analyzed)

obs.freq.pre <- table(sample.SC)/sum(table(sample.SC))
obs.freq.pre <- cbind.data.frame(names(obs.freq.pre), obs.freq.pre)
names(obs.freq.pre) <- c('','SC','obs.freq.pre')

obs.freq.post <- table(sample.SC.post)/sum(table(sample.SC.post))
obs.freq.post <- cbind.data.frame(names(obs.freq.post), obs.freq.post)
names(obs.freq.post) <- c('','SC','obs.freq.post')

pred.comp <- merge(pred.here, pred.orig, by='SC', all=T)
pred.comp <- merge(pred.comp, obs.freq.post, by='SC', all=T)
pred.comp <- merge(pred.comp, obs.freq.pre, by='SC', all=T)
pred.comp$pred.freq.colijn[is.na(pred.comp$pred.freq.colijn)] <-0
pred.comp$pred.freq.here[is.na(pred.comp$pred.freq.here)] <-0
pred.comp$obs.freq.post[is.na(pred.comp$obs.freq.post)] <-0
pred.comp$obs.freq.pre[is.na(pred.comp$obs.freq.pre)] <-0

par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(pred.comp$pred.freq.here, pred.comp$pred.freq.colijn, col='white')
text(pred.comp$pred.freq.here, pred.comp$pred.freq.colijn, pred.comp$SC)

plot(pred.comp$pred.freq.here, pred.comp$obs.freq.post, col='white')
text(pred.comp$pred.freq.here, pred.comp$obs.freq.post, pred.comp$SC)

plot(pred.comp$pred.freq.colijn, pred.comp$obs.freq.post, col='white')
text(pred.comp$pred.freq.colijn, pred.comp$obs.freq.post, pred.comp$SC)

plot(pred.comp$obs.freq.pre, pred.comp$obs.freq.post, col='white')
text(pred.comp$obs.freq.pre, pred.comp$obs.freq.post, pred.comp$SC)

#The simpler quad programming approach more accurately captures the SC frequencies than does the more complicated ODE model
#Best predictor of post GPSC is pre-GPSC
cor(pred.comp[,c('obs.freq.post','obs.freq.pre','pred.freq.colijn','pred.freq.here')])

