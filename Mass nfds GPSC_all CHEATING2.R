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
library(ggplot2)
source('plotNF.R')
source('easyNF_all.R')

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

#mass.pre <- x2[x2$`Year of Isolation`==2001,]
mass.pre <- x2

mass.post <- x2[x2$`Year of Isolation`==2007,]
gene.mat <- as.matrix(mass.pre[,grep('gene',names(mass.pre))])
gene.mat[is.na(gene.mat)] <- 0

gene.mat.pre.only <- gene.mat[x2$`Year of Isolation`==2001, ]


sample.serotypes <- mass.pre$st
sample.GPSCs <- mass.pre$GPSC
sample.gpsc.st <- paste(sample.serotypes, sample.GPSCs, sep='_')
sample.SC <- mass.pre$`Strain Cluster (SC)`
sample.SC.post <- mass.post$`Strain Cluster (SC)`
sample.GPSCs.post <- mass.post$GPSC
sample.year <-mass.pre$`Year of Isolation`
serotypes.post <- mass.post$st

df <-cbind.data.frame(x2$`Strain Cluster (SC)`, x2$GPSC, x2$`Consensus serotype`, 'year'=x2$`Year of Isolation`)
names(df) <-c('sc','gpsc','st','year')
test1 <-df[df$gpsc==1,]
table(test1$gpsc, test1$st)

nf.obj <- list('sample.GPSCs'=sample.gpsc.st,'sample.year'=sample.year, 'sample.serotypes'=sample.serotypes,'gene.mat.pre'=gene.mat.pre.only,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr, 'post.serotypes'=serotypes.post, country='US_MA')

#Run model
preds <- easyNF(nf.obj)

#Generate plots of observed and expected and calculate correlation for NVTs
plots <- plotNF(preds)

#Print plots
plots$ObsExpPlot +plots$PrePost

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


## Does this also hold when looking at GPSC?--doesn't seem to look as clean
gpsc.freq.post <- table(sample.GPSCs.post)/sum(table(sample.GPSCs.post))
gpsc.freq.post <- cbind.data.frame(names(gpsc.freq.post), gpsc.freq.post)
names(gpsc.freq.post) <- c('','gpsc','obs.freq.post')

gpsc.freq.pre <- table(sample.GPSCs)/sum(table(sample.GPSCs))
gpsc.freq.pre <- cbind.data.frame(names(gpsc.freq.pre), gpsc.freq.pre)
names(gpsc.freq.pre) <- c('','gpsc','obs.freq.pre')

pred.gps <- cbind.data.frame('gpsfreq'=preds$gpsfreq, 'gpsc.analyzed'=preds$gpsc.analyzed)
gpsc.freq.comp <- merge(gpsc.freq.pre,gpsc.freq.post, by='gpsc',all=T)
gpsc.freq.comp <- merge(gpsc.freq.comp,pred.gps, by.x='gpsc',by.y="gpsc.analyzed",all=T)

gpsc.freq.comp$obs.freq.pre[is.na(gpsc.freq.comp$obs.freq.pre)] <- 0
gpsc.freq.comp$obs.freq.post[is.na(gpsc.freq.comp$obs.freq.post)] <- 0
gpsc.freq.comp$gpsfreq[is.na(gpsc.freq.comp$gpsfreq)] <- 0

#
par(mfrow=c(1,2))
plot(gpsc.freq.comp$obs.freq.pre, gpsc.freq.comp$obs.freq.post, main='Pre vs Post', bty='l')
abline(a=0, b=1)

plot(gpsc.freq.comp$gpsfreq, gpsc.freq.comp$obs.freq.post, main='Expected vs Observed', bty='l')
abline(a=0, b=1)

#For GPSC, the NFDS model performs better than model assuming simple expansion
cor(gpsc.freq.comp[,c('obs.freq.post','obs.freq.pre',"gpsfreq")])



#Look at ST patterns by GPS pre and post
library(reshape2)
sample.combos.pre<-cbind.data.frame(sample.GPSCs, sample.serotypes)
sample.combos.pre$one <- 1
sample.combos.pre.m <- melt(sample.combos.pre,id.vars=c('sample.GPSCs','sample.serotypes'))
sample.combos.pre.c <- dcast(sample.combos.pre.m,sample.serotypes+ sample.GPSCs~., fun.aggregate = sum)
names(sample.combos.pre.c) <-c('st','gpsc','N_pre')

sample.combos.post<-cbind.data.frame(sample.GPSCs.post, sample.serotypes.post)
sample.combos.post$one <- 1
sample.combos.post.m <- melt(sample.combos.post,id.vars=c('sample.GPSCs.post','sample.serotypes.post'))
sample.combos.post.c <- dcast(sample.combos.post.m,sample.serotypes.post+ sample.GPSCs.post~., fun.aggregate = sum)
names(sample.combos.post.c) <-c('st','gpsc','N_post')

samples.combo.comp <- merge(sample.combos.pre.c,sample.combos.post.c, by=c('st','gpsc'), all=T)
samples.combo.comp$N_pre[is.na(samples.combo.comp$N_pre)] <- 0
samples.combo.comp$N_post[is.na(samples.combo.comp$N_post)] <- 0
samples.combo.comp$vt <- 0
samples.combo.comp$vt[samples.combo.comp$st %in% set.vts] <-1

samples.combo.comp.nvt <- samples.combo.comp[!(samples.combo.comp$st %in% set.vts),]
par(mfrow=c(1,1))
plot(samples.combo.comp.nvt$N_pre, samples.combo.comp.nvt$N_post)

library(lme4)
samples.combo.comp$gpsc <- as.factor(samples.combo.comp$gpsc)
samples.combo.comp$st <- as.factor(samples.combo.comp$st)

mod1 <- glmer(N_post ~ scale(log(N_pre+0.5)) + (1|gpsc) + (1|st), family='poisson', data=samples.combo.comp)
pred1 <- predict(mod1, type='response')
summary(mod1)

plot(pred1, samples.combo.comp$N_post )
abline(a=0,b=1)
  
#NNED TO CREATE 3 VARIABLES: PREVACCINE FREQUENCY OF THE GPSC/ST COMBP; TH PREVAX FREQUENCY OF THE GPSC AND PRE-VAX FREQUENCY OF THE SEROTYPE...USE THESE 3 IN A STANDARD GLM INSTEAD OF RANDOM EFFECTS
#Prev by GPSC
prev.gpsc <- dcast(sample.combos.pre.m,sample.GPSCs~., fun.aggregate = sum )
names(prev.gpsc) <- c('gpsc','N.gpsc')
prev.gpsc$prev.gpsc <- prev.gpsc$N.gpsc/length(sample.serotypes)
prev.st <- dcast(sample.combos.pre.m,sample.serotypes~., fun.aggregate = sum )
names(prev.st) <- c('st','N.st')
prev.st$prev.st <- prev.st$N.st/length(sample.serotypes)
prev.combo <- merge(samples.combo.comp,prev.st, by='st')
prev.combo <- merge(prev.combo,prev.gpsc, by='gpsc')
prev.combo$prev.st.gpsc.pre <- prev.combo$N_pre/length(sample.serotypes)

mod2 <- glm(N_post ~ prev.st + prev.gpsc +prev.st.gpsc.pre, data=prev.combo[!(prev.combo$st %in% set.vts),])
summary(mod2)
mod3 <- glm(N_post ~ log(prev.st+0.01) + log(prev.gpsc+0.01) +log(prev.st.gpsc.pre+0.01), data=prev.combo[!(prev.combo$st %in% set.vts),], family='poisson')
summary(mod3)
pred.mod3 <- predict(mod3, type='response')
plot(pred.mod3,prev.combo[!(prev.combo$st %in% set.vts),'N_post'] )
cor(pred.mod3,prev.combo[!(prev.combo$st %in% set.vts),'N_post'])


mod3b <- glm(N_post ~ log(prev.st+0.01) + log(prev.gpsc+0.01), data=prev.combo[!(prev.combo$st %in% set.vts),], family='poisson')
summary(mod3b)
pred.mod3b <- predict(mod3b, type='response')
plot(pred.mod3b,prev.combo[!(prev.combo$st %in% set.vts),'N_post'] )
cor(pred.mod3b,prev.combo[!(prev.combo$st %in% set.vts),'N_post'])

#Knowing the frequency of the GPSC/ST combo pre-PCV is a pretty correlate of post-PCV
mod3c <- glm(N_post ~ log(prev.st.gpsc.pre+0.01), data=prev.combo[!(prev.combo$st %in% set.vts),], family='poisson')
summary(mod3c)
pred.mod3c <- predict(mod3c, type='response')
plot(pred.mod3c,prev.combo[!(prev.combo$st %in% set.vts),'N_post'] )
cor(pred.mod3c,prev.combo[!(prev.combo$st %in% set.vts),'N_post'])

#how does it do predicting ST freq?
pred.reg.st <- aggregate(pred.mod3c, by=list(prev.combo$st[!(prev.combo$st %in% set.vts)]), FUN=sum)
names(pred.reg.st) <- c('st','pred.st.post')

obs.st.post <- table(serotypes.post)
obs.st.post <- cbind.data.frame('st'=names(obs.st.post),obs.st.post)
names(obs.st.post)[3] <- 'obs.freq.post'

comp.st.post <- merge(pred.reg.st, obs.st.post, by='st')

plot(comp.st.post$pred.st.post, comp.st.post$obs.freq.post)
cor(comp.st.post$pred.st.post, comp.st.post$obs.freq.post)
