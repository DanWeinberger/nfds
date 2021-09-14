#NOTE: would be ideal to use GPSCs instead of the BAPS

#This uses the data from azarian et al. PLOS pathogens and Biorxiv
#"For the present study, we focused on 937 pneumococcal carriage isolates collected during three epochs: pre-vaccine – population equilibrium (E1, 1998-2001); peri-vaccine – population perturbation (E2, 2006-2008); post-vaccine – population equilibration (E3, 2010-2012).  E1 preceded the introduction of PCV7, while E2 and E3 provided post-PCV7 snapshots 5 and 10 years, respectively, after the introduction of PCV7.  While E3 includes, in part, the introduction of PCV13, we have previously shown that the majority of the sample was obtained when the impact of PCV13 was minimal (4).  Additionally, we perform a sensitivity analysis to assess the effect of including the complete E3 sample."  

library(lme4)
library(fastDummies)
library(glmnet)
library(rmatio)
library(readxl)
library(ggplot2)
library(quadprog)
library(mgcv)
library(lme4)
library(patchwork)

#source('easyNF_migration.R')
source('easyNF.R')

source('plotNF.R')

#set.vts <- c('1','3','4','5','6A','6B','7F','9V','14','18C','19A','19F','23F')
set.vts <- c('4','6A','6B','9V','14','18C','19F','23F')

set.rr <-0.2

#ds <- read.csv('https://raw.githubusercontent.com/c2-d2/Projects/master/NFDS/data_southwestUS.csv')
#write.csv(ds,'./Data/azarian_navajo.csv')
d1 <- read.csv('./Data/azarian_navajo.csv')
d1$BAPS1 <- as.character(d1$BAPS1)
gps <-  read_excel('./Data/gladstone gps.xlsx', sheet='T21-HierBAPS')
gps <- gps[,c("GPSC...2","HierBAPs")]
names(gps) <-c('GPSC', 'BAPS1')

cogs <- readRDS('./Data/cog_annotations.rds')



gps <- unique(gps)

#d1 <- merge(d1, gps, by='BAPS1',all.x=T)
d1 <- d1[!is.na(d1$BAPS2),]

G <- d1[, 12:ncol(d1)]
freq.gene <- colMeans(G) #already seems to be filtered for intermediate-freq locui

names(G) <- paste0('gene', 1:ncol(G))

gene.mat <- as.matrix(G[d1$Epoch1=='E1',])
gene.mat[is.na(gene.mat)] <- 0

sample.serotypes <- d1$Sero[d1$Epoch1=='E1']
sample.GPSCs <- d1$BAPS1[d1$Epoch1=='E1']

serotypes.period2 <-d1$Sero[d1$Epoch1=='E2']
serotypes.period3 <-d1$Sero[d1$Epoch1=='E3']

nf.obj <- list('sample.GPSCs'=sample.GPSCs,'sample.serotypes'=sample.serotypes,'gene.mat'=gene.mat, 'set.vts'=set.vts,'set.rr'=set.rr,
               'post.serotypes'=serotypes.period3, country='Navajo')

#Run model
preds <- easyNF(nf.obj)

#
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


