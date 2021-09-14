##evaluate correlation between presence, absence of COGs and Ps components
library(pbapply)
library(readxl)
library(lme4)
library(stringr)
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

#Pull in ps structure data by Serotype
ps1 <- read.csv('./Data/ps_composition.csv')
ps1$Serotype <- as.character(ps1$Serotype)
ps1$Serotype[ps1$Serotype %in% c('15B','15C')] <-'15BC'
ps2 <- aggregate(ps1[,-1], by=list('st'=ps1$Serotype) ,FUN=mean)


d5 <- merge(d4, ps2, by.x='In_Silico_Serotype', by.y='st', all=T)
ac.sugs <- c("Ac","FucNAc" ,"GalNAc", "GlcNAc","ManNAc",'ManNAcA','PneNAc')

d5$NaC <- 0
d5$NaC[apply(d5[,ac.sugs],1,sum) >1] <- 1
cog.cols <-  grep('gene_',names(d5))

df.nac <- d5[,c('NaC', 'In_Silico_Serotype','GPSC.x',names(d5)[cog.cols]) ]

df.nac[,-c(1:3)][is.na(df.nac[,-c(1:3)])] <- 0
  

#first Just do a series of logistic regressions to screen

st.freq <- as.data.frame(table(df.nac$In_Silico_Serotype))
names(st.freq) <-c('st','st_freq')
df.nac <- merge(df.nac,st.freq, by.x='In_Silico_Serotype', by.y='st')
df.nac$In_Silico_Serotype[df.nac$st_freq<10] <- 'Other'

df.nac$In_Silico_Serotype <- as.factor(df.nac$In_Silico_Serotype)

mod1 <- glm(NaC ~ gene_lsrC + In_Silico_Serotype   , data=df.nac, family=binomial)

names(df.nac) <- gsub('-','_', names(df.nac), fixed=T)
names(df.nac) <- gsub('/','_', names(df.nac), fixed=T)
names(df.nac) <- gsub('(','', names(df.nac), fixed=T)
names(df.nac) <- gsub(')','', names(df.nac), fixed=T)
names(df.nac) <- gsub(' ','', names(df.nac), fixed=T)
names(df.nac) <- gsub('\\$','', names(df.nac))
names(df.nac) <-	gsub("as.factor(df.nac$",'',names(df.nac), fixed=T)
names(df.nac) <- gsub('-','_', names(df.nac))
names(df.nac) <- gsub(' ','_', names(df.nac))
df.nac$NaC <- as.factor(df.nac$NaC)
#form.fixed <- as.formula(paste0('NaC~ as.factor(In_Silico_Serotype) +', paste0(gen.names, collapse='+') ))
mod.df <- cbind.data.frame('NaC'=df.nac$NaC, X)

names(df.nac) <- sapply(names(df.nac), function(x) gsub('\\$','_',x))
names(df.nac) <- gsub('as.factor\\(','',names(df.nac))
names(df.nac) <- gsub('\\)','',names(df.nac))
names(df.nac) <- str_replace_all(names(df.nac), "[[:punct:]]", "")

gen.names <- names(df.nac)[grep('gene_',names(df.nac))]
# mod.fun <- pblapply(gen.names, function(x){
#               form1 <- as.formula(paste0('NaC ~ In_Silico_Serotype +',x))
#                   glm(form1   , data=df.nac, family=binomial)
#                   })
# coef1 <- sapply(mod.fun, '[[', 'coefficients')
# coef1 <- coef1[nrow(coef1),] #extract value for last coefficient, which is the one of interest
# 
# aics <- sapply(mod.fun, '[[', 'aic')
# aics <- cbind.data.frame(aics, coef1)
# aics <- cbind.data.frame(gen.names,aics)
# aics$Gene <- gsub('gene_','',aics$gen.names)
 cog.ann <- readRDS('./Data/cog_annotations.rds')
# aics$Gene <- gsub('group','group_',aics$Gene)
# 
# aics <- merge(aics, cog.ann, by='Gene', all.x=T)
# saveRDS(aics, './Results/univariate.glm.rds')


#Try to do same but adjust for GPSC

# gen.names <- names(df.nac)[grep('gene_',names(df.nac))]
# df.nac$GPSC.x <- as.factor(df.nac$GPSC.x)
# mod.fun2 <- pblapply(gen.names, function(x){
#   form1 <- as.formula(paste0('NaC ~ (1|In_Silico_Serotype) + (1|GPSC.x) + ',x))
#   glmer(form1   , data=df.nac, family=binomial,control=glmerControl(optimizer="bobyqa",
#                                                                    optCtrl=list(maxfun=2e5)))
# })


#Regular elastic net--accounts for serotype and GPSC as features
library(glmnet)
df.nac <- df.nac[!is.na(df.nac$GPSCx),]
st.mat <- model.matrix(~ df.nac$InSilicoSerotype)
gps.mat <- model.matrix(~ as.factor(df.nac$GPSCx))

X= cbind(as.matrix(df.nac[,gen.names]),st.mat[,-1], gps.mat[,-1])
Y= as.factor(df.nac$NaC)
fit = glmnet(x=X, y=Y, family = "binomial")
plot(fit, label = TRUE)
coefs2 <- as.data.frame(as.matrix(coef.glmnet(fit, s=0.0005)))
names(coefs2) <- 'beta'
coefs2$Gene <- gsub('gene_','',rownames(coefs2))
coefs2$Gene <- gsub('group','group_',coefs2$Gene)
cog.ann <- readRDS('./Data/cog_annotations.rds')
coefs2 <- merge(coefs2, cog.ann, by='Gene', all.x=T)
View(coefs2[coefs2$beta!=0 & 
              !grepl('Serotype',coefs2$Gene) &
              !grepl('GPSC', coefs2$Gene),])


#Not sure it is correct to adjust for serotype
##there is a direct correspondend between Serotype and Ps structure
##Serotype completely determines Ps structure (or is it the other way around)--
#If you are serotype X, you definitely have components A,B,C;
#if you have component, B, you are not necesarily serotype X
#adjusting for St doesn't make sense then, because it will be 100% yes or no
#but we still might want to adust correlation structure for non-independence
df.nac <- df.nac[!is.na(df.nac$GPSCx),]
gps.mat <- model.matrix(~ as.factor(df.nac$GPSCx))
gen.names <- names(df.nac)[grep('gene',names(df.nac))]
#As in pyseer elasticnet, add weight, which is 1/N_GPS

freq.gpsc <- cbind.data.frame('gps.freq'=table(df.nac$GPSCx))
names(freq.gpsc) <- c('GPSCx','gps.freq')

df.nac <- merge(df.nac,freq.gpsc, by='GPSCx' )
df.nac$gps.wgt <- 1/df.nac$gps.freq
X= cbind(as.matrix(df.nac[,gen.names])) #, gps.mat[,-1])
Y= as.factor(df.nac$NaC)

fit = glmnet(x=X, y=Y, family = "binomial", weights=df.nac$gps.wgt)

plot(fit, label = TRUE)
coefs2 <- as.data.frame(as.matrix(coef.glmnet(fit, s=0.001)))
names(coefs2) <- 'beta'
coefs2$Gene <- gsub('gene_','',rownames(coefs2))
coefs2$Gene <- gsub('group','group_',coefs2$Gene)
cog.ann <- readRDS('./Data/cog_annotations.rds')
cog.ann$Gene <- str_replace_all(cog.ann$Gene, "[[:punct:]]", "")
coefs2$Gene <- str_replace_all(coefs2$Gene, "[[:punct:]]", "")
cog.ann$Gene  <- paste0('gene',cog.ann$Gene )
coefs2 <- merge(coefs2, cog.ann, by='Gene', all.x=T)
coefs2$beta <- round(coefs2$beta,3)
View(coefs2[coefs2$beta!=0 & 
              !grepl('GPSC', coefs2$Gene),])

#FnlABC UDP=FucNAC (UDP-2 acetamio 26-d-deoxy...)



#################################

#install.packages('glmmLasso')
library(glmmLasso)#
library(stringr)
gen.names <- names(df.nac)[grep('gene_',names(df.nac))]
df.nac$NaC <- as.factor(df.nac$NaC)
#form.fixed <- as.formula(paste0('NaC~ as.factor(In_Silico_Serotype) +', paste0(gen.names, collapse='+') ))
mod.df <- cbind.data.frame('NaC'=df.nac$NaC, X)

names(df.nac) <- sapply(names(df.nac), function(x) gsub('\\$','_',x))
names(df.nac) <- gsub('as.factor\\(','',names(df.nac))
names(df.nac) <- gsub('\\)','',names(df.nac))
names(df.nac) <- str_replace_all(names(df.nac), "[[:punct:]]", "")


xnames <- names(df.nac)[grep('gene',names(df.nac))]
xvars <- paste0(xnames, collapse='+')
form.fixed <- as.formula(paste0('NaC~ ',xvars ))

#mod3 <- glmmLasso( fix=form.fixed, rnd = list(In_Silico_Serotype=~1,GPSC.x=~1 ),family = "binomial", data=df.nac)
#mod3 <- glmmLasso( fix=form.fixed, rnd = NULL,family = "binomial", data=df.nac)


#Alternative approach...randomly select 1 representative of each GPSC and repeatedly 



#Market basket
library(arules)
spl1 <- split(df.nac, df.nac$GPSCx)
samp1 <- lapply(spl1, function(x){
  nrows <- nrow(x)
  row.select <- sample(1:nrows, size=1)
  return(x[row.select,,drop=F])
})
samp1 <- do.call( 'rbind.data.frame', samp1)
gene.cols <- c(grep('gene',names(samp1)), which(names(samp1)=='NaC'))

#Convert dataset to transcation data type for market basket analysis
ps_mat <- as.matrix(samp1[,gene.cols])==1
ps_tran <- as(ps_mat,"transactions")
# rules_2 <- apriori(ps_tran) #,parameter = list(supp = 0.0001, conf = 0.0001,minlen=2, maxlen=8))
# inspect(head(rules_2, by = "lift"))
# # summary(ps_tran)  


#PCA first, then regression
#https://cran.r-project.org/web/packages/logisticPCA/vignettes/logisticPCA.html
library(logisticPCA)
gene.cols <- grep('gene',names(df.nac))
gen.mat <- df.nac[,gene.cols]
#Cross validation to determine minimum m value
#logpca_cv = cv.lpca(gen.mat, ks = 2, ms = 1:10)
plot(logpca_cv)
pca.mod1 <-  logisticPCA(gen.mat, k = 5, m = 5)

#Use minimm m value
library(ggplot2)
plot(pca.mod1, type = "scores") + geom_point(aes(colour = df.nac$InSilicoSerotype))

cps.genes <- rep(0, ncol(gen.mat))
cps.genes[c(grep('wz',colnames((gen.mat))) ,grep('cps',colnames((gen.mat))) ,grep('wc',colnames((gen.mat))), grep('cap',colnames((gen.mat))))] <- 1
plot(pca.mod1, type = "loadings") + geom_point(aes(colour = cps.genes))

extremes <- dimnames(gen.mat)[[2]][which(pca.mod1$U[,1]< (-0.05))]
View(extremes)
