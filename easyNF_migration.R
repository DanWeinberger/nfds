###easyNF
#Input: 
#gene.mat: an N*G matrix with N samples and G genes at intermediate frequency (5-95%); entries should be 0/1 for present/absent in that isolate
#sample.serotypes, a vector of length N giving the serotypes in gene.mat
#sample.GPSCs: a vector of length N giving the GPSCs in gene.mat
#set.vts: a vector listing which serotypes are in vaccine e.g. c('1','3','4','5', '6B')
#set.rr: population-level impact of vaccine against carriage (rr=0 is max; rr=1 is no effect)


easyNF <- function(nf.obj    ){
  x1b <-read_excel('./Data/gladstone gps.xlsx', sheet='T2-GPSC assignment dataset')
  
  meta1 <-read_excel('./Data/gladstone gps.xlsx')
  meta1.carr.pre <- merge(meta1,x1b, by='ERR' )
  #Restrict to GPSCs that were detected in pre-PCV samples
  meta1.pre <- meta1.carr.pre[ meta1.carr.pre$Vaccine_Period=='Pre-PCV', ] 
  gps.GPSC <- meta1.pre$GPSC.x
  gps.serotype <- meta1.pre$In_Silico_Serotype
  gps.st.gpsc.combos <- unique(cbind.data.frame(gps.GPSC, gps.serotype))
  
  sample.GPSCs=nf.obj$sample.GPSCs
  sample.serotypes=nf.obj$sample.serotypes
  #sample.serotypes.gpsc <- paste(sample.serotypes, sample.GPSCs, sep='_')
  gene.mat=nf.obj$gene.mat
  pcvsts=nf.obj$set.vts
  irr=nf.obj$set.rr
  
  #identify serotypes that are on a GPSC observed in sample.GPSCs that have discordant serotypes
  
  #GPSCs and serotypes that are detected in pre-vaccine carriage sample
  # sample.gps.st.gpsc.combos <- unique(cbind.data.frame(sample.GPSCs, sample.serotypes))
  # sample.gps.st.gpsc.combos.spl <- split(sample.gps.st.gpsc.combos, sample.gps.st.gpsc.combos$sample.GPSCs)
  # sample.gps.st.gpsc.combos.spl <- sapply(sample.gps.st.gpsc.combos.spl, function(x) x$sample.serotypes)
  # 
  #GPSCs that are detected in GPS collection, with various serotypes
  unobs.gps.st.gpsc.combos <- gps.st.gpsc.combos[gps.st.gpsc.combos$gps.GPSC %in% unique(sample.GPSCs),] #only keep GPSCs found in sample
  unobs.gps.st.gpsc.combos.spl <- split(unobs.gps.st.gpsc.combos, unobs.gps.st.gpsc.combos$gps.GPSC)
  unobs.gps.st.gpsc.combos.spl <- sapply(unobs.gps.st.gpsc.combos.spl, function(x) x$gps.serotype)
  unobs.gps.st.gpsc.combos.spl <- lapply(unobs.gps.st.gpsc.combos.spl, function(x) x <- x[!(x %in% pcvsts)])
  unobs.gps.st.gpsc.combos.spl <- sapply(unobs.gps.st.gpsc.combos.spl, function(x) paste(x, collapse=',') )
  unobs.gps.st.gpsc.combos.spl <- cbind.data.frame('gps.GPSC'=names(unobs.gps.st.gpsc.combos.spl), unobs.gps.st.gpsc.combos.spl)
  ##
 # Now need to take serotypes in unobs.gps.st.gpsc.combos.spl and randomly swap them into the corresponding GPSC in sample.serotypes (just make sure we don't replace same thing with itself)
  gpsc.st <- cbind.data.frame(sample.GPSCs, sample.serotypes)
  gpsc.st$orig.order <- 1:nrow(gpsc.st)
  set.seed(123)
  gpsc.st <- merge(gpsc.st,unobs.gps.st.gpsc.combos.spl,by.x='sample.GPSCs' ,by.y='gps.GPSC')
  
  rand.select <- t(apply(gpsc.st, 1, function(x){
    replace.st.vec <- strsplit(x['unobs.gps.st.gpsc.combos.spl'],',', fixed=T)[[1]]
    if(length(replace.st.vec)>0){
    x['replace.st'] <- sample(replace.st.vec,1, replace=T)
    }else{
      x['replace.st'] <-"MISS"
    }
    randnum <- runif(n=1)
    if(randnum<=0.1 &  x['replace.st'] != x["sample.serotypes"] & x['replace.st'] != 'MISS' ){
      x["sample.serotypes"] <- x['replace.st']
    }
      return(x)
  }))
  rand.select <- rand.select[order(rand.select[,'orig.order']),]
  sample.GPSCs.replace <- as.numeric(rand.select[,'sample.GPSCs'])
  sample.st.replace <- rand.select[,'sample.serotypes']
  
  ###
  
  #Note use OBSERVED serotypes pre, not the migrated STs
  st.freq.pre <- as.data.frame(table(sample.serotypes)/nrow(gene.mat))
  names(st.freq.pre) <- c('st', 'obs.pre.prev')
  
  pcvst.ind<- as.numeric(sample.st.replace %in% pcvsts)
  nvt.ind <- 1-pcvst.ind
  
  #1. we have a 'target', which is the frequency of each gene pre-PCV. This is'e_l'
  e_l <- apply(gene.mat, 2,mean, na.rm=T)
  
  #2. We have a set of (mostly) non-vaccine type strains, with info on the gene in each. 
  irr.carr.pcv <- rep(1, nrow(gene.mat))
  #irr.carr.pcv[pcvst.ind==1] <- irr
  irr.carr.pcv[pcvst.ind==1] <- irr
  
  
  ###
  
  #This matrix reduces the influence of genes in VTs; reduced based on IRR.carr.pcv
  gene.mat.vax.effect <- apply(gene.mat, 2, function(x) x*irr.carr.pcv )
  
  st.gps.all<- as.matrix(table(sample.st.replace, sample.GPSCs.replace))
  K0.vax <- aggregate(gene.mat.vax.effect, by=list('gpsc'=sample.GPSCs.replace), FUN=mean)
  K.vax <- as.matrix(K0.vax[,-1])
  
  st.gps.nvt<- as.matrix(table(sample.st.replace[pcvst.ind==0], sample.GPSCs.replace[pcvst.ind==0]))
  K0.nvt <- aggregate(gene.mat[pcvst.ind==0,], by=list('gpsc'=sample.GPSCs.replace[pcvst.ind==0]), FUN=mean)
  K.nvt <- as.matrix(K0.nvt[,-1])
  gpscs.analyzed <- sort(unique(sample.GPSCs.replace[pcvst.ind==0]))
  k_t <- t(K.nvt)
  
  
  #3. We want to determine the combination of K/vax that will produce e_l
  # Quadratic programming https://stats.stackexchange.com/questions/21565/how-do-i-fit-a-constrained-regression-in-r-so-that-coefficients-total-1
  #(this is basically the approach taken by Taj Azarian's paper)
  #Why quadratic program instead of linear regression?: it forces coefficients to add to 1 and be > 0, which 
  #is what we want because they represent prevalence
  
  ## this bit of code is for the quadratic programming. Output is a vector with coefficients, representing expected prevalence for each GPSC
  X<- as.matrix(k_t)
  X[X==0] <- 1e-5
  Y <- matrix(e_l, ncol=1)
  #x.nonzero <- X[,colSums(X)>0] #restrict to GPS that are >0 in K.vax

  Rinv <- solve(chol(t(X) %*% X))
  C <- cbind(rep(1,ncol(X)), diag(ncol(X)))
  b <- c(1,rep(0,ncol(X)))
  d <- t(Y) %*% X  
  res1 <-solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  prev.gps <- res1$solution #estimated prevalence of the GPS post-PCV

  # multiply coefficient for by the pre-vaccine proportion of each GPS caused by different serotypes
  #irr=0 means vaccine is 100% effective against carriage
  irr2 <- rep(1, nrow(st.gps.all))
  irr2[dimnames(st.gps.all)[[1]] %in% pcvsts] <- irr 
  
  st.prop.gps <- apply(st.gps.nvt,2, function(x) x/sum(x))
  
  pred_post_serotype_prev <- round(st.prop.gps %*% matrix(prev.gps, ncol=1),2)
  pred_post_serotype_prev <- cbind.data.frame(st=row.names(pred_post_serotype_prev),'pred.prev.post'=pred_post_serotype_prev[,1])
  
  ## Compare observed and predicted ST prevalence
  # see https://www.nature.com/articles/s41564-019-0651-y/figures/8?proof=trueIn%EF%BB%BF
  #Note the study from 
  val1 <- merge(pred_post_serotype_prev, st.freq.pre, by='st')
  val1$expected.rr <- val1$pred.prev.post / val1$obs.pre.prev
  out.obj <- list('val1'=val1, 'nf.obj'=nf.obj, 'gpsfreq'=prev.gps,'gpsc.analyzed'=gpscs.analyzed)
  #val1 <- merge(st.freq.pre, val1, by='st')
  #val1$rr.st <- val1$obs.post.prev/val1$obs.pre.prev
  return(out.obj)
  
}
