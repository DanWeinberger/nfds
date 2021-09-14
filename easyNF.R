###easyNF
#Input: 
#gene.mat: an N*G matrix with N samples and G genes at intermediate frequency (5-95%); entries should be 0/1 for present/absent in that isolate
#sample.serotypes, a vector of length N giving the serotypes in gene.mat
#sample.GPSCs: a vector of length N giving the GPSCs in gene.mat
#set.vts: a vector listing which serotypes are in vaccine e.g. c('1','3','4','5')
#set.rr: population-level impact of vaccine against carriage (rr=0 is max; rr=1 is no effect)


easyNF <- function(nf.obj    ){
  
  sample.GPSCs=nf.obj$sample.GPSCs
  sample.serotypes=nf.obj$sample.serotypes
  #sample.GPSCs <- paste(sample.serotypes, sample.GPSCs, sep='_')
  gene.mat=nf.obj$gene.mat
  pcvsts=nf.obj$set.vts
  irr=nf.obj$set.rr
  
  st.freq.pre <- as.data.frame(table(sample.serotypes)/nrow(gene.mat))
  names(st.freq.pre) <- c('st', 'obs.pre.prev')
  
  pcvst.ind<- as.numeric(sample.serotypes %in% pcvsts)
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
  
  st.gps.all<- as.matrix(table(sample.serotypes, sample.GPSCs))
  K0.vax <- aggregate(gene.mat.vax.effect, by=list('gpsc'=sample.GPSCs), FUN=mean)
  K.vax <- as.matrix(K0.vax[,-1])
  
  st.gps.nvt<- as.matrix(table(sample.serotypes[pcvst.ind==0], sample.GPSCs[pcvst.ind==0]))
  K0.nvt <- aggregate(gene.mat[pcvst.ind==0,], by=list('gpsc'=sample.GPSCs[pcvst.ind==0]), FUN=mean)
  K.nvt <- as.matrix(K0.nvt[,-1])
  
  gpscs.analyzed <- sort(unique(sample.GPSCs[pcvst.ind==0]))
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
  out.obj <- list('val1'=val1, 'nf.obj'=nf.obj,'gpsfreq'=prev.gps,'gpsc.analyzed'=gpscs.analyzed)
  #val1 <- merge(st.freq.pre, val1, by='st')
  #val1$rr.st <- val1$obs.post.prev/val1$obs.pre.prev
  return(out.obj)
  
}
