plotNF <- function(nf.out.obj){
  preds <- nf.out.obj
  st.freq.post2 <- as.data.frame(table(preds$nf.obj$post.serotypes) /sum(table(preds$nf.obj$post.serotypes)))
  names(st.freq.post2) <- c('st', 'obs.post.prev2')
  
  preds2 <- merge(nf.out.obj$val1, st.freq.post2, by='st', all=T)
  preds2$pred.prev.post[is.na(preds2$pred.prev.post)] <- 0
  preds2$obs.post.prev2[is.na( preds2$obs.post.prev2)] <- 0
  preds2$obs.pre.prev[is.na( preds2$obs.pre.prev)] <- 0
  
    
  preds2$rr.st2 <- preds2$obs.post.prev2/preds2$obs.pre.prev

  plot.nvt <- which(!(preds2$st %in% nf.out.obj$nf.obj$set.vts))
  
  yrange <-range(c( preds2$obs.post.prev2,preds2$pred.prev.post))
  preds.nvt <- preds2[plot.nvt,]
  preds.nvt$rescale.pre <- preds.nvt$obs.pre.prev/sum(preds.nvt$obs.pre.prev, na.rm=T)
  yrange2 <-range(c( preds.nvt$rescale.pre, preds.nvt$preds.nvt,preds2$pred.prev.post,preds2$obs.post.prev2,preds2$pred.prev.post))
  
  preds2$resid <- preds2$obs.post.prev2/preds2$pred.prev.post
  
  par(mfrow=c(1,1))
  p1 <- ggplot(preds2, aes(pred.prev.post, obs.post.prev2, label = st)) +
     geom_text() +
    xlim(yrange2) +
    ylim(yrange2) +
    labs(x = "Predicted prevalence post-PCV", y='Observed prevalence post-PCV', title=nf.out.obj$nf.out.obj$nf.obj$country ) +
    geom_abline(slope=1, intercept=0, color='gray', linetype='dashed') +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  
  p2 <- ggplot(preds.nvt, aes(rescale.pre, obs.post.prev2, label = st)) +
    geom_text() +
    xlim(yrange2) +
    ylim(yrange2) +
    labs(x = "Simple-expand-Pre", y='Observed prevalence post-PCV', title=nf.out.obj$nf.out.obj$nf.obj$country ) +
    geom_abline(slope=1, intercept=0, color='gray', linetype='dashed') +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  #Spearman Correlation for NVTs 
  #print('NFDS')
  nfds.corr1 <- cor(preds2$obs.post.prev2[plot.nvt],preds2$pred.prev.post[plot.nvt], method='spearman')
  #print('Simple expansion')
  simple.corr <- cor(preds2$obs.post.prev2[plot.nvt],preds2$obs.pre.prev[plot.nvt], method='spearman')
  
  compare.cors <- c(nfds.corr1,simple.corr)
  names(compare.cors) <- c('NFDS','Simple Expanasion')
  out.obj=list ('correlations'=compare.cors, 'ObsExpPlot'=p1,'PrePost'=p2, 'resid'=preds2[,c('st','resid')])
  
}