
###########################################################################################################
plot.summary.int.search<-function(summary.int){
  summary.int$info$log_lambda <- log(summary.int$info$lambda)
  breaks                        <- do.breaks(range(summary.int$info$deviance), 20)
  n_init                        <- 21 # number of initial alpha values at iteration zero
  summary.int$info$cols       <- level.colors(summary.int$info$deviance,at=breaks, col.regions = gray.colors)
  
  print(my.plot <- xyplot(log_lambda~alpha, data = summary.int$info, groups = cols, cex = 2, col = "black",  
               jitter.y=T, amount=0.01, ylab=expression(paste("log ",lambda)),xlab=expression(alpha),
               #             scales=list(x=list(log=T, equispaced.log = FALSE)), # x axis on log-scale
               panel = function(x, y, groups, ..., subscripts) { 
                 fill <- groups[subscripts] 
                 panel.grid(h = -1, v = -1)
                 panel.abline(h = log(summary.int$opt.lambda), v = summary.int$opt.alpha, col="red", lty = "dotted" )
                 panel.xyplot(x, y, pch = rep(c(22,21),c(n_init,nrow(summary.int$info)-n_init)), fill = fill, ...) ;
                 ltext(x=x, y=y, labels=n.features, pos=ifelse(y<0.1,3,4), offset=3, cex=1,col=1)},
               legend = list(top = list(fun = draw.colorkey, args = list(key = list(space = "top", col = gray.colors, at = breaks), draw = FALSE))),
               main="Cross-validated partial log likelihood deviance",
               #sub="number of selected features are printed next to symbol \n rectangles show initial alpha values"
  ))
}
