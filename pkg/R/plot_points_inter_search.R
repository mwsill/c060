
###########################################################################################################
plot.points.int.search<-function(summary.int, start.N=21,cex=1.5){
  # plot visited points vs. iteration steps
  niter<-nrow(summary.int$info) - start.N
  iter<-c(rep(0,start.N), c(1:niter))
  plot(summary.int$info$alpha, iter, xlab=expression(alpha), ylab="Iteration", pch=20,cex=cex)
  grid(NA, niter+1, lwd=2)
  abline(v=summary.int$opt.alpha, col="red")
}