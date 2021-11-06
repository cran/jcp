#' plot.jcp
#'
#' Plot method for class 'jcp'
#'
#' @param x object of class jcp
#' @param cex numeric, global sizes in plot
#' @param cex.main numeric, size of titles
#' @param ... additional arguments
#'
#' @return No return value, called for side effects
#'
#' @examples 
#' # Normal distributed sequence with 3 change points at
#' # c1=250 (change in expectation), 
#' # c2=500 (change in variance) and 
#' # c3=750 (change in expectation and variance) 
#' set.seed(0)
#' m      <- c(8,10,10,3);   s  <- c(4,4,10,5)
#' x      <- rnorm(1000, mean=rep(m,each=250), sd=rep(s,each=250))
#' result <- jcp(x)
#' summary(result)
#' plot(result)
#' 
#' # Set additional parameters (window set)
#' result2 <- jcp(x,H=c(80,160,240))
#' summary(result2)
#' plot(result2)
#' 
#' @seealso \code{\link{jcp}, \link{summary.jcp}}
#' @author Michael Messer
#' 
#' @references Michael Messer (2021) Bivariate change point detection - joint detection of changes in expectation and variance, Scandinavian Journal of Statistics, DOI 10.1111/sjos.12547. 
#' 
#' 
#' @rdname plot.jcp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


plot.jcp <- function(x,cex=1.0,cex.main=1.0,...)
{

####
#### Layout and parameter setting
####

oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))  
      
mat <- cbind(c(1,2,2,2),c(1,3,3,3),c(1,3,3,3))
layout(mat); par(cex=cex,cex.main=cex.main,mar=c(1.5,2,2,0.5));  

object     <-  x
x          <-  object$x
H          <-  object$H
q          <-  object$q
hat_cp     <-  object$CP_meta[,1]
hat_h      <-  object$CP_meta[,2]
hat_Ehc    <-  object$CP_meta[,3]
hat_Vhc    <-  object$CP_meta[,4]
hat_rho_cp <-  object$CP_meta[,5]
hat_mean   <-  object$mean_sd[,2]
hat_sd     <-  object$mean_sd[,3]
region     <-  object$region

col    <- rev(rainbow(length(H),start=0,end=0.35,s=1,v=0.75))
col_cp <- col[apply(as.matrix(hat_h),MARGIN=1,FUN=function(h,H=H){which(h==H)},H=H)]

####
#### Plot part 1 (of 3) - Data and estimated change points
####

ylm <- range(x)
plot(1:(length(x)+200),c(x,rep(-500,200)),bty="n",pch=19,cex=0.1,axes=FALSE,xlab="",ylab="",main=paste("Data and change point estimation (#CPs=",length(hat_cp),")",sep=""),ylim=ylm,xlim=c(0,length(x)),col=c(rep(1,length(x)),rep(0,200))); 
axis(1,at=c(1,hat_cp,length(x))); axis(2)
abline(v=hat_cp,col=col_cp,lwd=4)

####
#### Plot part 2 (of 3) - Estimated expectation and variance
####


par(mar=c(4,4,2,1)); # Change margins  
rngx <- range(hat_mean); rngy <- range(hat_sd)
rngxplot <- c(rngx[1]-0.15*diff(rngx),rngx[2]+0.15*diff(rngx))
rngyplot <- c(rngy[1]-0.15*diff(rngy),rngy[2]+0.15*diff(rngy))

plot(hat_mean,hat_sd,type="n",xlim=rngxplot,ylim=rngyplot,xlab="mean",ylab="standard deviation",axes=FALSE,main="",lwd=0.1,asp=1);
axis(1,col = 6,col.axis=6); axis(2,col = 4,col.axis=4)
par(xpd=TRUE); mtext("Parameter estimation",adj=0.5,font=2, cex=cex.main); par(xpd=FALSE)

points(hat_mean[1],hat_sd[1],col="darkgray",pch=3,lwd=3)
# plot rejection region
if(length(hat_cp)>0){  
  for(i in 1:length(hat_cp))
  {   
    # plot distribution 
    arrows(hat_mean[-length(hat_mean)],hat_sd[-length(hat_sd)],hat_mean[-1],hat_sd[-1],col="darkgray",length=0.1,lwd=2)
    text(hat_mean[-length(hat_mean)]+0.5*(diff(hat_mean)),hat_sd[-length(hat_sd)]+0.5*(diff(hat_sd)),hat_cp,cex=cex-0.4)
  }#end-for(i in 1:length(hat_cp))
}#end-if(length(hat_cp)>0)  


####
#### Plot part 3 (of 3) - Change point interpretation
####

xx <- seq(-q,q,length.out = 1000); yy <- sqrt(q^2-xx^2) # Rejection circle
rngx      <- c(min(xx,min(hat_Ehc,0)),max(xx,max(hat_Ehc,0)))
rngy      <- c(min(-yy,min(hat_Vhc,0)),max(yy,min(hat_Vhc,0)))
rngxplot  <- c(rngx[1]-0.25*diff(rngx),rngx[2]+0.65*diff(rngx))
rngyplot  <- c(rngy[1]-0.25*diff(rngy),rngy[2]+0.25*diff(rngy))

plot(0,0,type="n",xlim=rngxplot,ylim=rngyplot,xlab="E",ylab="V",axes=FALSE,asp=1,main="",lwd=0.1);
axis(1,col = 6,col.axis=6); axis(2,col = 4,col.axis=4)
lines(rep(0,2),rngyplot,lty=2,col=4)
lines(c(rngxplot[1],rngxplot[2]-5),rep(0,2),lty=2,col=6)
par(xpd=TRUE); mtext("Change point interpretation", cex=cex.main,font=2); par(xpd=FALSE)

# plot rejection region
xx <- seq(-q,q,length.out = 1000); yy <- sqrt(q^2-xx^2); 
if(region=="circle"){lines(xx,yy,col=2,lwd=2,lty=1); lines(xx,-yy,col=2,lwd=2,lty=1)}
if(region=="square"){lines(rep(-q,2),c(-q,q),col=2,lwd=2); lines(rep(q,2),c(-q,q),col=2,lwd=2); lines(c(-q,q),rep(-q,2),col=2,lwd=2); lines(c(-q,q),rep(q,2),col=2,lwd=2)}
if(region=="ellipse"){lines(xx,yy,col=2,lwd=2,lty=1); lines(xx,-yy,col=2,lwd=2,lty=1); lines(rep(-q,2),c(-q,q),col=2,lwd=2,lty=1); lines(rep(q,2),c(-q,q),col=2,lwd=2,lty=1); lines(c(-q,q),rep(-q,2),col=2,lwd=2,lty=1); lines(c(-q,q),rep(q,2),col=2,lwd=2,lty=1)}
lty <- 1

if(length(hat_cp)>0){  
  for(i in 1:length(hat_cp))
  {   
    # plot distribution 
    arrows(0,0,hat_Ehc[i],hat_Vhc[i],col=col_cp[i],code=3,length=0,lwd=2)
    points(hat_Ehc[i],hat_Vhc[i],pch=3,col=col_cp[i],lwd=3)
    CORmat <- cbind(c(1,hat_rho_cp[i]),c(hat_rho_cp[i],1)); #CORmat 
    root_CORmat <- 0.5 * cbind(c(sqrt(1+hat_rho_cp[i])+sqrt(1-hat_rho_cp[i]),sqrt(1+hat_rho_cp[i])-sqrt(1-hat_rho_cp[i])),c(sqrt(1+hat_rho_cp[i])-sqrt(1-hat_rho_cp[i]),sqrt(1+hat_rho_cp[i])+sqrt(1-hat_rho_cp[i])))
    quan1 <- sqrt(qchisq(0.95,df=2)); xxxx1 <- seq(-quan1,quan1,length.out = 1000); yyyy1 <- sqrt(quan1^2-xxxx1^2) 
    quan2 <- sqrt(qchisq(0.66,df=2)); xxxx2 <- seq(-quan2,quan2,length.out = 1000); yyyy2 <- sqrt(quan2^2-xxxx2^2) 
    ellipse_upper <- root_CORmat %*% rbind(xxxx1,yyyy1) + c(hat_Ehc[i],hat_Vhc[i]); ellipse_lower <- root_CORmat %*% rbind(xxxx1,-yyyy1) + c(hat_Ehc[i],hat_Vhc[i])
    lines(ellipse_upper[1,],ellipse_upper[2,],col=col_cp[i],lwd=0.5,lty=lty); lines(ellipse_lower[1,],ellipse_lower[2,],col=col_cp[i],lwd=0.5,lty=lty)
    ellipse_upper <- root_CORmat %*% rbind(xxxx2,yyyy2) + c(hat_Ehc[i],hat_Vhc[i]); ellipse_lower <- root_CORmat %*% rbind(xxxx2,-yyyy2) + c(hat_Ehc[i],hat_Vhc[i])
    lines(ellipse_upper[1,],ellipse_upper[2,],col=col_cp[i],lwd=0.5,lty=lty); lines(ellipse_lower[1,],ellipse_lower[2,],col=col_cp[i],lwd=0.5,lty=lty)  
  }#end-for(i in 1:length(hat_cp))
  
  hat_rho_cp; hat_Ehc; hat_Vhc
  omega <- ifelse(hat_Vhc>0,acos(hat_Ehc/sqrt(hat_Ehc^2+hat_Vhc^2)),2*pi-acos(hat_Ehc/sqrt(hat_Ehc^2+hat_Vhc^2)))
  psn <- ifelse(omega>=pi/4 & omega<(3/4)*pi,3,ifelse(omega>=(3/4)*pi & omega<(5/4)*pi,2,ifelse(omega>=(5/4)*pi & omega<(7/4)*pi,1,4)))
  if(length(hat_cp>0)){text(hat_Ehc,hat_Vhc,hat_cp,col=col_cp,cex=1,pos=psn)};
}#end-if(length(hat_cp)>0)  

cex.legend=1
legend(x="topright",inset=c(0,0.02),legend=as.character(H),pch=19,col=col,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h, " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.0)

}#end-plot.jcp




