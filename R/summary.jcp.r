#' summary.jcp
#'
#' Summary method for class 'jcp'
#'
#' @param object object of class jcp
#' @param ... additional arguments
#'
#' @return No return value, called for side effects
#'
#' @examples 
#' #' # Normal distributed sequence with 3 change points at
#' # c1=250 (change in expectation), 
#' # c2=500 (change in variance) and 
#' # c3=750 (change in expectation and variance) 
#' set.seed(0)
#' m      <- c(8,10,10,3);   s  <- c(4,4,10,5)
#' x      <- rnorm(1000, mean=rep(m,each=250), sd=rep(s,each=250))
#' result <- jcp(x)
#' plot(result)
#' summary(result)
#' 
#' # Set additional parameters (window set)
#' result2 <- jcp(x,H=c(80,160,240))
#' plot(result2)
#' summary(result2)
#' 
#' @seealso \code{\link{jcp}, \link{plot.jcp}}
#' @author Michael Messer
#' 
#' @references Michael Messer (2021) Bivariate change point detection - joint detection of changes in expectation and variance, Scandinavian Journal of Statistics, DOI 10.1111/sjos.12547.
#' 
#' 
#' @rdname summary.jcp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end

summary.jcp <- function(object,...)
{
  cat("","\n")
  cat("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM")   
  cat("","\n")
  cat("Joint change point detection - expectation and variance")   
  cat("","\n")
  cat("","\n")
  cat("Technical parameters: ");
  cat("","\n"); cat("         ")
  cat(paste("Tt = ",length(object$x)," and H = {",paste(object$H,collapse=", "),"}",sep="")); 
  cat("","\n"); cat("         ")
  cat(paste("Derivation of rejection boundary: ",object$method,sep=""));
  cat("","\n"); cat("         ")
  cat(paste("Significance level alpha = ",object$alpha,", Number of simulations = ",object$sim,sep=""));
  cat("","\n"); cat("","\n");
  
  cat(paste("Hypothesis test:",sep=""));
  cat("","\n"); cat("         ")
  if(object$M>object$q)
  {cat(paste("Stationarity was rejected: M = ",round(object$M,2)," > q = ", round(object$q,2),sep=""),sep="")}	
  else{cat(paste("Stationarity was not rejected: M = ",round(object$M,2)," < q = ", round(object$q,2),sep=""),sep="")}
  cat("","\n"); cat("         ")
  cat(paste("Rejection boundary: ",object$region,sep="")); 
  cat("","\n"); cat("","\n");
  
  cat("Joint change point detection and parameter estimation: "); 
  cat("","\n"); cat("         ")
  if(dim(object$CP_meta)[1]==0){cat("No change points detected")}
  if(dim(object$CP_meta)[1]==1){cat(paste(dim(object$CP_meta)[1],"change point detected: "))} 
  if(dim(object$CP_meta)[1]>1){cat(paste(dim(object$CP_meta)[1],"change points detected: "))}
  if(dim(object$CP_meta)[1]>0){cat(paste(object$CP_meta[,1]),sep=", ")}; 
  cat("","\n"); cat("         ")
  if(length(object$mean_sd[,2])==1){cat(paste(length(object$mean_sd[,2]),"section with"))}
  if(length(object$mean_sd[,2])>1){cat(paste(length(object$mean_sd[,2]),"sections with"))}
  cat("","\n"); cat("              ")
  if(length(object$mean_sd[,2])==1){cat(paste("estimated expectation: ",signif(object$mean_sd[,2],2)) )} else{cat("estimated expectations: ")
    cat(paste(signif(object$mean_sd[,2],2)),sep=", ")}; 
  cat("","\n"); cat("              ")
  if(length(object$mean_sd[,3])==1){cat(paste("estimated standard deviation: ",signif(object$mean_sd[,3],2)) )} else{cat("estimated standard deviations: "); cat(paste(signif(object$mean_sd[,3],2)),sep=", ")}
  cat("","\n");#cat("","\n");
  cat("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM")     
}# end-summary.jcp

