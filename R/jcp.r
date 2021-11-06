#' jcp
#'
#' Joint change point detection - expectation and variance - via bivariate moving sum statistics 
#'
#' @param x numeric vector. Input sequence of random variables.
#' @param H NA or numeric vector. Window set. If NA (default), then H is automatically set. If not NA, then H must an increasing vector of positive integers with maximum =< length(x)/2. 
#' @param q NA or numeric value. Rejection threshold. If NA (default), then the rejection boundary is derived in simulations (from Gaussian process limit) according to sim and alpha. If not NA, then q is considered predefined and must be set a postive real number.  
#' @param alpha numeric value. Significance level. Must be in (0,1), default = 0.05. In case of predefined q, alpha is set to NA.
#' @param sim numeric value. Number of simulations of limit process for approximation of q. Must be positive integer, default = 1000. In case of predefined q, sim is set to NA.
#' @param region character string. Defines rejection region, default = "square". Must be chosen either "square", "circle" or "ellipse". 

#' @return invisible list
#' \item{changepoints}{detected change points (increasingly ordered)}
#' \item{mean_sd}{matrix of estimated means and standard deviations}
#' \item{M}{test statistic}
#' \item{q}{rejection threshold}
#' \item{H}{window set}
#' \item{sim}{number of simulations of the limit process (approximation of q)}
#' \item{alpha}{significance level}
#' \item{region}{rejection region}
#' \item{method}{derivation of threshold q, either asymptotic or predefined}
#' \item{x}{input sequence}
#' \item{EVrho}{list containing the auxiliary processes E, V and correlation rho, for each element of H one list entry}
#' \item{CP_meta}{matrix containing meta information of estimation. Estimated change points (increasingly ordered), responsible window h, components E, V and rho of joint statistic at estimated change points (regarding responsible window)}
#' \item{SFA}{detected change points of single filter algorithms}
#' 
#'
#' @seealso \code{\link{plot.jcp}, \link{summary.jcp}}
#' @author Michael Messer
#' 
#' @references Michael Messer (2021) Bivariate change point detection - joint detection of changes in expectation and variance, Scandinavian Journal of Statistics, DOI 10.1111/sjos.12547.
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
#' 
#' @rdname jcp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


jcp <- function(x,H=NA,q=NA,alpha=0.05,sim=1000,region="square")
{

method <- ifelse(is.na(q),"asymptotic","predefined")

if(method == "predefined"){alpha <- NA; sim <- NA}  

if( (method == "predefined") & (!is.numeric(q) | q[1]<=0 | length(q) != 1)){stop("Invalid choice of q: must be either NA for asymptotic derivation, or set as predefined positive value")}

if( (method == "asymptotic") & (alpha*(1-alpha) <= 0)){stop("Invalid choice of significance level: alpha must be in (0,1)")}  

if( (method == "asymptotic") & (sim <= 0 | is.na(sim) | length(sim) != 1)){stop("Invalid choice of number of simulations: sim must be positve integer (~1000)")}

if( (region %in% c("square","circle","ellipse")) == FALSE){stop("Invalid choice of region: must be 'square', 'circle' or 'ellipse'")}

Tt <- length(x)    
  
if(is.na(H[1])){											                        
  if (Tt<=100) {stop("Automatic choice of windows failed. Time series is too short (need > 100 data points). Set H manually.")}
  H <- seq(50,min((Tt-1)/2,200),25)
}

if(!is.na(H[1])){
  if(is.null(H)){stop("If autoset.H is FALSE, the vector of window sizes H must be set")}
  if(!all(c(diff(c(2,H,(Tt/2))))>0) | !all(H%%1 == 0) | !is.numeric(H)){stop("Invalid choice of window set: H needs to be an increasing ordered vector of integers, with min(H) > 2 and max(H) <= T/2")}
}

#######
####### Simulate q
#######

if(method=="asymptotic"){
if( !sim%%1==0 | sim <= 0) {stop("Invalid choice of number of simulations: sim must be a positive integer")}
if(sim < 1000){cat("Warning: Number of simulations might be too low",sep="\n")}

  # Given a planar Brownian motion and window h, simulate maximal Euklidian distance of limit process
  max_h <- function(h,WE,WV)
  {  
    Tt <- length(WE)
    LE <- (WE[(2*h+1):(Tt)] - 2*WE[(1+h):(Tt-h)] + WE[1:(Tt-2*h)])/sqrt(2*h) 
    LV <- (WV[(2*h+1):(Tt)] - 2*WV[(1+h):(Tt-h)] + WV[1:(Tt-2*h)])/sqrt(2*h)
    max(sqrt(LE^2 + LV^2))
  }
  
  # Simulate planar Brownian motion and calculate global maximal Euklidian distance over all limit processes (for different h)    
  sim_global_max <- function(Tt,H)
  {
    WE      <- cumsum(rnorm(Tt))
    WV      <- cumsum(rnorm(Tt))
    maxL_H  <- vapply(as.matrix(H),FUN=max_h,WE=WE,WV=WV,numeric(length(1)))
    max(maxL_H)
  }
  
  q <- quantile(replicate(n=sim,sim_global_max(Tt=Tt,H=H)),1-alpha)
}#end-if-(method=="asymptotic")

#######
####### E, V and rho
####### 

# Given t and h, calculate E_ht, V_ht and rho_ht (three terms of the joint statistic).
EVrho__marginal__depending_on_h_and_t <- function(t,x,h) 
  {
    xr <- x[(t+1):(t+h)]; xl <- x[(t-h+1):t]
    EV_numerator <- c(mean(xr)-mean(xl) , var(xr)-var(xl))
    EV  <- sqrt(h) * EV_numerator / sqrt(c(var(xr) + var(xl), var((xr-mean(xr))^2)+var((xl-mean(xl))^2) )) 
    rho <- (mean((xr-mean(xr))^3) + mean((xl-mean(xl))^3)) / (sqrt(var(xr) + var(xl)) * sqrt(var((xr-mean(xr))^2)+var((xl-mean(xl))^2))  ) 
    E_V_rho <- c(EV,rho)
    E_V_rho
  }#End-EVrho__marginal__depending_on_h_and_t  

# Calculate E_ht, V_ht and rho_ht for fixed h and all t (matrix with three rows, column is time)
EVrho__matrix__fixed_h_depending_on_t <- function(x,h) 
  {
  tau <- h:(length(x)-h)
  apply(as.matrix(tau),MARGIN=1,FUN=EVrho__marginal__depending_on_h_and_t,x=x,h=h) 
  }#End-EVrho__matrix__fixed_h_depending_on_t

# Calculate E_ht, V_ht and rho_ht for all h and all t (list, i-th element is a matrix that corresponds to the i-th window)
EVrho_list <- lapply(as.matrix(H),FUN=EVrho__matrix__fixed_h_depending_on_t,x=x)

#######
####### Single Bivariate Moving Sum Algorithm
####### 

# Apply SFA for fixed h (depends on the rejection region)
SFA__vector <- function(EVrho_matrix_fixed_h,Tt=Tt,q=q,region=region)
{
  E  <- EVrho_matrix_fixed_h[1,]
  V  <- EVrho_matrix_fixed_h[2,]
  rho<- EVrho_matrix_fixed_h[3,]
  h  <- (Tt - length(E)+1)/2
  hat_cp <- c()           # Storage for estimated change points
  Enew   <- E; Vnew <- V  # Update statistics (cutout h-neighborhood, i.e., set to zero)
  total_deleted_timeg <- c(); total_deleted_time <- c()
  norm_euklid <- sqrt(Enew^2 + Vnew^2) # Change point estimation is always based on Euklidian distance...
  #... but decision to stop is based on type of the region
  # Euklidian distance <-> circle
  if(region=="circle"){nJ <- sqrt(Enew^2 + Vnew^2)} # Euklidian distances of E and V
  # Mahalanobis distance <-> time dependent ellipse
  if(region=="ellipse"){
    EVrho_new <- rbind(Enew,Vnew,rho)
    mahalanobis <- function(EVrho_new_column) # Mahalanobis distances of E and V (time dependent correlation)
      {
      rho_t <- EVrho_new_column[3]
      Gamma_t <- matrix(c(1,rho_t,rho_t,1),nrow=2)
      evp_t <- eigen(Gamma_t) 
      Gamma_t_inv <- evp_t$vectors %*% diag(1 / evp_t$values) %*% t(evp_t$vectors)
      sqrt(as.vector(EVrho_new_column[1:2] %*%  Gamma_t_inv %*% EVrho_new_column[1:2]))
      }
    nJ <- apply(EVrho_new,MARGIN=2,FUN=mahalanobis)
    }
  # Maximum distance <-> square
  if(region=="square"){nJ <- apply(rbind(Enew,Vnew),MARGIN=2,FUN=function(clm){max(abs(clm))})} # Maximal distances of E and V
  # Successive change point detection
  Mh <- max(nJ)
  if (is.nan(Mh)) {stop("Invalid input data. Statistics not well-defined.")}
  while(max(nJ) > q) 
  {
    rstrct <- which(nJ>q) # restrict to points where joint process lies in rejection region...
    hat_cg <- rstrct[which.max(abs(norm_euklid)[rstrct])]; hat_c <- hat_cg+(h-1); # ...from these choose maximizer wrt Euklid
    deleted_timeg <- max(1,(hat_cg-(h-1))):min(length(Enew),(hat_cg+(h-1)));  # time to delete regarding new CP
    # update upper quantities
    total_deleted_timeg <- c(total_deleted_timeg,deleted_timeg); # total time to be deleted
    Enew[deleted_timeg] <- 0; Vnew[deleted_timeg] <- 0; # update G1new, G2new
    norm_euklid <- sqrt(Enew^2 + Vnew^2)
    # Euklid
    if(region=="circle"){nJ <- sqrt(Enew^2 + Vnew^2)}
    # Mahalanobis
    if(region=="ellipse"){nJ <- sqrt(Enew^2 + Vnew^2)}
    # Maximum
    if(region=="square"){nJ <- apply(rbind(Enew,Vnew),MARGIN=2,FUN=function(clm){max(abs(clm))})}
    hat_cp <- c(hat_cp,hat_c)
  } # end-while
  hat_cp
  list(Mh=Mh,hat_cp=hat_cp)
} # end-SFA

# Apply SFA for all h in H
SFA_and_Mh_list <- lapply(EVrho_list,FUN=SFA__vector,Tt=Tt,q=q,region=region) 
M               <- max(unlist(lapply(SFA_and_Mh_list, FUN=function(x){x[[1]]})))
SFA_list        <- lapply(SFA_and_Mh_list, FUN=function(x){x[[2]]})
         
#######
####### Multi Bivariate Moving Sum Algorithm
#######

c_tmp		    <- SFA_list[[1]]				    # First, accept CPs detected with smallest h.
h_val_tmp 	<- rep(H[1],length(c_tmp)) 	# Corresponding h value.

if(length(H)>1){               
  for (k in 2:(length(H))){   # For window with index 2 and up ...
    if(!is.null(SFA_list[[k]][1])){ # Check if there are CPs detected with window hk
      ck <- SFA_list[[k]]
      for(candidate in 1:length(ck)){
      cps_in_neighborhood <- c_tmp[ c_tmp >= ck[candidate] - H[k] + 1  & c_tmp <= ck[candidate] + H[k] ]
      # If no CP in neighborhood of candidate, then add candidate
      if(length(cps_in_neighborhood)==0){c_tmp <- c(c_tmp,ck[candidate]); h_val_tmp <- c(h_val_tmp,H[k])}
      }# End-for-Candidate
    } # End-if
  } # End-for-k
} # End-if-length(H)

CP <- cbind(c_tmp,h_val_tmp) # Store CPs and associated windows
if(dim(CP)[1]>1){CP <- CP[order(CP[,1]),]}  # Sort CPs by their detection h value.
colnames(CP) <- c("changepoint","h")  
if(dim(CP)[1] > 0) {rownames(CP) <- as.character(1:dim(CP)[1])}

#######
####### Calculate estimated means and standard-deviations in detected sections
#######

# Means and standard deviations
 CPs <- c(1,as.vector(CP[,1]),length(x)) 
 section          <- 1:(length(CPs)-1)
 hat_mean <- rep(NA,(length(CPs)-1))        
 hat_sd   <- rep(NA,(length(CPs)-1))
 for(i in 1:(length(CPs)-1)){      
   hat_mean[i] <- mean(x[CPs[i]:CPs[i+1]])  
   hat_sd[i]   <-  sd(x[CPs[i]:CPs[i+1]])   
 }# end-for-i
 mean_sd <- cbind(section,hat_mean,hat_sd)
 
 # Estimated normal distribution at estimated change points
 hat_cpg <- CP[,1] - (CP[,2]-1)
 hat_rho_cp <- hat_Ehc <- hat_Vhc <- numeric(0)
 if(length(CP[,1])>0){
   hat_Ehc    <-  rep(NA,length(CP[,1]))
   hat_Vhc    <-  rep(NA,length(CP[,1]))
   hat_rho_cp <-  rep(NA,length(CP[,1]))
   for(i in 1:length(CP[,1]))
   {
     # calculate distribution
     cpi <- CP[i,1]; hi <- CP[i,2]
     hat_Ehc[i] <- EVrho_list[[which(H==hi)]][1,hat_cpg[i]]  
     hat_Vhc[i] <- EVrho_list[[which(H==hi)]][2,hat_cpg[i]]  
     xl <- x[(cpi-hi+1):(cpi)]; xr <- x[(cpi+1):(cpi+hi)]
     sigsqr    <- mean((xr - mean(xr))^2)
     mcenter3r <- mean((xr - mean(xr))^3)
     nusqr     <- mean((xr - mean(xr))^4) - sigsqr^2
     sigsql    <- mean((xl - mean(xl))^2)
     mcenter3l <- mean((xl - mean(xl))^3)
     nusql     <- mean((xl - mean(xl))^4) - sigsql^2
     hat_rho_cp[i] <- (mcenter3r + mcenter3l) / ( sqrt(sigsqr+sigsql) * sqrt(nusqr+nusql) )
     CORmat <- cbind(c(1,hat_rho_cp[i]),c(hat_rho_cp[i],1)); #CORmat 
     root_CORmat <- 0.5 * cbind(c(sqrt(1+hat_rho_cp[i])+sqrt(1-hat_rho_cp[i]),sqrt(1+hat_rho_cp[i])-sqrt(1-hat_rho_cp[i])),c(sqrt(1+hat_rho_cp[i])-sqrt(1-hat_rho_cp[i]),sqrt(1+hat_rho_cp[i])+sqrt(1-hat_rho_cp[i])))
   }#end-for(i in 1:length(CP[,1]))
 }#end-if(length(CP[,1])>0)  

CP_meta <- cbind(CP,hat_Ehc,hat_Vhc,hat_rho_cp) 
 
#######
####### Generate output 
#######

 lst <- list(changepoints=CP[,1],mean_sd=mean_sd,M=M,q=q,H=H,sim=sim,alpha=alpha,region=region,method=method,x=x,CP_meta=CP_meta,EVrho=EVrho_list,SFA=SFA_list) 
 class(lst) <- "jcp"
 invisible(lst)
}#end-jcp

