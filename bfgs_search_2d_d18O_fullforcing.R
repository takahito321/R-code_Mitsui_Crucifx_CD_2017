#################################################################################
###                                                                           ### 
###  Maximum likelihood estimation of the parameter of the oscillator model   ###
###                                                                           ###
###  << fit to d18O >>                                                        ### 
###                                                                           ###
###  Test R-code: bfgs_search_1d_test.R  (finally checked on 13 Dec 2019)     ###
###                                                                           ###
###  written in R by Takahito Mitsui (takahito321@gmail.com)                  ### 
###                                                                           ###
###  Ref. MC17: Mitsui and Crucifix, Clim Dyn, Vol. 48 pp. 2729-2749 (2017)   ###
#################################################################################

library(lhs) # do install.packages("lhs") in the R-terminal if you have not installed it

#################################################################################
###                                                                            
###  Proxy data and forcings to be used
###
###  !!! This part just prepares time series data (and might be just tedious). You can prepare the time series outside this R-script and read those data !!!  
###
###
ts<-26000              # newest data year BP
te<-90000              # oldest data year BP
L<-20                  # spacing yr of data (<-coming from Rasmmusen et al. QSR 2014/Seierstad et al. QSR 2014)
t<-seq(ts,te,by=L)     # observation times (yr)
n<-length(t)           # number of observation data points

###  NGRIP d18O from seierstad_20yr.dat 
tab2<-read.table("seierstad_20yr.dat") # data from Seierstad et al. QSR 2014 (some missing points were interpolated for Ca2+. See MC17) 
tab<-tab2[seq(1,nrow(tab2),by=2),] # 1,3,5,...,nrow(tab2)
for(i in 1:nrow(tab)){  # for-loop to find the indices for times ts and te
  if(tab[i,1]==ts) is<-i
  if(tab[i,1]==te) ie<-i
}
y<-tab[is:ie,5]
y<-y-(-41.15583) # subtract mean. see MC17
y<-y/2.235739    # divided by standard deviation. see MC17
y<-rev(y)        # change the direction of time series vector from the oldest point to the earliest point

###  Insolation data
ins<-numeric(n)  # vector for insolation values
tab2<-read.table("laskar_300_0_100yr.dat") # June-July monthely mean, 100-yr sampling, from Web: http://vo.imcce.fr/insola/earth/online/earth/online/index.php
tab2[,1]<--tab2[,1]*1000  # kyr to yr
for(i in 1:nrow(tab2)){   # for-loop to find the indices for times ts and te
  if(tab2[i,1]==te) ie<-i
  if(tab2[i,1]==ts) is<-i
}
tab<-tab2[ie:is,]
for(j in 1:(nrow(tab)-1)){
 for(i in 1:5){  # for-loop to interpolate LR04 with 20-yr spacing
     k=(j-1)*5+i
     ins[k]<-tab[j,2]+(tab[j+1,2]-tab[j,2])*(i-1)/5
   }
}
ins[length(ins)]<-tab[length(tab[,1]),2] # fill the edge of the vector
ins<-ins-474.9268  # subtract mean. see MC17
ins<-ins/15.08376  # divided by standard deviation. see MC17

### ice volume proxy data
LR04<-numeric(n)  # vector for ice volume values
tab2<-read.table("LR04-.dat") # LR04 data with time from past (negative) to present)
tab2[,1]<--tab2[,1]*1000  # kyr to yr and past -> positive 
for(i in 1:nrow(tab2)){   
  if(tab2[i,1]==te) ie<-i
  if(tab2[i,1]==ts) is<-i
}
tab<-tab2[ie:is,]
for(j in 1:(nrow(tab)-1)){  # for-loop to interpolate Laskar's insolation with 20-yr spacing
 for(i in 1:50){
     k=(j-1)*50+i
     LR04[k]<-tab[j,2]+(tab[j+1,2]-tab[j,2])*(i-1)/50
   }
}
LR04[length(LR04)]<-tab[length(tab[,1]),2] # fill the edge of the vector
LR04<-LR04-4.342722   # subtract mean. see MC17
LR04<-LR04/0.329534   # divided by standard deviation. see MC17
#################################################################################


##########################################################################################################################################
###
###  The function 'ukf', regardless of its name, returns the LOG-LIKELIHOOD of the above model when the model parameter values are given.
###    
ukf <-function(para){
na<-2                  # dimension of state space, here na=1 (Caution! na is written as m in Ref. MC17)
na2<-na*2              # the number of sigma points, here na2=na*2=4
H<-matrix(c(1,0),1,na) # observation matrix
h<-0.001               # time step size for propagating sigma points in UKF 
R<-para[7]*para[7]     # variance of the observation noise, epsilon^2 in Ref. MC17
Q<-diag(c(h*para[5]*para[5],h*para[6]*para[6]))  # covariance matrix, the equation below Eq.(6) in Ref. MC17
transient<-50          # ignore the first 50 data in the calculation of the log-likelihood 
beta<-para[4]/3        # calculate in advance (not important; just reflecting my coding history)
x.pre<-c(y[1],y[1])    # mean of the initial filtered density (this is a subjective choice but not critical)
p.pre<-diag(c(R,10))   # covariance of the initial filtered density (this is a subjective choice but not critical)
loc.lnlike<-numeric(n) # vector for the intermediate values used to calculate the log-likelihood finally

for(k in 2:n){
    
  # 1) Prediction step to get the predicted mean \hat x_k|k-1 (in Ref. MC17) and the predicted covarinance P_k|k-1 (in Ref. MC17) of the state  
  for(l in 1:L){ 
    A<-t(chol(na*p.pre)) # Choleski decomposition. This computation can make numerical divargence sometimes. Then take another initial guess of the parameters.
    sig.p<-cbind(matrix(x.pre,na,na)+A,matrix(x.pre,na,na)-A) # Form sigma points with Eq.(7) of Ref. MC17
    sig.pm<-sig.p        # just store the sigma points at k-1
    tmp1<-sig.pm[1,]-para[2]
    # One-step propagation of sigma points
    sig.p[1,]<-sig.pm[1,]+h*(para[1]*sig.pm[2,]-para[3]*sig.pm[1,]-beta*(tmp1*tmp1*tmp1+para[2]*para[2]*para[2]))
    sig.p[2,]<-sig.pm[2,]+h*(-sig.pm[1,]+para[8]+para[9]*ins[k-1]-para[10]*LR04[k-1])
    x.pre<-rowMeans(sig.p) # predicted mean \hat x_l (in Ref. MC17)
    # --------- Calculation of the predicted covariance P_l (in Ref. MC17) -------------#
    p.pre<-matrix(0,na,na)
    for(i in 1:na2){
      p.pre<-p.pre+(sig.p[,i]-x.pre) %*% t(sig.p[,i]-x.pre)
    } 
    p.pre<-p.pre/na2+Q  # predicted covariance
    # --------- Calculation of the predicted covariance P_l (in Ref. MC17) -------------#
  } # after this roop, you have \hat x_k|k-1 = \hat x_L = x.pre and P_k|k-1 = P_L = p.pre (in Ref. MC17)

  # 2) Filtering step to get the filtered mean \hat x_k|k (in Ref. MC17) and the filtered covariance P_k|k (in Ref. MC17)
  S<-H %*% p.pre %*% t(H)+R           # S_k in Eq.(8) of Ref. MC17
  K<-p.pre %*% t(H) %*% solve(S)      # K_k in Eq.(8) of Ref. MC17
  resi <-y[k]-H %*% x.pre             # \zeta_k in Eq.(8) of Ref. MC17 
  loc.lnlike[k]<--(log(det(S))+as.numeric(t(resi) %*% solve(S) %*% resi))/2 # ''local log-likelihood'' 
  x.pre<-x.pre + K %*% resi           # Filtered mean,  \hat x_k|k in Eq.(8) of Ref. MC17
  p.pre<-(diag(na)-K %*% H) %*% p.pre # Filtered covariance, P_k|k in Eq.(8) of Ref. MC17
}

return( -(n-transient)/2*log(2*pi)+sum(loc.lnlike[(transient+1):n]) ) # return the log-likelihood calculated by ignoring the first 50 local log-likelihoods  
}
##########################################################################################################################################


##########################################################################################################################################
###
###  Calculating the maximum likelihood and the corresponding parameter by L-BFGS-B method, see Ref. MC17
###
###  Latin hypercube sampling of the initial guesses of the parameters
ensemble<-10                                             # the number of initial guesses 
lo<-c(5, -0.5, -50, 5, 0, 0, 0.000001, -2, -2, -2)       # prescribed lower band of the parameters
up<-c(200,  0.5, 100, 200,  10, 3,     0.5,  2,  2,  2)  # prescribed upper band of the parameters
ps<-c(10, 0.1, 5, 10, 1, 1, 0.1, 0.5, 0.5, 0.5)          # scales of the parameters 
set<-maximinLHS(ensemble,length(ps))                     # An 'ensemble' by 'length(ps)' Latin Hypercube Sample matrix with values uniformly distributed on [0,1], see e.g., https://rdrr.io/cran/lhs/man/maximinLHS.html
set<-matrix(lo,ensemble,length(ps),byrow=T)+set %*% diag(up-lo) # the initial parameter guesses uniformly distributed between lo and up

###  maximizing the likelihood by the R-function 'optim' with the method L-BFGS-B. See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
#    - fnscale=-1 makes maximization enable
#    - parscale=ps is important to avoid numerical divergences, which would happen probably by prescribing very bad parameter values to the function 'ukf'
#    - trace = TRUE allows you to monitor the intermediate results of the optimization
fit<-optim(set[2,], fn=ukf, method=c("L-BFGS-B"), lower=lo, upper=up, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))
##########################################################################################################################################

##########################################################################################################################################
###
### To get the parameter values in Table 6: change FALSE to TRUE if you try
###
if(FALSE){ 
ensemble<-10
lo<-c(100, -0.65, -50, 50, 0, 0, 0.000001, -2, -2, -2)
up<-c(600,  0.5, 250, 200,  10, 3,     0.5,  2,  2,  2)
ps<-c(10, 0.1, 5, 10, 1, 1, 0.1, 0.5, 0.5, 0.5)
set<-maximinLHS(ensemble,length(ps))
set<-matrix(lo,ensemble,length(ps),byrow=T)+set %*% diag(up-lo)

set[1,]<-c(300.0000000,  -0.5752882, 114.1725173,  56.2275203,   5.1109058,   0.7862598, 0.0000010,  -0.1101134,   0.3338190,   0.5990490)
fit<-optim(set[1,], fn=ukf, method=c("L-BFGS-B"), lower=lo, upper=up, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))

hess<-optimHess(fit$par, fn=ukf, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))
#sqrt(diag(solve(-hess)))
}

##########################################################################################################################################
###
###  Uncertainty of the maximum likelihood estimate parameter
###
###  Calculating hessian matrix of the likelihood at the maximum likelihood estimate parameter
hess<-optimHess(fit$par, fn=ukf, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))
###  Standard error of the maximum likelihood estimate parameter
se<-sqrt(diag(solve(-hess)))
##########################################################################################################################################
