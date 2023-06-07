#################################################################################
###                                                                           ### 
###  Maximum likelihood estimation of the parameter of a 1D potential model   ###
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
###  Preparing a test time series with the following model:
###
###  - Stochastic dynamical model (Ito stochastic differential equation; SDE):
###    dz=-U'(z)*dt+a[6]*dW(t) with the potential U(z)=a[1]*z+a[2]*z^2+a[3]*z^3+exp(a[4])*z^4
###
###  - Observation model: y=z+(Gaussian white noise with stadnard deviation of a[5])
###
###  - a[1], ..., a[6] are the paremeters to be estimated
###  - The coefficient for z^4 is constrained to be positive such as exp(a[4]).
###    This is not necessary if we take a positive interval of the parameter search (see below)  

L<-20                     # spacing of observation times (yr), 20 yr
t<-seq(26000,90000,by=L)  # observation times (yr), 26000-90000 yr (0-64000 yr gives the same result!)
n<-length(t)              # number of observed data (length of observations), n=3201
h<-0.001                  # integration time step (_kyr_)  
y<-numeric(n)             # vector of length n to store the observations
z<-1                      # initial condition for the true state z

a<-c(0.2, -1.7,  0.02, -0.2 ,  0.06,  1.3) # True parameter values

### Time integration of the SDE by Euler-Maruyama method ### 
p<-a                      # input the true patameter values a to p 
set.seed(30)              # set a seed for noise with some integer, here 30
y[1]<-z+rnorm(1,0,p[5])   # initial observation at t=0        
for(k in 2:n){      
  for(l in 1:L){
    z<-z+h*(-4*exp(p[4])*z*z*z-3*p[3]*z*z-2*p[2]*z-p[1])+rnorm(1,0,p[6]*sqrt(h)) # one-step integration (h = 1yr = 0.001 kyr)
  }
  y[k]<-z+rnorm(1,0,p[5]) # observation at every L(=20) years
}
# plot(t,y,type="l") # for check
#################################################################################


##########################################################################################################################################
###
###  The function 'ukf', regardless of its name, returns the LOG-LIKELIHOOD of the above model when the model parameter values are given.
###                                        
ukf <-function(para){
na<-1                  # dimension of state space, here na=1 (Caution! na is written as m in Ref. MC17)
na2<-na*2              # the number of sigma points, here na2=na*2=2 
H<-matrix(1,1,na)      # observation matrix
h<-0.001               # time step same as the above
R<-para[5]*para[5]     # variance of the observation noise, epsilon^2 in Ref. MC17
Q<-h*para[6]*para[6]   # covariance matrix, the equation below Eq.(6) in Ref. MC17
transient<-50          # ignore the first 50 data in the calculation of the log-likelihood 
beta<-4*exp(para[4])   # calculate before the following steps
x.pre<-c(y[1])         # mean of the initial filtered density (this is a subjective choice but not critical)
p.pre<-matrix(R,na,na) # covariance of the initial filtered density (this is a subjective choice but not critical)
loc.lnlike<-numeric(n) # vector for the intermediate values used to calculate the log-likelihood finally

### main loop for the state estimation with the uncented Kalman filter (UKF)
for(k in 2:n){

  # 1) Prediction step to get the predicted mean \hat x_k|k-1 (in Ref. MC17) and the predicted covarinance P_k|k-1 (in Ref. MC17) of the state  
  for(l in 1:L){ 
    A<-t(chol(na*p.pre))    # Choleski decomposition. This computation can make numerical divargence sometimes. Then take another initial guess of the parameters.
    sig.p<-cbind(matrix(x.pre,na,na)+A,matrix(x.pre,na,na)-A) # Form sigma points with Eq.(7) of Ref. MC17
    # One-step propagation of sigma points
    sig.p[1,]<-sig.p[1,]+h*(-beta*sig.p[1,]*sig.p[1,]*sig.p[1,]-3.0*para[3]*sig.p[1,]*sig.p[1,]-2.0*para[2]*sig.p[1,]-para[1])     
    x.pre<-rowMeans(sig.p)  # predicted mean \hat x_l (in Ref. MC17)
    # --------- Calculation of the predicted covariance P_l (in Ref. MC17) -------------#
    p.pre<-matrix(0,na,na)
    for(i in 1:na2){
      p.pre<-p.pre+(sig.p[,i]-x.pre) %*% t(sig.p[,i]-x.pre)
    }
    p.pre<-p.pre/na2+Q      # predicted covariance
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
ensemble<-4                          # the number of initial guesses 
lo<-c(-3, -4, -2, -1, 0.000001, 0)   # prescribed lower band of the parameters, a[5]=S.D. of observation noise ,a[6]=S.D. of dynamical noise, see above
up<-c( 3,  3,  2,  1,      0.2, 3)   # prescribed upper band of the parameters
ps<-c(1, 1, 1, 1, 0.1, 1)            # scales of the parameters 
set<-maximinLHS(ensemble,length(ps)) # An 'ensemble' by 'length(ps)' Latin Hypercube Sample matrix with values uniformly distributed on [0,1], see e.g., https://rdrr.io/cran/lhs/man/maximinLHS.html
set<-matrix(lo,ensemble,length(ps),byrow=T)+set %*% diag(up-lo) # the initial parameter guesses uniformly distributed between lo and up

###  maximizing the likelihood by the R-function 'optim' with the method L-BFGS-B. See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
#    - fnscale=-1 makes maximization enable
#    - parscale=ps is important to avoid numerical divergences, which would happen probably by prescribing very bad parameter values to the function 'ukf'
#    - trace = TRUE allows you to monitor the intermediate results of the optimization
fit<-optim(set[1,], fn=ukf, method=c("L-BFGS-B"), lower=lo, upper=up, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))
##########################################################################################################################################


##########################################################################################################################################
###
###  Uncertainty of the maximum likelihood estimate parameter
###
###  Calculating hessian matrix of the likelihood at the maximum likelihood estimate parameter
hess<-optimHess(fit$par, fn=ukf, control=list(fnscale=-1, parscale=ps, maxit=200000, trace=TRUE))
###  Standard error of the maximum likelihood estimate parameter
se<-sqrt(diag(solve(-hess)))
##########################################################################################################################################


##########################################################################################################################################
###
###  Results: note that the ture paremeter values are within the 95% confidence intervals as follows.
###
a            # true paremeter values (a[1],a[2],a[3],a[4],a[5],a[6])
fit$par      # maximum likelihood estimate parameter
fit$par-2*se # lower bound of 95% confidence interval (when Gaussian assumed)
fit$par+2*se # upper bound of 95% confidence interval (when Gaussian assumed)

###  Comparison between the given time series and the time series simulated with the maximum likelihood estimate parameters (noises are the same!) 
p<-fit$par
w<-numeric(n)
z<-1
set.seed(30) # the same seed for the noise is taken again (in order to check how the estimated paremeter values are well)
w[1]<-z+rnorm(1, 0, p[5]) 
for(k in 2:n){
  for(l in 1:L){
    z<-z+h*(-4*exp(p[4])*z*z*z-3*p[3]*z*z-2*p[2]*z-p[1])+rnorm(1,0,p[6]*sqrt(h))
  }
  w[k]<-z+rnorm(1,0,p[5]) 
}

### plotting two time series
#postscript(file="bfgs_search_1d_test.eps", horizontal=TRUE, encoding="WinAnsi.enc") 
plot(t,y,type="l",main="Given time series (black) and time series simulatied with maximum likelihood estimate (green)")
lines(t,w,col=3)
#dev.off()
##########################################################################################################################################
