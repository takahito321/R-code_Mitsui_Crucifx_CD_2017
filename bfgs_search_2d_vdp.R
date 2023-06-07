##################################################################################
###                                                                            ### 
###  Maximum likelihood estimation of the parameters of the van der Pol model  ###
###                                                                            ###
###  Test R-code: bfgs_search_2d_vdp.R  (finally checked on 7 June 2023)       ###
###                                                                            ###
###  written in R by Takahito Mitsui (takahito321@gmail.com)                   ### 
###                                                                            ###
###  Ref. MC17: Mitsui and Crucifix, Clim Dyn, Vol. 48 pp. 2729-2749 (2017)    ###
##################################################################################

library(lhs) # do install.packages("lhs") in the R-terminal if you have not installed it


#################################################################################
###                                                                            
###  Preparing a test time series with the following model:
###
###  - Stochastic van der Pol oscillator (Ito stochastic differential equation; SDE):
###    dz=(p[1]*z-p[2]*z^3/3+p[3]*v)*dt+p1[5]*dW1 
###    dv=                 (-z+p[4])*dt+p1[6]*dW2 (W1, W2 are standard Wiener processes)
###
###  - Observation model: y=z+(Gaussian white noise with stadnard deviation of p[7])
###
###  - p[1], ..., p[7] are the paremeters to be estimated 

L<-20                     # spacing of observation times (yr), 20 yr
t<-seq(1,20000,by=L)      # observation times (yr), 26000-90000 yr (0-64000 yr gives the same result!)
n<-length(t)              # number of observed data (length of observations), n=3201
h<-0.001                  # integration time step (_kyr_)  
y<-numeric(n)             # vector of length n to store the observations
z<-1                      # initial condition for the true state z
v<-0                      # initial condition for the true state v

p<-c(10, 10, 10, -0.5, 0.4, 0.2, 0.1)  # True parameter values

### Time integration of the SDE by Euler-Maruyama method ### 
set.seed(30)              # set a seed for noise with some integer, here 30
y[1]<-z+rnorm(1,0,p[7])   # initial observation at t=0        
for(k in 2:n){      
    for(l in 1:L){
        zn<-z+(p[1]*z-p[2]*z^3/3+p[3]*v)*h+rnorm(1,0,p[5]*sqrt(h))
        vn<-v+(-z+p[4])*h+rnorm(1,0,p[6]*sqrt(h))    # one-step integration (h = 1yr = 0.001 kyr)
        z<-zn
        v<-vn
  }
  y[k]<-z+rnorm(1,0,p[7]) # observation at every L(=20) years
}
plot(t,y,type="l") # for check
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
Q<-diag(c(h*para[5]*para[5],h*para[6]*para[6]))  # covariance matrix
transient<-50          # ignore the first 50 data in the calculation of the log-likelihood 
x.pre<-c(y[1],y[1])    # mean of the initial filtered density (this is a subjective choice but not critical)
p.pre<-diag(c(R,10))   # covariance of the initial filtered density (this is a subjective choice but not critical)
loc.lnlike<-numeric(n) # vector for the intermediate values used to calculate the log-likelihood finally

for(k in 2:n){
    
  # 1) Prediction step to get the predicted mean \hat x_k|k-1 (in Ref. MC17) and the predicted covarinance P_k|k-1 (in Ref. MC17) of the state  
  for(l in 1:L){ 
    A<-t(chol(na*p.pre)) # Choleski decomposition. This computation can make numerical divargence sometimes. Then take another initial guess of the parameters.
    sig.p<-cbind(matrix(x.pre,na,na)+A,matrix(x.pre,na,na)-A) # Form sigma points with Eq.(7) of Ref. MC17
    sig.pm<-sig.p        # just store the sigma points at k-1
    # One-step propagation of sigma points
    sig.p[1,]<-sig.pm[1,]+h*(para[1]*sig.pm[1,]-para[2]*sig.pm[1,]^3/3+para[3]*sig.pm[2,])
    sig.p[2,]<-sig.pm[2,]+h*(-sig.pm[1,]+para[4])
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
ensemble<-10                                         # the number of initial guesses 
lo<-c( 1,   1,   1,   -1,  0.01, 0.01, 0.01)         # prescribed lower band of the parameters
up<-c(50,  50,  50,    1,     1,    1,    1)         # prescribed upper band of the parameters
ps<-c(10,  10,  10,  0.5,   0.1,  0.1,  0.1)         # scales of the parameters 
set<-maximinLHS(ensemble,length(ps))                 # An 'ensemble' by 'length(ps)' Latin Hypercube Sample matrix with values uniformly distributed on [0,1], see e.g., https://rdrr.io/cran/lhs/man/maximinLHS.html
set<-matrix(lo,ensemble,length(ps),byrow=T)+set %*% diag(up-lo) # the initial parameter guesses uniformly distributed between lo and up

###  maximizing the likelihood by the R-function 'optim' with the method L-BFGS-B. See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
#    - fnscale=-1 makes maximization enable
#    - parscale=ps is important to avoid numerical divergences, which would happen probably by prescribing very bad parameter values to the function 'ukf'
#    - trace = TRUE allows you to monitor the intermediate results of the optimization
fit<-optim(set[1,], fn=ukf, method=c("L-BFGS-B"), lower=lo, upper=up, control=list(fnscale=-1, parscale=ps, maxit=1000, trace=TRUE))  # MLE computed from initial parameter guess set[1,]
#fit<-optim(set[2,], fn=ukf, method=c("L-BFGS-B"), lower=lo, upper=up, control=list(fnscale=-1, parscale=ps, maxit=1000, trace=TRUE))  # MLE computed from initial parameter guess set[2,]

##########################################################################################################################################

##########################################################################################################################################
###
###  Uncertainty of the maximum likelihood estimate parameter
###
###  Calculating hessian matrix of the likelihood at the maximum likelihood estimate parameter
hess<-optimHess(fit$par, fn=ukf, control=list(fnscale=-1, parscale=ps, maxit=1000, trace=TRUE))
###  Standard error of the maximum likelihood estimate parameter
se<-sqrt(diag(solve(-hess)))
##########################################################################################################################################
