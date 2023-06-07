# R-code_Mitsui_Crucifx_CD_2017

## R-codes for the maximum likelihood estimates:

### bfgs_search_2d_vdp.R 

- to calculate the maximum likelihood estimate of a 2D van der Pol oscillator model:
- dz=(p[1]*z-p[2]*z^3/3+p[3]*v)*dt+p[5]*dW1, dv=(-z+p[4])*dt+p[6]*dW2 (see R-code)
- p<-c(10, 10, 10, -0.5, 0.4, 0.2, 0.1) : true parameter values
- Maximal likelihood estimate fit$par : 9.81218211  9.99239906  9.79013847 -0.48895473  0.42612936  0.16403794  0.09005438 (sample)
- initial parameter guess : set[1,] set[2,] ... set[ensemble, ]
- maxit=1000 : maximal iteration steps in likelihood optimization by "L-BFGS-B"

### bfgs_search_1d_test.R 

- to calculate the maximum likelihood estimate of a 1D potential model

### bfgs_search_2d_ca_fullforcing.R 

- to calculate the maximum likelihood estimate of the 2D oscillator model with full forcing (ice and insolation) for NGRIP -log10(Ca2+)

### bfgs_search_2d_d18O_fullforcing.R 

- to calculate the maximum likelihood estimate of the oscillator model with full forcing (ice and insolation) for NGRIP d18O

- Note that the fit of the oscillator model to -log10(Ca2+) is significantly better than that of NGRIP d18O. Thus, most of the conclusions of our paper are based on the results with -log10(Ca2+).


---
## Data time series:

- seierstad_20yr.dat : NGRIP Ca2+ and d18O (see MC17 and code)

- laskar_300_0_100yr.dat : Laskar summer insolation at 65N (see MC17 and code)

- LR04-.dat : LR04 ice volume proxy (see MC17 and code)


---
### Ref. MC17: Mitsui and Crucifix, Clim Dyn, Vol. 48 pp. 2729-2749 (2017)  

