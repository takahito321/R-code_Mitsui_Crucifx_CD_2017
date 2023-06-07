# R-code_Mitsui_Crucifx_CD_2017

R-code for calculating the maximum likelihood estimates of the models:

- bfgs_search_1d_test.R :

to calculate the maximum likelihood estimate of the 1D potential model for a sample time series

- bfgs_search_2d_ca_fullforcing.R :

to calculate the maximum likelihood estimate of the 2D oscillator model with full forcing (ice and insolation) for NGRIP -log10(Ca2+)

- bfgs_search_2d_d18O_fullforcing.R :

to calculate the maximum likelihood estimate of the oscillator model with full forcing (ice and insolation) for NGRIP d18O

# Note that the fit of the oscillator model to -log10(Ca2+) is significantly better than that of NGRIP d18O. Thus, most of the conclusions of our paper are based on the results with -log10(Ca2+).


---
Data time series:

- seierstad_20yr.dat : NGRIP Ca2+ and d18O (see MC17 and code)

- laskar_300_0_100yr.dat : Laskar summer insolation at 65N (see MC17 and code)

- LR04-.dat : LR04 ice volume proxy (see MC17 and code)


---
Ref. MC17: Mitsui and Crucifix, Clim Dyn, Vol. 48 pp. 2729-2749 (2017)  

