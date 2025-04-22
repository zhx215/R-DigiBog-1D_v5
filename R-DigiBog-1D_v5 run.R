### The R code for 1-D DigiBog model (in Fortran) presented in Morris et al. (2015), modified to include seasonality, surface pond, snowmelt, and variable density ###
### Author: Zhengyu Xia ###

## start ##
setwd("C:/Users/zhyxi/Dropbox/digibog/LLK")
rm(list=ls())
set.seed(888)
library(pracma)

## parameter ##
t_extent <- 11573 # yr
lateral_extent <- 40000 # cm, LLK peatland mean radius

oxic_decay_base <- 0.014 # yr-1, used for tuning
anoxic_decay_base <- 0.00018 # yr-1, used for tuning
base_temp <- 8.3 # degree C, Apr-Oct mean temperature at LLK peatland
Q10_oxic <- 2 # (from Liu et al., 2024)
Q10_anoxic <- 2 # (from Liu et al., 2024)

porosity <- 0.3 # (from Morris et al., 2015)
k_param_a <- 15800 # cm yr-1, hydraulic conductivity model parameter (from LLK data)
k_param_b <- 5.45 # hydraulic conductivity model parameter (from LLK data)
drain_effi <- 2 # draining efficiency, circular = 2, ellipse = 1 (Laolike peatland is like circular)
x_factor <- 0.5 # recalcitrance (from Morris et al., 2015)
ddf <- 0.15 # cm degree C-1 day-1, degree day factor to calculate snowmelt (from Ramirez et al., 2023)
pond_size <- 0 # cm

density_mode <- "invariable"
initial_density <- 0.15 # g cm-3, it is the initial density in variable mode or the constant density in constant mode (from Dong et al., 2021)

z1_soil <- 30 # cm, soil evaporation threshold 1 (from Swinnen et al., 2021)
z2_soil <- 50 # cm, soil evaporation threshold 2  (from Swinnen et al., 2021)
z_plant <- 80 # cm, plant transpiration threshold (including capillary uptake, 60+20) (from Ge et al., 2023)
Ellenberg <- 9 # waterlogging sensitivity (from Swinnen et al., 2021)
LAI <- 2 # leaf area index (from Mejdova et al., 2021)
prod_scaling <- 0.26 # reduce production

## climate forcing data ##
source("R-DigiBog-1D_v5 climate input.R")

## function ##
# get aet and productivity annually
aet_prod_func <- function(wtd,pet,z1_soil,z2_soil,z_plant,Ellenberg,LAI,prod_scaling){ # wtd cm, annual potential evapotranspiration mm, cm, cm, cm, unitless, unitless, unitless
  pevap <- pet*exp(-0.4*LAI)
  ptran <- pet-pevap
  at_func <- pchipfun(c(0,z_plant/2,z_plant),c(ptran*Ellenberg/12,ptran,0))
  if (wtd<z1_soil){
    aevap <- pevap
  } else if (wtd>z2_soil){
    aevap <- 0
  } else {
    aevap <- pevap*(z2_soil-wtd)/(z2_soil-z1_soil)
  }
  if (wtd<=z_plant){
    atran <- at_func(wtd)
  } else {
    atran <- 0
  }
  aet <- aevap+atran
  row <- 3000*(1-exp(-0.0009695*(aet-20)))/1000/10*prod_scaling # g cm-2 yr-1
  if (row<=0) {row <- 0} # no production
  return(c(aet,row)) # annual actual evapotranspiration mm, annual production g cm-2 yr-1
}

# water level change #
dHdt <- function(p,et,melt,theta_d,K_sum,H,L,childs_factor){ #precip cm/yr, ET cm/yr, snowmelt cm/yr, porosity, hydraulic conductivity sum cm yr-1, water table height cm, radius cm, draining efficiency 
  return(((p+melt-et)-(childs_factor*K_sum*H)/(L^2))/theta_d)
}

# calculate density based on hydraulic conductivity #
derived_density <- function(layer_k){
  if (density_mode=="variable"){
    initial_k <- k_param_a*exp(k_param_b*1)/3153600000 # m s-1
    current_k <- layer_k/3153600000 #m s-1
    density_change <- (log10(current_k)-log10(initial_k))*-20 # kg m-3
    return(initial_density+density_change/1000) # g cm-3
  } else {
    return(initial_density)
  }
}

## run model ##
source("R-DigiBog-1D_v5 algorithm.R")

## run model ##
source("R-DigiBog-1D_v5 validation.R")

# ## plot results ##
# source("R-DigiBog-1D_v5 plot.R")