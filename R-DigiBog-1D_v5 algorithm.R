## initiate model parameters ##
# cohorts
layer_mass <- rep(NA,t_extent); layer_mass[1] <- 0.01 # g cm-2, gives an initial peat mass to start the model
layer_initial_mass <- rep(NA,t_extent); layer_initial_mass[1] <- layer_mass[1] # g cm-2
layer_remaining_mass <- rep(NA,t_extent); layer_remaining_mass[1] <- 1
layer_density <- rep(NA,t_extent); layer_density[1] <- initial_density # g cm-3
layer_thickness <- rep(NA,t_extent); layer_thickness[1] <- layer_mass[1]/layer_density[1] # cm
layer_elevation <- rep(NA,t_extent); layer_elevation[1] <- layer_thickness[1] # cm
wet_proportion <- rep(NA,t_extent); wet_proportion[1] <- 1
layer_hydro_k <- rep(NA,t_extent); layer_hydro_k[1] <- k_param_a*(exp(k_param_b*layer_remaining_mass[1])) # cm yr-1
layer_transmissivity <- rep(NA,t_extent); layer_transmissivity[1] <- layer_thickness[1]*layer_hydro_k[1] # cm2 yr-1
# parameters
peat_height <- layer_elevation[1] # cm
wt_height <- layer_elevation[1] # cm, required to calculate wt change
peat_transmissivity <- layer_transmissivity[1] # cm2 yr-1, required to calculate wt change
snow <- 0 # cm
pond <- 0 # cm
wtd <- 0 # cm

## initiate model outputs ##
# final annual time series
wtd_output <- gswtd_output <- prod_output <- decay_output <- height_output <- mass_output <- rep(NA,t_extent)
# intermediate seasonal time series
daily_wtd_output <- daily_prod_output <- daily_decay_output <- daily_melt_output <- daily_pond_output <- vector()
daily_temp_output <- daily_precip_output <- daily_pet_output <- daily_aet_output <- vector()

## algorithm ##
time_counter <- 0
timestep <- 1/365 # yr
for (i in 1:t_extent){
  # start a new year
  print(paste("year",i))
  daily_wtd <- vector() # cm, for calculating annual average wtd
  daily_prod <- vector() # g cm-2 day-1, for calculating annual production
  daily_decay <- vector() # g cm-2 day-1, for calculating annual decay
  daily_melt <- vector() # cm
  daily_pond <- vector() # cm
  daily_aet <- vector() # cm
  # daily climate inputs
  if (i>=11508){
    climate_data <- generate_station_daily(1950-(11500-i+1))
  } else {
    climate_data <- generate_trace_daily(i)
  }
  daily_t <- climate_data[[3]] # degree C
  daily_p <- climate_data[[4]]/10 # cm
  reg <- monthly_pet_equation(climate_data[[1]]) 
  daily_pet <- calculate_pet(daily_t,reg[1],reg[2]) # mm
  # daily step
  repeat{
    time_counter <- time_counter+1
    if (time_counter > 365) {break}
    # non-freezing
    if (daily_t[time_counter] > 0){
      # decomposition
      oxic_decay <- oxic_decay_base*Q10_oxic^((daily_t[time_counter]-base_temp)/10) # yr-1
      anoxic_decay <- anoxic_decay_base*Q10_anoxic^((daily_t[time_counter]-base_temp)/10) # yr-1
      temporary_mass <- sum(layer_mass[!is.na(layer_mass)]) # g cm-2, mass before daily decay, for calculating daily decay
      for (j in 1:i){ # mass/thickness/elevation for each cohort
        layer_mass[j] <- layer_mass[j]*((1-wet_proportion[j])*exp(-timestep*layer_remaining_mass[j]^x_factor*oxic_decay)+
                                          wet_proportion[j]*exp(-timestep*layer_remaining_mass[j]^x_factor*anoxic_decay)) # g cm-2
        layer_remaining_mass[j] <- layer_mass[j]/layer_initial_mass[j]
        layer_thickness[j] <- layer_mass[j]/layer_density[j] # cm
        if (j==1){
          layer_elevation[j] <- layer_thickness[j] # cm
        } else {
          layer_elevation[j] <- layer_thickness[j]+layer_elevation[j-1] # cm
        }
      }
      peat_height <- layer_elevation[i] # cm, peat height at top
      daily_decay <- c(daily_decay,temporary_mass-sum(layer_mass[!is.na(layer_mass)])) # g cm-2
      # snowmelt
      if (snow > 0){
        snowmelt <- daily_t[time_counter]*ddf
        if (snowmelt > snow){
          snowmelt <- snow # no enough snow
        }
        snow <- snow-snowmelt
      } else {
        snowmelt <- 0
      }
      # aet and production calculation
      aet_prod <- aet_prod_func(wtd,daily_pet[time_counter]*365,z1_soil,z2_soil,z_plant,Ellenberg,LAI,prod_scaling)
      aet <- aet_prod[1]/365/10 # cm, daily aet
      prod <- aet_prod[2]/365 # g cm-2 yr-1, daily production
      # wt and pond
      wt_change_flux <- timestep*dHdt(daily_p[time_counter]*365,aet*365,snowmelt*365,porosity,peat_transmissivity,wt_height,lateral_extent,
                                      drain_effi)*porosity # cm, wt change assuming no porosity for pond level calculation
      pond_test <- pond # remember this pond value
      if ((pond_test>0)&(wt_change_flux>=0)){ # w/ pond, wetter
        pond <- pond+wt_change_flux
        wt_height <- peat_height
        if (pond>pond_size) {pond <- pond_size}
      }
      if ((pond_test>0)&(wt_change_flux<0)){ # w/ pond, drier
        if ((pond+wt_change_flux)>=0){
          pond <- pond+wt_change_flux
          wt_height <- peat_height
        } else {
          wt_height <- wt_height+(pond+wt_change_flux)/porosity
          pond <- 0
        }
      }
      if ((pond_test==0)&(wt_change_flux>0)){ # w/o pond, wetter
        if ((wt_height+wt_change_flux/porosity)<=peat_height){
          pond <- 0
          wt_height <- wt_height+wt_change_flux/porosity
        } else {
          pond <- ((wt_height+wt_change_flux/porosity)-peat_height)*porosity
          wt_height <- peat_height
          if (pond>pond_size) {pond <- pond_size}
        }
      }
      if ((pond_test==0)&(wt_change_flux<=0)){ # w/o pond, drier
        pond <- 0
        wt_height <- wt_height+wt_change_flux/porosity
      }
      if (wt_height <= 0) {wt_height <- 0.0000001} # force no negative wt
      wtd <- peat_height-wt_height-pond # cm, daily wtd
      # peat property
      for (j in 1:i){ # wet proportion/conductivity/transmissivity for each cohort
        if (j==1){
          if (wt_height >= layer_elevation[j]) {wet_proportion[1] <- 1} else {wet_proportion[1] <- wt_height/layer_thickness[1]}
        } else {
          if (wt_height >= layer_elevation[j]) {wet_proportion[j] <- 1}
          else if (wt_height <= layer_elevation[j-1]) {wet_proportion[j] <- 0}
          else {wet_proportion[j] <- (wt_height-layer_elevation[j-1])/layer_thickness[j]}
        }
        layer_hydro_k[j] <- k_param_a*(exp(k_param_b*layer_remaining_mass[j])) # cm yr-1
        layer_transmissivity[j] <- layer_hydro_k[j]*layer_thickness[j]*wet_proportion[j] # cm2 yr-1
        layer_density[j] <- derived_density(layer_hydro_k[j]) # g cm-3
      }
      peat_transmissivity <- sum(layer_transmissivity[!is.na(layer_transmissivity)]) # cm2 yr-1
      # production
      daily_prod <- c(daily_prod,prod) # g cm-2 day-1
    } else {
      daily_decay <- c(daily_decay,0) # no decay
      snow <- snow+daily_p[time_counter] # snow accumulation
      snowmelt <- 0 # no snowmelt
      wtd <- peat_height-wt_height-pond # cm, daily wtd
      aet <- 0 # cm, daily wtd
      prod <- 0 # no production
      daily_prod <- c(daily_prod,prod) # should be no production
    }
    daily_wtd <- c(daily_wtd,wtd)
    daily_melt <- c(daily_melt,snowmelt)
    daily_pond <- c(daily_pond,pond)
    daily_aet <- c(daily_aet,aet)
  }
  # prepare the next year
  layer_mass[i+1] <- sum(daily_prod) # g cm-2, new production by summing up all daily
  layer_initial_mass[i+1] <- layer_mass[i+1] # g cm-2
  layer_remaining_mass[i+1] <- 1
  layer_density[i+1] <- initial_density # g cm-3
  layer_thickness[i+1] <- layer_mass[i+1]/layer_density[i+1] # cm
  layer_elevation[i+1] <- peat_height+layer_thickness[i+1] # cm
  wet_proportion[i+1] <- 0
  # final annual time series output
  wtd_output[i] <- mean(daily_wtd) # cm, annual average wtd
  gswtd_output[i] <- mean(daily_wtd[daily_t>0]) # cm, growing season average wtd
  prod_output[i] <- layer_mass[i+1] # g cm-2
  decay_output[i] <- sum(daily_decay) # g cm-2
  height_output[i] <- layer_elevation[i+1] # cm, after new layer added
  mass_output[i] <- sum(layer_mass[!is.na(layer_mass)]) # g cm-2, after new layer added
  # final intermediate seasonal time series output
  daily_wtd_output <- rbind(daily_wtd_output,daily_wtd)
  daily_prod_output <- rbind(daily_prod_output,daily_prod)
  daily_decay_output <- rbind(daily_decay_output,daily_decay)
  daily_melt_output <- rbind(daily_melt_output,daily_melt)
  daily_pond_output <- rbind(daily_pond_output,daily_pond)
  daily_temp_output <- rbind(daily_temp_output,daily_t)
  daily_precip_output <- rbind(daily_precip_output,daily_p)
  daily_pet_output <- rbind(daily_pet_output,daily_pet/10) # cm
  daily_aet_output <- rbind(daily_aet_output,daily_aet)
  # end
  time_counter <- 0
}