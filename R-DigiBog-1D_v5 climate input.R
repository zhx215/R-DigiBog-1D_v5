## read data ##
donggang_data <- read.csv("Donggang station.csv") # year, month, day, temp degree C, precip mm
trace_data <- read.csv("Trace11500.csv") # 2-13 temp degree C, 14-25 precip mm

## calculate monthly pet and establish a pet-temp relationship ##
k_value <- c(80,89,99,110,120,125,123,115,104,93,83,78)/100 # for Thornthwaite PET calculation at 40N (from Ponce, 1958)
monthly_pet_equation <- function(monthly_temp){
  temp <- monthly_temp
  temp[temp<0] <- 0
  I <- (temp/5)^1.514
  J <- sum(I)
  c <- 0.000000675*J^3-0.0000771*J^2+0.01792*J+0.49239
  pet <- k_value*1.6*(10*temp/J)^c # cm
  regres <- summary(lm(pet[temp>0]~temp[temp>0]))$coefficients
  intercept <- regres[1]
  slope <- regres[2]
  return(c(slope,intercept)) # predict monthly pet (cm/month) from monthly temp, degree C
}

## calculate daily pet ##
divide_day <- c(rep(31,31),rep(28,28),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31))
calculate_pet <- function(daily_temp,a,b){
  daily_pet <- (daily_temp*a+b)*10/divide_day # mm/day
  daily_pet[daily_pet<0] <- 0
  return(daily_pet) # unit same as input
}

## select a random year and return Donggang station data ##
donggang_select_year <- function(which_year){
  if (which_year=="random"){
    year <- sample(seq(1957,2022,1),1)
  } else {
    year <- which_year
  }
  id <- which(donggang_data[,1]==year)
  selected_data <- donggang_data[id,]
  monthly_temp <- monthly_precip <- daily_temp <- daily_precip <- vector()
  for (mn in 1:12){
    id <- which(selected_data[,2]==mn)
    if (mn==2){
      id <- id[1:28] # give up Feb 29
    }
    monthly_temp <- c(monthly_temp,mean(selected_data[id,4]))
    monthly_precip <- c(monthly_precip,sum(selected_data[id,5]))
    daily_temp <- c(daily_temp,selected_data[id,4])
    daily_precip <- c(daily_precip,selected_data[id,5])
  }
  return(list(monthly_temp,monthly_precip,daily_temp,daily_precip))
}

## generate trace daily temp and precip data for a year ##
first_day <- c(1,32,60,91,121,152,182,213,244,274,305,335,366) # last day - 1 = 365
temp_correction <- 4+1.5 # degree C
generate_trace_daily <- function(yr){
  monthly_temp <- as.vector(t(trace_data[yr,2:13])) # trace
  monthly_precip <- as.vector(t(trace_data[yr,14:25])) # trace
  select_year_output <- donggang_select_year("random")
  select_year_monthly_temp <- select_year_output[[1]] # station
  select_year_monthly_precip <- select_year_output[[2]] # station
  select_year_daily_temp <- select_year_output[[3]] # station
  select_year_daily_precip <- select_year_output[[4]] # station
  monthly_temp_diff <- monthly_temp-select_year_monthly_temp # trace-station temp
  monthly_precip_ratio <- monthly_precip/select_year_monthly_precip # trace/station precip
  daily_temp <- daily_precip <- rep(NA,365)
  for (i in 1:12){
    daily_temp[first_day[i]:(first_day[i+1]-1)] <- select_year_daily_temp[first_day[i]:(first_day[i+1]-1)]+monthly_temp_diff[i]
    daily_precip[first_day[i]:(first_day[i+1]-1)] <- select_year_daily_precip[first_day[i]:(first_day[i+1]-1)]*monthly_precip_ratio[i]
  }
  return(list(monthly_temp-temp_correction,monthly_precip,daily_temp-temp_correction,daily_precip))
}

## collect Donggang station daily and precip data for a year ##
temp_correction2 <- 4 # degree C
generate_station_daily <- function(yr){
  get <- donggang_select_year(yr)
  o1 <- get[[1]]-temp_correction2
  o2 <- get[[2]]
  o3 <- get[[3]]-temp_correction2
  o4 <- get[[4]]
  return(list(o1,o2,o3,o4))
}