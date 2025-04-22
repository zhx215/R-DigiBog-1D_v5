## moving average ##
x_scale <- c(0,t_extent)
window <- 50
yr <- seq(1+window/2,t_extent-window/2,1)
rm <- layer_remaining_mass[1:t_extent]
mb <- (prod_output-decay_output)*10000 # g m-2 yr-1, mass balance
mar <- layer_mass[1:t_extent]*10000 # g m-2 yr-1, mass accumulation rate
ms <- mass_output*10000
wtd_window <- gswtd_window <- rm_window <- mb_window <- mar_window <- vector()
for (i in 1:length(yr)){
  wtd_window[i] <- mean(wtd_output[(yr[i]-window/2):(yr[i]+window/2)])
  gswtd_window[i] <- mean(gswtd_output[(yr[i]-window/2):(yr[i]+window/2)])
  rm_window[i] <- mean(layer_remaining_mass[(yr[i]-window/2):(yr[i]+window/2)])
  mb_window[i] <- mean(mb[(yr[i]-window/2):(yr[i]+window/2)])
  mar_window[i] <- mean(mar[(yr[i]-window/2):(yr[i]+window/2)])
}

## plotting ##
# annual time series
par(mfrow=c(5,1),mar=c(4.5,5,1,2))
plot(NA,NA,xlim=x_scale,ylim=c(max(wtd_output[!is.na(wtd_output)]),min(wtd_output[!is.na(wtd_output)])),xlab="time yr",ylab="WTD cm")
lines(seq(1,t_extent,1),wtd_output,col="gray")
lines(yr,wtd_window,lwd=2)

plot(NA,NA,xlim=x_scale,ylim=c(min(rm[!is.na(rm)]),max(rm[!is.na(rm)])),xlab="time yr",ylab="remained mass")
lines(seq(1,t_extent,1),rm,col="gray")
lines(yr,rm_window,lwd=2)

plot(NA,NA,xlim=x_scale,ylim=c(min(mb[!is.na(mb)]),max(mb[!is.na(mb)])),xlab="time yr",ylab="mass balance g m-2 yr-1")
lines(seq(1,t_extent,1),mb,col="gray")
lines(yr,mb_window,lwd=2)

plot(NA,NA,xlim=x_scale,ylim=c(min(mar[!is.na(mar)]),max(mar[!is.na(mar)])),xlab="time yr",ylab="mass accmulation rate g m-2 yr-1")
lines(seq(1,t_extent,1),mar,col="gray")
lines(yr,mar_window,lwd=2)

plot(NA,NA,xlim=x_scale,ylim=c(0,max(ms)),xlab="time yr",ylab="total mass g m-2")
lines(seq(1,t_extent,1),ms,lwd=1)

# seasonal pattern
par(mfrow=c(4,1),mar=c(4.5,5,1,2))
daily_wtd_50 <- daily_wtd_plus68 <- daily_wtd_minus68 <- daily_wtd_plus95 <- daily_wtd_minus95 <- vector()
for (i in 1:365){
  rank <- sort(daily_wtd_output[,i])[t_extent*c(2.5,16,50,84,97.5)/100] # 2.5%, 16%, 50%, 84%, 97.5%
  daily_wtd_50 <- c(daily_wtd_50,rank[3])
  daily_wtd_plus95 <- c(daily_wtd_plus95,rank[5])
  daily_wtd_minus95 <- c(daily_wtd_minus95,rank[1])
  daily_wtd_plus68 <- c(daily_wtd_plus68,rank[4])
  daily_wtd_minus68 <- c(daily_wtd_minus68,rank[2])
}
plot(NA,NA,xlim=c(1,365),ylim=c(max(daily_wtd_plus95),min(daily_wtd_minus95)),xlab="doy",ylab="WTD cm")
polygon(c(1:365,365:1),c(daily_wtd_plus95,daily_wtd_minus95[365:1]),col="gray",border=NA)
polygon(c(1:365,365:1),c(daily_wtd_plus68,daily_wtd_minus68[365:1]),col="darkgray",border=NA)
lines(1:365,daily_wtd_50)

daily_mb_50 <- daily_mb_plus68 <- daily_mb_minus68 <- daily_mb_plus95 <- daily_mb_minus95 <- vector()
for (i in 1:365){
  rank <- sort(((daily_prod_output-daily_decay_output)*10000)[,i])[t_extent*c(2.5,16,50,84,97.5)/100] # 2.5%, 16%, 50%, 84%, 97.5%
  daily_mb_50 <- c(daily_mb_50,rank[3])
  daily_mb_plus95 <- c(daily_mb_plus95,rank[5])
  daily_mb_minus95 <- c(daily_mb_minus95,rank[1])
  daily_mb_plus68 <- c(daily_mb_plus68,rank[4])
  daily_mb_minus68 <- c(daily_mb_minus68,rank[2])
}
plot(NA,NA,xlim=c(1,365),ylim=c(max(daily_mb_plus95),min(daily_mb_minus95)),xlab="doy",ylab="mass balance g m-2 day-1")
polygon(c(1:365,365:1),c(daily_mb_plus95,daily_mb_minus95[365:1]),col="gray",border=NA)
polygon(c(1:365,365:1),c(daily_mb_plus68,daily_mb_minus68[365:1]),col="darkgray",border=NA)
lines(1:365,daily_mb_50)

daily_melt_50 <- daily_melt_plus68 <- daily_melt_minus68 <- daily_melt_plus95 <- daily_melt_minus95 <- vector()
for (i in 1:365){
  rank <- sort(daily_melt_output[,i])[t_extent*c(2.5,16,50,84,97.5)/100] # 2.5%, 16%, 50%, 84%, 97.5%
  daily_melt_50 <- c(daily_melt_50,rank[3])
  daily_melt_plus95 <- c(daily_melt_plus95,rank[5])
  daily_melt_minus95 <- c(daily_melt_minus95,rank[1])
  daily_melt_plus68 <- c(daily_melt_plus68,rank[4])
  daily_melt_minus68 <- c(daily_melt_minus68,rank[2])
}
plot(NA,NA,xlim=c(1,365),ylim=c(min(daily_melt_plus95),max(daily_melt_minus95)),xlab="doy",ylab="snowmelt cm day-1")
polygon(c(1:365,365:1),c(daily_melt_plus95,daily_melt_minus95[365:1]),col="gray",border=NA)
polygon(c(1:365,365:1),c(daily_melt_plus68,daily_melt_minus68[365:1]),col="darkgray",border=NA)
lines(1:365,daily_melt_50)

daily_pond_50 <- daily_pond_plus68 <- daily_pond_minus68 <- daily_pond_plus95 <- daily_pond_minus95 <- vector()
for (i in 1:365){
  rank <- sort(daily_pond_output[,i])[t_extent*c(2.5,16,50,84,97.5)/100] # 2.5%, 16%, 50%, 84%, 97.5%
  daily_pond_50 <- c(daily_pond_50,rank[3])
  daily_pond_plus95 <- c(daily_pond_plus95,rank[5])
  daily_pond_minus95 <- c(daily_pond_minus95,rank[1])
  daily_pond_plus68 <- c(daily_pond_plus68,rank[4])
  daily_pond_minus68 <- c(daily_pond_minus68,rank[2])
}
plot(NA,NA,xlim=c(1,365),ylim=c(min(daily_pond_minus95),max(daily_pond_plus95)),xlab="doy",ylab="snowpond cm day-1")
polygon(c(1:365,365:1),c(daily_pond_plus95,daily_pond_minus95[365:1]),col="gray",border=NA)
polygon(c(1:365,365:1),c(daily_pond_plus68,daily_pond_minus68[365:1]),col="darkgray",border=NA)
lines(1:365,daily_pond_50)

# correlation
par(mfrow=c(3,2),mar=c(4.5,5,1,2))
plot(rowSums(daily_precip_output[(t_extent-1000):t_extent,])*10,mb[(t_extent-1000):t_extent],xlab="precip mm yr-1",ylab="mass balance g m-2 yr-1")
plot(rowMeans(daily_temp_output[(t_extent-1000):t_extent,]),mb[(t_extent-1000):t_extent],xlab="annual temp degree C",ylab="mass balance g m-2 yr-1")