library("readxl")
data1 <- read_excel("validation.xlsx",sheet=1)
data2 <- read_excel("validation.xlsx",sheet=2)
data3 <- read_excel("validation.xlsx",sheet=3)

par(mfrow=c(2,2))
depth <- as.vector(t(data1[,1]))
age <- as.vector(t(data1[,2]))
plot(age,depth,ylim=c(700,0),type="l")
lines(length(layer_elevation):1,layer_elevation[length(layer_elevation)]-layer_elevation,col="red")

core <- as.vector(t(data2[,1]))
depth <- as.vector(t(data2[,2]))
ksat <- as.vector(t(data2[,3]))
plot(ksat*31560000*100,depth,ylim=c(80,0),xlim=c(1e4,1e7),log="x",bg=core,pch=21,cex=1.5)
lines(layer_hydro_k,layer_elevation[length(layer_elevation)]-layer_elevation[1:(length(layer_elevation)-1)],col="red")

doy <- as.vector(t(data3[,1]))
wt_2022 <- as.vector(t(data3[,2]))
plot(doy,-wt_2022-mean(-wt_2022),type="l",ylim=c(20,-20))
lines(1:365,daily_wtd_output[nrow(daily_wtd_output),]-mean(daily_wtd_output[nrow(daily_wtd_output),]),col="red")

plot(1:573,wtd_output[11001:11573],type="l",ylim=c(50,0),col="red",lwd=0.5)
text(500,30,round(mean(wtd_output[(length(wtd_output))-50:length(wtd_output)])))

plot(1822:2022,rowSums(daily_precip_output[11373:11573,]),type="l")