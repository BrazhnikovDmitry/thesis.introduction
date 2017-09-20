require("ncdf")
require("oce")
require("fields")
require("misc3d")
require("akima")

### ETOPO1
get_etopo <- function()
{
	path <- "/home/dmitry/Work/datasets/TSUNAMI/ETOPO1_Ice_g_gmt4.grd"

	nc <- open.ncdf(path)

	x <- get.var.ncdf(nc, "x") ### Longitude
	y <- get.var.ncdf(nc, "y") ### Latitude

	### Define cut-region for Koko guyot
	###  Koko's coords: 35.25N 171.58E
	x0 <- 171.58
	y0 <- 35.25

	lim_x <- c(x0-1.5, x0+1)
	lim_y <- c(y0-1.5, y0+1)

	### Cut that place, find the nearest indices
	i1 <- which.min (abs(lim_x[1] - x))
	i2 <- which.min (abs(lim_x[2] - x))

	j1 <- which.min (abs(lim_y[1] - y))
	j2 <- which.min (abs(lim_y[2] - y))

	koko_1 <- get.var.ncdf(nc, "z", start=c(i1, j1), count=c(i2-i1+1, j2-j1+1))
	
	return (list(koko=koko_1, lon=x[i1:i2], lat=y[j1:j2]))
}

### TOPO30
get_topo30 <- function(dx1=1, dx2=1, dy1=1, dy2=1)
{
	path <- "/home/dmitry/Work/datasets/TSUNAMI/topo30.grd"

	nc <- open.ncdf(path)

	x <- get.var.ncdf(nc, "lon") ### Longitude
	y <- get.var.ncdf(nc, "lat") ### Latitude

	### Define cut-region for Koko guyot
	###  Koko's coords: 35.25N 171.58E
	x0 <- 171.58
	y0 <- 35.25

	lim_x <- c(x0-dx1, x0+dx2)
	lim_y <- c(y0-dy1, y0+dy2)

	### Cut that place, find the nearest indices
	i1 <- which.min (abs(lim_x[1] - x))
	i2 <- which.min (abs(lim_x[2] - x))

	j1 <- which.min (abs(lim_y[1] - y))
	j2 <- which.min (abs(lim_y[2] - y))

	koko_30 <- get.var.ncdf(nc, "z", start=c(i1, j1), count=c(i2-i1+1, j2-j1+1))
	print(paste(i1, i2, j1, j2))
	return (list(koko=koko_30, lon=x[i1:i2], lat=y[j1:j2]))
}

### Plot
koko_1 <- get_etopo()
koko_30 <- get_topo30()
gr <- expand.grid(x=koko_1$lon, y=koko_1$lat)
out <- interp(gr[,1], gr[,2], as.vector(koko_1$koko), xo=koko_30$lon, yo=koko_30$lat)

#pdf("emp_seam.pdf", width=4, height=10)
tiff("emp_seam.tiff", width=400, height=1000)
image.plot(koko_30$lon, koko_30$lat, koko_30$koko, col=two.colors(n=100, start="black", middle="grey10", end="white"), xlab="Longitude", ylab="Latitude", legend.lab="Depth, m")
dev.off()
two.colors(n=100, start="white", middle="salmon", end="orangered3 ")

pdf("koko.pdf", width=8.5, height=10)
par(cex.lab=1.7)
par(cex.axis=1.7)
image.plot(koko_30$lon, koko_30$lat, koko_30$koko, xlab="Longitude", ylab="Latitude", horizontal=T)
contour(koko_30$lon, koko_30$lat, koko_30$koko, add=T, labcex=2, col="white")
dev.off()
#svg("koko.svg", width=7, height=7)
par(mfrow=c(1,3))

image.plot(koko_1$lon, koko_1$lat, koko_1$koko, xlab="Longitude", ylab="Latitude", horizontal=T)
contour(koko_1$lon, koko_1$lat, koko_1$koko, add=T, labcex=1, col="white")
title(main="Koko guyot, ETOPO1")

image.plot(koko_30$lon, koko_30$lat, koko_30$koko, xlab="Longitude", ylab="Latitude", horizontal=T)
contour(koko_30$lon, koko_30$lat, koko_30$koko, add=T, labcex=1, col="white")
title(main="Koko guyot, TOPO30")

image.plot(koko_30$lon, koko_30$lat, out$z-koko_30$koko, horizontal=T)
contour(koko_30$lon, koko_30$lat, out$z-koko_30$koko, add=T, labcex=1, col="white")
title(main="Diff two bathys, ETOPO1 interpolated to 30' grid \n - TOPO30")

dev.off()

##### Make plot of transects in EW direction
pdf("transects.pdf", width=15, height=5)
xx <- geodDist(rep(koko_1$lat[1], length(koko_1$lon)), koko_1$lon, alongPath=T)

#lims_x <- c(min(koko_1$lon), max(koko_1$lon))
lims_x <- c(min(xx), max(xx))
lims_y <- c(min(koko_1$koko), max(koko_1$koko))
par(mfrow=c(1,3))
plot(lims_x, lims_y, type="n")
for (i in 60:120) {#1:length(koko_1$lat)) {
	xx <- geodDist(rep(koko_1$lat[i], length(koko_1$lon)), koko_1$lon, alongPath=T)
	lines(xx, koko_1$koko[,i])	# koko_1$lon
}
lim2aver <- 1:dim(koko_1$koko)[1]#60:120
bb <- apply(koko_1$koko[lim2aver,], 1, mean)
lines(xx, bb, col="red", lwd=2)	# koko_1$lon

lim2aver <- 60:120
bb <- apply(koko_1$koko[lim2aver,], 1, mean)
lines(xx, bb, col="blue", lwd=2)	# koko_1$lon
title(main="EW transects")
#image.plot(1:length(koko_1$lon), 1:length(koko_1$lat), koko_1$koko)

##### Make plot of transects in SN direction
#X11()
xx <- geodDist(koko_1$lat, rep(koko_1$lon[1], length(koko_1$lat)), alongPath=T)

#lims_x <- c(min(koko_1$lon), max(koko_1$lon))
lims_x <- c(min(xx), max(xx))
lims_y <- c(min(koko_1$koko), max(koko_1$koko))
#par(mfrow=c(1,2))
plot(lims_x, lims_y, type="n")
for (j in 1:length(koko_1$lon)) {
	xx <- geodDist(koko_1$lat, rep(koko_1$lon[j], length(koko_1$lat)), alongPath=T)
	lines(xx, koko_1$koko[j,])	# koko_1$lon
}
lim2aver <- 1:dim(koko_1$koko)[2]#60:120
bb <- apply(koko_1$koko[,lim2aver], 2, mean)
lines(xx, bb, col="red", lwd=2)	# koko_1$lon

lim2aver <- 60:120
bb <- apply(koko_1$koko[,lim2aver], 2, mean)
lines(xx, bb, col="blue", lwd=2)	# koko_1$lon
title(main="SN transects")

xx1 <- geodDist(rep(koko_1$lat[1], length(koko_1$lon)), koko_1$lon, alongPath=T)
xx2 <- geodDist(koko_1$lat, rep(koko_1$lon[1], length(koko_1$lat)), alongPath=T)
#image.plot(1:length(koko_1$lon), 1:length(koko_1$lat), koko_1$koko)
image.plot(xx1, xx2, koko_1$koko)
dev.off()
### Ugly so far
surface3d(i1:i2, j1:j2, koko/100)



########### Make ellipse out of Koko guoyt
require("ncdf")
require("oce")
require("fields")
require("misc3d")
require("akima")
require("plotrix")

koko_30 <- get_topo30(dx1=0.75, dx2=0.75, dy1=0.75, dy2=0.75)

# koko_30$lon, koko_30$lat, 
pdf("koko_ell.pdf")
image.plot(1:length(koko_30$lon), 1:length(koko_30$lat), koko_30$koko, col=two.colors(n=100, start="black", middle="grey10", end="white"), xlab="Longitude", ylab="Latitude", legend.lab="Depth, m")
contour(1:length(koko_30$lon), 1:length(koko_30$lat), koko_30$koko, levels=c(-1000, -1500, -800, seq(-400, -600, by=-50)), add=T)
p_lon <- c(25, 140, 80, 120)
p_lat <- c(150, 50, 80, 120)
points(p_lon, p_lat, col="red")

maj <- sqrt( (p_lat[1] - p_lat[2])^2 + (p_lon[1] - p_lon[2])^2 )/2
min <- sqrt( (p_lat[3] - p_lat[4])^2 + (p_lon[3] - p_lon[4])^2 )/2
angle <- atan2(p_lat[1] - c_lat, p_lon[1] - c_lon)*180/pi - 10
c_lat <- 103; c_lon <- 87
draw.ellipse(c_lon, c_lat, maj, min, angle, lwd=2, border="navy")
dev.off()

maj <- geodDist(koko_30$lat[p_lat[1]], koko_30$lon[p_lon[1]], koko_30$lat[p_lat[2]], koko_30$lon[p_lon[2]])/2
min <- geodDist(koko_30$lat[p_lat[3]], koko_30$lon[p_lon[3]], koko_30$lat[p_lat[4]], koko_30$lon[p_lon[4]])/2

points(c_lon, c_lat, col="navy")

