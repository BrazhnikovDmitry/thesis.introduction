require("ncdf")
require("oce")
require("fields")
require("misc3d")
require("akima")
require("maptools")

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

koko_1 <- get_etopo()
koko_30 <- get_topo30()

### Tsunami spectra
g <- 9.81

solve_disp_relation <- function(om, H)
{
	disp_relation <- function(x) return (om - sqrt(g*x*tanh(x*H)))
	
	k <- uniroot(disp_relation, c(0, 1e3), tol=1e-15)$root
	
	return (k)
}

T <- seq(1/60, 15, len=1e3)*60

om <- 2*pi/T; H <- c(350, 1000, 2e3, 3e3, 5.5e3)

k <- array(0, c(length(T), length(H)))

for (j in 1:length(H)) {
	for (i in 1:length(T)) {
		k[i,j] <- solve_disp_relation(2*pi/T[i], H[j])
	}
}

lambda <- 2*pi/k
cr <- t(apply(lambda, c(1), "/", H))

COLS <- c("black", "darkred", "navy", "lightsalmon4", "saddlebrown")

pdf("limit_tsunami.pdf", width=10, height=6)
par(cex.axis=1.5)
par(cex.lab=1.5)
par(mar=c(5,5,4,2)+0.1)
plot(c(min(T/60), max(T/60)+2), c(min(cr), max(cr)), type="n", xlab="Wave period, mins", ylab=expression(lambda/H))#, xlim=c(0, 15), ylim=c(0, 200))
for (j in 1:length(H)) {
	lines(T/60, cr[,j], col=COLS[j], lwd=2.5)
}
axis(1, at=8, labels=8)
#pointLabel(rep(1.01, length(H)), cr[1e3,], labels=paste(H, "m"), col=COLS, cex=1.5, srt=90)
text(rep(16, length(H)), cr[1e3,], labels=paste(H, "m"), col=COLS, cex=1.5)
abline(h=20, lwd=2, lty=2)
abline(v=8, lwd=2, lty=2)
text(16, 20, cex=1.5, labels=expression(lambda/H==20), pos=1)
dev.off()

l <- 2*pi/seq(10, 100e3)
H <- 1e3
om <-  sqrt(g*2*pi/l*tanh(2*pi/l*H))
c <- om/2*pi/l
