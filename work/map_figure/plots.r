require("ncdf")
require("oce")
require("fields")
require("misc3d")
require("akima")

### Plot Koko and Tassie shelf together
path <- "/home/dmitry/Work/datasets/TSUNAMI/ETOPO1_Ice_g_gmt4.grd"

nc <- open.ncdf(path)

x <- get.var.ncdf(nc, "x") ### Longitude
y <- get.var.ncdf(nc, "y") ### Latitude

### Define cut-region
lim_x <- c(130, 180)
lim_y <- c(-60, 60)

### Cut that place, find the nearest indices
i1 <- which.min (abs(lim_x[1] - x))
i2 <- which.min (abs(lim_x[2] - x))

j1 <- which.min (abs(lim_y[1] - y))
j2 <- which.min (abs(lim_y[2] - y))

bp <- get.var.ncdf(nc, "z", start=c(i1, j1), count=c(i2-i1+1, j2-j1+1))

idx <- i1:i2
jdx <- j1:j2

del <- 10
sub_x <- seq(1, dim(bp)[1], by=del)
sub_y <- seq(1, dim(bp)[2], by=del)

bp[which (bp > 0)] <- 0
image.plot(x[idx[sub_x]], y[jdx[sub_y]], bp[sub_x, sub_y], col=oceColorsGebco(500), horizontal=TRUE, xlab="", ylab="")

##############
lim_x <- c(170, 173)
lim_y <- c(33.75, 36)

i1 <- which.min (abs(lim_x[1] - x))
i2 <- which.min (abs(lim_x[2] - x))

j1 <- which.min (abs(lim_y[1] - y))
j2 <- which.min (abs(lim_y[2] - y))

bp <- get.var.ncdf(nc, "z", start=c(i1, j1), count=c(i2-i1+1, j2-j1+1))

idx <- i1:i2
jdx <- j1:j2

del <- 1
sub_x <- seq(1, dim(bp)[1], by=del)
sub_y <- seq(1, dim(bp)[2], by=del)

bp[which (bp > 0)] <- 0
image.plot(x[idx[sub_x]], y[jdx[sub_y]], bp[sub_x, sub_y], col=oceColorsGebco(500), horizontal=TRUE, xlab="", ylab="")

##############
path <- "/home/dmitry/Work/datasets/TSUNAMI/topo30.grd"

nc <- open.ncdf(path)

x <- get.var.ncdf(nc, "lon") ### Longitude
y <- get.var.ncdf(nc, "lat") ### Latitude

### Define cut-region
lim_x <- c(146, 155)
lim_y <- c(-45, -40.5)

i1 <- which.min (abs(lim_x[1] - x))
i2 <- which.min (abs(lim_x[2] - x))

j1 <- which.min (abs(lim_y[1] - y))
j2 <- which.min (abs(lim_y[2] - y))

bp <- get.var.ncdf(nc, "z", start=c(i1, j1), count=c(i2-i1+1, j2-j1+1))

idx <- i1:i2
jdx <- j1:j2

del <- 1
sub_x <- seq(1, dim(bp)[1], by=del)
sub_y <- seq(1, dim(bp)[2], by=del)

bp[which (bp > 0)] <- 0
image.plot(x[idx[sub_x]], y[jdx[sub_y]], bp[sub_x, sub_y], col=oceColorsGebco(500), horizontal=TRUE, xlab="", ylab="")
