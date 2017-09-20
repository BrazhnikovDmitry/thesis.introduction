source("/home/dmitry/Work/toolkit_it/misc.r")

#######################################
### Error, measured as mode 2 amplitude
#mesh0 <- mesh
load("final_w_sub_1.dat")
mesh0 <- mesh

source("sol_wunsch.r")
#mesh <- mesh0

alp <- (mesh0$H[1] - mesh0$H[mesh0$nx])/(mesh0$xx[1] - mesh0$xx[mesh0$nx])
c <- sqrt(om^2/(N2 - om^2))
gamma <- atan(alp)/c[1]

sol_wunsch_d2 <- analyt_sol(mesh0, n=1, N2=N2[1], gamma=gamma)

tm <- 1
num_ts <- 3
i0 <- which.min(abs(mesh$H - 0.5e3))
i1 <- which.min(abs(mesh$H - 3e3))
seq_x <- i0:i1

seq_ts <- seq(i0+30, i1-30, len=num_ts)
dx <- diff(mesh$xx[seq_ts])[1]
dx <- dx*0.5
scale <- dx/max(abs(sol_wunsch_d2$uu[seq_ts, 1, , tm]), na.rm=T)

#for (tm in 1:mesh$nt) {
#	png(paste0(tm, ".png"), width=600, height=250)
	pdf(paste0("gamma_d2", ".pdf"), width=10, height=4)
	plot(mesh0$xx[seq_x]/1e3, -mesh$H[seq_x], type="l", xlab="x coordinate, km", ylab="depth, m", xlim=c(mesh$xx[i0], mesh$xx[i1])/1e3, ylim=c(-max(mesh$H[seq_x]), 0))

	for (i in seq_ts) {
		j0 <- which.min(abs(mesh$H[i] + mesh$zz))
		lines((mesh$xx[i] + sol_wunsch_d2$uu[i,1,1:j0,tm]*scale)/1e3, mesh$zz[1:j0], lwd=2)
		lines(rep(mesh$xx[i], j0)/1e3, mesh$zz[1:j0], lty=2)
		lines(c(mesh$xx[i]/1e3, (mesh$xx[i] + sol_wunsch_d2$uu[i,1,1,tm]*scale)/1e3), c(mesh$zz[1], mesh$zz[1]), lty=2)
		lines(c(mesh$xx[i]/1e3, (mesh$xx[i] + sol_wunsch_d2$uu[i,1,j0,tm]*scale)/1e3), c(mesh$zz[j0], mesh$zz[j0]), lty=2)
	}
		
	tmp <- which( sol_wunsch_d2$uu[seq_ts, 1, , tm] == max(sol_wunsch_d2$uu[seq_ts, 1, , tm], na.rm=T), arr.ind=T)
	i_a <- seq_ts[tmp[1]]; j_a <- tmp[2]
	j0 <- which.min(abs(mesh$H[i_a] + mesh$zz))
	arrows( mesh$xx[i_a]/1e3, mesh$zz[j_a], (mesh$xx[i_a] + sol_wunsch_d2$uu[i_a,1,j_a,tm]*scale)/1e3, mesh$zz[j_a], lty=1, lwd=2, angle=20, len=0.15)
	text((mesh$xx[i_a] + sol_wunsch_d2$uu[i_a,1,j_a,tm]*scale)/1e3, mesh$zz[j_a], pos=4, labels=sprintf("u = %.1f cm/s", sol_wunsch_d2$uu[i_a,1,j_a,tm]*1e3), lwd=2, cex=1.4, col="darkred", font=2)

	title(main=expression(gamma==0.2), cex.main=1.4)

	dev.off()
#}

### gamma = 0.7
gamma <- 0.7
c <- sqrt(om^2/(N2 - om^2))
alp <- tan(gamma*c)[1]
mesh$H[,1] <- (mesh$xx - mesh$xx[1])*alp
mesh$nz <- 1e3
mesh$zz <- seq(0, -max(mesh$H[,1]), len=mesh$nz)
mesh$mask_rho <- array(0, c(mesh$nx, 1, mesh$nz))

for (i in 1:mesh$nx) {
	j0 <- which.min (abs(mesh$H[i,1] + mesh$zz))
	mesh$mask_rho[i,1,1:j0] <- 1
}

sol_wunsch <- analyt_sol(mesh, n=1, N2=N2[1], gamma=gamma)

tm <- 1
num_ts <- 3
seq_x <- i0:i1

seq_ts <- seq(i0+30, i1-30, len=num_ts)
dx <- diff(mesh$xx[seq_ts])[1]
dx <- dx*0.5
scale <- dx/max(abs(sol_wunsch$uu[seq_ts, 1, , tm]), na.rm=T)

#for (tm in 1:mesh$nt) {
#	png(paste0(tm, ".png"), width=600, height=250)
	pdf(paste0("gamma_d7", ".pdf"), width=10, height=4)
	plot(mesh0$xx[seq_x]/1e3, -mesh$H[seq_x], type="l", xlab="x coordinate, km", ylab="depth, m", xlim=c(mesh$xx[i0], mesh$xx[i1])/1e3, ylim=c(-max(mesh$H[seq_x]), 0))

	for (i in seq_ts) {
		j0 <- which.min(abs(mesh$H[i] + mesh$zz))
		lines((mesh$xx[i] + sol_wunsch$uu[i,1,1:j0,tm]*scale)/1e3, mesh$zz[1:j0], lwd=2)
		lines(rep(mesh$xx[i], j0)/1e3, mesh$zz[1:j0], lty=2)
		lines(c(mesh$xx[i]/1e3, (mesh$xx[i] + sol_wunsch$uu[i,1,1,tm]*scale)/1e3), c(mesh$zz[1], mesh$zz[1]), lty=2)
		lines(c(mesh$xx[i]/1e3, (mesh$xx[i] + sol_wunsch$uu[i,1,j0,tm]*scale)/1e3), c(mesh$zz[j0], mesh$zz[j0]), lty=2)
	}
		
	tmp <- which( sol_wunsch$uu[seq_ts, 1, , tm] == max(sol_wunsch$uu[seq_ts, 1, , tm], na.rm=T), arr.ind=T)
	i_a <- seq_ts[tmp[1]]; j_a <- tmp[2]
	j0 <- which.min(abs(mesh$H[i_a] + mesh$zz))
	arrows( mesh$xx[i_a]/1e3, mesh$zz[j_a], (mesh$xx[i_a] + sol_wunsch$uu[i_a,1,j_a,tm]*scale)/1e3, mesh$zz[j_a], lty=1, lwd=2, angle=20, len=0.15)
	text((mesh$xx[i_a] + sol_wunsch$uu[i_a,1,j_a,tm]*scale)/1e3, mesh$zz[j_a], pos=3, labels=sprintf("u = %.1f cm/s", sol_wunsch$uu[i_a,1,j_a,tm]*1e3), lwd=2, cex=1.4, col="darkred", font=2)
	title(main=expression(gamma==0.7), cex.main=1.4)
	dev.off()

sol_wunsch_d2 <- analyt_sol(mesh0, n=1, N2=N2[1], gamma=0.2)
sol_wunsch_d7 <- analyt_sol(mesh, n=1, N2=N2[1], gamma=0.7)

j_d2 <- which.min( abs(mesh0$H[i0,1]/2 + mesh0$zz) )
j_d7 <- which.min( abs(mesh$H[i0,1]/2 + mesh$zz) )

pdf(paste0("middepth_vel.pdf"), width=8, height=4)
two_at_one(mesh$tM2, sol_wunsch_d2$uu[i1, 1, j_d2, ], mesh$tM2, sol_wunsch_d7$uu[i1, 1, j_d7, ], xlab=expression(M[2]~period), ylab="middepth velocity, m/s")
abline(h=0, lty=2)
legend("bottomleft", legend=c(expression(gamma==0.2), expression(gamma==0.7)), col=c("black", "darkred"), lwd=2, bg="white")
dev.off()
