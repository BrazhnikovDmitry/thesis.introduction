source("/home/dmitry/Work/toolkit_it/misc.r")

#######################################
### Error, measured as mode 2 amplitude
#mesh0 <- mesh
load("final_w_sub_1.dat")
mesh0 <- mesh

source("sol_wunsch.r")
mesh <- mesh0

gamma <- -seq(0, 1, len=101)
err <- array(0, dim=c(dim(lmodes$pmodes)[4]+1, length(gamma)))
amp_md <- array(0, dim=c(mesh$nt, length(gamma)))

kin_en <- double(length(gamma))
kin_en_m <- array(0, dim=c(length(gamma), dim(lmodes$pmodes)[4]+1))

ii <- 200

mesh$xx <- mesh0$xx[ii]
mesh$nx <- 1
mesh$mask_rho[1,1,] <- mesh0$mask_rho[ii,1,]
mesh$H[,1] <- mesh$H[ii,1]

#lmodes$pmodes[1,1,,] <- lmodes$pmodes[ii,1,,]
#lmodes$norms[1,1,] <- lmodes$norms[ii,1,]

for (l in 2:(length(gamma)-1)) {
	print(l)

	if (l != 1) {
		mesh$t <- seq(0, 2*pi, len=mesh$nt)*12.42*3600
		c <- sqrt(om^2/(N2 - om^2))
		alp <- tan(gamma[l]*c)[1]
		mesh$H[1,1] <- -(mesh$xx[1] - mesh0$xx[1])*alp
#		mesh$zz[] <- 1
		mesh$zz <- seq(0, -mesh$H[1,1], len=mesh$nz)
		mesh$mask_rho[1,1,] <- 1
		
		lmodes <- eigen_modes_loc(mesh, N2, smooth=F, derivs=F)
		lmodes$pmodes[1,1,,] <- aperm(apply(lmodes$pmodes[1,1,,], 1, "/", apply(lmodes$pmodes[1,1,,], 2, max)), c(2,1))
	}

	sol_wunsch <- analyt_sol(mesh, n=1, N2=N2[1], gamma=gamma[l])
	t <- mesh$t; nt <- mesh$nt

	A <- array(0, dim=c(length(t), 4))
	A[,1] <- 1
	A[,2] <- t
	A[,3] <- cos(om*t)
	A[,4] <- sin(om*t)
	
	### Perform fit
	u_a <- apply(sol_wunsch$uu, c(1,2,3), smp_comp_regr, A=A)

	u_an <- array(0, c(mesh$nx, mesh$ny, dim(lmodes$pmodes)[4]))
#v_an <- array(0, c(mesh$nx, mesh$ny, dim(lmodes$pmodes)[4]))

	t_f$u_f <- u_a

	ii <- complex(real=0, imaginary=1)
	#t <- mesh$t[nt_seq]
	for (i in 1:mesh$nx) {
		print(i)
		for (j in 1:mesh$ny) {
			tmp <- (Re(t_f$u_f[i,j,]) + Im(t_f$u_f[i,j,])*ii)/max(sol_wunsch$uu[1,1,,])
			u_an[i,j,] <- fit_em_all(Re(tmp), lmodes$pmodes[i,j,,])$coefs + ii*fit_em_all(Im(tmp), lmodes$pmodes[i,j,,])$coefs#/max(abs(tmp))
	#		tmp <- Re(t_f$v_f[i,j,]) + Im(t_f$v_f[i,j,])*ii
	#		v_mod[i,j,] <- fit_em_all(Re(tmp), lmodes$pmodes[i,j,,])$coefs + ii*fit_em_all(Im(tmp), lmodes$pmodes[i,j,,])$coefs#/max
		}
	}

	
#	pdf(paste0("cmp_fits_wunsch_",l,".pdf"))
#	par(mfrow=c(2,2))
#	plot(Re(t_f$u_f)[1,1,], mesh$zz, type="l", ylab="depth, m", xlab="amp of cosine")
#	title(main=paste("gamma=", gamma[l]))
#	plot(Im(t_f$u_f)[1,1,], mesh$zz, type="l", ylab="depth, m", xlab="amp of sine")
#	image.plot(mesh$t, -mesh$zz, t(Re(sol_wunsch$uu)[1,1,,]), xlab="time", ylab="depth, m")
#	title(main="vertical structure")
#	plot(mesh$t, sol_wunsch$uu[1,1,200,]/max(sol_wunsch$uu[1,1,,]), type="l", xlab="time", ylab="amp")
#	title(main="normalized velocity at mid-depth")
#	dev.off()

#	ii <- 1
#	png(paste0("figs_err/decomp_wunsch_",l,".png"))
#	plot(abs((t_f$u_f[i,j,]))/max( abs((t_f$u_f[i,j,])), na.rm=T ), mesh$zz, type="l", xlab="Velocity, m/s", ylab="depth, m", lwd=2)
#	for (m in 2:5) lines(Re(u_an[ii,1,m])*lmodes$pmodes[ii,1,,m], mesh$zz, col="navy", lwd=2)
#	legend("topright", legend=c("Wunsch solution", "Flat modes"), col=c("black", "navy"), lwd=2)
#	dev.off()

#	pdf(paste0("fit_",l,".pdf"))
#	plot(Re(tmp), mesh$zz, type="l", lwd=2, xlab="Amplitude of semidiurnal fit", ylab="depth, m")
#	cc <- 0
#	for (k in 2:6) {
#		cc <- cc+lmodes$pmodes[1,1,,k]*Re(u_an[1,1,k])
#		lines(lmodes$pmodes[1,1,,k]*Re(u_an[1,1,k]), mesh$zz)
#	}
#	lines(cc, mesh$zz, col="red", lwd=2)
#	dev.off()

#	err[,l] <- u_an[1,1,]#*lmodes$norms[1,1,]
#	amp_md[,l] <- sol_wunsch$uu[1,1,200,]/max(sol_wunsch$uu[1,1,,])
	kin_en[l] <- sum( abs(u_a)^2)/max(sol_wunsch$uu[1,1,,])^2
	kin_en_m[l,1] <- abs(u_an[1,1,2])^2*sum(lmodes$pmodes[1,1,,2]^2)
	for (k in 3:6) kin_en_m[l,k-1] <- kin_en_m[l,k-2] + sum( (abs(u_an[1,1,k])*lmodes$pmodes[1,1,,k])^2 )
}

seq <- c(21, 71)
COLS <- tim.colors(6)

for (l in seq) {
	c <- sqrt(om^2/(N2 - om^2))
	alp <- tan(gamma[l]*c)[1]
	mesh$H[1,1] <- -(mesh$xx[1] - mesh0$xx[1])*alp
	mesh$zz <- seq(0, -mesh$H[1,1], len=mesh$nz)
	mesh$mask_rho[1,1,] <- 1
		
	lmodes <- eigen_modes_loc(mesh, N2, smooth=F, derivs=F)
	lmodes$pmodes[1,1,,] <- aperm(apply(lmodes$pmodes[1,1,,], 1, "/", apply(lmodes$pmodes[1,1,,], 2, max)), c(2,1))
	sol_wunsch <- analyt_sol(mesh, n=1, N2=N2[1], gamma=gamma[l])
	t <- mesh$t; nt <- mesh$nt

	A <- array(0, dim=c(length(t), 4))
	A[,1] <- 1
	A[,2] <- t
	A[,3] <- cos(om*t)
	A[,4] <- sin(om*t)
	
	### Perform fit
	u_a <- apply(sol_wunsch$uu, c(1,2,3), smp_comp_regr, A=A)
	u_an <- array(0, c(mesh$nx, mesh$ny, dim(lmodes$pmodes)[4]))

	t_f$u_f <- u_a

	ii <- complex(real=0, imaginary=1)

	for (i in 1:mesh$nx) {
		print(i)
		for (j in 1:mesh$ny) {
			tmp <- (Re(t_f$u_f[i,j,]) + Im(t_f$u_f[i,j,])*ii)/max(sol_wunsch$uu[1,1,,])
			u_an[i,j,] <- fit_em_all(Re(tmp), lmodes$pmodes[i,j,,])$coefs + ii*fit_em_all(Im(tmp), lmodes$pmodes[i,j,,])$coefs#/max(abs(tmp))
		}
	}
		
	pdf(paste0("fit_",l,".pdf"), height=7, width=5)
	plot(Re(tmp), mesh$zz, type="l", lwd=3, xlab="Amplitude of semidiurnal fit", ylab="depth, m")
	cc <- 0
	for (k in 2:6) {
		cc <- cc+lmodes$pmodes[1,1,,k]*Re(u_an[1,1,k])
		lines(lmodes$pmodes[1,1,,k]*Re(u_an[1,1,k]), mesh$zz, col=COLS[k], lwd=2, lty=2)
	}
	lines(cc, mesh$zz, col="orange", lwd=2)
	if (l != 71) legend("bottomleft", lwd=2, col=c("black", "orange", COLS[2:6]), legend=c("analytical", "total fit", paste("mode ", 1:5)))
	abline(v=0, lty=2)
	if (l == 21) title(main=expression(gamma==0.2)) else title(main=expression(gamma==0.7))
	dev.off()
}

seq_g <- 2:100
pdf("errors1.pdf", width=8, height=4)
plot(-gamma[seq_g], kin_en[seq_g], type="l", xlab=expression(list(slope~angle,gamma)), ylab="kinetic energy", lwd=2)
lines(-gamma[seq_g], kin_en_m[seq_g,1], type="l", col="navy", lwd=2)
lines(-gamma[seq_g], kin_en_m[seq_g,2], type="l", col="darkred", lwd=2)
legend("bottomleft", lwd=2, col=c("black", "navy", "darkred"), legend=c("analytical", "mode 1", "mode 1 and 2"))
dev.off()

pdf("errors2.pdf", width=8, height=4)
plot(-gamma[seq_g], abs(kin_en[seq_g] - kin_en_m[seq_g,1])/kin_en[seq_g]*100, type="l", xlab=expression(list(slope~angle,gamma)), ylab="error, %", col="navy", lwd=2, ylim=c(0, 0.3)*100)
lines(-gamma[seq_g], abs(kin_en[seq_g] - kin_en_m[seq_g,2])/kin_en[seq_g]*100, type="l", col="darkred", lwd=2)
axis(1, at=0.35, label=0.35)
legend("topleft", lwd=2, col=c("navy", "darkred"), legend=c("mode 1", "mode 1 and 2"))
abline(v=c(0.35), lty=2)
abline(h=5, lty=2)
dev.off()

