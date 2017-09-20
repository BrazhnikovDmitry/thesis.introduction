analyt_sol <- function(mesh, n, N2, om=2*pi/12.4206/3600, gamma=0.5, x0=1000)
{
	c <- sqrt(om^2/(N2 - om^2))
	gamma <- tan(gamma*c)
	###
	# n - mode number
	delta <- (c + gamma)/(c - gamma)
#	print(delta)
	q <- 2*pi*n/log(delta)
#	x0 <- mesh$xx[1]
	tt <- mesh$t

	### Gerkema comment, eq. 8
	u_c <- array(0, dim=c(mesh$nx, mesh$nz))
	u_s <- array(0, dim=c(mesh$nx, mesh$nz))

	uu <- array(0, dim=c(mesh$nx, 1, mesh$nz, mesh$nt))

	ii <- complex(real=0, imaginary=1)
	
	for (i in 1:mesh$nx) {
		if (mesh$H[i] == 0) next
		for (j in 1:mesh$nz) {
			if (mesh$mask_rho[i,1,j] == 0) next

			xi_m <- (mesh$xx[i]-x0)*c - mesh$zz[j]
			xi_p <- (mesh$xx[i]-x0)*c + mesh$zz[j]

			u_c[i,j] <-  q/xi_m * cos(q*log(xi_m)) + q/xi_p * cos(q*log(xi_p))
			u_s[i,j] <- -q/xi_m * sin(q*log(xi_m)) - q/xi_p * sin(q*log(xi_p))
			print(u_c[i,j])
			uu[i,1,j,] <- u_c[i,j]*cos(om*tt) + u_s[i,j]*sin(om*tt)
		}
	}
	
#	return (list(u_c=u_c, u_s=u_s))
	return (list(uu=uu, u_c=u_c, u_s=u_s))
}

analyt_super <- function(mesh, n, N2, om=2*pi/12.4206/3600, gamma=0.5, x0=1000)
{
	c <- sqrt(om^2/(N2 - om^2))
	gamma <- tan(gamma*c)
	###
	# n - mode number
	delta <- (c + gamma)/(c - gamma)
	delta_p <- abs( (c - gamma)/(c + gamma) )

	eta <- 2*n*pi*log(delta_p)/(log(delta_p)^2 + pi^2)
	del <- -2*n*pi^2/(log(delta_p)^2 + pi^2)
	print(paste(eta, delta_p, del))
	tt <- mesh$t

	### Gerkema comment, eq. 8
#	u_c <- array(0, dim=c(mesh$nx, mesh$nz))
#	u_s <- array(0, dim=c(mesh$nx, mesh$nz))

	uu <- array(0, dim=c(mesh$nx, 1, mesh$nz, mesh$nt))

	ii <- complex(real=0, imaginary=1)
	
	for (i in 1:mesh$nx) {
		if (mesh$H[i] == 0) next
		for (j in 1:mesh$nz) {
			if (mesh$mask_rho[i,1,j] == 0) next

			xi_m <- (mesh$xx[i]-x0)*c - mesh$zz[j]
			xi_p <- (mesh$xx[i]-x0)*c + mesh$zz[j]
			if (xi_p < 0) {
				ll <- log(abs(xi_p)) + ii*pi
				pp <- ii*abs(xi_p)^(del-1)
			} else {
				ll <- log(xi_p)
				pp <- (xi_p)^(del-1)
			}
#print(xi_p)
			u1 <- ii*(eta+ii*del)*(xi_m)^(del-1)*exp(-ii*eta*log(xi_m))
			u2 <- ii*(eta+ii*del)*pp*exp(-ii*eta*ll)

#			print (exp(-ii*eta*log(xi_p)))
#			print (xi_p)
#			print(u_c[i,j])
			uu[i,1,j,] <- -u1-u2
		}
	}
	
#	return (list(u_c=u_c, u_s=u_s))
	return (list(uu=uu))
}
#sol_wunsch <- analyt_sol(mesh, n=1, N2=N2, gamma=-0.16)
