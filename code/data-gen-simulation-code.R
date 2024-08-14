# Generate X data 
rm(list = ls())
library(mvtnorm)

error.fn = function(i, n.size, err.type)
{
	if(err.type == 1)
	{
		err.vec = rnorm(n.size, sd = 1)
	}
	else
	{
		err.vec = (rchisq(n.size, df = 1)- 1)/sqrt(2)
	}
	return(err.vec) 
}

yx.gen = function(n, beta.0.true, rho, M)
{
	p = length(beta.0.true)
	# covariance matrix from which x are generated #
	Sigma = matrix(0, nrow = p, ncol = p)
	for(i in 1:p)
	{
		Sigma[i, ] = (rho)^(abs(1:p - i))	
	}

	## X matrix generation # initial matrix is changed to get column.mean = 0, column.norm = 1 #
	X.prelim.1 = rmvnorm(n = n, mean = rep(0, p), sigma = Sigma) 
	X.prelim.2 = apply(X.prelim.1, 2, function(u){u-mean(u)}) # X.fix is not changed over simulation # col. means = 0 #
	X.fix = apply(X.prelim.2, 2, function(u){u/sqrt(sum(u^2))}) # columns have norm = 1 #
	
	signal.true = X.fix%*%beta.0.true	
	
	# snr calc #
	snr.norm = sqrt(sum(signal.true^2))/(sqrt(n)*1)
	snr.gamma = sqrt(sum(signal.true^2))/(sqrt(n)*(1))
	print(c(snr.norm, snr.gamma))
	
	y.array = array(0, dim = c(2, n, M))
	for(j in 1:2)
	{
		y.array[j, , ] = as.vector(signal.true) + sapply(1:M, error.fn, n.size = n, err.type = j) 
	}
	yx.list = list()
	yx.list[[1]] <- y.array
	yx.list[[2]] <- X.fix
	yx.list[[3]] <- beta.0.true
	yx.list[[4]] <- c(snr.norm, snr.gamma) 
	
	my.file.name = paste("yx-all-n-",n,"-p-",p,".Rdata", sep = "")
	save(yx.list, file = my.file.name)
}
#beta.0.true = c(rep(5,3), rep(0,87)) 
#yx.gen(n = 25, beta.0.true, rho = 0.6, M = 400)
