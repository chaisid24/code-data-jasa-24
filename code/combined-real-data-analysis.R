rm(list = ls())
library(parallel)
library(glmnet)
source("alasso-real-data.R")
source("post-lasso-ols-real-data.R")

# Running CV once # finding one "good" set of variables #
good.set.m = function(m, y, X, CV.k)
{
	n = nrow(X)
	p = ncol(X)
	    
    A = cv.glmnet(X, y, nfolds = CV.k, intercept = F, standardize = F)
	CV.lam = A$lambda.min
    
    beta.init.lasso = as.numeric(coef(A, s = CV.lam)[-1])
    A.hat = (1:p)[beta.init.lasso!=0]
    beta.plols = rep(0, p)
    beta.plols[A.hat] = as.vector(lm(y ~ -1+X[ ,A.hat])$coefficients)
    
    # alasso estimate # computed at the same CV.lam #
	wt.vec = ifelse(beta.init.lasso != 0, abs(beta.init.lasso), 0)
    B = glmnet(X, y, intercept = F, standardize = F, penalty.factor = 1/(wt.vec), lambda = c(2,1)*CV.lam)
	beta.alasso = as.numeric(B$beta[,2])
	A.hat.alasso = (1:p)[beta.alasso!=0]
    
    # good.set 
    A.good = intersect(A.hat, A.hat.alasso) # sends back an index.set which appears in both ALASSO and Lasso selected variables #
    return(list(A.good, CV.lam))	
}

good.set = function(rho.thresh, CV.k, M)
{
    # Initial processing # # Read in data
    x.main <- t(read.table("Ro131.csv", header=TRUE, sep=","))
    # Standardise x.main
    x.new <- apply(x.main, 2, scale)
   
    # Remove predictor gene, which is tightly related to Ro131 as done in Segal et al
    x.new <- x.new[,c(1:6077,6079:6320)]
    # Read in the response variable, entered manually
	y.main <- c(143,84,98,83,153,141,191,130,744,381,1047,806,621,849,475,966,708,487,1447,1693,1731,1025,376,126,102,149,91,153,235,68)
    y.new = (y.main - mean(y.main))/sd(y.main)

    # correlation based screening 
    cor.vals = cor(x.new, y.new)
    index.high = (1:ncol(x.new))[abs(cor.vals) >= rho.thresh]
    x.new.sub = x.new[ , index.high]
    datayx.list = list(y.new, x.new.sub)
    
    # running CV method M times #
    out.good.list = mclapply(1:M, good.set.m, mc.cores = 8, y.new, x.new.sub, CV.k)
    p = ncol(x.new.sub)
    count.vec = numeric(p)
    CV.lam = numeric(M)
    for(m in 1:M)
    {
        u.j = out.good.list[[m]][[1]]
        CV.lam[m] = out.good.list[[m]][[2]]
        count.vec[u.j] = count.vec[u.j] + 1
    }
    count.mat = cbind(1:p, count.vec, count.vec/M)
    good.set = (1:p)[count.vec/M >= 0.9]
    return(list(good.set, mean(CV.lam), count.mat, y.new, x.new.sub))
}

#U = good.set(rho.thresh = 0.6, CV.k = 5, M = 200)
alasso.plols.fn = function(y, X.scale, lam.fix)
{
	n = nrow(X.scale)
	p = ncol(X.scale)

	A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = c(2,1)*lam.fix)
	beta.init.lasso = as.numeric(A1$beta[,2])

    A.hat = (1:p)[beta.init.lasso!=0]
    beta.plols = rep(0, p)
    beta.plols[A.hat] = as.vector(lm(y ~ -1+X.scale[ ,A.hat])$coefficients)
   
	wt.vec = ifelse(beta.init.lasso != 0, abs(beta.init.lasso), 0) 
	B1 = glmnet(X.scale, y, intercept = F, standardize = F, penalty.factor = 1/wt.vec, lambda = c(2,1)*lam.fix)
	beta.alasso = as.numeric(B1$beta[,2])
	beta.all = cbind(beta.init.lasso, beta.alasso, beta.plols)	
	return(beta.all)
}

ci.all = function(rho.thresh, CV.k, M, B, alpha)
{
    # -------------------- #
    U.main = good.set(rho.thresh, CV.k, M)
    A.good = U.main[[1]]
    lam.fix = U.main[[2]]
    count.mat = U.main[[3]]
    y = U.main[[4]]
    X = U.main[[5]]
    cnames.vec = colnames(X)
    # -------------------------------------- #

	n = nrow(X)
	p = ncol(X)
	
	C = alasso.plols.fn(y, X, lam.fix)
	beta.init.lasso = C[ ,1]
	beta.alasso = C[ ,2]
	beta.plols = C[ ,3]
	
	A.hat.alasso = intersect((1:p)[beta.alasso!=0], A.good)
	A.hat.plols = intersect((1:p)[beta.plols!=0], A.good)
    
    ci.mat.alasso <- matrix(0, nrow = length(A.hat.alasso), ncol = (1+6+1))
    ci.mat.plols <- matrix(0, nrow = length(A.hat.plols), ncol = (1+6+1))
    colnames(ci.mat.alasso) <- colnames(ci.mat.plols) <- c("Var (j)","RB-L","RB-U","PB-L","PB-U","Or-L","Or-U","Est.")
	ci.mat.alasso[ ,1] <- A.hat.alasso
	ci.mat.alasso[ ,8] <- beta.alasso[A.hat.alasso]
	ci.mat.plols[ ,1] <- A.hat.plols
	ci.mat.plols[ ,8] <- beta.plols[A.hat.plols]
	
    # CI's for ALASSO case # RB, PB and Oracle based two sided CI's 
	for(j in 1:length(A.hat.alasso))
	{
	    d.vec = numeric(p)
		d.vec[A.hat.alasso[j]] = 1
		theta.hat = sum(d.vec*beta.alasso)
		init.info.list = list(beta.init.lasso, beta.alasso, lam.fix)
		
	    # oracle ci #
		C11 = t(X[ , A.hat.alasso])%*%(X[ , A.hat.alasso])/n
		C11.inv = solve(C11)
		rho.sqr = as.numeric(matrix(d.vec[A.hat.alasso], nrow = 1)%*%C11.inv%*%matrix(d.vec[A.hat.alasso], ncol = 1))
		res.vec = y - X%*%beta.alasso
		sigma.check = sqrt(mean(res.vec^2))
		ci.mat.alasso[j, 6:7] = signif(rep(theta.hat, 2) - qnorm(c(1-alpha/2, alpha/2))*sigma.check*sqrt(rho.sqr)/sqrt(n), 3)
	    	
		# RB #
		ci.mat.alasso[j, 2:3] = signif(rbci.fn(y, X, CV.k, A.hat.alasso[j], B, alpha, init.info.list), 3) #vector with two elements 
		# PB #
		ci.mat.alasso[j, 4:5] = signif(pbci.fn(y, X, CV.k, A.hat.alasso[j], B, alpha, init.info.list), 3)
	}
	
	# CI's for Post-Lasso OLS case # RB, PB and Oracle based two sided CI's 
	for(j in 1:length(A.hat.plols))
	{
	    d.vec = numeric(p)
		d.vec[A.hat.plols[j]] = 1
		theta.hat = sum(d.vec*beta.plols)
		init.info.list = list(beta.plols, lam.fix, 1) # k.fac = 1 #
		
	    # oracle ci #
		C11 = t(X[ , A.hat.plols])%*%(X[ , A.hat.plols])/n
		C11.inv = solve(C11)
		rho.sqr = as.numeric(matrix(d.vec[A.hat.plols], nrow = 1)%*%C11.inv%*%matrix(d.vec[A.hat.plols], ncol = 1))
		res.vec = y - X%*%beta.plols
		sigma.check = sqrt(mean(res.vec^2))
		ci.mat.plols[j, 6:7] = signif(rep(theta.hat, 2) - qnorm(c(1-alpha/2, alpha/2))*sigma.check*sqrt(rho.sqr)/sqrt(n), 3)
		
		# RB #
		ci.mat.plols[j, 2:3] = signif(rbci.fn.rd(y, X, CV.k, A.hat.plols[j], B, alpha, init.info.list), 3)
        # PB #
		ci.mat.plols[j, 4:5] = signif(pbci.fn.rd(y, X, CV.k, A.hat.plols[j], B, alpha, init.info.list), 3)
		
	}
	cat("ALASSO Based CIs","\n")
	print(ci.mat.alasso)
	print(cnames.vec[A.hat.alasso])
	
	cat("PLOLS Based CIs","\n")
	print(ci.mat.plols)
	print(cnames.vec[A.hat.plols])
	
	return(list(ci.mat.alasso, ci.mat.plols, lam.fix, A.good, count.mat, cnames.vec[A.hat.alasso], cnames.vec[A.hat.plols]))
}
A.main = ci.all(rho.thresh = 0.6, CV.k = 5, M = 200, B = 500, alpha = 0.1)
save(A.main, file = "real-data-output.Rdata")
