rm(list = ls())
library(parallel)
library(glmnet)

postlasso.ols.fn = function(y, X.scale, CV.list)
{
	n = nrow(X.scale)
	p = ncol(X.scale)
	lam.search.logic = CV.list[[1]] # = 1 emphasises that CV based search is required #
	CV.k = CV.list[[2]] # the k-fold parameter #

	if(lam.search.logic == 1) # do a CV based search for lambda for the initial Lasso estimator #
	{
		# we are using scaled and centered columns in design matrix
		# no intercept term is included in the model fitting stage as the same is not done in data-generation #
		# the output from A is equal to lars(..., normalize = F, intercept = F) with lambda.lars = n*lambda.glmnet #

		A = cv.glmnet(X.scale, y, nfolds = CV.k, intercept = F, standardize = F)
		CV.lam = A$lambda.min

		# by default coef.est.init contains an intercept term (=0) which is dropped #
		beta.init.lasso = as.numeric(coef(A, s = CV.lam)[-1])

        A.hat = (1:p)[beta.init.lasso!=0]
       
        beta.hat = rep(0, p)
        beta.hat[A.hat] = solve(t(X.scale[ , A.hat])%*%X.scale[ ,A.hat])%*%t(X.scale[ ,A.hat])%*%y
		out.all = list(beta.hat, CV.lam)	
	}
	else # if CV.lam is already supplied #
	{
		CV.lam = CV.list[[3]]

		A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = c(2,1)*CV.lam)
		beta.init.lasso = as.numeric(A1$beta[,2])

        A.hat = (1:p)[beta.init.lasso!=0]
        
        beta.hat = rep(0, p)
        beta.hat[A.hat] = solve(t(X.scale[ , A.hat])%*%X.scale[ ,A.hat])%*%t(X.scale[ ,A.hat])%*%y
    	out.all = list(beta.hat, CV.lam)	
	}
	return(out.all)
}

postlasso.ols.within.boot <- function(b, y.boot.mat, X, CV.list)
{
    beta.star = as.vector(postlasso.ols.fn(y.boot.mat[ ,b], X, CV.list)[[1]]) # only return beta.star
	return(beta.star) 
}

c1.fn = function(b, beta.s.mat, X.mat, d.vec, y.star.mat, beta.hat)
{
	n = nrow(X.mat)
	p = ncol(X.mat)

    beta.star = beta.s.mat[ ,b]
    y.star = y.star.mat[ ,b]
    e.star = y.star - X.mat%*%beta.star
    sigma.star.sqr = mean(e.star^2)
    
    T.star = sqrt(n)*sum(d.vec*(beta.star - beta.hat))
    R.star = T.star/sqrt(sigma.star.sqr)
    return(abs(R.star))    
}

c2.fn = function(b, beta.ss.mat, X, y, CV.lam, d.vec, G.mat, mu.G, beta.hat)
{
    n = nrow(X)
	p = ncol(X)

	beta.ss = beta.ss.mat[ ,b]
	G.vec = G.mat[ ,b]
	wt.vec.mod = ((G.vec - mu.G)^2)/(mu.G^2)

	Tn.ss = sqrt(n)*sum((d.vec)*as.numeric(beta.ss - beta.hat))
    res.ss = y - X%*%beta.ss
    sigma.ss.sqr = mean(wt.vec.mod*(res.ss^2))

    return(Tn.ss/sqrt(sigma.ss.sqr))	
}
rbci.fn = function(y, X, CV.k, coef.index, B, alpha, beta.0, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]]
	CV.lam = init.info.list[[2]]
	
	theta.hat = sum(d.vec*beta.hat)
	theta.0 = sum(d.vec*beta.0)

	res.vec = y - X%*%beta.hat
    sigma.check.sqr = mean(res.vec^2)
	c.res.vec = res.vec - mean(res.vec)

	# RB part #
	rboot.imat = sapply(1:B, function(i, z){sample(z, size = length(z), replace = T)}, z = 1:n) # n x B matrix of sampled indices
	y.boot.rb = as.vector(X%*%beta.hat) + sapply(1:B, function(i, b1.mat, c.res){c.res[b1.mat[ ,i]]}, b1.mat = rboot.imat, c.res = c.res.vec)
    
	beta.star.mat = sapply(1:B, postlasso.ols.within.boot, y.boot.rb, X, list(0, CV.k, CV.lam)) # pxB matrix 
    R.star.dot.2.abs = sapply(1:B, c1.fn, beta.star.mat, X, d.vec, y.boot.rb, beta.hat)
	quant.rb = quantile(R.star.dot.2.abs, probs = 1-alpha)

	ci.rb = rep(theta.hat, 2) + c(-1, 1)*quant.rb*sqrt(sigma.check.sqr)/sqrt(n)
	logic.rb = ifelse(ci.rb[1]<= theta.0 & theta.0 <= ci.rb[2], 1, 0)
	length.rb = abs(ci.rb[2] - ci.rb[1])

	rb.ci.out = cbind(c(logic.rb, length.rb), ci.rb)
	return(rb.ci.out)
}

pbci.fn = function(y, X, CV.k, coef.index, B, alpha, beta.0, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]]
	CV.lam = init.info.list[[2]]
	
    theta.hat = sum(d.vec*beta.hat)
	theta.0 = sum(d.vec*beta.0)

	res.vec = y - X%*%beta.hat
	A.hat = (1:p)[beta.hat != 0]
	d.hat = d.vec[A.hat]

	C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
    sigma.check.sqr = mean(res.vec^2)
	Sigma.hat = matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%matrix(d.hat, ncol = 1)
	X.hat.mod = t(sapply(1:nrow(X), function(i, a, v){a[i, ]*v[i]}, a = X[ ,A.hat], v = as.numeric(res.vec))) # n x p0.hat mat
	Sigma.tilde = as.numeric(matrix(d.hat,nrow=1)%*%solve(C11.hat)%*%(t(X.hat.mod)%*%X.hat.mod/n)%*%solve(C11.hat)%*%matrix(d.hat,ncol=1))
	
	term.1 = sqrt(sigma.check.sqr)*(Sigma.tilde^(-1/2))  # to be used in the pb.pivot construction #
    term.3 = sqrt(sigma.check.sqr)*sqrt(Sigma.hat)/sqrt(n) # to be used in PB CI construction #

	# computing w2 and w4 #
	u1 = as.numeric(matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%t(X[ ,A.hat])) # n-vector
	w2 = -(1/sigma.check.sqr)*(1/Sigma.tilde)*mean((u1^2)*(res.vec^4)) + (1/sigma.check.sqr^2)*mean(res.vec^4)
	w4 = (2/(Sigma.tilde^2))*mean((u1^4)*(res.vec^4)) + 4*(1/sigma.check.sqr)*(1/Sigma.tilde)*mean((u1^2)*(res.vec^4)) - 3*(1/(sigma.check.sqr^2))*mean(res.vec^4) + 1
	z.alpha = qnorm(1 - alpha/2)
	C.np = -(z.alpha/n)*(w2/2 + (w4/24)*(z.alpha^2 - 3))

	# PB part #
	mu.star = (1/2)/(1/2 + 3/2)
	G.mat = sapply(1:B, function(i, n.size){rbeta(n.size, 1/2, 3/2)}, n.size = n) # n x B matrix 
	pboot.imat = apply(G.mat, 2, function(v2, c1){(1/c1)*(v2-c1)}, c1 = mu.star) # n x B transformed from G.mat#
	z.mat = as.vector(X%*%beta.hat) + sapply(1:B, function(b, u.mat, r.vec){r.vec*u.mat[ ,b]}, u.mat = pboot.imat, r.vec = res.vec) # nxB matrix
	
	beta.ss.mat = sapply(1:B, postlasso.ols.within.boot, z.mat, X, list(0, CV.k, CV.lam)) # pxB matrix of pboot.mcp estimates #
    term2.ss = sapply(1:B, c2.fn, beta.ss.mat, X, y, CV.lam, d.vec, G.mat, mu.star, beta.hat) # B-vector
	Rn.ss = as.numeric(term2.ss)*term.1 # PB pivot # 
	quant.absRn.ss = quantile(abs(Rn.ss), probs = 1-alpha)
	hn.tilde.alpha = quant.absRn.ss + C.np

	pb.ci = rep(theta.hat, 2) + c(-1,1)*as.numeric(term.3*hn.tilde.alpha)
	len.pb.ci = abs(pb.ci[2] - pb.ci[1])
	logic.pb.ci = ifelse(pb.ci[1]<= theta.0 & theta.0 <= pb.ci[2], 1, 0)
	pb.ci.out = cbind(c(logic.pb.ci, len.pb.ci), pb.ci)
    return(pb.ci.out)
}

data.set.m = function(m, y.mat, X.scale, CV.k, coef.index, B, alpha, beta.0)
{
    y = y.mat[ ,m]
	n = nrow(X.scale)
	p = ncol(X.scale)
	p0 = length(beta.0[beta.0!=0])
	
	d.vec = numeric(p)
	d.vec[coef.index] <- 1

    C = postlasso.ols.fn(y, X.scale, list(1, CV.k, -1)) # CV based search required # third argument is dummy value # 
	beta.hat = as.vector(C[[1]])
    CV.lam = C[[2]]

    beta.ols = as.vector(solve(t(X.scale[ ,1:p0])%*%X.scale[ ,1:p0])%*%t(X.scale[ ,1:p0])%*%y)
    theta.hat = sum(d.vec*beta.hat)
    
	theta.0 = sum(d.vec*beta.0)
	
    beta.ols <- c(beta.ols, rep(0, p-p0))
    
    out.list = list()
	if(theta.hat!=0)
	{
		init.info.list = list(beta.hat, CV.lam)
        
		# RB #
        rb.ci.out = rbci.fn(y, X.scale, CV.k, coef.index, B, alpha, beta.0, init.info.list)
		
		# oracle ci #
    	C11 = t(X.scale[ ,1:p0])%*%(X.scale[ ,1:p0])/n
		C11.inv = solve(C11)
		d1.vec = d.vec[1:p0] # p0 dimensional vector #
		rho.sqr = as.numeric(matrix(d1.vec, nrow = 1)%*%C11.inv%*%matrix(d1.vec, ncol = 1))
		res.vec = y - X.scale%*%beta.hat
		sigma.check = sqrt(mean(res.vec^2))
		oracle.ci = rep(theta.hat, 2) - qnorm(c(1-alpha/2, alpha/2))*sigma.check*sqrt(rho.sqr)/sqrt(n)
		
		logic.oracle = ifelse(oracle.ci[1]<= theta.0 & theta.0 <= oracle.ci[2], 1, 0)
		length.oracle = abs(oracle.ci[2] - oracle.ci[1])
		oracle.ci.out = cbind(c(logic.oracle, length.oracle), oracle.ci)

		# PB #
        pb.ci.out = pbci.fn(y, X.scale, CV.k, coef.index, B, alpha, beta.0, init.info.list)
        # -------------------------------------------------------------- #
		ci.mat <- cbind(rb.ci.out, pb.ci.out, oracle.ci.out) #2 x 6 matrix 
		        
        if(m%%100 == 0)
		{
            cat("m = ",m,"\n")
		    #print(ci.mat)
        }
	}
	else
	{
		ci.mat <- -1
	}
    #print(c(theta.hat, CV.lam))
    #return(cbind(beta.hat, beta.ols))
	return(list(ci.mat, theta.hat, CV.lam, cbind(beta.hat, beta.ols)))
}


postlasso.ols.all <- function(n, p, CV.k, coef.index, alpha, B, err.type)
{
	my.file.name = paste("yx-all-n-",n,"-p-",p,".Rdata", sep = "")
	load(my.file.name)
	
	y.mat = yx.list[[1]][err.type, , ] # loads n x M response matrix 
	X.fix = yx.list[[2]]
	beta.0 = yx.list[[3]]
	snr.vec = yx.list[[4]]
    M = ncol(y.mat)
    
	output.super <- mclapply(1:M, data.set.m, mc.cores = 8, y.mat, X.fix, CV.k, coef.index, B, alpha, beta.0)
    beta.array = array(0, dim = c(M, p, 2))
	cov.count <- av.len <- rep(0, 3)
    theta.hat <- CV.lam <- numeric(M)
    int.count = 0
    
	for(m in 1:M)
	{
	    u.m = output.super[[m]][[1]]
        theta.hat[m] = output.super[[m]][[2]]
        CV.lam[m] = output.super[[m]][[3]]
		if(is.matrix(u.m))
		{
			int.count = int.count + 1
			cov.count[1] = cov.count[1] + u.m[1, 1] #RB
			cov.count[2] = cov.count[2] + u.m[1, 3] #PB
			cov.count[3] = cov.count[3] + u.m[1, 5] #Oracle

			av.len[1] = av.len[1] + u.m[2, 1] # RB
			av.len[2] = av.len[2] + u.m[2, 3] # PB
			av.len[3] = av.len[3] + u.m[2, 5] # Oracle
 		}

        beta.array[m, , ] = output.super[[m]][[4]]
	}
    #par(mfrow = c(2,5))
    #plot(beta.array[,1,1], beta.array[,1,2], xlab = "pl.ols", ylab = "ols", type = "p", pch = .75, col = 4)
    #abline(a = 0, b = 1, col = 1)
    #for(j in 2:10)
    #{
    #    plot(beta.array[,j,1], beta.array[,j,2], xlab = "pl.ols", ylab = "ols", type = "p", pch = .75, col = 4)
    #    abline(a = 0, b = 1, col = 1)
    #}
    #prop.total = numeric(490)
    #for(j in 1:490)
    #{
    #    prop.total[j] = length((1:M)[beta.array[ ,10+j, 1]==0])/M
    #}
    #print(summary(prop.total))
    cat("emp.cov = ", signif(cov.count/int.count,3),"\n")
    cat("av.len =", signif(av.len/int.count,3),"\n")
    cat("int.count = ", int.count,"\n")
    cat("err.type = ", err.type,"\n")
    cat("coef.index = ", coef.index, "\n")
    cat("n=",n," p=",p,"\n")   
	my.file.name.new = paste("out-plols-n-",n,"-p-",p,"-coef-beta-",coef.index,"-errtype-",err.type,"-95.Rdata", sep = "")
	save(output.super, file = my.file.name.new)

	#output.super <- lapply(1:M, data.set.m, y.mat, X.fix, CV.k, coef.index, B, alpha, beta.0)
	#return(beta.array)	
}

# uncomment as per requirement
# for n = 25 case, CV.k = 3 can be tried
# for (n,p) = (25,90), we used p_0 = 3
# for (n,p) = (50,90), we used p_0 = 4
# for p = 500 cases, we used p_0 = 10

#postlasso.ols.all(n = 25, p = 90, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 25, p = 90, CV.k = 5, coef.index = 3, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 25, p = 90, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 2)
#postlasso.ols.all(n = 25, p = 90, CV.k = 5, coef.index = 3, alpha = 0.1, B = 500, err.type = 2)

#postlasso.ols.all(n = 50, p = 90, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 50, p = 90, CV.k = 5, coef.index = 4, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 50, p = 90, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 2)
#postlasso.ols.all(n = 50, p = 90, CV.k = 5, coef.index = 4, alpha = 0.1, B = 500, err.type = 2)

#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 2)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B = 500, err.type = 2)

#postlasso.ols.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B = 500, err.type = 1)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B = 500, err.type = 2)
#postlasso.ols.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B = 500, err.type = 2)
