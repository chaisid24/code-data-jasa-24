#rm(list = ls())
#library(parallel)
#library(glmnet)

postlasso.ols.fn.rd = function(y, X.scale, CV.list)
{
	n = nrow(X.scale)
	p = ncol(X.scale)
	lam.search.logic = CV.list[[1]] # = 1 emphasises that CV based search is required #
	CV.k = CV.list[[2]] # the k-fold parameter #
    k.fac = CV.list[[4]]

	if(lam.search.logic == 1) # do a CV based search for lambda for the initial Lasso estimator #
	{
		# we are using scaled and centered columns in design matrix
		# no intercept term is included in the model fitting stage as the same is not done in data-generation #
		# the output from A is equal to lars(..., normalize = F, intercept = F) with lambda.lars = n*lambda.glmnet #

		A = cv.glmnet(X.scale, y, nfolds = CV.k, intercept = F, standardize = F)
		lam = k.fac*A$lambda.min

		# by default coef.est.init contains an intercept term (=0) which is dropped #
		beta.init.lasso = as.numeric(coef(A, s = lam)[-1])

        A.hat = (1:p)[beta.init.lasso!=0]
        beta.hat = rep(0, p)
        beta.hat[A.hat] = solve(t(X.scale[ , A.hat])%*%X.scale[ ,A.hat])%*%t(X.scale[ ,A.hat])%*%y
		out.all = list(beta.hat, lam)	
	}
	else # if lam is already supplied #
	{
		lam = CV.list[[3]]

		A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = c(2,1)*lam)
		beta.init.lasso = as.numeric(A1$beta[,2])
        
        A.hat = (1:p)[beta.init.lasso!=0]
        beta.hat = rep(0, p)
        beta.hat[A.hat] = as.vector(lm(y~-1+X.scale[ ,A.hat])$coefficients)
        out.all = list(beta.hat, lam)	
	}
	return(out.all)
}

postlasso.ols.within.boot.rd <- function(b, y.boot.mat, X, CV.list)
{
    beta.star = as.vector(postlasso.ols.fn.rd(y.boot.mat[ ,b], X, CV.list)[[1]]) # only return beta.star
	return(beta.star) 
}

c1.fn.rd = function(b, beta.s.mat, X.mat, d.vec, y.star.mat, beta.hat)
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

c2.fn.rd = function(b, beta.ss.mat, X, y, d.vec, G.mat, mu.G, beta.hat)
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

rbci.fn.rd = function(y, X, CV.k, coef.index, B, alpha, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]]
	lam = init.info.list[[2]]
	k.fac = init.info.list[[3]]
	
	theta.hat = sum(d.vec*beta.hat)
	
	res.vec = y - X%*%beta.hat
    sigma.check.sqr = mean(res.vec^2)
	c.res.vec = res.vec - mean(res.vec)

	# RB part #
	rboot.imat = sapply(1:B, function(i, z){sample(z, size = length(z), replace = T)}, z = 1:n) # n x B matrix of sampled indices
	y.boot.rb = as.vector(X%*%beta.hat) + sapply(1:B, function(i, b1.mat, c.res){c.res[b1.mat[ ,i]]}, b1.mat = rboot.imat, c.res = c.res.vec)
    
	beta.star.mat = sapply(1:B, postlasso.ols.within.boot.rd, y.boot.rb, X, list(0, CV.k, lam, k.fac)) # pxB matrix 
    R.star.dot.2.abs = sapply(1:B, c1.fn.rd, beta.star.mat, X, d.vec, y.boot.rb, beta.hat)
    #print(R.star.dot.2.abs)
	quant.rb = quantile(R.star.dot.2.abs, na.rm = T, probs = 1-alpha)

	ci.rb = rep(theta.hat, 2) + c(-1, 1)*quant.rb*sqrt(sigma.check.sqr)/sqrt(n)
	return(ci.rb)
}

pbci.fn.rd = function(y, X, CV.k, coef.index, B, alpha, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]]
	lam = init.info.list[[2]]
	k.fac = init.info.list[[3]]
	
    theta.hat = sum(d.vec*beta.hat)

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
	
	beta.ss.mat = sapply(1:B, postlasso.ols.within.boot.rd, z.mat, X, list(0, CV.k, lam, k.fac)) # pxB matrix
    term2.ss = sapply(1:B, c2.fn.rd, beta.ss.mat, X, y, d.vec, G.mat, mu.star, beta.hat) # B-vector
	Rn.ss = as.numeric(term2.ss)*term.1 # PB pivot # 
	quant.absRn.ss = quantile(abs(Rn.ss), na.rm = T, probs = 1-alpha)
	hn.tilde.alpha = quant.absRn.ss + C.np

	pb.ci = rep(theta.hat, 2) + c(-1,1)*as.numeric(term.3*hn.tilde.alpha)
	return(pb.ci)
}


