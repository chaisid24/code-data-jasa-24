# this code computes alasso based symmetric confidence intervals for observed non-zero coefficients
# using both residual bootstrap and perturbation bootstrap based approaches.
# it also computes alasso based confidence intervals using the oracle limit law

#rm(list = ls())
#library(parallel)
#library(glmnet)


# computes alasso estimates at an optimal lambda, which is obtained by k-fold CV
# if (previously found) lambda is supplied and alasso estimate at that lambda value is computed
alasso.fn = function(y, X.scale, CV.list)
{
	n = nrow(X.scale)
	p = ncol(X.scale)
	lam.search.logic = CV.list[[1]] 	## = 1 emphasises that CV based search is required 
	CV.k = CV.list[[2]] 			## the k-fold CV parameter 

	if(lam.search.logic == 1) # do a CV based search for lambda for the initial Lasso estimator #
	{
		# we are using scaled and centered columns in design matrix
		# no intercept term is included in the model fitting stage as the same is not done in data-generation #
		# the output from A is equal to lars(..., normalize = F, intercept = F) with lambda.lars = n*lambda.glmnet #

		A = cv.glmnet(X.scale, y, nfolds = CV.k, intercept = F, standardize = F)
		CV.lam = A$lambda.min

		# by default coef.est.init contains an intercept term (=0) which is dropped #
		beta.init.lasso = as.numeric(coef(A, s = CV.lam)[-1])
		
		# alasso estimate # computed at the same CV.lam #
		delta = 1/sqrt(n)
		wt.vec = ifelse(beta.init.lasso != 0, abs(beta.init.lasso), 0)

		B = glmnet(X.scale, y, intercept = F, standardize = F, penalty.factor = 1/(wt.vec), lambda = c(2,1)*CV.lam)
		beta.hat = as.numeric(B$beta[,2])
	
		beta.all = cbind(beta.init.lasso, beta.hat)
		out.all = list(beta.all, CV.lam)	
	}
	else # if CV.lam is already supplied #
	{
		CV.lam = CV.list[[3]]

		A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = c(2,1)*CV.lam)
		beta.init.lasso = as.numeric(A1$beta[,2])
		wt.vec = ifelse(beta.init.lasso != 0, abs(beta.init.lasso), 0) 
		
		# alasso estimate #
		B1 = glmnet(X.scale, y, intercept = F, standardize = F, penalty.factor = 1/wt.vec, lambda = c(2,1)*CV.lam)
		beta.hat = as.numeric(B1$beta[,2])
		
		beta.all = cbind(beta.init.lasso, beta.hat)
		out.all = list(beta.all, CV.lam)	
	}
	return(out.all)
}


# alasso estimates for bootstrapped data set
# b : index for the b-th bootstrap data set, b = 1,...,B
alasso.within.boot <- function(b, y.boot.mat, X, CV.related.list)
{
	beta.star.all = alasso.fn(y.boot.mat[ ,b], X, CV.related.list)[[1]] # a px2 matrix #
	return(beta.star.all) 
}

# intermediate step for computing symmetric bias adjusted bootstrap statistic
c1.fn = function(b, beta.s.init.mat, beta.s.mat, X.mat, lam, d.vec, y.star.mat, beta.hat)
{
	n = nrow(X.mat)
	p = ncol(X.mat)

	v1 = (1:p)[beta.s.mat[, b]!=0]	
	s.breve.star.vec = (lam/(1*1))*sign(beta.s.mat[v1, b])/abs(beta.s.init.mat[v1, b])
	C11.star.hat.inv = solve((t(X.mat[ ,v1])%*%X.mat[ ,v1])/n)
	b.breve.star = -as.numeric(matrix(d.vec[v1], nrow = 1)%*%(C11.star.hat.inv)%*%matrix(s.breve.star.vec, ncol = 1))

	term.1.star = (C11.star.hat.inv)%*%matrix(s.breve.star.vec, ncol = 1)
	sigma.star.dot.2.sqr = mean((y.star.mat[ ,b] - X.mat[ ,v1]%*%(beta.s.mat[v1, b] - term.1.star))^2)
	T.star = sqrt(n)*matrix(d.vec, nrow = 1)%*%(beta.s.mat[ ,b] - beta.hat)
	
	R.star.dot.2 = (T.star - sqrt(n)*b.breve.star)/sqrt(sigma.star.dot.2.sqr)
	return(abs(R.star.dot.2))
}

# intermediate step
c2.fn = function(b, beta.ss.mat, beta.ss.init.mat, X, y, CV.lam, d.vec, G.mat, mu.G, beta.hat)
{
	beta.ss = beta.ss.mat[ ,b]
	beta.ss.init = beta.ss.init.mat[ ,b]
	G.vec = G.mat[ ,b]
	wt.vec.mod = ((G.vec - mu.G)^2)/(mu.G^2)

	n = nrow(X)
	p = ncol(X)	
	Tn.ss = sqrt(n)*sum((d.vec)*as.numeric(beta.ss - beta.hat))

	A.ss = (1:p)[beta.ss!=0]
	
	C11.ss = t(X[ ,A.ss])%*%(X[ ,A.ss])/n
	s.breve.ss = (CV.lam/(1*1))*sign(beta.ss[A.ss])/abs(beta.ss.init[A.ss])
	beta.ss.dot2 = beta.ss[A.ss] - solve(C11.ss)%*%s.breve.ss
	res.ss = y - X[ ,A.ss]%*%beta.ss.dot2
	b.breve.ss = as.numeric(-matrix(d.vec[A.ss], nrow = 1)%*%solve(C11.ss)%*%matrix(s.breve.ss, ncol = 1))
	
	sigma.ss.dot2.sqr = mean((res.ss^2)*wt.vec.mod)

	term.2 = (Tn.ss - sqrt(n)*b.breve.ss)/sqrt(sigma.ss.dot2.sqr)
	return(term.2)	
}

# computing residual bootstrap based confidence interval
rbci.fn = function(y, X, CV.k, coef.index, B, alpha, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.init = init.info.list[[1]]
	beta.hat = init.info.list[[2]]
	CV.lam = init.info.list[[3]]
	
	theta.hat = sum(d.vec*beta.hat)
	
	res.vec = y - X%*%beta.hat
	c.res.vec = res.vec - mean(res.vec)
	A.hat = (1:p)[beta.hat != 0]
	d.hat = d.vec[A.hat]

	C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
	s.breve = (CV.lam/(1*1))*sign(beta.hat[A.hat])/abs(beta.init[A.hat])
	beta.hat.dot2 = beta.hat[A.hat] - solve(C11.hat)%*%s.breve
	res.temp = y - X[ ,A.hat]%*%beta.hat.dot2
	sigma.dot2.sqr = mean((res.temp)^2)
	b.breve = as.numeric(-matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%matrix(s.breve, ncol = 1))
	
	term.1 = sqrt(sigma.dot2.sqr)/sqrt(n) # needed for RB CI's

	# RB part #
	rboot.imat = sapply(1:B, function(i, z){sample(z, size = length(z), replace = T)}, z = 1:n) # n x B matrix of sampled indices
	y.boot.rb = as.vector(X%*%beta.hat) + sapply(1:B, function(i, b1.mat, c.res){c.res[b1.mat[ ,i]]}, b1.mat = rboot.imat, c.res = c.res.vec) #
	E.rb = lapply(1:B, alasso.within.boot, y.boot.rb, X, list(0, CV.k, CV.lam)) # B list: each is a px2 matrix # 
	beta.star.init.mat <- sapply(1:B, function(i, a1){a1[[i]][ ,1]}, a1 = E.rb) # pxB matrix of rboot.lasso initial estimates #
	beta.star.mat <- sapply(1:B, function(i, a1){a1[[i]][ ,2]}, a1 = E.rb) # pxB matrix of rboot.alasso estimates #

	R.star.dot.2.abs = sapply(1:B, c1.fn, beta.star.init.mat, beta.star.mat, X, CV.lam, d.vec, y.boot.rb, beta.hat)
	quant.rb = quantile(R.star.dot.2.abs, probs = 1-alpha)
	ci.rb = rep(theta.hat - b.breve, 2) + c(-1, 1)*term.1*quant.rb
	#length.rb = abs(ci.rb[2] - ci.rb[1])
    	return(ci.rb)
}

# computing perturbation bootstrap based confidence interval
pbci.fn = function(y, X, CV.k, coef.index, B, alpha, init.info.list)
{
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.init = init.info.list[[1]]
	beta.hat = init.info.list[[2]]
	CV.lam = init.info.list[[3]]
	
	theta.hat = sum(d.vec*beta.hat)

	res.vec = y - X%*%beta.hat
	A.hat = (1:p)[beta.hat != 0]
	d.hat = d.vec[A.hat]

	C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
	s.breve = (CV.lam/(1*1))*sign(beta.hat[A.hat])/abs(beta.init[A.hat])
	beta.hat.dot2 = beta.hat[A.hat] - solve(C11.hat)%*%s.breve
	res.temp = y - X[ ,A.hat]%*%beta.hat.dot2
	sigma.dot2.sqr = mean((res.temp)^2)
	b.breve = as.numeric(-matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%matrix(s.breve, ncol = 1))
	Sigma.hat = matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%matrix(d.hat, ncol = 1)
	X.hat.mod = t(sapply(1:nrow(X), function(i, a, v){a[i, ]*v[i]}, a = X[ ,A.hat], v = as.numeric(res.temp))) # n x p0.hat mat
	Sigma.dot2 = as.numeric(matrix(d.hat,nrow=1)%*%solve(C11.hat)%*%(t(X.hat.mod)%*%X.hat.mod/n)%*%solve(C11.hat)%*%matrix(d.hat,ncol=1))
	
	term.1 = sqrt(sigma.dot2.sqr)*sqrt(Sigma.hat)/sqrt(n) # to be used in the CI construction #

	# computing w2 and w4 #
	u1 = as.numeric(matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%t(X[ ,A.hat])) # n-vector
	w2 = -(1/sigma.dot2.sqr)*(1/Sigma.dot2)*mean((u1^2)*(res.temp^4)) + (1/sigma.dot2.sqr^2)*mean(res.temp^4)
	w4 = 2/(Sigma.dot2^2)*mean((u1^4)*(res.temp^4)) + 4*(1/sigma.dot2.sqr)*(1/Sigma.dot2)*mean((u1^2)*(res.temp^4)) - 3*(1/(sigma.dot2.sqr^2))*mean(res.temp^4) + 1
	z.alpha = qnorm(1 - alpha/2)
	C.np = -(z.alpha/n)*(w2/2 + (w4/24)*(z.alpha^2 - 3))

	# PB part #
	mu.star = (1/2)/(1/2 + 3/2)
	G.mat = sapply(1:B, function(i, n.size){rbeta(n.size, 1/2, 3/2)}, n.size = n) # n x B matrix 
	pboot.imat = apply(G.mat, 2, function(v2, c1){(1/c1)*(v2-c1)}, c1 = mu.star) # n x B transformed from G.mat#
	z.mat = as.vector(X%*%beta.hat) + sapply(1:B, function(b, u.mat, r.vec){r.vec*u.mat[ ,b]}, u.mat = pboot.imat, r.vec = res.vec) # nxB matrix
	
	E.pb = lapply(1:B, alasso.within.boot, z.mat, X, list(0, CV.k, CV.lam)) # B list: each is a px2 matrix #
	beta.ss.init.mat <- sapply(1:B, function(i, a1){a1[[i]][ ,1]}, a1 = E.pb) # pxB matrix of pboot.lasso initial est#
	beta.ss.mat <- sapply(1:B, function(i, a1){a1[[i]][ ,2]}, a1 = E.pb) # pxB matrix of pboot.alasso estimates #

	term2.ss = sapply(1:B, c2.fn, beta.ss.mat, beta.ss.init.mat, X, y, CV.lam, d.vec, G.mat, mu.star, beta.hat) # B-vector
	Rn.ss = as.numeric(term2.ss)*sqrt(sigma.dot2.sqr)*(Sigma.dot2^(-1/2)) # PB pivot # 
	quant.absRn.ss = quantile(abs(Rn.ss), probs = 1-alpha)
	hn.tilde.alpha = quant.absRn.ss + C.np

	pb.ci = rep(theta.hat - b.breve, 2) + c(-1,1)*as.numeric(term.1*hn.tilde.alpha)
	#len.pb.ci = abs(pb.ci[2] - pb.ci[1])
	return(pb.ci)
}


