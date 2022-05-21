


### [1] probability plot

plotprob <- function(probmatrix, label, main)
{
   o <- order(label)

   y <- label[o]

   probmatrix <- probmatrix[o,]

   n <- nrow(probmatrix)

   nc <- ncol(probmatrix)

   plot(row(probmatrix), probmatrix, type = "n", xlab = "observation",
        ylab = "probabilities", ylim = c(0, 1.2), axes = FALSE, main = main)

   axis(1)
   axis(2, at = seq(0, 1.2, by = 0.2), labels = c("0.0", "0.2",
        "0.4", "0.6", "0.8", "1.0", ""))

   axis(4)

    for (j in 1:nc)
    {
        points(1:n, probmatrix[,j], col = j + 1, pch=16)
    }

    for (j in 1:(nc - 1))
    {
        abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
    }

    abline(h = 1)

    for(j in 1:nc)
      {
        text(cumsum(table(y))[j] - table(y)[j]/2,
            1.1, labels = as.character(j-1), col=j+1,
            cex=0.9)
      }
}

### [2] for stratified cross-validation

roundvector <- function(x, maxint){
   fx <- floor(x)
   aftercomma <- x-fx
   roundorder <- order(aftercomma, decreasing=TRUE)
   i <- 1
   while(sum(fx) < maxint){ 
    fx[roundorder[i]] <- ceiling(x[roundorder[i]])
    i <- i+1
   }
 return(fx)
}

### [3] for stratified cross-validation -- rowswaps

rowswaps <- function(blocklist){

cols <- length(blocklist)
fold <- nrow(blocklist[[1]])
learnmatrix <- blocklist[[1]]
for(i in 2:cols) learnmatrix <- cbind(learnmatrix, blocklist[[i]])
rs <- rowSums( learnmatrix == 0)
allowedzeros <- ceiling(sum(rs)/fold)
indmatrix <-  matrix(rep(1:fold, each=cols), nrow=fold, byrow=TRUE) 
while(any(rs > allowedzeros)){
  indmatrix <- replicate(cols, sample(1:fold))
  temp2list <- blocklist
  for(i in 1:cols) temp2list[[i]] <- blocklist[[i]][indmatrix[,i], ]
  learnmatrix <- temp2list[[1]]
  for(i in 2:cols) learnmatrix <- cbind(learnmatrix, temp2list[[i]])
  rs <- rowSums( learnmatrix == 0)
}
return(indmatrix)
}

### [4] Receiver Operator characteristic

ROCinternal <- function(test, resp, plot, ...)
{
    dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          if(!hasArg(xlab)) ll$xlab <- "Threshold for assignment to class 1"
          if(!hasArg(ylab)) ll$ylab <- "specificity for class 0"
          if(!hasArg(main)) ll$main <- "Receiver Operator Characteristic"
          if(!hasArg(lwd)) ll$lwd <- 2
    m <- as.matrix(table(test, resp))
    fv <- sort(unique(test))
    nr <- dim(m)[1]
    a <- apply(m, 2, sum)
    m <- addmargins(m, 2)
    m <- apply(m[nr:1, ], 2, cumsum)[nr:1, ]
    sn <- c(m[, 2]/a[2], 0)
    sp <- c((a[1] - m[, 1])/a[1], 1)
    pvp <- c(m[, 2]/m[, 3], 1)
    pvn <- (a[1] - m[, 1])/(sum(a) - m[, 3])
    pvn <- c(pvn, rev(pvn)[1])
    res <- data.frame(cbind(sn, sp, pvp, pvn, c(NA, fv)))
    auc <- sum((res[-1, 1] + res[-nr, 1])/2 * diff(res[, 2]))
    #xl <- range(test)
    #ll$x <- xl
    #ll$y <- 0:1
    #ll$xlim <- xl
    #ll$ylim <- 0:1
    #ll$type <- "n"
    ll$x <- 1-res[,2]
    ll$y <- res[,1]
    ll$xlim <- 0:1
    ll$xlab <- "1-specificity"
    ll$ylim <- 0:1
    ll$ylab <- "Sensitivity"
    ll$type <- "n"
    if(plot){
    do.call("plot", args=ll)
            #plot(xl, 0:1, xlim = xl, xlab = paste(deparse(substitute(test)), 
            #    "(grid at deciles)"), ylim = 0:1, ylab = " ", 
            #    type = "n")
            
     #plot(1 - res[, 2], res[, 1], xlim = 0:1, xlab = "1-Specificity",
     #       ylim = 0:1, ylab = "Sensitivity", type = "n", ...)
     #   if (grid)
     #       abline(h = 0:10/10, v = 0:10/10, col = gray(0.9))
     #   abline(0, 1, col = gray(0.4))
        box()
     ll$type <- "l"
     do.call("lines", args = ll)

            
    #box()
    #ll$xlim <- NULL
    #ll$ylim <- NULL
    #ll$type <- "l"
    #ll$x <- fv
    #ll$y <- res[, 2]
    #ll$left <- TRUE
    #ll$order <- TRUE
    #do.call("steplines", args=ll)
    plot(function(x) x, from=0, to=1, lty="dashed", add=TRUE)
    text(0.8, 0.1, cex=2, label=paste("AUC=", round(auc,3), sep=""))
    }
    names(auc) <- "auc"
    return(invisible(auc))
}

### [6] penalized logistic regression

bklr <- function(y, Ka, Kp, lambda, w = 1, eps = 0.001, maxit = 20)
{
  N <- length(y)
  if(is.vector(Ka)) {
    p <- 1
    D <- rbind(0, cbind(0, 1))
    d <- Kp
    U <- 1/sqrt(d)
    bigU <- diag(c(1, U))
    Ka <- cbind(1, Ka * U)
  }
  else {
    p <- ncol(Ka)
    Kpeigen <- eigen(Kp, symmetric = TRUE)
    D <- rbind(0, cbind(0, diag(p)))
    d <- abs(Kpeigen$values)
    d[d < .Machine$double.eps] <- .Machine$double.eps
    U <- Kpeigen$vectors %*% diag(1/sqrt(d))
    bigU <- cbind(c(1, rep(0, p)), rbind(0, U))
    Ka <- cbind(1, Ka %*% U)
  }
  ybar <- mean(y)
  mu <- rep(ybar, N)
  beta0 <- binomial()$linkfun(ybar)
  beta <- c(beta0, rep(0, p))
  eta <- rep(beta0, N)
  wt <- sqrt((binomial()$mu.eta(mu)^2)/binomial()$variance(mu))
  z1 <- eta + (y - mu)/(binomial()$mu.eta(mu))
  tmpWZ <- wt * z1
  z2 <- t(Ka) %*% tmpWZ
  normd <- 1
  iter <- 0
  while((normd > eps) && (iter < maxit)) {
    beta.old <- beta
    tmpRes <- bkreg(z2, wt, lambda, Ka, D)
    beta <- tmpRes$beta
    eta <- tmpRes$fit
    mu <- binomial()$linkinv(eta)
    wt <- sqrt((binomial()$mu.eta(mu)^2)/binomial()$variance(mu))
    z1 <- eta + (y - mu)/(binomial()$mu.eta(mu))
    tmpWZ <- wt * z1
    z2 <- t(Ka) %*% tmpWZ
    normd <- sum((beta.old - beta)^2)/sum(beta.old^2)
    iter <- iter + 1
  }
  alpha <- bigU %*% beta
  predict <- (sign(eta) + 1)/2
  list(fit = eta, mu = mu, alpha = alpha, predict = predict, wt = wt)
}


bkreg <- function(z2, wt, lambda, Ka, D)
{
  K <- t(Ka) %*% (Ka * as.vector(wt)) + lambda * D
  svdK <- svd(K)
  if(any(svdK$d < .Machine$doubl.eps))
  stop("Numerical problems occured in plrCMA: kernel matrix is computationally singualar. Please check the input X. \n")
  invK <- svdK$v %*% diag(1/svdK$d, nrow(K)) %*% t(svdK$u)
  beta <- invK %*% z2
  fit <- Ka %*% beta
  list(beta = beta, fit = fit, invK = invK)
}

bklr.predict <- function(alpha, kernel, y = NULL)
{
  eta <- cbind(1, kernel) %*% alpha
  junk <- care.exp(eta)
  mu <- junk/(1 + junk)
  if (is.null(y))
    return( list(mu = mu, fit = eta, dev = NULL) )
  dev <- care.dev(mu, y)
  list(mu = mu, fit = eta, dev = dev)
}

care.exp <- function(x, thresh = 100) {
  about36 <-  - log(.Machine$double.eps)
  thresh <- max(c(thresh, about36))
  if( any(abs(x) > thresh) )
    warning("Fitted values over thresh")
  x[x > thresh] <- thresh
  x[x < ( - thresh)] <-  - thresh
  exp(x)
}

care.dev <- function(mu, y) {
  dev <- y * log(mu) + (1 - y) * log(1 - mu)
  if(any(small <- mu * (1 - mu) < .Machine$double.eps)) {
    warning("Fitted values close to 0 or 1")
    smu <- mu[small]
    sy <- y[small]
    smu <- ifelse(smu < .Machine$double.eps, .Machine$double.eps,
                  smu)
    onemsmu <- ifelse((1 - smu) < .Machine$double.eps,
                      .Machine$double.eps, 1 - smu)
    dev[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
  }
  -sum(dev)
}

mklr <- function(y, Ka, lambda, eps=0.001, maxit=30)
{
  N <- dim(y)[1]
  M <- dim(y)[2]
  bigU <- vector("list", length=M)
  bigKa <- vector("list", length=M)
  lend <- NULL
  beta <- vector("list", length=M)
  beta0 <- NULL
  z1 <- vector("list", length=M)
  z2 <- vector("list", length=M)
  z10 <- NULL
  z20 <- NULL
  for (i in 1:M) {
    Kpeigen <- eigen(Ka[,,i], symmetric=TRUE)
    d <- Kpeigen$values
    d <- d[d >= .Machine$double.eps]
    lend[i] <- length(d)
    bigU[[i]] <- Kpeigen$vectors[,seq(lend[i])] %*% diag(1/sqrt(d))
    bigKa[[i]] <- Kpeigen$vectors[,seq(lend[i])] %*% diag(sqrt(d))
    beta[[i]] <- rep(0, lend[i])
    z1[[i]] <- rep(0, lend[i])
    z2[[i]] <- rep(0, lend[i])
  }
  
  ybar <- apply(y, 2, mean)
  mu <- matrix(ybar, nrow=N, ncol=M, byrow=TRUE)
  beta0 <- log(ybar) - mean(log(ybar))
  eta <- matrix(beta0, nrow=N, ncol=M, byrow=TRUE)

  wteta <- (1-mu)*mu*eta
  for (i in 1:M) {
    z1[[i]] <- t(bigKa[[i]]) %*% wteta[,i]
    z2[[i]] <- t(bigKa[[i]]) %*% (y[,i] - mu[,i])
    z10[i] <- sum(wteta[,i])
    z20[i] <- sum(y[,i] - mu[,i])
  }
  dev <- -2 * (sum(eta * y, na.rm = TRUE) - sum(log(apply(my.care.exp(eta,
                              thresh=300), 1, sum))))
  pdev <- dev
  
  crit1 <- 1
  crit2 <- 1
  crit3 <- 1
  iter <- 0
  while((crit1 > eps || crit2 > eps/10 || crit3 > eps) & (iter < maxit)) {
    beta.old <- beta
    pdev.old <- pdev
    tmpRes <- mkreg(z1, mu, lambda, bigKa, lend, y, beta.old, pdev.old,
                    z2, z10, z20)
    beta <- tmpRes$beta
    beta0 <- tmpRes$beta0
    eta <- tmpRes$fit

    if ( any(abs(eta) > 37) )
      warning("eta in mklr > 37\n")
    junk <- my.care.exp(eta, thresh=300)
    mu <- junk/apply(junk, 1, sum)

    wteta <- (1-mu)*mu*eta
    for (i in 1:M) {
      z1[[i]] <- t(bigKa[[i]]) %*% wteta[,i] + lambda * beta[[i]]
      z2[[i]] <- t(bigKa[[i]]) %*% (y[,i] - mu[,i]) - lambda * beta[[i]]
      z10[i] <- sum(wteta[,i])
      z20[i] <- sum(y[,i] - mu[,i])
    }

    bold <- 0
    bnew <- 0
    znew <- 0
    for (i in 1:M) {
      bold <- bold + sum(beta.old[[i]]^2)
      bnew <- bnew + sum((beta.old[[i]] - beta[[i]])^2)
      znew <- znew + sum(z2[[i]]^2) + z20[i]^2
    }
    crit3 <- znew/(sum(lend)+M)
    pdev <- tmpRes$pdev
    crit2 <- abs((pdev.old - pdev)/pdev.old)
    crit1 <- bnew/bold
    iter <- iter + 1
  }
  
  alpha <- matrix(0, nrow=N+1, ncol=M)
  alpha[1,] <- beta0
  for (i in 1:M)
    alpha[2:(N+1),i] <- bigU[[i]] %*% beta[[i]]
  
  predict <- matrix(apply(mu, 1, order), nrow=M)[M,]
  wt <- mu*(1-mu)

  df <- 0
  for (i in 1:M) {
    KK <- t(bigKa[[i]]) %*% (bigKa[[i]] * wt[,i])
    df <- df + sum(diag(solve(KK + lambda * diag(lend[i]), KK)))
  }
  df <- df + M
  
  list(kclass=M, alpha=alpha, mu=mu, fit=eta, predict=predict, wt=wt,
       lend=lend, df=df)
}


mkreg <- function(z1, mu, lambda, bigKa, lend, y, beta.old, pdev.old, z2,
                  z10, z20) {
  N <- dim(mu)[1]
  M <- dim(mu)[2]

  j <- 0
  repeat {
    beta <- beta.old
    beta0 <- NULL
    tmppenal <- 0
    fit <- matrix(0, nrow=N, ncol=M)
    #cat("j", j, "\n")
    for (i in 1:M) {
      wt <- mu[,i] * (1-mu[,i])
      tmpW <- t(bigKa[[i]]) %*% wt
      tmpZ <- z1[[i]] + z2[[i]]/2^j
      tmpZ0 <- z10[i] + z20[i]/2^j
      KK <- t(bigKa[[i]]) %*% (bigKa[[i]]*wt) + lambda * diag(lend[i])
      tcoef <- solve(KK, cbind(tmpW, tmpZ))
      zcoef <- bigKa[[i]] %*% tcoef
      beta0[i] <- (tmpZ0 - sum(wt * zcoef[,2]))/sum(wt * (1 - zcoef[,1]))
      beta[[i]] <- tcoef[,2] - tcoef[,1] * beta0[i]
      fit[,i] <- zcoef[,2] - zcoef[,1] * beta0[i] + beta0[i]
      tmppenal <- tmppenal + sum(beta[[i]]^2)
    }
    beta0 <- beta0 - mean(beta0)

    ##if ( any(abs(fit) > 37) )
      #warning("fit in mkreg > 37\n")
    tmpDev <- -2 * (sum(fit * y, na.rm = TRUE) - sum(log(apply(my.care.exp(fit, thresh=300),1, sum))))
    tmpPdev <- tmpDev + lambda*tmppenal
    
    j <- j+1
    if (tmpPdev < pdev.old || j > 10)
      break
  }
  
  list(beta = beta, fit = fit, pdev=tmpPdev, beta0 = beta0)
}

mklr.predict <- function(klrfit, kernel, y = NULL)
{
  kclass <- klrfit$kclass
  N <- dim(y)[1]
  tmpind <- !is.na(y[,1])
  eta <- matrix(0, nrow=N, ncol=kclass)
  for (k in 1:kclass)
    eta[,k] <- cbind(1, kernel[,,k]) %*% klrfit$alpha[,k]
  junk1 <- my.care.exp(eta)
  junk2 <- apply(junk1, 1, sum)
  mu <- junk1/junk2
  predict <- matrix(apply(mu, 1, order), nrow=kclass)[kclass,]
  dev <- -2 * (sum(eta*y, na.rm=TRUE) - sum(log(junk2[tmpind]))) / sum(tmpind)
  
  list(mu = mu, predict = predict, dev = dev)
}

my.care.exp <- function(x, thresh = about36)
{
  about36 <-  - log(.Machine$double.eps)
  thresh <- max(c(thresh, about36))
  x[x > thresh] <- thresh
  x[x < ( - thresh)] <-  - thresh
  exp(x)
}

### [7] helper function for variable importance plot 

characterplot <- function(char, x, y, deltax, deltay, cex=1){
spltchar <- unlist(strsplit(char, ""))
for(s in seq(along=spltchar)) points(x+(s-1)*deltax, y+deltay, pch=spltchar[s], cex=cex)
}

### [8] for probability calculations in discriminant analysis for high dimensional
### data

safeexp <- function (x)
{
    xx = sign(x) * pmin(abs(x), 500)
    return(exp(xx))
}

### [9] L_2 penalized logistic regression routine for pls_lr

penlogitfit <- function(Z, y, lambda){
 ### preparations
 fam <- binomial()
 invlink <- fam$linkinv
 pp <- ncol(Z)-1
 eps <- 1e-10
 converged <- FALSE
 iter <- 0
 maxiter <- 25
 ###
 beta <- c(glm.fit(Z[,1],y)$coef, rep(0, pp))
  while((iter <= maxiter) & !converged){
    eta <- Z %*% beta
    mu <- invlink(eta)
    res <- y - mu
    grad <- t(Z) %*% res - c(0, rep(lambda, pp))
    d <- fam$mu.eta(eta)
    sigma <- fam$variance(mu)
    w <- drop(sqrt(d*d/sigma)) 
    Ztw <- t(w*Z)
    H <- -tcrossprod(Ztw)
    diag(H)[-1] <- diag(H)[-1] - lambda
    dir <- qr.solve(-H, grad)
    betaold <- beta
    beta <- beta + dir
    if(sqrt(sum((beta - betaold)^2)/sum(betaold * betaold)) < eps)
    converged <- TRUE
    iter <- iter+1
  }
  if(!converged) stop("Convergence failure in penalized logistic regression \n")
  else return(beta)
}



### [10] mean using rm=T (needed in compare.r) 

meanrm<-function(x){
mean(x,na.rm=T)
}







