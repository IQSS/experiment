ATEtrunc <- function(Y, Z, grp = NULL, data = parent.frame(),
                     match = NULL, weights = NULL, fpc = TRUE,
                     monotone = "negative") {

  trunc.internal <- function(Y, gamma, lower = TRUE) {
    n <- length(Y)
    m <- round(gamma*n)
    Ys <- sort(Y)
    if (lower)
      return(mean(Ys[1:(n-m)]))
    else
      return(mean(Ys[(m+1):n]))
  }
  
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  W <- is.na(Y)
  grp <- eval(call$grp, envir = data)
  match <- eval(call$match, envir = data)
  weights <- eval(call$weights, envir = data)

  if (is.null(grp)) {
    stop("this option is not available yet")
  } else {
    res <- list(call = call, Y = Y, Z = Z, grp = grp,
                match = match, weights = weights) 
    if (is.null(match))
      stop("This option is not yet available.")
    else {
      res$m <- m <- length(unique(match))
      ## ignoring the truncation
      res$Y1bar <- Y1bar <- tapply(Y[Z==1 & W==0], match[Z==1 & W==0], mean)
      res$Y0bar <- Y0bar <- tapply(Y[Z==0 & W==0], match[Z==0 & W==0], mean)
      res$diff <- diff <- Y1bar-Y0bar
      res$n1 <- n1 <- tapply(rep(1, sum(Z==1)), match[Z==1], sum)
      res$n0 <- n0 <- tapply(rep(1, sum(Z==0)), match[Z==0], sum)
      res$n1.w <- n1.w <- tapply(rep(1, sum(Z==1 & W==0)), match[Z==1 & W==0], sum)
      res$n0.w <- n0.w <- tapply(rep(1, sum(Z==0 & W==0)), match[Z==0 & W==0], sum)
      n <- sum(n1+n0)
      n.w <- sum(n1.w + n0.w)
      ## with monotonicity assumption
      p0 <- mean(W[Z==0])
      p1 <- mean(W[Z==1])
      if (monotone == "positive") {
        if (p0 > p1)
          stop("Pr(W=1|Z=1) must be greater than Pr(W=1|Z=0)")
        else {
          res$Y0bar.lb <- Y0bar.lb <-
            trunc.internal(Y[Z == 0 & W == 0], 1-(1-p1)/(1-p0), lower = TRUE) 
          res$Y0bar.ub <- Y0bar.ub <-
            trunc.internal(Y[Z == 0 & W == 0], 1-(1-p1)/(1-p0), lower = FALSE) 
        }
      } else if (monotone == "negative") {
        if (p1 > p0)
          stop("Pr(W=1|Z=0) must be greater than Pr(W=1|Z=1)")
        else {
          res$Y1bar.lb <- Y1bar.lb <-
            trunc.internal(Y[Z == 1 & W == 0], 1-(1-p0)/(1-p1), lower = TRUE) 
          res$Y1bar.ub <- Y1bar.ub <-
            trunc.internal(Y[Z == 1 & W == 0], 1-(1-p0)/(1-p1), lower = FALSE) 
        }
      } else if (is.null(monotone)) {
        stop("this option is not yet available.")
      } else {
        stop("invalid input for `monotone'")
      }
    }
    
    if (is.null(weights)) {
      ## variance for PATE1 (sampling of clusters)
      N1 <- w1 <- n1
      N0 <- w0 <- n0
      N1.w <- w1.w <- n1.w
      N0.w <- w0.w <- n0.w
    } else {
      ## variance for PATE2 (double sampling)
      w1 <- w1.w <- N1 <- tapply(weights[Z==1], match[Z==1], mean)
      w0 <- w0.w <- N0 <- tapply(weights[Z==0], match[Z==0], mean)
    }
    w <- w1 + w0
    w <- n*w/sum(w)
    w.w <- w1.w + w0.w
    w.w <- n.w*w.w/sum(w.w)
    ## estimates ignoring the truncation
    ATE.est <- weighted.mean(diff, w.w)
    ATE.var <- m*sum((w.w*diff-n.w*ATE.est/m)^2)/((m-1)*(n.w^2))
    ## bounds with truncation
    if (monotone == "negative") {
      res$ATE.ub <- Y1bar.ub - mean(Y[Z==0], na.rm = TRUE)
      res$ATE.lb <- Y1bar.lb - mean(Y[Z==0], na.rm = TRUE)
    } else {
      res$ATE.ub <- mean(Y[Z==1], na.rm = TRUE) - Y0bar.lb
      res$ATE.lb <- mean(Y[Z==1], na.rm = TRUE) - Y0bar.ub
    }
    ## return the resutls
    res$est <- ATE.est
    res$var <- ATE.var
    res$w <- w
    if (!is.null(match)) {
      res$w1 <- w1
      res$w0 <- w0
      res$N0 <- N0
      res$N1 <- N1
    }
    class(res) <- "ATEtrunc"
    return(res)
  }
}
