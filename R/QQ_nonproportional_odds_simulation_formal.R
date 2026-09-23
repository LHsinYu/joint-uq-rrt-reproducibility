# Non-proportional-odds misspecification simulation for the UQ-RRT joint model
#
# Purpose:
#   Generate ordinal outcomes with cutpoint-specific covariate slopes, then fit
#   the existing proportional-odds EM model without changing its assumptions.
#   This directly assesses robustness to violations of proportional odds.
#
# The original scripts are not modified. This program does nothing unless
# QQ_RUN=1 is supplied. Use QQ_QUICK_TEST=1 for a small validation run.

rm(list = ls())
gc()

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(MASS))

safe_inverse <- function(x) {
  ans <- tryCatch(solve(x), error = function(e) NULL)
  if (is.null(ans) || any(!is.finite(ans))) return(NULL)
  ans
}

ordinal_prob_npo <- function(x, cutpoints, slopes, tolerance = 1e-12) {
  x <- as.numeric(x)
  cutpoints <- as.numeric(cutpoints)
  slopes <- as.numeric(slopes)
  if (length(cutpoints) != length(slopes)) {
    stop("Each ordinal cutpoint must have one data-generating slope.")
  }
  cumulative <- vapply(
    seq_along(cutpoints),
    function(k) plogis(cutpoints[k] + slopes[k] * x),
    numeric(length(x))
  )
  if (is.null(dim(cumulative))) cumulative <- matrix(cumulative, ncol = 1L)
  if (ncol(cumulative) > 1L && any(t(apply(cumulative, 1L, diff)) < -tolerance)) {
    stop("Cumulative probabilities cross; adjust cutpoints or slopes.")
  }
  category <- cbind(
    cumulative[, 1L],
    if (ncol(cumulative) > 1L) {
      cumulative[, 2:ncol(cumulative), drop = FALSE] -
        cumulative[, 1:(ncol(cumulative) - 1L), drop = FALSE]
    },
    1 - cumulative[, ncol(cumulative)]
  )
  category[abs(category) < tolerance] <- 0
  if (any(!is.finite(category)) || any(category < 0) ||
      any(abs(rowSums(category) - 1) > 1e-10)) {
    stop("Invalid ordinal category probabilities.")
  }
  category
}

draw_categories <- function(probability_matrix) {
  vapply(
    seq_len(nrow(probability_matrix)),
    function(i) sample.int(ncol(probability_matrix), 1L,
                           prob = probability_matrix[i, ]),
    integer(1L)
  )
}

Greenberg.RR.data.NPO <- function(nn, pp, cc, beta,
                                  cut0, slopes0, cut1, slopes1) {
  x <- rbinom(nn, 1L, 0.5)
  chi.1 <- cbind(`(Intercept)` = 1, x = x)
  latent_probability <- drop(plogis(chi.1 %*% beta))
  Y1 <- rbinom(nn, 1L, latent_probability)

  prob_z0 <- ordinal_prob_npo(x, cut0, slopes0)
  prob_z1 <- ordinal_prob_npo(x, cut1, slopes1)
  Z0 <- draw_categories(prob_z0)
  Z1 <- draw_categories(prob_z1)
  ZZ <- ifelse(Y1 == 1L, Z1, Z0)

  directed_to_sensitive <- rbinom(nn, 1L, pp)
  unrelated_answer <- rbinom(nn, 1L, cc)
  Y1.0 <- directed_to_sensitive * Y1 +
    (1L - directed_to_sensitive) * unrelated_answer

  list(
    data = cbind(Y1 = Y1, Y1.0 = Y1.0, ZZ = ZZ, chi.1),
    latent_probability = latent_probability,
    prob_z0 = prob_z0,
    prob_z1 = prob_z1
  )
}

true_posterior_npo <- function(Y1.0, ZZ, chi.1, pp, cc, beta,
                               prob_z0, prob_z1) {
  prior1 <- drop(plogis(chi.1 %*% beta))
  q1 <- pp + (1 - pp) * cc
  q0 <- (1 - pp) * cc
  rr1 <- ifelse(Y1.0 == 1L, q1, 1 - q1)
  rr0 <- ifelse(Y1.0 == 1L, q0, 1 - q0)
  idx <- cbind(seq_along(ZZ), ZZ)
  numerator <- prior1 * rr1 * prob_z1[idx]
  denominator <- numerator + (1 - prior1) * rr0 * prob_z0[idx]
  numerator / denominator
}

po_category_prob <- function(chi.1, ordinal_parameters, BB) {
  KK <- ncol(chi.1)
  x <- as.matrix(chi.1[, 2:KK, drop = FALSE])
  cut <- ordinal_parameters[seq_len(BB - 1L)]
  slope <- ordinal_parameters[BB:(BB + KK - 2L)]
  cumulative <- vapply(
    cut,
    function(a) plogis(a + drop(x %*% slope)),
    numeric(nrow(chi.1))
  )
  cbind(
    cumulative[, 1L],
    if (BB > 2L) cumulative[, 2:(BB - 1L), drop = FALSE] -
      cumulative[, 1:(BB - 2L), drop = FALSE],
    1 - cumulative[, BB - 1L]
  )
}

# Existing proportional-odds estimator core
CEY=function(Y1.0,ZZ,chi.1,pp,cc, Theta.hat,BB)
{
  KK=ncol(chi.1)
  xz.mat=chi.1[,2:KK]
  beta.h=Theta.hat[1:KK]
  alpha.00h=Theta.hat[(KK+1):length(Theta.hat)]
  alpha.0h=alpha.00h[1:(BB-1+(KK-1))]
  alpha.1h=alpha.00h[(BB+KK-1):(2*(BB-1+(KK-1)))]
  b1.h=chi.1%*%beta.h
  H1=1/(1+exp(-b1.h))
  
  
  gz0.pom=as.matrix(xz.mat)%*%alpha.0h[BB:(BB+KK-2)]
  gz1.pom=as.matrix(xz.mat)%*%alpha.1h[BB:(BB+KK-2)]
  HHZ.0=matrix(0,nrow(chi.1),BB)
  HHZ.1=matrix(0,nrow(chi.1),BB)
  HHZ.0[,1]=1/(1+exp(-alpha.0h[1]-gz0.pom))
  HHZ.1[,1]=1/(1+exp(-alpha.1h[1]-gz1.pom))
  
  for(kk in 2:(BB-1))
  {
    HHZ.0[,kk]=1/(1+exp(-alpha.0h[kk]-gz0.pom))-1/(1+exp(-alpha.0h[kk-1]-gz0.pom))
    HHZ.1[,kk]=1/(1+exp(-alpha.1h[kk]-gz1.pom))-1/(1+exp(-alpha.1h[kk-1]-gz1.pom))
  }
  HHZ.0[,BB]=1-1/(1+exp(-alpha.0h[BB-1]-gz0.pom))
  HHZ.1[,BB]=1-1/(1+exp(-alpha.1h[BB-1]-gz1.pom))
  NZ=cbind(c(1:nrow(chi.1)),ZZ)
  
  P1.Y10.ZZ=(pp+(1-pp)*cc)*H1*HHZ.1[NZ]+(1-pp)*cc*(1-H1)*HHZ.0[NZ]
  P0.Y10.ZZ=(1-pp)*(1-cc)*H1*HHZ.1[NZ]+(pp+(1-pp)*(1-cc))*(1-H1)*HHZ.0[NZ]
  ff.11.ZZ=(pp+(1-pp)*cc)*H1*HHZ.1[NZ]/P1.Y10.ZZ
  ff.10.ZZ=(1-pp)*(1-cc)*H1*HHZ.1[NZ]/P0.Y10.ZZ
  
  EE.Y=Y1.0*ff.11.ZZ+(1-Y1.0)*ff.10.ZZ
  return(EE.Y)
}

## Three part of Q-function
#
Greenberg.and.Z.reg1=function(pp,cc,beta.1,E.Y,ZZ,chi.1)
{
  KK=ncol(chi.1)
  ## beta.1=Theta.0[1:KK]
  b1=chi.1%*%beta.1
  H1=1/(1+exp(-b1))
  log.like.GZ=sum(E.Y*log(H1)+(1-E.Y)*log(1-H1))
  return(-log.like.GZ)
}
##
Greenberg.and.Z.reg11=function(pp,cc,alpha.1,E.Y,ZZ,chi.1,BB)
{
  KK=ncol(chi.1)
  xz.mat=chi.1[,2:KK]
  HHZ.1=matrix(0,nrow(chi.1),BB)
  gz1.pom=as.matrix(xz.mat)%*%alpha.1[BB:(BB+KK-2)]
  
  HHZ.1[,1]=1/(1+exp(-alpha.1[1]-gz1.pom))
  
  for(kk in 2:(BB-1))
  {
    HHZ.1[,kk]=1/(1+exp(-alpha.1[kk]-gz1.pom))-1/(1+exp(-alpha.1[kk-1]-gz1.pom))
  }
  HHZ.1[,BB]=1-1/(1+exp(-alpha.1[BB-1]-gz1.pom))
  if(any(!is.finite(HHZ.1)) || any(HHZ.1 <= 0)) return(.Machine$double.xmax)
  NZ=cbind(c(1:nrow(chi.1)),ZZ)
  log.like.GZ1=sum(E.Y*log(HHZ.1[NZ]))
  return(-log.like.GZ1)
}
Greenberg.and.Z.reg10=function(pp,cc,alpha.0,E.Y,ZZ,chi.1,BB)
{
  KK=ncol(chi.1)
  xz.mat=chi.1[,2:KK]
  gz0.pom=as.matrix(xz.mat)%*%alpha.0[BB:(BB+KK-2)]
  HHZ.0=matrix(0,nrow(chi.1),BB)
  HHZ.0[,1]=1/(1+exp(-alpha.0[1]-gz0.pom))
  
  for(kk in 2:(BB-1))
  {
    HHZ.0[,kk]=1/(1+exp(-alpha.0[kk]-gz0.pom))-1/(1+exp(-alpha.0[kk-1]-gz0.pom))
  }
  HHZ.0[,BB]=1-1/(1+exp(-alpha.0[BB-1]-gz0.pom))
  if(any(!is.finite(HHZ.0)) || any(HHZ.0 <= 0)) return(.Machine$double.xmax)
  NZ=cbind(c(1:nrow(chi.1)),ZZ)
  log.like.GZ0=sum((1-E.Y)*log(HHZ.0[NZ]))
  return(-log.like.GZ0)
}
EM.alogoritmH1=function(Y1.0,ZZ,chi.1,pp,cc, Theta.hat,BB)
{
  KK=ncol(chi.1)
  KS=0
  Err=1
  index.convergence=1
  Est.Theta=rep(-99999,length(Theta.hat))
  while(KS<=500 & Err >=0.00001)
  {
    E.Y=CEY(Y1.0,ZZ,chi.1,pp,cc, Theta.hat,BB)
    beta.hat=Theta.hat[1:KK]
    alpha.0.hat=Theta.hat[(KK+1):(KK+(BB-1)+(KK-1))]
    alpha.1.hat=Theta.hat[(2*KK+BB-1):length(Theta.hat)]
    beta.1.est=nlminb(beta.hat, Greenberg.and.Z.reg1,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE) #, lower = -4, upper =4)
    alpha.0.est=nlminb(alpha.0.hat, Greenberg.and.Z.reg10,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,BB=BB,hessian=TRUE) #, lower = -4, upper =4)
    alpha.1.est=nlminb(alpha.1.hat, Greenberg.and.Z.reg11,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,BB=BB,hessian=TRUE) #, lower = -4, upper =4)
    Theta.new.hat=c(beta.1.est$par,alpha.0.est$par,alpha.1.est$par)
    Err=sum(abs(Theta.hat-Theta.new.hat))/length(Theta.hat)
    if(Err<=0.00001) {Est.Theta=Theta.hat; index.convergence=0 }
    Theta.hat=Theta.new.hat
    KS=KS+1
    #cat("KS=",KS," ","Err=",Err," ")
  }
  return(c(Est.Theta,Err,index.convergence))
}
### Q-function
QQ.f=function(pp,cc,Y1.0,ZZ,chi.1, Theta.hat.EY, Theta.H)
{
  KK=ncol(chi.1)
  beta.1=Theta.H[1:KK];
  alpha.0=Theta.H[(KK+1):(KK+(BB-1)+(KK-1))];  alpha.1=Theta.H[(2*KK+BB-1):length(Theta.H)]
  E1.Y=CEY(Y1.0,ZZ,chi.1,pp,cc, Theta.hat.EY,BB)
  QQ=Greenberg.and.Z.reg1(pp,cc,beta.1,E1.Y,ZZ,chi.1,BB)+Greenberg.and.Z.reg11(pp,cc,alpha.1,E1.Y,ZZ,chi.1,BB)+Greenberg.and.Z.reg10(pp,cc,alpha.0,E1.Y,ZZ,chi.1,BB)
  return(-QQ)
}
CP.95=function(est1,est1.se,Theta.t)
{
  NP=length(Theta.t)
  A.025=matrix(0,nrow(est1),NP)
  A.975=matrix(0,nrow(est1),NP)
  A.95CP=matrix(0,nrow(est1),NP)
  for(k in 1:nrow(est1))
  {
    A.025[k,]=est1[k,]-1.96*est1.se[k,]
    A.975[k,]=est1[k,]+1.96*est1.se[k,]
    for(j in 1:NP)
    {
      if(Theta.t[j]>=A.025[k,j] & Theta.t[j]<=A.975[k,j]) {A.95CP[k,j]=1}
    } #end for j
  } ## end for k
  return(A.95CP)
}
##Orginal Log likelihood for (Y1.0,ZZ)
Log.like.GRY0Z=function(pp,cc,Y1.0,ZZ,chi.1,Theta.hat,BB)
{
  KK=ncol(chi.1)
  xz.mat=chi.1[,2:KK]
  beta.h=Theta.hat[1:KK]
  alpha.00h=Theta.hat[(KK+1):length(Theta.hat)]
  alpha.0h=alpha.00h[1:(BB-1+(KK-1))]
  alpha.1h=alpha.00h[(BB+KK-1):(2*(BB-1+(KK-1)))]
  b1.h=chi.1%*%beta.h
  H1=1/(1+exp(-b1.h))
  
  
  
  gz0.pom=as.matrix(xz.mat)%*%alpha.0h[BB:(BB+KK-2)]
  gz1.pom=as.matrix(xz.mat)%*%alpha.1h[BB:(BB+KK-2)]
  HHZ.0=matrix(0,nrow(chi.1),BB)
  HHZ.1=matrix(0,nrow(chi.1),BB)
  HHZ.0[,1]=1/(1+exp(-alpha.0h[1]-gz0.pom))
  HHZ.1[,1]=1/(1+exp(-alpha.1h[1]-gz1.pom))
  
  for(kk in 2:(BB-1))
  {
    HHZ.0[,kk]=1/(1+exp(-alpha.0h[kk]-gz0.pom))-1/(1+exp(-alpha.0h[kk-1]-gz0.pom))
    HHZ.1[,kk]=1/(1+exp(-alpha.1h[kk]-gz1.pom))-1/(1+exp(-alpha.1h[kk-1]-gz1.pom))
  }
  HHZ.0[,BB]=1-1/(1+exp(-alpha.0h[BB-1]-gz0.pom))
  HHZ.1[,BB]=1-1/(1+exp(-alpha.1h[BB-1]-gz1.pom))
  NZ=cbind(c(1:nrow(chi.1)),ZZ)
  P1.Y10.ZZ=(pp+(1-pp)*cc)*H1*HHZ.1[NZ]+(1-pp)*cc*(1-H1)*HHZ.0[NZ]
  P0.Y10.ZZ=(1-pp)*(1-cc)*H1*HHZ.1[NZ]+(pp+(1-pp)*(1-cc))*(1-H1)*HHZ.0[NZ]
  
  log.like.Y0Z=sum(Y1.0*log(P1.Y10.ZZ)+(1-Y1.0)*log(P0.Y10.ZZ))
  return(log.like.Y0Z)
}
#partial derivate of Q in Theta.hat part,  to let Q function be a vectors
Df1.LogLike.Y0Z=function(pp,cc,Y1.0,ZZ,chi.1, Theta.hat,BB)
{
  delta=0.00005
  k=length(Theta.hat)
  d1=diag(rep(delta,k))                        
  df.1=matrix(0,1,k)                            
  for(i in 1:k)                                 
  {       
    df.1[1,i]=(Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.hat+d1[i,],BB)-Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.hat-d1[i,],BB))/(2*delta)    
  }       
  return(t(df.1)) 
}
##
Df2.LogLike.Y0Z=function(pp,cc,Y1.0,ZZ,chi.1, Theta.hat,BB)
{
  delta=0.00005
  k=length(Theta.hat)                               
  d1=diag(rep(delta,k))                        
  df.1=matrix(0,k,k)                            
  for(i in 1:k)                                 
  {       
    df.1[,i]=(Df1.LogLike.Y0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.hat+d1[i,],BB)-Df1.LogLike.Y0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.hat-d1[i,],BB))/(2*delta) 
  }       
  return(df.1) 
}

parameter_names <- function(BB, KK) {
  c(
    paste0("beta", 0:(KK - 1L)),
    paste0("cut0_", seq_len(BB - 1L)),
    paste0("slope0_", seq_len(KK - 1L)),
    paste0("cut1_", seq_len(BB - 1L)),
    paste0("slope1_", seq_len(KK - 1L))
  )
}

make_start <- function(beta, cut0, slopes0, cut1, slopes1,
                       random = FALSE, sd = 0.20) {
  common0 <- mean(slopes0)
  common1 <- mean(slopes1)
  ans <- c(beta, cut0, common0, cut1, common1)
  if (random) {
    ans <- ans + rnorm(length(ans), 0, sd)
    KK <- length(beta)
    BB <- length(cut0) + 1L
    idx_cut0 <- (KK + 1L):(KK + BB - 1L)
    idx_cut1 <- (2L * KK + BB - 1L):(2L * KK + 2L * BB - 3L)
    ans[idx_cut0] <- sort(ans[idx_cut0])
    ans[idx_cut1] <- sort(ans[idx_cut1])
  }
  ans
}

fit_po_multistart <- function(Y1.0, ZZ, chi.1, pp, cc, base_start, BB,
                              n_starts = 5L) {
  starts <- vector("list", n_starts)
  starts[[1L]] <- base_start
  if (n_starts > 1L) {
    for (s in 2:n_starts) starts[[s]] <- make_start_from_vector(base_start, BB)
  }
  fits <- lapply(starts, function(start) {
    ans <- tryCatch(
      EM.alogoritmH1(Y1.0, ZZ, chi.1, pp, cc, start, BB),
      error = function(e) e
    )
    if (inherits(ans, "error")) {
      return(list(ok = FALSE, reason = paste0("EM_error: ", conditionMessage(ans))))
    }
    p <- length(start)
    theta <- ans[seq_len(p)]
    err <- ans[p + 1L]
    flag <- ans[p + 2L]
    if (flag != 0 || !is.finite(err) || err > 1e-5 || any(!is.finite(theta))) {
      return(list(ok = FALSE, reason = "EM_not_converged"))
    }
    loglik <- tryCatch(
      Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, theta, BB),
      error = function(e) NA_real_
    )
    list(ok = is.finite(loglik), theta = theta, error = err,
         loglik = loglik, reason = if (is.finite(loglik)) "" else "invalid_loglik")
  })
  good <- which(vapply(fits, function(z) isTRUE(z$ok), logical(1L)))
  if (!length(good)) {
    reasons <- paste(unique(vapply(fits, `[[`, character(1L), "reason")), collapse = "; ")
    return(list(ok = FALSE, reason = reasons))
  }
  chosen <- good[which.max(vapply(fits[good], `[[`, numeric(1L), "loglik"))]
  best <- fits[[chosen]]
  best$n_converged_starts <- length(good)
  best$start_agreement <- if (length(good) > 1L) {
    max(vapply(fits[good], function(z) abs(z$loglik - best$loglik), numeric(1L))) < 1e-5
  } else NA
  best
}

fit_reduced_homogeneous <- function(Y1.0, ZZ, chi.1, pp, cc, full_theta, BB) {
  KK <- ncol(chi.1)
  ordinal_length <- (BB - 1L) + (KK - 1L)
  a0 <- full_theta[(KK + 1L):(KK + ordinal_length)]
  a1 <- full_theta[(KK + ordinal_length + 1L):length(full_theta)]
  start <- c(full_theta[seq_len(KK)], (a0 + a1) / 2)
  objective <- function(reduced) {
    cutpoints <- reduced[(KK + 1L):(KK + BB - 1L)]
    if (any(!is.finite(reduced)) || any(diff(cutpoints) <= 1e-6)) return(1e100)
    common <- reduced[(KK + 1L):length(reduced)]
    expanded <- c(reduced[seq_len(KK)], common, common)
    loglik <- tryCatch(
      Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, expanded, BB),
      error = function(e) NA_real_
    )
    if (!is.finite(loglik)) 1e100 else -loglik
  }
  fit <- tryCatch(
    nlminb(start, objective, control = list(iter.max = 2000L, eval.max = 4000L,
                                            rel.tol = 1e-10)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0L || !is.finite(fit$objective)) {
    return(list(ok = FALSE, loglik = NA_real_))
  }
  list(ok = TRUE, theta = fit$par, loglik = -fit$objective)
}

make_start_from_vector <- function(base, BB, sd = 0.20) {
  ans <- base + rnorm(length(base), 0, sd)
  KK <- (length(base) - 2L * (BB - 1L)) / 3L
  KK <- as.integer(KK)
  idx_cut0 <- (KK + 1L):(KK + BB - 1L)
  start_alpha1 <- KK + (BB - 1L) + (KK - 1L) + 1L
  idx_cut1 <- start_alpha1:(start_alpha1 + BB - 2L)
  ans[idx_cut0] <- sort(ans[idx_cut0])
  ans[idx_cut1] <- sort(ans[idx_cut1])
  ans
}

empty_failure <- function(scenario, n, replication, attempt, reason) {
  data.frame(
    scenario = scenario, n = n, replication = replication, attempt = attempt,
    status = "failed", failure_reason = reason,
    convergence_error = NA_real_, loglik = NA_real_,
    n_converged_starts = 0L, start_agreement = NA,
    hessian_valid = FALSE, extreme_solution = NA,
    prevalence_true = NA_real_, prevalence_est = NA_real_,
    prevalence_error = NA_real_, posterior_mae = NA_real_,
    posterior_rmse = NA_real_, brier_fitted = NA_real_, brier_true = NA_real_,
    ordinal_prob_mae = NA_real_, ordinal_prob_rmse = NA_real_,
    lr_statistic = NA_real_, lr_p_value = NA_real_,
    wald_joint_statistic = NA_real_, wald_joint_p_value = NA_real_,
    stringsAsFactors = FALSE
  )
}

simulate_one <- function(scenario_row, n, replication, attempt, pp, cc, beta,
                         cut0, cut1, n_starts) {
  scenario <- scenario_row$scenario
  slopes0 <- unlist(scenario_row[c("slope0_k1", "slope0_k2")], use.names = FALSE)
  slopes1 <- unlist(scenario_row[c("slope1_k1", "slope1_k2")], use.names = FALSE)
  generated <- tryCatch(
    Greenberg.RR.data.NPO(n, pp, cc, beta, cut0, slopes0, cut1, slopes1),
    error = function(e) e
  )
  if (inherits(generated, "error")) {
    return(empty_failure(scenario, n, replication, attempt,
                         paste0("data_generation: ", conditionMessage(generated))))
  }
  dat <- generated$data
  Y1 <- dat[, "Y1"]
  Y1.0 <- dat[, "Y1.0"]
  ZZ <- dat[, "ZZ"]
  chi.1 <- dat[, 4:ncol(dat), drop = FALSE]
  BB <- ncol(generated$prob_z0)
  KK <- ncol(chi.1)
  base_start <- make_start(beta, cut0, slopes0, cut1, slopes1)
  fit <- fit_po_multistart(Y1.0, ZZ, chi.1, pp, cc, base_start, BB, n_starts)
  if (!isTRUE(fit$ok)) {
    return(empty_failure(scenario, n, replication, attempt, fit$reason))
  }
  theta <- fit$theta
  names(theta) <- parameter_names(BB, KK)

  hessian <- tryCatch(
    Df2.LogLike.Y0Z(pp, cc, Y1.0, ZZ, chi.1, theta, BB),
    error = function(e) NULL
  )
  covariance <- if (is.null(hessian) || any(!is.finite(hessian))) NULL else
    safe_inverse(-hessian)
  hessian_valid <- !is.null(covariance) && all(diag(covariance) > 0)
  se <- rep(NA_real_, length(theta))
  names(se) <- names(theta)
  if (hessian_valid) {
    se <- sqrt(diag(covariance))
    names(se) <- names(theta)
  }

  fitted_posterior <- CEY(Y1.0, ZZ, chi.1, pp, cc, theta, BB)
  true_posterior <- true_posterior_npo(
    Y1.0, ZZ, chi.1, pp, cc, beta,
    generated$prob_z0, generated$prob_z1
  )
  a0_idx <- (KK + 1L):(KK + BB + KK - 2L)
  a1_idx <- (2L * KK + BB - 1L):length(theta)
  fitted_z0 <- po_category_prob(chi.1, theta[a0_idx], BB)
  fitted_z1 <- po_category_prob(chi.1, theta[a1_idx], BB)
  prob_error <- c(fitted_z0 - generated$prob_z0,
                  fitted_z1 - generated$prob_z1)
  prevalence_est <- mean(plogis(drop(chi.1 %*% theta[seq_len(KK)])))
  extreme <- any(abs(theta) > 10) || any(se > 20, na.rm = TRUE)

  reduced_fit <- fit_reduced_homogeneous(Y1.0, ZZ, chi.1, pp, cc,
                                         theta, BB)
  lr_statistic <- if (isTRUE(reduced_fit$ok))
    max(0, 2 * (fit$loglik - reduced_fit$loglik)) else NA_real_
  homogeneity_df <- (BB - 1L) + (KK - 1L)
  lr_p_value <- if (is.finite(lr_statistic))
    pchisq(lr_statistic, df = homogeneity_df, lower.tail = FALSE) else NA_real_

  wald_joint_statistic <- NA_real_
  wald_joint_p_value <- NA_real_
  if (hessian_valid) {
    ordinal_length <- (BB - 1L) + (KK - 1L)
    contrast <- cbind(matrix(0, ordinal_length, KK),
                      diag(ordinal_length), -diag(ordinal_length))
    contrast_covariance <- contrast %*% covariance %*% t(contrast)
    inverse_contrast_covariance <- safe_inverse(contrast_covariance)
    difference <- contrast %*% theta
    if (!is.null(inverse_contrast_covariance)) {
      wald_joint_statistic <- drop(t(difference) %*%
                                     inverse_contrast_covariance %*% difference)
      wald_joint_p_value <- pchisq(wald_joint_statistic,
                                   df = homogeneity_df,
                                   lower.tail = FALSE)
    }
  }

  base <- data.frame(
    scenario = scenario, n = n, replication = replication, attempt = attempt,
    status = "converged", failure_reason = "",
    convergence_error = fit$error, loglik = fit$loglik,
    n_converged_starts = fit$n_converged_starts,
    start_agreement = fit$start_agreement,
    hessian_valid = hessian_valid, extreme_solution = extreme,
    prevalence_true = mean(generated$latent_probability),
    prevalence_est = prevalence_est,
    prevalence_error = prevalence_est - mean(generated$latent_probability),
    posterior_mae = mean(abs(fitted_posterior - true_posterior)),
    posterior_rmse = sqrt(mean((fitted_posterior - true_posterior)^2)),
    brier_fitted = mean((fitted_posterior - Y1)^2),
    brier_true = mean((true_posterior - Y1)^2),
    ordinal_prob_mae = mean(abs(prob_error)),
    ordinal_prob_rmse = sqrt(mean(prob_error^2)),
    lr_statistic = lr_statistic, lr_p_value = lr_p_value,
    wald_joint_statistic = wald_joint_statistic,
    wald_joint_p_value = wald_joint_p_value,
    stringsAsFactors = FALSE
  )
  for (nm in names(theta)) {
    base[[paste0("estimate_", nm)]] <- theta[[nm]]
    base[[paste0("se_", nm)]] <- se[[nm]]
  }
  base
}

increment_run_dir <- function(root, stem) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  i <- 1L
  repeat {
    candidate <- file.path(root, paste0(stem, "__run_", i))
    if (!dir.exists(candidate)) {
      dir.create(candidate, recursive = TRUE)
      return(candidate)
    }
    i <- i + 1L
  }
}

rbind_fill <- function(frames) {
  all_names <- unique(unlist(lapply(frames, names), use.names = FALSE))
  aligned <- lapply(frames, function(frame) {
    missing <- setdiff(all_names, names(frame))
    for (name in missing) frame[[name]] <- NA
    frame[all_names]
  })
  do.call(rbind, aligned)
}

summarize_results <- function(results, beta, po_truth) {
  converged <- results[results$status == "converged", , drop = FALSE]
  scenarios <- unique(results[c("scenario", "n")])
  scenario_summary <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    key <- scenarios[i, ]
    all_i <- results[results$scenario == key$scenario & results$n == key$n, , drop = FALSE]
    ok <- all_i[all_i$status == "converged", , drop = FALSE]
    data.frame(
      scenario = key$scenario, n = key$n,
      attempts = nrow(all_i), converged = nrow(ok),
      convergence_rate = nrow(ok) / nrow(all_i),
      valid_hessian_rate = if (nrow(ok)) mean(ok$hessian_valid) else NA_real_,
      extreme_rate = if (nrow(ok)) mean(ok$extreme_solution, na.rm = TRUE) else NA_real_,
      mean_prevalence_error = if (nrow(ok)) mean(ok$prevalence_error) else NA_real_,
      prevalence_rmse = if (nrow(ok)) sqrt(mean(ok$prevalence_error^2)) else NA_real_,
      mean_posterior_mae = if (nrow(ok)) mean(ok$posterior_mae) else NA_real_,
      mean_posterior_rmse = if (nrow(ok)) mean(ok$posterior_rmse) else NA_real_,
      mean_ordinal_prob_mae = if (nrow(ok)) mean(ok$ordinal_prob_mae) else NA_real_,
      mean_ordinal_prob_rmse = if (nrow(ok)) mean(ok$ordinal_prob_rmse) else NA_real_,
      lr_valid = if (nrow(ok)) sum(is.finite(ok$lr_p_value)) else 0L,
      lr_rejection_rate_005 = if (nrow(ok) && any(is.finite(ok$lr_p_value)))
        mean(ok$lr_p_value[is.finite(ok$lr_p_value)] < 0.05) else NA_real_,
      lr_mcse = if (nrow(ok) && any(is.finite(ok$lr_p_value))) {
        rate <- mean(ok$lr_p_value[is.finite(ok$lr_p_value)] < 0.05)
        sqrt(rate * (1 - rate) / sum(is.finite(ok$lr_p_value)))
      } else NA_real_,
      wald_valid = if (nrow(ok)) sum(is.finite(ok$wald_joint_p_value)) else 0L,
      wald_rejection_rate_005 = if (nrow(ok) && any(is.finite(ok$wald_joint_p_value)))
        mean(ok$wald_joint_p_value[is.finite(ok$wald_joint_p_value)] < 0.05) else NA_real_,
      wald_mcse = if (nrow(ok) && any(is.finite(ok$wald_joint_p_value))) {
        rate <- mean(ok$wald_joint_p_value[is.finite(ok$wald_joint_p_value)] < 0.05)
        sqrt(rate * (1 - rate) / sum(is.finite(ok$wald_joint_p_value)))
      } else NA_real_
    )
  }))

  estimate_cols <- grep("^estimate_", names(converged), value = TRUE)
  parameter_summary <- do.call(rbind, lapply(split(converged, list(converged$scenario,
                                                                   converged$n), drop = TRUE),
    function(d) {
      do.call(rbind, lapply(estimate_cols, function(ec) {
        parameter <- sub("^estimate_", "", ec)
        sc <- paste0("se_", parameter)
        est <- d[[ec]]
        se <- d[[sc]]
        truth <- NA_real_
        if (parameter == "beta0") truth <- beta[1L]
        if (parameter == "beta1") truth <- beta[2L]
        if (d$scenario[1L] == "PO") truth <- po_truth[[parameter]]
        valid_se <- is.finite(se)
        data.frame(
          scenario = d$scenario[1L], n = d$n[1L], parameter = parameter,
          truth = truth, n_estimate = sum(is.finite(est)),
          mean_estimate = mean(est, na.rm = TRUE),
          bias = if (is.finite(truth)) mean(est - truth, na.rm = TRUE) else NA_real_,
          empirical_sd = sd(est, na.rm = TRUE),
          mean_se = if (any(valid_se)) mean(se[valid_se]) else NA_real_,
          rmse = if (is.finite(truth)) sqrt(mean((est - truth)^2, na.rm = TRUE)) else NA_real_,
          coverage_95 = if (is.finite(truth) && any(valid_se))
            mean(est[valid_se] - 1.96 * se[valid_se] <= truth &
                   est[valid_se] + 1.96 * se[valid_se] >= truth) else NA_real_
        )
      }))
    }))
  list(scenario = scenario_summary, parameter = parameter_summary)
}

run_requested <- identical(Sys.getenv("QQ_RUN", "0"), "1")
if (!run_requested) {
  cat("Program created successfully; no simulation was started.\n")
  cat("Quick test: QQ_RUN=1 QQ_QUICK_TEST=1 Rscript R/QQ_nonproportional_odds_simulation_formal.R\n")
  cat("Formal run: QQ_RUN=1 Rscript R/QQ_nonproportional_odds_simulation_formal.R\n")
} else {
  quick <- identical(Sys.getenv("QQ_QUICK_TEST", "0"), "1")
  seed <- as.integer(Sys.getenv("QQ_SEED", "20260917"))
  set.seed(seed)
  pp <- as.numeric(Sys.getenv("QQ_P", "0.5"))
  cc <- as.numeric(Sys.getenv("QQ_C", "0.25"))
  beta <- c(-1, -1)
  cut0 <- c(-2, -0.5)
  cut1 <- c(-2, -0.5)
  sample_sizes <- if (quick) 120L else c(500L, 1000L)
  replications <- if (quick) 2L else as.integer(Sys.getenv("QQ_KT", "1000"))
  max_extra_attempts <- if (quick) 2L else as.integer(Sys.getenv("QQ_MAX_EXTRA", "500"))
  n_starts <- if (quick) 2L else as.integer(Sys.getenv("QQ_N_STARTS", "5"))
  workers <- if (quick) 1L else as.integer(Sys.getenv(
    "QQ_WORKERS", as.character(max(1L, parallel::detectCores() - 1L))))

  scenarios <- data.frame(
    scenario = c("PO", "Mild_NPO", "Moderate_NPO"),
    slope0_k1 = c(1.00, 0.75, 0.50),
    slope0_k2 = c(1.00, 1.25, 1.50),
    slope1_k1 = c(1.00, 0.75, 0.50),
    slope1_k2 = c(1.00, 1.25, 1.50),
    stringsAsFactors = FALSE
  )
  design <- merge(scenarios, data.frame(n = sample_sizes), all = TRUE)
  root <- Sys.getenv("QQ_OUTPUT_DIR", file.path(getwd(), "qq_npo_results"))
  stem <- sprintf("NPO_p%s_c%s_KT%s_seed%s",
                  format(pp, trim = TRUE), format(cc, trim = TRUE),
                  replications, seed)
  requested_run_dir <- Sys.getenv("QQ_RUN_DIR", "")
  if (nzchar(requested_run_dir)) {
    run_dir <- normalizePath(requested_run_dir, mustWork = FALSE)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    run_dir <- increment_run_dir(root, stem)
  }
  write.csv(scenarios, file.path(run_dir, "scenarios.csv"), row.names = FALSE)
  settings <- data.frame(
    setting = c("seed", "p", "c", "sample_sizes", "replications",
                "max_extra_attempts", "n_starts", "workers", "quick_test"),
    value = c(seed, pp, cc, paste(sample_sizes, collapse = ";"), replications,
              max_extra_attempts, n_starts, workers, quick)
  )
  write.csv(settings, file.path(run_dir, "settings.csv"), row.names = FALSE)

  tasks <- lapply(seq_len(nrow(design)), function(i) design[i, , drop = FALSE])
  run_design_row <- function(row) {
    target <- replications
    max_attempts <- target + max_extra_attempts
    checkpoint_file <- file.path(
      run_dir,
      sprintf("checkpoint_%s_n%s.rds", row$scenario, row$n)
    )
    checkpoint <- if (file.exists(checkpoint_file))
      tryCatch(readRDS(checkpoint_file), error = function(e) NULL) else NULL
    if (is.list(checkpoint) && identical(checkpoint$target, target)) {
      accepted <- checkpoint$accepted
      attempt <- checkpoint$attempt
      out <- checkpoint$out
    } else {
      accepted <- 0L
      attempt <- 0L
      out <- list()
    }
    next_report <- max(10L, 10L * (floor(accepted / 10L) + 1L))
    while (accepted < target && attempt < max_attempts) {
      attempt <- attempt + 1L
      one <- simulate_one(row, row$n, accepted + 1L, attempt, pp, cc, beta,
                          cut0, cut1, n_starts)
      out[[length(out) + 1L]] <- one
      if (one$status == "converged") accepted <- accepted + 1L
      if (accepted %% 25L == 0L || accepted == target) {
        saveRDS(list(target = target, accepted = accepted, attempt = attempt,
                     out = out), checkpoint_file)
      }
      pct <- floor(100 * accepted / target)
      if (pct >= next_report) {
        cat(sprintf("[%s, N=%d] %d%% complete (%d attempts)\n",
                    row$scenario, row$n, pct, attempt))
        next_report <- next_report + 10L
      }
    }
    rbind_fill(out)
  }
  `%||%` <- function(x, y) if (is.null(x)) y else x
  if (workers > 1L && .Platform$OS.type != "windows") {
    result_list <- parallel::mclapply(tasks, run_design_row, mc.cores = workers,
                                      mc.set.seed = TRUE)
  } else {
    result_list <- lapply(tasks, run_design_row)
  }
  results <- do.call(rbind, result_list)
  rownames(results) <- NULL
  write.csv(results, file.path(run_dir, "replications_full.csv"), row.names = FALSE)

  po_truth <- c(
    beta0 = beta[1L], beta1 = beta[2L],
    cut0_1 = cut0[1L], cut0_2 = cut0[2L], slope0_1 = 1,
    cut1_1 = cut1[1L], cut1_2 = cut1[2L], slope1_1 = 1
  )
  summaries <- summarize_results(results, beta, po_truth)
  write.csv(summaries$scenario, file.path(run_dir, "scenario_summary.csv"), row.names = FALSE)
  write.csv(summaries$parameter, file.path(run_dir, "parameter_summary.csv"), row.names = FALSE)
  saveRDS(list(settings = settings, scenarios = scenarios, results = results,
               scenario_summary = summaries$scenario,
               parameter_summary = summaries$parameter),
          file.path(run_dir, "simulation_results.rds"))
  cat("Results written to:", normalizePath(run_dir), "\n")
}
