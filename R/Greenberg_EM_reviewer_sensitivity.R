rm(list = ls())
gc()

# ============================================================
# 審查意見：固定模擬次數，完整保留收斂、失敗與極端解。
# 1. 不以真值篩選估計值，也不在失敗後補跑。
# 2. EM 收斂：平均絕對參數變化 <= EM_TOL，且 M-step 均正常結束。
# 3. 參數估計成功與 Hessian/SE 成功分開記錄，各統計量使用其實際可用分母。
# 4. 已收斂結果若 max(abs(theta_hat)) > 10 或 max(SE) > 20，只標記為
#    extreme solution，不從主要結果刪除。
# 5. 同時輸出「保留全部結果」與「排除 extreme solution」的參數摘要，
#    並另存極端解比例、判定門檻及各參數觸發次數。
# ============================================================
EM_TOL <- as.numeric(Sys.getenv("EM_TOL", "1e-6"))
EM_MAX_ITER <- as.integer(Sys.getenv("EM_MAX_ITER", "2000"))
HESSIAN_RCOND_TOL <- as.numeric(Sys.getenv("EM_HESSIAN_RCOND_TOL", "1e-12"))
EXTREME_COEF_CUTOFF <- as.numeric(Sys.getenv("EM_EXTREME_COEF", "10"))
EXTREME_SE_CUTOFF <- as.numeric(Sys.getenv("EM_EXTREME_SE", "20"))
N_STARTS <- as.integer(Sys.getenv("EM_N_STARTS", "10"))
START_LOGLIK_TOL <- as.numeric(Sys.getenv("EM_START_LOGLIK_TOL", "1e-6"))
START_PARAMETER_TOL <- as.numeric(Sys.getenv("EM_START_PARAMETER_TOL", "1e-3"))
detected.cores <- suppressWarnings(parallel::detectCores(logical=TRUE))
if(!is.finite(detected.cores)) detected.cores <- 2L
N_CORES <- as.integer(Sys.getenv("EM_N_CORES", as.character(max(1L,min(4L,detected.cores-1L)))))
CHECKPOINT_EVERY <- as.integer(Sys.getenv("EM_CHECKPOINT_EVERY", "10"))
RESUME_RUN <- tolower(Sys.getenv("EM_RESUME", "true")) %in% c("true","1","yes")
CREATE_PLOTS <- tolower(Sys.getenv("EM_CREATE_PLOTS", "false")) %in%
  c("true","1","yes")

stopifnot(
  is.finite(EM_TOL), EM_TOL > 0,
  is.finite(EM_MAX_ITER), EM_MAX_ITER >= 1,
  is.finite(HESSIAN_RCOND_TOL), HESSIAN_RCOND_TOL > 0,
  is.finite(EXTREME_COEF_CUTOFF), EXTREME_COEF_CUTOFF > 0,
  is.finite(EXTREME_SE_CUTOFF), EXTREME_SE_CUTOFF > 0,
  is.finite(N_STARTS), N_STARTS >= 1,
  is.finite(START_LOGLIK_TOL), START_LOGLIK_TOL > 0,
  is.finite(START_PARAMETER_TOL), START_PARAMETER_TOL > 0,
  is.finite(N_CORES), N_CORES >= 1,
  is.finite(CHECKPOINT_EVERY), CHECKPOINT_EVERY >= 1
)
############################################
library(MASS)
Greenberg.RR.data=function(nn,pp,cc,beta.1, alpha.00,BB)
{
  ##xx=runif(nn,-1,1)
  #xx=sample(c(0,1,2),nn,replace=T, prob=c(0.5,0.3,0.2))
  xx=rbinom(nn,1, 0.5)
  ##zz=sample(c(0,1,2),nn,replace=T, prob=c(0.5,0.3,0.2))
  #zz=rbinom(nn,1, 0.5)
  #zz=sample(c(-1,0,1),nn,replace=T, prob=c(0.3,0.3,0.4))
  chi.1=cbind(1,xx)
  xz.mat=cbind(xx)
  alpha.0=alpha.00[c(1:(BB))]             #1,2,...,BB-1; alpha01,alpha02
  alpha.1=alpha.00[c((BB+1):(2*BB))]  # 1,2....,BB-1; alpha11, alpha12
  b1=chi.1%*%beta.1
  
  gz0.pom=xz.mat%*%alpha.0[BB]
  gz1.pom=xz.mat%*%alpha.1[BB]
  
  H1=1/(1+exp(-b1))
  
  n.obs=nrow(chi.1)
  HHZ.0=matrix(0,n.obs,BB)
  HHZ.1=matrix(0,n.obs,BB)
  HHZ.0[,1]=1/(1+exp(-alpha.0[1]-gz0.pom))
  HHZ.1[,1]=1/(1+exp(-alpha.1[1]-gz1.pom))
  
  for(kk in 2:(BB-1))
  {
    HHZ.0[,kk]=1/(1+exp(-alpha.0[kk]-gz0.pom))-1/(1+exp(-alpha.0[kk-1]-gz0.pom))
    HHZ.1[,kk]=1/(1+exp(-alpha.1[kk]-gz1.pom))-1/(1+exp(-alpha.1[kk-1]-gz1.pom))
  }
  HHZ.0[,BB]=1-1/(1+exp(-alpha.0[BB-1]-gz0.pom))
  HHZ.1[,BB]=1-1/(1+exp(-alpha.1[BB-1]-gz1.pom))
  
  Y1=rbinom(nn,1,H1)
  
  QQ.Z0=matrix(0,nn,BB)
  QQ.Z1=matrix(0,nn,BB)
  Z0=rep(0, nn)
  Z1=rep(0,nn)
  for(i in 1:nn)
  {
    Z0[i]=sample(c(1:BB), 1, replace=TRUE,HHZ.0[i,])
    QQ.Z0[i,Z0[i]]=1
    Z1[i]=sample(c(1:BB), 1, replace=TRUE,HHZ.1[i,])
    QQ.Z1[i,Z1[i]]=1
  }
  ZZ=Y1*Z1+(1-Y1)*Z0
  QQ.ZZ=Y1*QQ.Z1+(1-Y1)*QQ.Z0
  TT=rbinom(nn,1,pp)
  DD=rbinom(nn,1,cc)
  Y1.0=TT*Y1+(1-TT)*DD
  return(cbind(Y1,Y1.0,ZZ,chi.1)) #,QQ.ZZ))
}

###Conditional Expection of Y[i] given Y0[i] and ZZ[i] and X
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
  n.obs=nrow(chi.1)
  HHZ.0=matrix(0,n.obs,BB)
  HHZ.1=matrix(0,n.obs,BB)
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
  mstep.convergence=rep(NA_integer_,3)
  while(KS<EM_MAX_ITER && Err >= EM_TOL)
  {
    E.Y=CEY(Y1.0,ZZ,chi.1,pp,cc, Theta.hat,BB)
    beta.hat=Theta.hat[1:KK]
    alpha.0.hat=Theta.hat[(KK+1):(KK+(BB-1)+(KK-1))]
    alpha.1.hat=Theta.hat[(2*KK+BB-1):length(Theta.hat)]
    beta.1.est=nlminb(beta.hat, Greenberg.and.Z.reg1,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE) #, lower = -4, upper =4)
    alpha.0.est=nlminb(alpha.0.hat, Greenberg.and.Z.reg10,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,BB=BB,hessian=TRUE) #, lower = -4, upper =4)
    alpha.1.est=nlminb(alpha.1.hat, Greenberg.and.Z.reg11,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,BB=BB,hessian=TRUE) #, lower = -4, upper =4)
    mstep.convergence=c(beta.1.est$convergence,alpha.0.est$convergence,alpha.1.est$convergence)
    Theta.new.hat=c(beta.1.est$par,alpha.0.est$par,alpha.1.est$par)
    Err=sum(abs(Theta.hat-Theta.new.hat))/length(Theta.hat)
    if(Err<=EM_TOL && all(mstep.convergence==0)) {
      Est.Theta=Theta.new.hat
      index.convergence=0
    }
    Theta.hat=Theta.new.hat
    KS=KS+1
    #cat("KS=",KS," ","Err=",Err," ")
  }
  ans=c(Est.Theta,Err,index.convergence)
  attr(ans,"iterations")=KS
  attr(ans,"mstep.convergence")=mstep.convergence
  return(ans)
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
  truth=matrix(Theta.t,nrow=nrow(est1),ncol=length(Theta.t),byrow=TRUE)
  valid=is.finite(est1) & is.finite(est1.se) & est1.se>=0
  covered=matrix(NA_real_,nrow=nrow(est1),ncol=length(Theta.t))
  lower=est1-1.96*est1.se
  upper=est1+1.96*est1.se
  covered[valid]=(truth[valid]>=lower[valid] & truth[valid]<=upper[valid])*1
  return(covered)
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
  n.obs=nrow(chi.1)
  HHZ.0=matrix(0,n.obs,BB)
  HHZ.1=matrix(0,n.obs,BB)
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
########H0: independent model 
### Log.like.Y0
Log.like.Y0.H0=function(pp,cc,Y1.0,ZZ,chi.1,beta.H0)
{
  KK=ncol(chi.1)
  b1=chi.1%*%beta.H0
  H1=1/(1+exp(-b1))
  P.Y0.1=(pp+(1-pp)*cc)*H1+(1-pp)*cc*(1-H1)
  P.Y0.0=(1-pp)*(1-cc)*H1+(pp+(1-pp)*(1-cc))*(1-H1)
  Log.like.Y0=sum(Y1.0*log(P.Y0.1)+(1-Y1.0)*log(P.Y0.0))
  return(-Log.like.Y0)
}
## Log.like.ZZ
Log.like.ZZ.H0=function(Y1.0,ZZ,chi.1,alpha.H0,BB)
{
  KK=ncol(chi.1)
  xz.mat=chi.1[,2:KK]
  gz0.pom=as.matrix(xz.mat)%*%alpha.H0[BB:(BB+KK-2)]
  
  HHZ=matrix(0,nrow(chi.1),BB)
  
  HHZ[,1]=1/(1+exp(-alpha.H0[1]-gz0.pom))
  
  
  for(kk in 2:(BB-1))
  {
    HHZ[,kk]=1/(1+exp(-alpha.H0[kk]-gz0.pom))-1/(1+exp(-alpha.H0[kk-1]-gz0.pom))
  }
  HHZ[,BB]=1-1/(1+exp(-alpha.H0[BB-1]-gz0.pom))
  
  NZ=cbind(c(1:nrow(chi.1)),ZZ)
  Log.like.ZZH0=sum(log(HHZ[NZ]))
  return(-Log.like.ZZH0)
}
Estimation.H0=function(pp,cc,Y1.0,ZZ,chi.1,Theta.H0,BB)
{
  KK=ncol(chi.1)
  beta.H0=Theta.H0[1:KK]
  alpha.H0=Theta.H0[(KK+1):(BB-1+KK+(KK-1))]
  beta.1.H0=nlminb(beta.H0, Log.like.Y0.H0,gr=NULL,Y1.0=Y1.0, ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
  alpha.0.H0=nlminb(alpha.H0, Log.like.ZZ.H0,gr=NULL, Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1,BB=BB,hessian=TRUE)
  
  Theta.hat.H0=c(beta.1.H0$par,alpha.0.H0$par)
  return(Theta.hat.H0)
}

# 對每次 full-model EM 結果分開判定「參數估計成功」與「SE 成功」。
check.EM.result=function(Theta.est,hessian,em.error,em.flag)
{
  out=list(point.ok=FALSE,se.ok=FALSE,reason="",variance=NULL,
           min.variance.diag=NA_real_,hessian.rcond=NA_real_,
           min.information.eigenvalue=NA_real_)

  if(!isTRUE(em.flag==0) || !is.finite(em.error) || em.error>EM_TOL) {
    out$reason="EM_not_converged"; return(out)
  }
  if(any(!is.finite(Theta.est))) {
    out$reason="nonfinite_estimate"; return(out)
  }
  out$point.ok=TRUE
  if(any(!is.finite(hessian))) {
    out$reason="nonfinite_hessian"; return(out)
  }

  observed.info=-(hessian+t(hessian))/2
  info.eigen=tryCatch(eigen(observed.info,symmetric=TRUE,only.values=TRUE)$values,
                      error=function(e) NA_real_)
  out$min.information.eigenvalue=min(info.eigen)
  if(any(!is.finite(info.eigen)) || out$min.information.eigenvalue<=0) {
    out$reason="nonpositive_definite_information"; return(out)
  }
  out$hessian.rcond=rcond(observed.info)
  if(!is.finite(out$hessian.rcond) || out$hessian.rcond<HESSIAN_RCOND_TOL) {
    out$reason="singular_hessian"; return(out)
  }
  variance.par=tryCatch(solve(observed.info),error=function(e) NULL)
  if(is.null(variance.par)) {
    out$reason="hessian_inverse_failed"; return(out)
  }
  variance.diag=diag(variance.par)
  out$min.variance.diag=min(variance.diag)
  if(any(!is.finite(variance.diag)) || !all(variance.diag>0)) {
    out$reason="nonpositive_variance_diagonal"; return(out)
  }

  out$se.ok=TRUE
  out$reason="valid_SE"
  out$variance=variance.par
  out
}

################################
# The data-generating and fitted randomized-response mechanisms can differ.
# Equal true/fitted values evaluate a correctly specified design; unequal values
# evaluate sensitivity to device-probability misspecification/noncompliance.
pp.true=as.numeric(Sys.getenv("EM_P_TRUE", "0.50"))
cc.true=as.numeric(Sys.getenv("EM_C_TRUE", "0.25"))
pp=as.numeric(Sys.getenv("EM_P_FIT", as.character(pp.true)))
cc=as.numeric(Sys.getenv("EM_C_FIT", as.character(cc.true)))
if(any(!is.finite(c(pp.true,cc.true,pp,cc))) ||
   any(c(pp.true,cc.true,pp,cc)<0) || any(c(pp.true,cc.true,pp,cc)>1) ||
   pp.true<=0 || pp<=0) {
  stop("EM_P_TRUE, EM_C_TRUE, EM_P_FIT, and EM_C_FIT must be probabilities; p must be > 0.")
}
nn=as.integer(Sys.getenv("EM_SIM_N", "500"))
# 文章 Section 4 / Table 1 的 simulation settings。
# 注意：文章將 latent logistic 參數記為 alpha=(alpha_c,alpha_x)，
# 但本程式沿用舊變數名稱 beta.t。
PREVALENCE <- tolower(Sys.getenv("EM_PREVALENCE", "low"))
CASE_ID <- toupper(Sys.getenv("EM_CASE", "C1"))

beta.t <- switch(
  PREVALENCE,
  high=c(0,-1),    # Pr(Y=1) 約 0.4
  low=c(-1,-1),    # Pr(Y=1) 約 0.2；文章正文圖表採此設定
  stop("EM_PREVALENCE must be 'high' or 'low'.")
)

case.parameters <- list(
  C1=list(alpha.t0=c(-2,-0.5,1), alpha.t1=c(-2,-0.5,1),
          model="homogeneous"),
  C2=list(alpha.t0=c(-1,0.5,1), alpha.t1=c(-1,0.5,1),
          model="homogeneous"),
  C3=list(alpha.t0=c(-1,0.5,2), alpha.t1=c(-1,0.5,2),
          model="homogeneous"),
  C4=list(alpha.t0=c(-2,-0.5,1), alpha.t1=c(-1,0,2),
          model="heterogeneous"),
  C5=list(alpha.t0=c(-1,0.5,1), alpha.t1=c(-2,0,1.5),
          model="heterogeneous"),
  C6=list(alpha.t0=c(-1,0.5,2), alpha.t1=c(-0.5,1,2.5),
          model="heterogeneous")
)
if(!CASE_ID %in% names(case.parameters)) {
  stop("EM_CASE must be one of C1, C2, C3, C4, C5, C6.")
}
selected.case <- case.parameters[[CASE_ID]]
alpha.t0 <- selected.case$alpha.t0
alpha.t1 <- selected.case$alpha.t1
CASE_MODEL <- selected.case$model
cat("Article setting:",CASE_ID,"/",CASE_MODEL,"/",PREVALENCE,
    "prevalence / n=",nn,
    "/ true(p,c)=",paste0("(",pp.true,",",cc.true,")"),
    "/ fitted(p,c)=",paste0("(",pp,",",cc,")"),"\n")
BB=length(alpha.t0)-(length(beta.t)-1)+1   # 一般化的設定
alpha.t=c(alpha.t0,alpha.t1)
Theta.t=c(beta.t,alpha.t0,alpha.t1)
Theta.t.H0=c(beta.t,alpha.t0)
KT=as.integer(Sys.getenv("EM_SIM_KT", "1000"))

# 每筆資料使用1組資料驅動起始值與 N_STARTS-1 組隨機擾動起始值。
# 不再使用模擬真值作為 EM 起始值。
make.EM.starts=function(Y1.0,ZZ,chi.1,pp,cc,BB,n.starts=N_STARTS)
{
  KK=ncol(chi.1)
  latent.prevalence=(mean(Y1.0)-(1-pp)*cc)/pp
  latent.prevalence=min(max(latent.prevalence,0.05),0.95)
  beta.start=c(qlogis(latent.prevalence),rep(0,KK-1))
  cumulative.z=sapply(seq_len(BB-1),function(k) mean(ZZ<=k))
  cumulative.z=pmin(pmax(cumulative.z,0.05),0.95)
  ordinal.start=c(sort(qlogis(cumulative.z)),rep(0,KK-1))
  fixed.start=c(beta.start,ordinal.start,ordinal.start)
  starts=list(fixed.start)
  if(n.starts>1) {
    for(s in 2:n.starts) {
      beta.random=beta.start+rnorm(length(beta.start),0,0.5)
      ordinal0=ordinal.start+rnorm(length(ordinal.start),0,0.5)
      ordinal1=ordinal.start+rnorm(length(ordinal.start),0,0.5)
      ordinal0[seq_len(BB-1)]=sort(ordinal0[seq_len(BB-1)])
      ordinal1[seq_len(BB-1)]=sort(ordinal1[seq_len(BB-1)])
      starts[[s]]=c(beta.random,ordinal0,ordinal1)
    }
  }
  starts
}

fit.EM.multistart=function(Y1.0,ZZ,chi.1,pp,cc,BB)
{
  starts=make.EM.starts(Y1.0,ZZ,chi.1,pp,cc,BB)
  n.par=length(starts[[1]])
  fits=lapply(starts,function(start) {
    tryCatch(EM.alogoritmH1(Y1.0,ZZ,chi.1,pp,cc,start,BB),
             error=function(e) {
               ans=c(rep(NA_real_,n.par),NA_real_,1)
               attr(ans,"iterations")=NA_integer_
               ans
             })
  })
  loglik=sapply(fits,function(fit) {
    theta=fit[seq_len(n.par)]
    flag=fit[n.par+2]
    if(!isTRUE(flag==0) || any(!is.finite(theta))) return(NA_real_)
    tryCatch(Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,theta,BB),error=function(e) NA_real_)
  })
  converged=which(is.finite(loglik))
  if(length(converged)) {
    best=converged[which.max(loglik[converged])]
  } else {
    errors=sapply(fits,function(fit) fit[n.par+1])
    errors[!is.finite(errors)]=Inf
    best=if(all(is.infinite(errors))) 1L else which.min(errors)
  }
  coincide=NA
  max.loglik.diff=NA_real_
  max.parameter.diff=NA_real_
  if(length(converged)>=2) {
    best.theta=fits[[best]][seq_len(n.par)]
    loglik.diff=abs(loglik[converged]-loglik[best])
    parameter.diff=sapply(converged,function(s)
      max(abs(fits[[s]][seq_len(n.par)]-best.theta)))
    max.loglik.diff=max(loglik.diff)
    max.parameter.diff=max(parameter.diff)
    same.loglik=loglik.diff<=START_LOGLIK_TOL
    same.parameters=parameter.diff<=START_PARAMETER_TOL
    coincide=all(same.loglik & same.parameters)
  }
  list(fit=fits[[best]],n.starts=length(starts),n.converged=length(converged),
       solutions.coincide=coincide,best.loglik=loglik[best],
       max.loglik.diff=max.loglik.diff,max.parameter.diff=max.parameter.diff)
}

new.diagnostic.row=function(replication)
{
  data.frame(
    replication=replication,completed=FALSE,em_converged=FALSE,
    point_estimate_available=FALSE,se_available=FALSE,accepted=FALSE,
    reason=NA_character_,error.message=NA_character_,n.starts=N_STARTS,
    n.converged.starts=NA_integer_,start.solutions.coincide=NA,
    max.start.loglik.diff=NA_real_,max.start.parameter.diff=NA_real_,
    observed.loglik=NA_real_,em.iterations=NA_integer_,em.error=NA_real_,
    max.abs.estimate=NA_real_,max.SE=NA_real_,extreme.coefficient=FALSE,
    extreme.SE=FALSE,extreme.solution=FALSE,extreme.reason="not_available",
    min.variance.diag=NA_real_,hessian.rcond=NA_real_,
    min.information.eigenvalue=NA_real_,runtime.seconds=NA_real_,
    worker.pid=NA_integer_,stringsAsFactors=FALSE
  )
}

# 每個 worker 只負責一次 replication，不寫檔；主程式收回後才寫 checkpoint。
run.one.replication=function(replication,seed.state)
{
  started=proc.time()[["elapsed"]]
  diagnostic=new.diagnostic.row(replication)
  diagnostic$worker.pid=Sys.getpid()
  result=list(
    diagnostic=diagnostic,
    est1=rep(NA_real_,length(Theta.t)),se1=rep(NA_real_,length(Theta.t)),
    est0=rep(NA_real_,length(Theta.t.H0)),se0=rep(NA_real_,length(Theta.t.H0)),
    LR=rep(NA_real_,2),wald=rep(NA_real_,6),last.dataset=NULL
  )

  tryCatch({
    assign(".Random.seed",seed.state,envir=.GlobalEnv)
    GZ.data=Greenberg.RR.data(nn,pp.true,cc.true,beta.t,alpha.t,BB)
    Y1=GZ.data[,1]
    Y1.0=GZ.data[,2]
    ZZ=GZ.data[,3]
    chi.1=GZ.data[,4:ncol(GZ.data)]
    KK=ncol(chi.1)
    n.par=length(Theta.t)

    multistart.fit=fit.EM.multistart(Y1.0,ZZ,chi.1,pp,cc,BB)
    Est.Theta.all=multistart.fit$fit
    Theta.H=Est.Theta.all[seq_len(n.par)]
    em.error=Est.Theta.all[n.par+1]
    em.flag=Est.Theta.all[n.par+2]
    diagnostic$n.converged.starts=multistart.fit$n.converged
    diagnostic$start.solutions.coincide=multistart.fit$solutions.coincide
    diagnostic$max.start.loglik.diff=multistart.fit$max.loglik.diff
    diagnostic$max.start.parameter.diff=multistart.fit$max.parameter.diff
    diagnostic$observed.loglik=multistart.fit$best.loglik
    diagnostic$em.iterations=attr(Est.Theta.all,"iterations")
    diagnostic$em.error=em.error

    if(isTRUE(em.flag==0) && is.finite(em.error) && em.error<=EM_TOL &&
       all(is.finite(Theta.H))) {
      full.hessian=tryCatch(Df2.LogLike.Y0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H,BB),
                            error=function(e) matrix(NA_real_,n.par,n.par))
    } else {
      full.hessian=matrix(NA_real_,n.par,n.par)
    }
    convergence.check=check.EM.result(Theta.H,full.hessian,em.error,em.flag)
    diagnostic$em_converged=convergence.check$point.ok
    diagnostic$point_estimate_available=convergence.check$point.ok
    diagnostic$se_available=convergence.check$se.ok
    diagnostic$accepted=convergence.check$se.ok
    diagnostic$reason=convergence.check$reason
    diagnostic$min.variance.diag=convergence.check$min.variance.diag
    diagnostic$hessian.rcond=convergence.check$hessian.rcond
    diagnostic$min.information.eigenvalue=convergence.check$min.information.eigenvalue

    if(convergence.check$point.ok) {
      result$est1=Theta.H
      diagnostic$max.abs.estimate=max(abs(Theta.H))
      diagnostic$extreme.coefficient=diagnostic$max.abs.estimate>EXTREME_COEF_CUTOFF
    }

    if(convergence.check$se.ok) {
      variance.par=convergence.check$variance
      result$se1=sqrt(diag(variance.par))
      diagnostic$max.SE=max(result$se1)
      diagnostic$extreme.SE=diagnostic$max.SE>EXTREME_SE_CUTOFF

      reduced.result=tryCatch({
        Theta.H0=c(beta.t,alpha.t[seq_len(length(alpha.t0))])
        Est.Theta.H0=Estimation.H0(pp,cc,Y1.0,ZZ,chi.1,Theta.H0,BB)
        if(any(!is.finite(Est.Theta.H0))) stop("nonfinite reduced-model estimate")
        Theta.H0.1=c(Est.Theta.H0,Est.Theta.H0[(KK+1):(2*KK+(BB-2))])
        LR.stat=2*(Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H,BB)-
          Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H0.1,BB))
        beta.2.H0=optim(Est.Theta.H0[seq_len(length(beta.t))],Log.like.Y0.H0,gr=NULL,
          Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1,pp=pp,cc=cc,hessian=TRUE)
        alpha0.2.H0=optim(Est.Theta.H0[(length(beta.t)+1):length(Theta.H0)],
          Log.like.ZZ.H0,gr=NULL,Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1,BB=BB,hessian=TRUE)
        reduced.variance.diag=c(diag(solve(beta.2.H0$hessian)),
                                diag(solve(alpha0.2.H0$hessian)))
        if(any(!is.finite(reduced.variance.diag)) || any(reduced.variance.diag<=0))
          stop("invalid reduced-model variance")
        list(theta=Est.Theta.H0,se=sqrt(reduced.variance.diag),LR=LR.stat)
      },error=function(e) NULL)
      if(!is.null(reduced.result)) {
        result$est0=reduced.result$theta
        result$se0=reduced.result$se
        result$LR=c(reduced.result$LR,1-pchisq(reduced.result$LR,(KK-1)+(BB-1)))
      }

      A.matrix=cbind(matrix(0,BB-1,KK),diag(rep(1,BB-1)),
                     matrix(0,BB-1,KK-1),diag(rep(-1,BB-1)),matrix(0,BB-1,KK-1))
      B.matrix=cbind(matrix(0,KK-1,KK+BB-1),diag(rep(1,KK-1)),
                     matrix(0,KK-1,BB-1),diag(rep(-1,KK-1),nrow=KK-1))
      wald.int=as.numeric(t(A.matrix%*%Theta.H)%*%
        ginv(A.matrix%*%variance.par%*%t(A.matrix))%*%(A.matrix%*%Theta.H))
      wald.coef=as.numeric(t(B.matrix%*%Theta.H)%*%
        ginv(B.matrix%*%variance.par%*%t(B.matrix))%*%(B.matrix%*%Theta.H))
      Full.matrix=rbind(A.matrix,B.matrix)
      wald.joint=as.numeric(t(Full.matrix%*%Theta.H)%*%
        ginv(Full.matrix%*%variance.par%*%t(Full.matrix))%*%
        (Full.matrix%*%Theta.H))
      result$wald=c(wald.joint,1-pchisq(wald.joint,nrow(Full.matrix)),
                    wald.int,1-pchisq(wald.int,nrow(A.matrix)),
                    wald.coef,1-pchisq(wald.coef,nrow(B.matrix)))
    }

    if(convergence.check$point.ok) {
      diagnostic$extreme.solution=diagnostic$extreme.coefficient || diagnostic$extreme.SE
      diagnostic$extreme.reason=
        if(diagnostic$extreme.coefficient && diagnostic$extreme.SE) {
          "large_coefficient_and_SE"
        } else if(diagnostic$extreme.coefficient) {
          "large_coefficient"
        } else if(diagnostic$extreme.SE) {
          "large_SE"
        } else {
          "not_extreme"
        }
    }

    result$last.dataset=data.frame(Y1.0=Y1.0,Y1=Y1,X=GZ.data[,5],Z=ZZ)
    diagnostic$completed=TRUE
    result$diagnostic=diagnostic
    result
  },error=function(e) {
    diagnostic$completed=TRUE
    diagnostic$reason="worker_error"
    diagnostic$error.message=conditionMessage(e)
    result$diagnostic=diagnostic
    result
  }) -> result

  result$diagnostic$runtime.seconds=proc.time()[["elapsed"]]-started
  result
}

# 預設存在程式所在資料夾的 Greenberg_EM_results，macOS/Windows 皆可用。
source.path=tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
rscript.argument=grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if((is.null(source.path) || !length(source.path)) && length(rscript.argument))
  source.path=sub("^--file=","",rscript.argument[1])
script.path=tryCatch(normalizePath(source.path,winslash="/",mustWork=TRUE),
                     error=function(e) NA_character_)
script.dir=if(is.character(script.path) && length(script.path)==1 && !is.na(script.path))
  dirname(script.path) else getwd()
DEFAULT_OUTPUT_ROOT <- file.path(script.dir,"Greenberg_EM_results")
OUTPUT_ROOT <- Sys.getenv("EM_OUTPUT_DIR", DEFAULT_OUTPUT_ROOT)
SCENARIO_ID <- paste0(CASE_ID,"_",PREVALENCE,"_n",nn,"_KT",KT)
make.unique.run.dir <- function(root,scenario)
{
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  candidate=file.path(root,scenario)
  if(!dir.exists(candidate)) {
    dir.create(candidate,recursive=TRUE)
    return(candidate)
  }
  run.number=2L
  repeat {
    candidate=file.path(root,paste0(scenario,"__run_",run.number))
    if(!dir.exists(candidate)) {
      dir.create(candidate,recursive=TRUE)
      return(candidate)
    }
    run.number=run.number+1L
  }
}
requested.run.id=Sys.getenv("EM_RUN_ID",SCENARIO_ID)
resume.candidate=file.path(OUTPUT_ROOT,requested.run.id)
resume.checkpoint=file.path(resume.candidate,"checkpoint.rds")
resume.previous=file.path(resume.candidate,"checkpoint_previous.rds")
if(RESUME_RUN && dir.exists(resume.candidate) &&
   (file.exists(resume.checkpoint) || file.exists(resume.previous))) {
  RUN_OUTPUT_DIR=resume.candidate
} else {
  RUN_OUTPUT_DIR=make.unique.run.dir(OUTPUT_ROOT,requested.run.id)
}
RUN_OUTPUT_DIR <- normalizePath(RUN_OUTPUT_DIR,winslash="/",mustWork=TRUE)
CHECKPOINT_FILE=file.path(RUN_OUTPUT_DIR,"checkpoint.rds")
CHECKPOINT_PREVIOUS_FILE=file.path(RUN_OUTPUT_DIR,"checkpoint_previous.rds")
cat("Output directory:",RUN_OUTPUT_DIR,"\n")

est1=matrix(NA_real_,KT,length(Theta.t))
est1.se=matrix(NA_real_,KT,length(Theta.t))
est1.cp=matrix(NA_real_,KT,length(Theta.t))

est0=matrix(NA_real_,KT,length(Theta.t.H0))
est0.se=matrix(NA_real_,KT,length(Theta.t.H0))
est0.cp=matrix(NA_real_,KT,length(Theta.t.H0))

LR.test=matrix(NA_real_,KT,2)
#for(K1 in 1:KTR)
wald.test.result <- matrix(NA_real_, nrow=KT, ncol=6)
colnames(wald.test.result) <- c(
  "joint_stat", "joint_pval",
  "intercept_stat", "intercept_pval",
  "coefX_stat", "coefX_pval"
)

jj=0
KTR=KT # 固定嘗試 KT 次；失敗不補跑
KTR.0=0
convergence.log=do.call(rbind,lapply(seq_len(KT),new.diagnostic.row))
last.dataset=NULL
last.accepted.theta=NULL
last.accepted.hessian=NULL
last.accepted.variance=NULL

SIM_SEED=as.integer(Sys.getenv("EM_SIM_SEED","20260916"))
RNGkind("L'Ecuyer-CMRG")
set.seed(SIM_SEED)
replication.seeds=vector("list",KT)
stream.seed=.Random.seed
for(i in seq_len(KT)) {
  replication.seeds[[i]]=stream.seed
  stream.seed=parallel::nextRNGStream(stream.seed)
}

settings.signature=list(
  schema_version=3L,case=CASE_ID,case_model=CASE_MODEL,prevalence=PREVALENCE,
  beta=beta.t,alpha0=alpha.t0,alpha1=alpha.t1,p=pp,c=cc,n=nn,KT=KT,
  seed=SIM_SEED,EM_TOL=EM_TOL,EM_MAX_ITER=EM_MAX_ITER,
  HESSIAN_RCOND_TOL=HESSIAN_RCOND_TOL,EXTREME_COEF_CUTOFF=EXTREME_COEF_CUTOFF,
  EXTREME_SE_CUTOFF=EXTREME_SE_CUTOFF,N_STARTS=N_STARTS,
  START_LOGLIK_TOL=START_LOGLIK_TOL,START_PARAMETER_TOL=START_PARAMETER_TOL
)

read.checkpoint=function(primary,previous)
{
  for(path in c(primary,previous)) {
    if(file.exists(path)) {
      value=tryCatch(readRDS(path),error=function(e) NULL)
      if(!is.null(value)) return(value)
    }
  }
  NULL
}

save.checkpoint=function(status="running")
{
  checkpoint=list(
    signature=settings.signature,status=status,saved_at=as.character(Sys.time()),
    completed=which(convergence.log$completed),convergence_log=convergence.log,
    est1=est1,est1.se=est1.se,est0=est0,est0.se=est0.se,
    LR.test=LR.test,wald.test.result=wald.test.result,last.dataset=last.dataset
  )
  temporary=paste0(CHECKPOINT_FILE,".tmp")
  saveRDS(checkpoint,temporary)
  verified=tryCatch(readRDS(temporary),error=function(e) NULL)
  if(is.null(verified)) stop("Temporary checkpoint verification failed: ",temporary)
  if(file.exists(CHECKPOINT_FILE))
    file.copy(CHECKPOINT_FILE,CHECKPOINT_PREVIOUS_FILE,overwrite=TRUE)
  if(file.exists(CHECKPOINT_FILE)) file.remove(CHECKPOINT_FILE)
  if(!file.rename(temporary,CHECKPOINT_FILE))
    stop("Could not promote temporary checkpoint: ",temporary)
  invisible(NULL)
}

checkpoint=if(RESUME_RUN) read.checkpoint(CHECKPOINT_FILE,CHECKPOINT_PREVIOUS_FILE) else NULL
if(!is.null(checkpoint)) {
  if(!identical(checkpoint$signature,settings.signature)) {
    stop("Checkpoint settings do not match the current run. Use a different EM_RUN_ID or set EM_RESUME=false.")
  }
  convergence.log=checkpoint$convergence_log
  est1=checkpoint$est1
  est1.se=checkpoint$est1.se
  est0=checkpoint$est0
  est0.se=checkpoint$est0.se
  LR.test=checkpoint$LR.test
  wald.test.result=checkpoint$wald.test.result
  last.dataset=checkpoint$last.dataset
  cat("Resuming checkpoint: ",sum(convergence.log$completed),"/",KT,
      " replications already completed.\n",sep="")
}

# 進度條以已完成的固定模擬次數為基準，失敗也計入進度。
SHOW_PROGRESS <- tolower(Sys.getenv("EM_PROGRESS","true")) %in% c("true","1","yes")
show.progress <- function(done,target,n.point,n.se,last.reason,width=30L)
{
  if(!SHOW_PROGRESS) return(invisible(NULL))
  proportion=min(1,done/target)
  filled=floor(width*proportion)
  bar=paste0(strrep("=",filled),if(filled<width) ">" else "",
             strrep(" ",max(0,width-filled-1L)))
  cat(sprintf("\rProgress [%-30s] %6.2f%% | runs %d/%d | estimates %d | valid SE %d | last: %-35s",
              bar,100*proportion,done,target,n.point,n.se,last.reason))
  flush.console()
  invisible(NULL)
}
if(FALSE) { # 舊的單核迴圈保留作為參考，實際改由下方批次平行引擎執行。
show.progress(0,KT,0,0,"starting")
while(KTR.0<KTR)
{
  KTR.0=KTR.0+1
  GZ.data=Greenberg.RR.data(nn,pp.true,cc.true,beta.t, alpha.t,BB)
  Y1.0=GZ.data[,2]
  ZZ=GZ.data[,3]
  chi.1=GZ.data[,4:ncol(GZ.data)]
  KK=ncol(chi.1)
  b1=chi.1%*%beta.t
  gz0.pom=as.matrix(chi.1[,2,drop=FALSE])%*%alpha.t0[BB:(BB+KK-2)]
  gz1.pom=as.matrix(chi.1[,2,drop=FALSE])%*%alpha.t1[BB:(BB+KK-2)]
  H1=1/(1+exp(-b1))
  HZZ.0=matrix(0,nn,BB)
  HZZ.1=matrix(0,nn,BB)
  HZZ.0[,1]=1/(1+exp(-alpha.t0[1]-gz0.pom))
  HZZ.1[,1]=1/(1+exp(-alpha.t1[1]-gz1.pom))
  for(kk in 2:(BB-1))
  {
    HZZ.0[,kk]=1/(1+exp(-alpha.t0[kk]-gz0.pom))-1/(1+exp(-alpha.t0[kk-1]-gz0.pom))
    HZZ.1[,kk]=1/(1+exp(-alpha.t1[kk]-gz1.pom))-1/(1+exp(-alpha.t1[kk-1]-gz1.pom))
  }
  HZZ.0[,BB]=1-1/(1+exp(-alpha.t0[BB-1]-gz0.pom))
  HZZ.1[,BB]=1-1/(1+exp(-alpha.t1[BB-1]-gz1.pom))
  mean(H1)
  apply(HZZ.0,2,mean)
  apply(HZZ.1,2,mean)
  n.par=length(Theta.t)
  multistart.fit=fit.EM.multistart(Y1.0,ZZ,chi.1,pp,cc,BB)
  Est.Theta.all=multistart.fit$fit
  convergence.log$n.converged.starts[KTR.0]=multistart.fit$n.converged
  convergence.log$start.solutions.coincide[KTR.0]=multistart.fit$solutions.coincide
  convergence.log$max.start.loglik.diff[KTR.0]=multistart.fit$max.loglik.diff
  convergence.log$max.start.parameter.diff[KTR.0]=multistart.fit$max.parameter.diff
  Theta.H=Est.Theta.all[seq_len(n.par)]
  em.error=Est.Theta.all[n.par+1]
  em.flag=Est.Theta.all[n.par+2]
  convergence.log$em.iterations[KTR.0]=attr(Est.Theta.all,"iterations")
  convergence.log$em.error[KTR.0]=em.error

  if(isTRUE(em.flag==0) && is.finite(em.error) && em.error<=EM_TOL && all(is.finite(Theta.H))) {
    AA.11=tryCatch(Df2.LogLike.Y0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H,BB),
                   error=function(e) matrix(NA_real_,n.par,n.par))
  } else {
    AA.11=matrix(NA_real_,n.par,n.par)
  }
  convergence.check=check.EM.result(Theta.H,AA.11,em.error,em.flag)
  convergence.log$em_converged[KTR.0]=convergence.check$point.ok
  convergence.log$point_estimate_available[KTR.0]=convergence.check$point.ok
  convergence.log$se_available[KTR.0]=convergence.check$se.ok
  convergence.log$reason[KTR.0]=convergence.check$reason
  convergence.log$min.variance.diag[KTR.0]=convergence.check$min.variance.diag
  convergence.log$hessian.rcond[KTR.0]=convergence.check$hessian.rcond
  convergence.log$min.information.eigenvalue[KTR.0]=convergence.check$min.information.eigenvalue

  if(convergence.check$point.ok) {
    est1[KTR.0,]=Theta.H
    convergence.log$max.abs.estimate[KTR.0]=max(abs(Theta.H))
    convergence.log$extreme.coefficient[KTR.0]=
      convergence.log$max.abs.estimate[KTR.0]>EXTREME_COEF_CUTOFF
  }

  if(convergence.check$se.ok) {
    convergence.log$accepted[KTR.0]=TRUE
    variance.par=convergence.check$variance
    last.accepted.theta=Theta.H
    last.accepted.hessian=AA.11
    last.accepted.variance=variance.par
    jj=jj+1
    ##
    est1.se[KTR.0,]=sqrt(diag(variance.par))
    convergence.log$max.SE[KTR.0]=max(est1.se[KTR.0,])
    convergence.log$extreme.SE[KTR.0]=convergence.log$max.SE[KTR.0]>EXTREME_SE_CUTOFF
    reduced.result=tryCatch({
      Theta.H0=c(beta.t,alpha.t[seq_len(length(alpha.t0))])
      Est.Theta.H0=Estimation.H0(pp,cc,Y1.0,ZZ,chi.1,Theta.H0,BB)
      if(any(!is.finite(Est.Theta.H0))) stop("nonfinite reduced-model estimate")
      Theta.H0.1=c(Est.Theta.H0,Est.Theta.H0[(KK+1):(2*KK+(BB-2))])
      LR.stat=2*(Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H,BB)-
        Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H0.1,BB))
      beta.2.H0=optim(Est.Theta.H0[seq_len(length(beta.t))],Log.like.Y0.H0,gr=NULL,
        Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1,pp=pp,cc=cc,hessian=TRUE)
      alpha0.2.H0=optim(Est.Theta.H0[(length(beta.t)+1):length(Theta.H0)],
        Log.like.ZZ.H0,gr=NULL,Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1,BB=BB,hessian=TRUE)
      reduced.variance.diag=c(diag(solve(beta.2.H0$hessian)),
                              diag(solve(alpha0.2.H0$hessian)))
      if(any(!is.finite(reduced.variance.diag)) || any(reduced.variance.diag<=0))
        stop("invalid reduced-model variance")
      list(theta=Est.Theta.H0,se=sqrt(reduced.variance.diag),LR=LR.stat)
    },error=function(e) NULL)
    if(!is.null(reduced.result)) {
      est0[KTR.0,]=reduced.result$theta
      est0.se[KTR.0,]=reduced.result$se
      LR.test[KTR.0,]=c(reduced.result$LR,
                         1-pchisq(reduced.result$LR,(KK-1)+(BB-1)))
    }
    # === Wald Test 計算 ===
    A.matrix <- cbind(matrix(0, BB - 1, KK), diag(rep(1, BB - 1)), matrix(0, BB - 1, KK - 1), diag(rep(-1, BB - 1)), matrix(0, BB - 1, KK - 1))
    wald.int <- t(A.matrix %*% Theta.H) %*% ginv(A.matrix %*% variance.par %*% t(A.matrix)) %*% (A.matrix %*% Theta.H)
    wald.int.p <- 1 - pchisq(wald.int, nrow(A.matrix))
    
    B.matrix <- cbind(matrix(0, KK - 1, KK + BB - 1), diag(rep(1, KK - 1)), matrix(0, KK - 1, BB - 1), diag(rep(-1, KK - 1), nrow = KK - 1))
    wald.coef <- t(B.matrix %*% Theta.H) %*% ginv(B.matrix %*% variance.par %*% t(B.matrix)) %*% (B.matrix %*% Theta.H)
    wald.coef.p <- 1 - pchisq(wald.coef, nrow(B.matrix))

    # Joint Wald test for H0: cutpoints and slopes are equal across latent groups.
    # This is the same joint null hypothesis tested by the LR statistic above.
    Full.matrix <- rbind(A.matrix, B.matrix)
    wald.joint <- t(Full.matrix %*% Theta.H) %*%
      ginv(Full.matrix %*% variance.par %*% t(Full.matrix)) %*%
      (Full.matrix %*% Theta.H)
    wald.joint.p <- 1 - pchisq(wald.joint, nrow(Full.matrix))
    
    wald.test.result[KTR.0, ] <- c(
      as.numeric(wald.joint), as.numeric(wald.joint.p),
      as.numeric(wald.int), as.numeric(wald.int.p),
      as.numeric(wald.coef), as.numeric(wald.coef.p)
    )
  }
  if(convergence.check$point.ok) {
    convergence.log$extreme.solution[KTR.0]=
      convergence.log$extreme.coefficient[KTR.0] || convergence.log$extreme.SE[KTR.0]
    convergence.log$extreme.reason[KTR.0]=
      if(convergence.log$extreme.coefficient[KTR.0] && convergence.log$extreme.SE[KTR.0]) {
        "large_coefficient_and_SE"
      } else if(convergence.log$extreme.coefficient[KTR.0]) {
        "large_coefficient"
      } else if(convergence.log$extreme.SE[KTR.0]) {
        "large_SE"
      } else {
        "not_extreme"
      }
  }
  show.progress(KTR.0,KT,sum(convergence.log$point_estimate_available[seq_len(KTR.0)]),
                sum(convergence.log$se_available[seq_len(KTR.0)]),
                convergence.log$reason[KTR.0])
}
if(SHOW_PROGRESS) cat("\n")
}

apply.replication.result=function(value)
{
  i=value$diagnostic$replication
  convergence.log[i,] <<- value$diagnostic
  est1[i,] <<- value$est1
  est1.se[i,] <<- value$se1
  est0[i,] <<- value$est0
  est0.se[i,] <<- value$se0
  LR.test[i,] <<- value$LR
  wald.value=value$wald
  if(length(wald.value)!=ncol(wald.test.result)) {
    warning("Replication ",i," returned ",length(wald.value),
            " Wald values; expected ",ncol(wald.test.result),
            ". The Wald row was recorded as undefined.",call.=FALSE)
    wald.value=rep(NA_real_,ncol(wald.test.result))
  }
  wald.test.result[i,] <<- wald.value
  if(!is.null(value$last.dataset)) last.dataset <<- value$last.dataset
  invisible(NULL)
}

remaining=which(!convergence.log$completed)
completed.before=KT-length(remaining)
show.progress(completed.before,KT,sum(convergence.log$point_estimate_available),
              sum(convergence.log$se_available),
              if(completed.before) "resumed" else "starting")

cluster=NULL
if(length(remaining) && N_CORES>1) {
  worker.count=min(N_CORES,length(remaining))
  cluster=parallel::makeCluster(worker.count,type="PSOCK")
  on.exit(if(!is.null(cluster)) parallel::stopCluster(cluster),add=TRUE)
  worker.exports=c(
    "Greenberg.RR.data","CEY","Greenberg.and.Z.reg1","Greenberg.and.Z.reg10",
    "Greenberg.and.Z.reg11","EM.alogoritmH1","Log.like.GRY0Z","Df1.LogLike.Y0Z",
    "Df2.LogLike.Y0Z","Log.like.Y0.H0","Log.like.ZZ.H0","Estimation.H0",
    "check.EM.result","make.EM.starts","fit.EM.multistart","new.diagnostic.row",
    "run.one.replication","EM_TOL","EM_MAX_ITER","HESSIAN_RCOND_TOL",
    "EXTREME_COEF_CUTOFF","EXTREME_SE_CUTOFF","N_STARTS","START_LOGLIK_TOL",
    "START_PARAMETER_TOL","Theta.t","Theta.t.H0","beta.t","alpha.t","alpha.t0",
    "alpha.t1","pp.true","cc.true","pp","cc","nn","BB"
  )
  parallel::clusterExport(cluster,worker.exports,envir=environment())
  parallel::clusterEvalQ(cluster,library(MASS))
  cat("Parallel workers:",worker.count,"(PSOCK)\n")
} else if(length(remaining)) {
  cat("Parallel workers: 1 (serial fallback)\n")
}

if(length(remaining)) {
  batches=split(remaining,ceiling(seq_along(remaining)/CHECKPOINT_EVERY))
  for(batch in batches) {
    if(is.null(cluster)) {
      values=lapply(batch,function(i) run.one.replication(i,replication.seeds[[i]]))
    } else {
      values=parallel::parLapplyLB(cluster,batch,function(i,seeds)
        run.one.replication(i,seeds[[i]]),seeds=replication.seeds)
    }
    invisible(lapply(values,apply.replication.result))
    save.checkpoint("running")
    done=sum(convergence.log$completed)
    last.reason=tail(convergence.log$reason[convergence.log$completed],1)
    show.progress(done,KT,sum(convergence.log$point_estimate_available),
                  sum(convergence.log$se_available),last.reason)
  }
}
if(SHOW_PROGRESS) cat("\n")
if(!is.null(cluster)) {
  parallel::stopCluster(cluster)
  cluster=NULL
}
KTR.0=sum(convergence.log$completed)
jj=sum(convergence.log$se_available)
convergence.log=convergence.log[seq_len(KTR.0),,drop=FALSE]
write.csv(convergence.log,file.path(RUN_OUTPUT_DIR,"convergence_diagnostics.csv"),row.names=FALSE)
failure.reason.summary=as.data.frame(table(convergence.log$reason),stringsAsFactors=FALSE)
names(failure.reason.summary)=c("reason","count")
failure.reason.summary$proportion=failure.reason.summary$count/nrow(convergence.log)
convergence.summary=data.frame(
  metric=c("total_replications","EM_converged_with_point_estimate","valid_standard_error",
           "point_estimation_failure","SE_failure_given_point_estimate",
           "extreme_solutions_all_runs","extreme_solutions_among_converged",
           "multi_start_solutions_coincide"),
  count=c(KT,sum(convergence.log$point_estimate_available),sum(convergence.log$se_available),
          sum(!convergence.log$point_estimate_available),
          sum(convergence.log$point_estimate_available & !convergence.log$se_available),
          sum(convergence.log$extreme.solution),sum(convergence.log$extreme.solution),
          sum(convergence.log$start.solutions.coincide %in% TRUE)),
  denominator=c(KT,KT,KT,KT,sum(convergence.log$point_estimate_available),KT,
                sum(convergence.log$point_estimate_available),
                sum(!is.na(convergence.log$start.solutions.coincide))),stringsAsFactors=FALSE
)
convergence.summary$proportion=ifelse(convergence.summary$denominator>0,
                                      convergence.summary$count/convergence.summary$denominator,NA_real_)
write.csv(convergence.summary,file.path(RUN_OUTPUT_DIR,"convergence_summary.csv"),row.names=FALSE)
write.csv(failure.reason.summary,file.path(RUN_OUTPUT_DIR,"failure_reason_summary.csv"),row.names=FALSE)
simulation.settings=data.frame(
  setting=c("article_case","case_model","prevalence","logistic_alpha",
            "ordinal_Y0_beta_gamma","ordinal_Y1_beta_gamma",
            "p_true","c_true","p_fitted","c_fitted",
            "seed","n","fixed_replications_no_replacement","output_directory",
            "EM_tolerance","EM_max_iterations","Hessian_rcond_tolerance",
            "extreme_coefficient_cutoff","extreme_SE_cutoff",
            "extreme_solutions_removed_from_primary_summary",
            "extreme_solution_reporting","number_of_EM_starts",
            "multi_start_loglik_tolerance","multi_start_parameter_tolerance",
            "parallel_workers","checkpoint_every","resume_enabled",
            "create_plots_after_simulation","run_id"),
  value=c(CASE_ID,CASE_MODEL,PREVALENCE,paste(beta.t,collapse=","),
          paste(alpha.t0,collapse=","),paste(alpha.t1,collapse=","),
          pp.true,cc.true,pp,cc,
          Sys.getenv("EM_SIM_SEED","20260916"),nn,KT,RUN_OUTPUT_DIR,EM_TOL,EM_MAX_ITER,
          HESSIAN_RCOND_TOL,EXTREME_COEF_CUTOFF,EXTREME_SE_CUTOFF,"FALSE",
          "including_and_excluding_for_full_and_reduced_models",N_STARTS,
          START_LOGLIK_TOL,START_PARAMETER_TOL,N_CORES,CHECKPOINT_EVERY,RESUME_RUN,
          CREATE_PLOTS,
          basename(RUN_OUTPUT_DIR))
)
write.csv(simulation.settings,file.path(RUN_OUTPUT_DIR,"simulation_settings.csv"),row.names=FALSE)
if(jj==0) warning("No replication produced a valid standard error; all attempts remain in the output files.")
write.csv(wald.test.result,file.path(RUN_OUTPUT_DIR,"wald_replications.csv"),row.names=FALSE)
# 顯著水準 alpha = 0.05
alpha <- 0.05
sig_counts <- colSums(wald.test.result[, c(2, 4, 6),drop=FALSE] < alpha, na.rm = TRUE)
names(sig_counts) <- c(
  "Joint sig. count", "Intercept sig. count", "X coefficient sig. count"
)
est1.cp=CP.95(est1,est1.se,Theta.t)
est0.cp=CP.95(est0,est0.se,Theta.t.H0)
finite.mean=function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
finite.sd=function(x) if(sum(is.finite(x))>=2) sd(x[is.finite(x)]) else NA_real_
parameter.summary=function(est,se,coverage,truth,parameter.names)
{
  result=lapply(seq_along(truth),function(k) {
    estimate.ok=is.finite(est[,k])
    se.ok=estimate.ok & is.finite(se[,k])
    errors=est[estimate.ok,k]-truth[k]
    data.frame(
      parameter=parameter.names[k],true_parameter=truth[k],n_estimate=sum(estimate.ok),
      mean_estimate=if(length(errors)) mean(est[estimate.ok,k]) else NA_real_,
      bias=if(length(errors)) mean(errors) else NA_real_,
      empirical_SD=finite.sd(est[,k]),
      RMSE=if(length(errors)) sqrt(mean(errors^2)) else NA_real_,
      n_valid_SE=sum(se.ok),mean_ASE=finite.mean(se[se.ok,k]),
      coverage_95=finite.mean(coverage[se.ok,k]),stringsAsFactors=FALSE
    )
  })
  out=do.call(rbind,result)
  metrics=c("mean_estimate","bias","empirical_SD","RMSE","mean_ASE","coverage_95")
  out[metrics]=lapply(out[metrics],function(x) round(x,4))
  out
}

full.parameter.names=c("alpha_c","alpha_x","beta_01","beta_02","gamma_0",
                       "beta_11","beta_12","gamma_1")
reduced.parameter.names=c("alpha_c","alpha_x","beta_1","beta_2","gamma")
output.1=parameter.summary(est1,est1.se,est1.cp,Theta.t,full.parameter.names)
output.2=parameter.summary(est0,est0.se,est0.cp,Theta.t.H0,reduced.parameter.names)

# 老師要求：主要摘要保留極端解，另提供排除極端解後的敏感性分析。
# 極端解定義沿用上方預先設定的門檻；不以真值或估計誤差篩選。
nonextreme.rows=which(
  convergence.log$point_estimate_available %in% TRUE &
  !(convergence.log$extreme.solution %in% TRUE)
)
output.1.nonextreme=parameter.summary(
  est1[nonextreme.rows,,drop=FALSE],
  est1.se[nonextreme.rows,,drop=FALSE],
  est1.cp[nonextreme.rows,,drop=FALSE],
  Theta.t,full.parameter.names
)
output.1.comparison=rbind(
  data.frame(analysis="all_available_including_extreme",output.1,
             stringsAsFactors=FALSE),
  data.frame(analysis="excluding_extreme_solutions",output.1.nonextreme,
             stringsAsFactors=FALSE)
)

# Reduced model 也獨立套用相同的極端解門檻，並同時輸出保留與排除極端解的結果。
# 此處不沿用 full model 的 extreme.solution，避免因 full model 極端而誤刪
# reduced model 本身穩定的估計結果。
reduced.point.available=apply(est0,1,function(x) all(is.finite(x)))
reduced.se.available=apply(est0.se,1,function(x) all(is.finite(x) & x>0))
reduced.extreme.coefficient=reduced.point.available &
  apply(abs(est0)>EXTREME_COEF_CUTOFF,1,any,na.rm=TRUE)
reduced.extreme.SE=reduced.se.available &
  apply(est0.se>EXTREME_SE_CUTOFF,1,any,na.rm=TRUE)
reduced.extreme.solution=reduced.extreme.coefficient | reduced.extreme.SE
reduced.nonextreme.rows=which(reduced.point.available & !reduced.extreme.solution)
output.2.nonextreme=parameter.summary(
  est0[reduced.nonextreme.rows,,drop=FALSE],
  est0.se[reduced.nonextreme.rows,,drop=FALSE],
  est0.cp[reduced.nonextreme.rows,,drop=FALSE],
  Theta.t.H0,reduced.parameter.names
)
output.2.comparison=rbind(
  data.frame(analysis="all_available_including_extreme",output.2,
             stringsAsFactors=FALSE),
  data.frame(analysis="excluding_extreme_solutions",output.2.nonextreme,
             stringsAsFactors=FALSE)
)

reduced.extreme.solution.summary=data.frame(
  metric=c(
    "total_replications","point_estimate_available","valid_standard_error",
    "extreme_solution","nonextreme_point_estimate",
    "valid_SE_and_nonextreme","extreme_coefficient","extreme_SE"
  ),
  count=c(
    KT,sum(reduced.point.available),sum(reduced.se.available),
    sum(reduced.extreme.solution),length(reduced.nonextreme.rows),
    sum(reduced.se.available & !reduced.extreme.solution),
    sum(reduced.extreme.coefficient),sum(reduced.extreme.SE)
  ),
  denominator=rep(KT,8),
  stringsAsFactors=FALSE
)
reduced.extreme.solution.summary$proportion=
  reduced.extreme.solution.summary$count/reduced.extreme.solution.summary$denominator
reduced.extreme.solution.summary$definition=c(
  "fixed Monte Carlo replications; failures are not replaced",
  "finite converged reduced-model point estimate",
  "finite positive reduced-model SE for every parameter",
  paste0("max(abs(estimate)) > ",EXTREME_COEF_CUTOFF,
         " or max(SE) > ",EXTREME_SE_CUTOFF),
  "reduced-model point estimate available and not extreme",
  "valid reduced-model SE and not extreme",
  paste0("max(abs(estimate)) > ",EXTREME_COEF_CUTOFF),
  paste0("max(SE) > ",EXTREME_SE_CUTOFF)
)

reduced.extreme.by.parameter=do.call(rbind,lapply(seq_along(reduced.parameter.names),function(k) {
  estimate.extreme=is.finite(est0[,k]) & abs(est0[,k])>EXTREME_COEF_CUTOFF
  SE.extreme=is.finite(est0.se[,k]) & est0.se[,k]>EXTREME_SE_CUTOFF
  either.extreme=estimate.extreme | SE.extreme
  data.frame(
    parameter=reduced.parameter.names[k],
    estimate_cutoff=EXTREME_COEF_CUTOFF,
    n_extreme_estimate=sum(estimate.extreme),
    SE_cutoff=EXTREME_SE_CUTOFF,
    n_extreme_SE=sum(SE.extreme),
    n_extreme_either=sum(either.extreme),
    proportion_among_all_replications=sum(either.extreme)/KT,
    stringsAsFactors=FALSE
  )
}))

extreme.solution.summary=data.frame(
  metric=c(
    "total_replications","point_estimate_available","valid_standard_error",
    "extreme_solution","nonextreme_point_estimate",
    "valid_SE_and_nonextreme","extreme_coefficient","extreme_SE"
  ),
  count=c(
    KT,
    sum(convergence.log$point_estimate_available %in% TRUE),
    sum(convergence.log$se_available %in% TRUE),
    sum(convergence.log$extreme.solution %in% TRUE),
    length(nonextreme.rows),
    sum(convergence.log$se_available %in% TRUE &
          !(convergence.log$extreme.solution %in% TRUE)),
    sum(convergence.log$extreme.coefficient %in% TRUE),
    sum(convergence.log$extreme.SE %in% TRUE)
  ),
  denominator=rep(KT,8),
  stringsAsFactors=FALSE
)
extreme.solution.summary$proportion=
  extreme.solution.summary$count/extreme.solution.summary$denominator
extreme.solution.summary$definition=c(
  "fixed Monte Carlo replications; failures are not replaced",
  "finite converged full-model point estimate",
  "positive-definite observed information and finite positive SE",
  paste0("max(abs(estimate)) > ",EXTREME_COEF_CUTOFF,
         " or max(SE) > ",EXTREME_SE_CUTOFF),
  "point estimate available and not an extreme solution",
  "valid full-model SE and not an extreme solution",
  paste0("max(abs(estimate)) > ",EXTREME_COEF_CUTOFF),
  paste0("max(SE) > ",EXTREME_SE_CUTOFF)
)

extreme.by.parameter=do.call(rbind,lapply(seq_along(full.parameter.names),function(k) {
  estimate.extreme=is.finite(est1[,k]) & abs(est1[,k])>EXTREME_COEF_CUTOFF
  SE.extreme=is.finite(est1.se[,k]) & est1.se[,k]>EXTREME_SE_CUTOFF
  either.extreme=estimate.extreme | SE.extreme
  data.frame(
    parameter=full.parameter.names[k],
    estimate_cutoff=EXTREME_COEF_CUTOFF,
    n_extreme_estimate=sum(estimate.extreme),
    SE_cutoff=EXTREME_SE_CUTOFF,
    n_extreme_SE=sum(SE.extreme),
    n_extreme_either=sum(either.extreme),
    proportion_among_all_replications=sum(either.extreme)/KT,
    stringsAsFactors=FALSE
  )
}))

# 逐次估計與參數摘要輸出。
full.replications=data.frame(convergence.log,est1,est1.se,est1.cp,
                             check.names=FALSE)
names(full.replications)=c(names(convergence.log),paste0("estimate_",full.parameter.names),
                           paste0("SE_",full.parameter.names),
                           paste0("covered_",full.parameter.names))
reduced.replications=data.frame(
                                replication=seq_len(nrow(est0)),
                                point_estimate_available=reduced.point.available,
                                se_available=reduced.se.available,
                                extreme_coefficient=reduced.extreme.coefficient,
                                extreme_SE=reduced.extreme.SE,
                                extreme_solution=reduced.extreme.solution,
                                est0,est0.se,est0.cp,
                                check.names=FALSE)
names(reduced.replications)=c("replication","point_estimate_available","se_available",
                              "extreme_coefficient","extreme_SE","extreme_solution",
                              paste0("estimate_",reduced.parameter.names),
                              paste0("SE_",reduced.parameter.names),
                              paste0("covered_",reduced.parameter.names))
LR.output=data.frame(replication=seq_len(nrow(LR.test)),LR_statistic=LR.test[,1],
                     LR_p_value=LR.test[,2])
write.csv(output.1,file.path(RUN_OUTPUT_DIR,"parameter_summary_full.csv"),row.names=FALSE)
write.csv(output.1.nonextreme,
          file.path(RUN_OUTPUT_DIR,"parameter_summary_full_nonextreme.csv"),row.names=FALSE)
write.csv(output.1.comparison,
          file.path(RUN_OUTPUT_DIR,"parameter_summary_full_comparison.csv"),row.names=FALSE)
write.csv(output.2,file.path(RUN_OUTPUT_DIR,"parameter_summary_reduced.csv"),row.names=FALSE)
write.csv(output.2.nonextreme,
          file.path(RUN_OUTPUT_DIR,"parameter_summary_reduced_nonextreme.csv"),row.names=FALSE)
write.csv(output.2.comparison,
          file.path(RUN_OUTPUT_DIR,"parameter_summary_reduced_comparison.csv"),row.names=FALSE)
write.csv(extreme.solution.summary,
          file.path(RUN_OUTPUT_DIR,"extreme_solution_summary.csv"),row.names=FALSE)
write.csv(extreme.by.parameter,
          file.path(RUN_OUTPUT_DIR,"extreme_by_parameter.csv"),row.names=FALSE)
write.csv(reduced.extreme.solution.summary,
          file.path(RUN_OUTPUT_DIR,"extreme_solution_summary_reduced.csv"),row.names=FALSE)
write.csv(reduced.extreme.by.parameter,
          file.path(RUN_OUTPUT_DIR,"extreme_by_parameter_reduced.csv"),row.names=FALSE)
write.csv(full.replications,file.path(RUN_OUTPUT_DIR,"replications_full.csv"),row.names=FALSE)
write.csv(reduced.replications,file.path(RUN_OUTPUT_DIR,"replications_reduced.csv"),row.names=FALSE)
write.csv(LR.output,file.path(RUN_OUTPUT_DIR,"likelihood_ratio_replications.csv"),row.names=FALSE)
if(!is.null(last.dataset)) {
  DD=last.dataset
  DD$Y1.0=as.factor(DD$Y1.0)
  DD$Y1=as.factor(DD$Y1)
  DD$X=as.factor(DD$X)
  if(interactive() && requireNamespace("gtsummary",quietly=TRUE) &&
     requireNamespace("dplyr",quietly=TRUE)) {
    print(gtsummary::tbl_summary(DD,by="Z"))
  }
  write.csv(DD,file.path(RUN_OUTPUT_DIR,"last_simulated_dataset.csv"),row.names=FALSE)
}

rate.summary=function(p.values,label)
{
  valid=is.finite(p.values)
  rate=if(any(valid)) mean(p.values[valid]<=0.05) else NA_real_
  data.frame(metric=label,value=rate,n_valid=sum(valid),n_total=length(p.values),
             MCSE=if(is.finite(rate)) sqrt(rate*(1-rate)/sum(valid)) else NA_real_)
}
test.performance=rbind(
  rate.summary(LR.test[,2],if(CASE_MODEL=="homogeneous") "LR_type_I_error" else "LR_power"),
  rate.summary(wald.test.result[,"joint_pval"],
               if(CASE_MODEL=="homogeneous") "Wald_joint_type_I_error" else "Wald_joint_power"),
  rate.summary(wald.test.result[,"intercept_pval"],"Wald_intercept_rejection_rate"),
  rate.summary(wald.test.result[,"coefX_pval"],"Wald_slope_rejection_rate")
)
test.performance$undefined_rate=1-test.performance$n_valid/test.performance$n_total
write.csv(test.performance,file.path(RUN_OUTPUT_DIR,"test_performance_summary.csv"),row.names=FALSE)
valid.wald.rows=which(apply(wald.test.result,1,function(x) all(is.finite(x))))
if(length(valid.wald.rows)) {
  last.wald=data.frame(test=colnames(wald.test.result),
                       value=as.numeric(wald.test.result[tail(valid.wald.rows,1),]))
} else {
  last.wald=data.frame(test="none",value=NA_real_)
}
write.csv(last.wald,file.path(RUN_OUTPUT_DIR,"last_accepted_wald_tests.csv"),row.names=FALSE)

# Bias 圖：使用已輸出的摘要表後處理，不會重新執行模擬。
# 圖中同時呈現保留與排除極端解的 Bias，誤差線為 Bias 的 95% Monte Carlo interval。
bias.plot.script=file.path(script.dir,"Greenberg_EM_bias_plots.R")
if(CREATE_PLOTS) {
  if(file.exists(bias.plot.script)) {
    tryCatch({
      source(bias.plot.script,local=TRUE)
      create_bias_plots(RUN_OUTPUT_DIR)
    },error=function(e) {
      warning("Plots were not created: ",conditionMessage(e),call.=FALSE)
    })
  } else {
    warning("Plotting script not found: ",bias.plot.script,call.=FALSE)
  }
} else {
  cat("Plot generation is disabled (EM_CREATE_PLOTS=false). ",
      "Existing results can be plotted later with Greenberg_EM_bias_plots.R.\n",
      sep="")
}

saveRDS(list(
  settings=simulation.settings,
  convergence_log=convergence.log,
  convergence_summary=convergence.summary,
  failure_reason_summary=failure.reason.summary,
  full_estimates=est1,full_standard_errors=est1.se,full_coverage=est1.cp,
  reduced_estimates=est0,reduced_standard_errors=est0.se,reduced_coverage=est0.cp,
  likelihood_ratio=LR.test,wald_replications=wald.test.result,
  parameter_summary_full=output.1,
  parameter_summary_full_nonextreme=output.1.nonextreme,
  parameter_summary_full_comparison=output.1.comparison,
  parameter_summary_reduced=output.2,
  parameter_summary_reduced_nonextreme=output.2.nonextreme,
  parameter_summary_reduced_comparison=output.2.comparison,
  extreme_solution_summary=extreme.solution.summary,
  extreme_by_parameter=extreme.by.parameter,
  extreme_solution_summary_reduced=reduced.extreme.solution.summary,
  extreme_by_parameter_reduced=reduced.extreme.by.parameter,
  test_performance=test.performance,last_accepted_wald_tests=last.wald
),file.path(RUN_OUTPUT_DIR,"complete_results.rds"))
save.checkpoint("complete")

cat("Completed. Fixed replications: ",KT,
    "; point estimates: ",sum(convergence.log$point_estimate_available),
    "; valid SE: ",sum(convergence.log$se_available),
    "; extreme solutions: ",sum(convergence.log$extreme.solution),
    "; nonextreme point estimates: ",length(nonextreme.rows),
    "; reduced-model extreme solutions: ",sum(reduced.extreme.solution),
    "; reduced-model nonextreme point estimates: ",length(reduced.nonextreme.rows),
    "; workers: ",min(N_CORES,KT),
    ".\nResults saved to:\n",RUN_OUTPUT_DIR,"\n",sep="")
