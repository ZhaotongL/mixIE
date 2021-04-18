#' Fixed effect IVW model
#'
#' @import stats
#' @keywords internal
#'
mr_ivw_fe <- function (b_exp, b_out, se_exp, se_out)
{
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
          !is.na(se_out)) < 2)
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  ivwlm <- stats::lm(b_out ~ -1 + b_exp, weights = 1/se_out^2)
  ivw.res <- summary(ivwlm)
  b <- ivw.res$coef["b_exp", "Estimate"]
  se <- ivw.res$coef["b_exp", "Std. Error"]/ivw.res$sigma
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail = FALSE)
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), sig=ivw.res$sigma))
}

#' Egger regression
#'
#' @import stats
#' @keywords internal
#'
mr_egger_ll <-function (b_exp, b_out, se_exp, se_out,flip=0)
{
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  nulllist <- list(b = NA, se = NA, pval = NA, nsnp = NA, b_i = NA,
                   se_i = NA, pval_i = NA,sig=NA)
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) &
          !is.na(se_out)) < 3) {
    return(nulllist)
  }
  if(flip==1){
    sign0 <- function(x) {
      x[x == 0] <- 1
      return(sign(x))
    }
    to_flip <- sign0(b_exp) == -1
    b_out = b_out * sign0(b_exp)
    b_exp = abs(b_exp)
  }
  mod <- stats::lm(b_out ~ b_exp, weights = 1/se_out^2)
  smod <- summary(mod)
  if (nrow(stats::coefficients(smod)) > 1) {
    b <- stats::coefficients(smod)[2, 1]
    se <- stats::coefficients(smod)[2, 2]/min(1, smod$sigma)
    pval <- 2 * stats::pt(abs(b/se), length(b_exp) - 2, lower.tail = FALSE)
    b_i <- stats::coefficients(smod)[1, 1]
    se_i <- stats::coefficients(smod)[1, 2]/min(1, smod$sigma)
    pval_i <- 2 * stats::pt(abs(b_i/se_i), length(b_exp) - 2, lower.tail = FALSE)
  }
  else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), sig=smod$sigma,
              b_i = b_i, se_i = se_i, pval_i = pval_i))
}

#' sum finite
#'
#' sum over only finite number
#'
#' @keywords internal
#'
sum.finite <- function(x) {
  sum(x[is.finite(x)])
}

#' Likelihood of valid IVs
#'
#' calculate the likelihood of valid IVs up to a constant
#'
#' @keywords internal
#'
L0_i <- function(b_exp,b_out,se_exp,se_out,
                 theta){
  L0 =     1/sqrt(2*pi*se_out^2) * exp(-(b_out-theta*b_exp)^2/(2*se_out^2))
  return(L0)
}

#' Likelihood of invalid IVs
#'
#' calculate the likelihood of invalid IVs up to a constant
#'
#' @keywords internal
#'

L1_i <- function(b_exp,b_out,se_exp,se_out,
                 theta,r,c){
  L1 =     1/sqrt(2*pi*c*se_out^2) * exp(-(b_out-theta*b_exp-r)^2/(2*(c*se_out^2)))
  return(L1)
}

#' posterior probability
#'
#' calculate the posterior probability in the E step
#'
#' @keywords internal
#'
tau_i <- function(b_exp,b_out,se_exp,se_out,
                  theta,r,c,p){
  L0 = L0_i(b_exp,b_out,se_exp,se_out,theta) * (1-p)
  L1 = L1_i(b_exp,b_out,se_exp,se_out,theta,r,c)  * p
  tau_i0 = L0/(L0+L1)
  tau_i1 = L1/(L0+L1)
  return(list(tau_i0=tau_i0,tau_i1=tau_i1))
}

#' information matrix
#'
#' calculate the informationmatrix based on the complete data log-likelihood for standard error
#'
#' @keywords internal
#'
B_Complete <- function(b_exp,b_out,se_exp,se_out,
                       theta,r,c,p,tau_0,tau_1){
  B = matrix(0,nrow = 4,ncol = 4)
  B[1,1] = sum(tau_1*b_exp^2/(c*se_out^2) +
                 tau_0*b_exp^2/se_out^2)
  B[1,2] = sum(tau_1*b_exp/(c*se_out^2))
  B[1,3] = sum(tau_1*b_exp*(b_out-theta*b_exp-r)/(c^2*se_out^2))
  B[2,2] = sum(tau_1/(c*se_out^2))
  B[2,3] = sum(tau_1*(b_out-theta*b_exp-r)/(c^2*se_out^2))
  B[3,3] = sum(tau_1*((b_out-theta*b_exp-r)^2/(c^3*se_out^2)-1/(2*c^2)))
  B[4,4] = sum(tau_1/p^2+tau_0/(1-p)^2)
  B[lower.tri(B)] = t(B)[lower.tri(B)]
  return(B)
}

#' gradient
#'
#' calculate the gradient based on the complete data log-likelihood for standard error
#'
#' @keywords internal
#'

S_Complete <- function(b_exp,b_out,se_exp,se_out,
                       theta,r,c,p,tau_0,tau_1){
  partial_theta = sum( tau_0 * b_exp*(b_out-theta*b_exp)/se_out^2 +
                         tau_1 * b_exp*(b_out-theta*b_exp-r)/(c*se_out^2) )
  partial_r = sum( tau_1 * (b_out-theta*b_exp-r)/(c*se_out^2) )
  partial_c = sum(1/2 * tau_1 * ( (b_out-theta*b_exp-r)^2/(c^2*se_out^2) - 1/c))
  partial_p = sum(tau_1/p - tau_0/(1-p))
  Sc = matrix(c(partial_theta,partial_r,partial_c,partial_p),ncol=1)
  return( Sc %*% t(Sc))

}

#' EM algorithm
#'
#' For starting the CEM algorithm
#'
#' @import HelpersMG stats
#' @keywords internal
#'

EM <-function(b_exp,b_out,se_exp,se_out,
              initial_theta,initial_r,initial_c,initial_p,
              maxit,n){
  m = length(b_exp)
  theta = initial_theta
  theta_old = theta-1
  r = initial_r
  c = initial_c
  p = initial_p
  niter = 0
  while(abs(theta-theta_old)>=1e-6 & niter<maxit){
    niter = niter + 1
    theta_old = theta
    ## E-step
    tau_vec = tau_i(b_exp,b_out,se_exp,se_out,theta,r,c,p)
    tau_0_vec = tau_vec$tau_i0
    tau_1_vec = tau_vec$tau_i1

    ## M-step
    p = sum.finite(tau_1_vec)/m
    r = sum.finite(tau_1_vec*(b_out-theta*b_exp)/se_out^2) /
      sum.finite(tau_1_vec/se_out^2)
    c = max(1,sum.finite(tau_1_vec*(b_out-theta*b_exp-r)^2 / se_out^2) /
              sum.finite(tau_1_vec))
    theta =
      sum.finite(tau_1_vec*(b_out - r)*b_exp / (c*se_out^2) +
                   tau_0_vec*b_out*b_exp/se_out^2) /
      sum.finite(tau_1_vec*b_exp^2 / (c*se_out^2) +
                   tau_0_vec*b_exp^2/se_out^2)
    if(is.na(theta)) break
  }
  if(is.na(theta) || any(is.na(tau_1_vec))){
    warning('please choose another starting value or increase maxit!')
    return(list(theta=NaN,se=NaN,pval=NaN,
                r=NaN,se_r=NaN,pval_r=NaN,
                c=NaN,p=NaN,
                lobs=NaN,lcom=NaN,AIC=NaN,BIC=NaN,
                niter=niter,
                tau_1_vec=tau_1_vec))
  }
  if(p>1) {
    p=1
    tau_1_vec = rep(1,m)
    tau_0_vec = rep(0,m)
  }

  B_complete = B_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,p,tau_0_vec,tau_1_vec)
  B_missing = S_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,p,tau_0_vec,tau_1_vec)
  FIM = B_complete - B_missing
  if(is.finite(FIM[4,4])){
    se = HelpersMG::SEfromHessian(FIM)
  }else{
    se = HelpersMG::SEfromHessian(FIM[-4,-4])
  }
  se_theta = se[1]
  se_r = se[2]
  lobs =  sum.finite(log((1-p)*L0_i(b_exp,b_out,se_exp,se_out,theta) +
                           p*L1_i(b_exp,b_out,se_exp,se_out,theta,r,c)))
  pval = 2*stats::pnorm(abs(theta/se_theta),lower.tail=FALSE)
  pval_r = 2*stats::pnorm(abs(r/se_r),lower.tail=FALSE)
  BIC =  -2*lobs + log(n)*(2+sum(c(r!=0,c>1)))

  return(list(theta=theta,se=se_theta,pval=pval,
              r=r,se_r=se_r,pval_r=pval_r,
              c=c,p=p,
              BIC=BIC,
              niter=niter,
              tau_1_vec=tau_1_vec))
}

#' CEM algorithm
#'
#' Run CEM algorithm with one starting value
#'
#' @import HelpersMG stats
#' @keywords internal
#'

CEM <-function(b_exp,b_out,se_exp,se_out,
               initial_theta,initial_r,initial_c,initial_p,
               EM_start,EM_maxit,
               maxit,n){
  m = length(b_exp)
  theta = initial_theta
  theta_old = theta -1
  r = initial_r
  c = initial_c
  p = initial_p
  niter = 0
  ## EM-start
  if(EM_start==T){
    EM_first = EM(b_exp = b_exp,b_out=b_out,se_exp = se_exp,se_out = se_out,
                  initial_theta = initial_theta,initial_p=initial_p,
                  initial_r=initial_r,initial_c=initial_c,
                  maxit = EM_maxit,n=n)
    theta = EM_first$theta
    theta_old = theta -1
    r = EM_first$r
    c = EM_first$c
    p = EM_first$p
    if(is.na(theta)){
      warning('please choose another starting value or increase maxit for EM start!')
      out=EM_first
      out$BIC = NaN
      return(out)
    }
  }

  while(abs(theta-theta_old)>=1e-6 & niter<maxit){
    niter = niter + 1
    theta_old = theta
    ## E-step
    tau_vec = tau_i(b_exp,b_out,se_exp,se_out,theta,r,c,p)
    tau_0_vec = tau_vec$tau_i0
    tau_1_vec = tau_vec$tau_i1

    ## C-step
    invalid_ind = which(tau_1_vec>0.5)
    valid_ind = setdiff(1:m,invalid_ind)
    ## M-step
    if(length(invalid_ind)>0){
      p = length(invalid_ind)/m
      r_vec = rep(0,m)
      r = sum((b_out[invalid_ind]-theta*b_exp[invalid_ind])/se_out[invalid_ind]^2) /
        sum(1/se_out[invalid_ind]^2)
      r_vec[invalid_ind] = r
      c_vec = rep(1,m)
      c = max(1,sum((b_out[invalid_ind]-theta*b_exp[invalid_ind]-r_vec[invalid_ind])^2 /
                      se_out[invalid_ind]^2)/(length(invalid_ind)))
      c_vec[invalid_ind] = c
    } else{
      invalid_ind = integer(0)
      p = 0
      r = 0
      r_vec = rep(0,m)
      c = 1
      c_vec = rep(1,m)
    }

    theta =
      sum((b_out - r_vec)*b_exp / (c_vec*se_out^2)) /
      sum(b_exp^2 / (c_vec*se_out^2))
    if(is.na(theta)) break
  }
  if(is.na(theta) || any(is.na(tau_1_vec))){
    warning('please choose another starting value or increase maxit!')
    return(list(theta=NaN,se=NaN,pval=NaN,
                r=NaN,se_r=NaN,pval_r=NaN,
                c=NaN,p=NaN,
                BIC=NaN,
                niter=niter,
                tau_1_vec=tau_1_vec))
  }

  tau_1_full = rep(0,m)
  tau_1_full[which(tau_1_vec>0.5)] = 1
  tau_0_full = 1 - tau_1_full
  B_complete = B_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,p,tau_0_vec,tau_1_vec)
  B_missing = S_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,p,tau_0_vec,tau_1_vec)
  FIM = B_complete - B_missing
  if(is.finite(FIM[4,4])){
    se = suppressWarnings(sqrt(diag(solve(FIM))))
    if(any(is.na(se))){
      B_complete = B_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,mean(tau_1_vec),tau_0_vec,tau_1_vec)
      B_missing = S_Complete(b_exp,b_out,se_exp,se_out,theta,r,c,mean(tau_1_vec),tau_0_vec,tau_1_vec)
      se = HelpersMG::SEfromHessian(B_complete-B_missing)}
  }else{
    se = HelpersMG::SEfromHessian(FIM[-4,-4])
  }

  se_theta = se[1]
  se_r = se[2]
  lcom = sum(log(tau_0_full*L0_i(b_exp,b_out,se_exp,se_out,theta) +
                   tau_1_full*L1_i(b_exp,b_out,se_exp,se_out,theta,r,c)))
  pval = 2*stats::pnorm(abs(theta/se_theta),lower.tail=FALSE)
  pval_r = 2*stats::pnorm(abs(r/se_r),lower.tail=FALSE)
  BIC = -2*lcom + log(n)*(2+sum(c(r!=0,c>1)))
  return(list(theta=theta,se=se_theta,pval=pval,
              r=r,se_r=se_r,pval_r=pval_r,
              c=c,p=p,
              BIC=BIC,
              niter=niter,
              tau_1_vec=tau_1_vec))
}

