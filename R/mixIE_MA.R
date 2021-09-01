#' mixIE_multiple_start
#'
#' Perform CEM algorithm with multiple starting values
#'
#' @param b_exp A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param b_out A vector of SNP effects on the outcome variable, usually obtained from a GWAS.
#' @param se_exp A vector of standard errors of \code{b_exp}.
#' @param se_out A vector of standard errors of \code{b_out}.
#' @param n Sample size of either one of the GWAS dataset.
#' @param initial_theta_n Number of different thetas generated, default is 50.
#' @param initial_r Initial value of r, the averaged pleiotropic effect, default is 0.
#' @param initial_c Initial value of c, the overdispersion parameter for invalid IVs, default is 1.
#' @param initial_p Initial value of the proportion of invalid IVs, default is 0.2.
#' @param flip Whether to reorient the SNPs like Egger regression?
#' @param EM_start Whether to use EM algorithm to start the CEM? Default is TRUE.
#' @param EM_maxit Number of iterations used in EM algorithm to start the CEM, default is 2.
#' @param maxit Maximum number of iterations for each optimization, default is 200.
#' @param ivw Whether to add the fixed effect IVW in the model candidates explictly? Default is TRUE.
#' @param egger Whether to add the Egger regression in the model candidates explictly? Default is TRUE.
#' @param lb.theta Lower bound of theta from which starting values are generated.
#' @param ub.theta Upper bound of theta from which starting values are generated.
#'
#' @import tidyr
#'
#' @return A list
#' \describe{
#' \item{theta}{A vector of estimated causal effect for different starting values ordered by BIC}
#' \item{se}{A vector of corresponding standard error of \code{theta}}
#' \item{pval}{A vector of corresponding two-sided p-value of \code{theta}}
#' \item{c}{A vector of estimated overdispersion parameter for different starting values ordered by BIC}
#' \item{r}{A vector of estimated pleiotropic effect for different starting values ordered by BIC}
#' \item{ser}{A vector of corresponding standard error of \code{r}}
#' \item{p}{A vector of estimated proportion of invalid IVs for different starting values ordered by BIC}
#' \item{BIC}{A vector of BIC for different starting values sorting increasingly}
#' \item{niter}{A vector of number of iterations for different starting values ordered by BIC}
#' \item{tau_1_mat}{A matrix of posterior probabilities of each IVs being invalid for different starting values ordered by BIC}
#' }
#'
#' @export
#'
mixIE_multiple_start <- function(b_exp,b_out,se_exp,se_out,n,
                                 initial_theta_n=50,initial_r=0,initial_c=1,initial_p=0.2,
                                 flip=1, EM_start=T,EM_maxit=2,
                                 maxit=200,ivw=T,egger=T,lb.theta=NULL,ub.theta=NULL){
  if(flip==1){
    sign0 <- function(x) {
      x[x == 0] <- 1
      return(sign(x))
    }
    to_flip <- sign0(b_exp) == -1
    b_out = b_out * sign0(b_exp)
    b_exp = abs(b_exp)
  }

  m = length(b_exp)
  bound.theta <- max(abs(b_out/b_exp))
  mrivw = mr_ivw_fe(b_exp,b_out,se_exp,se_out)
  mregger = mr_egger_ll(b_exp,b_out,se_exp,se_out)
  if(is.null(lb.theta)) {lb.theta = -bound.theta}
  if(is.null(ub.theta)) {ub.theta = bound.theta}
  theta_seq = runif(initial_theta_n,min=lb.theta,max=ub.theta)
  theta_seq = c(0,theta_seq,min(b_out/b_exp),max(b_out/b_exp),mrivw$b,mregger$b)
  starting_values = crossing(theta_seq,initial_r)
  r_seq = starting_values$initial_r
  theta_seq = starting_values$theta_seq


  theta_vec = se_vec =ser_vec=pval_vec =  p_vec = c_vec = r_vec = BIC_vec = rep(NaN,length(theta_seq)+2)
  tau_1_mat = matrix(NaN,nrow=length(theta_seq)+2,ncol=m)
  niter_vec = rep(maxit,length(theta_seq)+2)
  for(i in 1:length(theta_seq)){
    CEM_i = CEM(b_exp,b_out,se_exp,se_out,
                initial_theta = theta_seq[i],
                initial_r = r_seq[i], initial_c=initial_c, initial_p=initial_p,
                EM_start = EM_start,EM_maxit = EM_maxit,
                maxit = maxit,n=n)
    if(is.na(CEM_i$theta)) next
    if(CEM_i$niter==maxit) next
    theta_vec[i]=CEM_i$theta
    se_vec[i] = CEM_i$se
    ser_vec[i] = CEM_i$se_r
    pval_vec[i] = CEM_i$pval
    BIC_vec[i] = CEM_i$BIC
    p_vec[i] = CEM_i$p
    r_vec[i] = CEM_i$r
    c_vec[i] = CEM_i$c
    niter_vec[i] = CEM_i$niter
    tau_1_mat[i,] = CEM_i$tau_1_vec
  }
  if(ivw==T){
    mivw = mr_ivw_fe(b_exp,b_out,se_exp,se_out)
    theta_vec[i+1]=mivw$b
    se_vec[i+1]=mivw$se
    pval_vec[i+1]=mivw$pval
    lcom_ivw=sum(log(L0_i(b_exp,b_out,se_exp,se_out,mivw$b)))
    BIC_vec[i+1]=-2*lcom_ivw+log(n)*(1)
    p_vec[i+1]=0
    r_vec[i+1]=0
    c_vec[i+1]=1
    tau_1_mat[i+1,]=rep(0,m)
    niter_vec[i+1] = 1

  }
  if(egger==T){
    megger = mr_egger_ll(b_exp,b_out,se_exp,se_out,flip=flip)
    theta_vec[i+2]=megger$b
    se_vec[i+2]=megger$se
    ser_vec[i+2]=megger$se_i
    pval_vec[i+2]=megger$pval
    lcom_egger=sum(log(L1_i(b_exp,b_out,se_exp,se_out,megger$b,megger$b_i,(sqrt((m-2)/m)*megger$sig)^2)))
    BIC_vec[i+2]=-2*lcom_egger+log(n)*3
    p_vec[i+2]=1
    r_vec[i+2]=megger$b_i
    c_vec[i+2]=(sqrt((m-2)/m)*megger$sig)^2
    tau_1_mat[i+2,]=rep(1,m)
    niter_vec[i+2] = 1

  }

  BIC_order = order(BIC_vec,decreasing = F)
  theta_BIC = theta_vec[BIC_order]
  r_BIC = r_vec[BIC_order]
  c_BIC = c_vec[BIC_order]
  p_BIC = p_vec[BIC_order]
  se_BIC = se_vec[BIC_order]
  ser_BIC = ser_vec[BIC_order]
  pval_BIC = pval_vec[BIC_order]
  BIC_BIC = BIC_vec[BIC_order]
  theta_seq_BIC = theta_seq[BIC_order]
  niter_BIC = niter_vec[BIC_order]
  tau_1_BIC = tau_1_mat[BIC_order,]

  if(is.na(theta_BIC[1])){
    warning('please use another set of starting values or increase maxit!')
    out = CEM(b_exp,b_out,se_exp,se_out,
              initial_theta = 0,
              initial_r = 0, initial_c=initial_c, initial_p=initial_p,
              EM_start = EM_start,EM_maxit = EM_maxit,
              maxit = maxit,n=n)
    return(list(BIC_model=out))
  }

  out = list(theta=theta_BIC,se=se_BIC,pval=pval_BIC,r=r_BIC,ser=ser_BIC,p=p_BIC,c=c_BIC,niter=niter_BIC,
             BIC=BIC_BIC,
             tau_1_mat=tau_1_BIC)

  return(out)

}

#' mixIE with model averaging
#'
#' This is the main function of mixIE-MA
#'
#' @param b_exp A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param b_out A vector of SNP effects on the outcome variable, usually obtained from a GWAS.
#' @param se_exp A vector of standard errors of \code{b_exp}.
#' @param se_out A vector of standard errors of \code{b_out}.
#' @param n Sample size of either one of the GWAS dataset.
#' @param n_model Number of top models in the final candidate list, default is 5.
#' @param ... Arguments to be passed to \code{mixIE_multiple_start}
#'
#' @return A list
#' \describe{
#' \item{theta_BIC_MA}{Estimated causal effect from mixIE-MA}
#' \item{se_BIC_MA}{Estimate standard error for \code{theta_BIC_MA}}
#' \item{pval_BIC_MA}{Two-sided p-value of \code{theta_BIC_MA}}
#' \item{c_BIC_MA}{Estimated overdispersion parameter from mixIE-MA}
#' \item{r_BIC_MA}{Estimated pleiotropic effect from mixIE-MA}
#' \item{p_BIC_MA}{Estimated proportion of invalid IVs from mixIE-MA}
#' \item{tau_BIC_MA}{A vector of (averaged) posterior probabilities of each IVs being invalid from mixIE-MA}
#' \item{theta_vec}{A vector of estimated causal effects from \code{n_model} models}
#' }
#'
#' @export
#'
mixIE_MA <- function(b_exp,b_out,se_exp,se_out,n,
                     n_model=5,...){
  if(length(b_exp)==2){
    ivw_res = mr_ivw_fe(b_exp,b_out,se_exp,se_out)
    theta_BIC_MA=ivw_res$b
    se_BIC_MA=ivw_res$se
    pval_BIC_MA=ivw_res$pval
    r_BIC_MA=0
    p_BIC_MA=0
    c_BIC_MA=1
    tau_BIC_MA=rep(0,2)
  } else{
    mixIEms_BIC = mixIE_multiple_start(b_exp,b_out,se_exp,se_out,n=n,...)
    if(TRUE){
      keep_ind = is.finite(mixIEms_BIC$theta) & is.finite(mixIEms_BIC$BIC)
      theta_vec = mixIEms_BIC$theta[keep_ind]
      p_vec = mixIEms_BIC$p[keep_ind]
      r_vec = mixIEms_BIC$r[keep_ind]
      c_vec = mixIEms_BIC$c[keep_ind]
      se_vec = mixIEms_BIC$se[keep_ind]
      ser_vec = mixIEms_BIC$ser[keep_ind]
      BIC_vec = mixIEms_BIC$BIC[keep_ind]
      BIC_vec = BIC_vec - min(BIC_vec)
      tau_1_mat = mixIEms_BIC$tau_1_mat[keep_ind,]

      BICmodel_ind =  c()
      i = j = 1
      while(i<length(theta_vec)){
        BICmodel_ind = c(BICmodel_ind,i)
        for(j in (i+1):length(theta_vec)){
          if(abs(theta_vec[j]-theta_vec[i])<1e-04) next
          if(abs(theta_vec[j]-theta_vec[i])>1e-04) break
        }
        i = j
      }
      if(!is.null(BICmodel_ind)){
        if(abs(theta_vec[i]-theta_vec[BICmodel_ind[length(BICmodel_ind)]]) > 1e-04) BICmodel_ind=c(BICmodel_ind,i)}
      BICmodel_ind = BICmodel_ind[1:min(length(BICmodel_ind),n_model)]
      if(is.null(BICmodel_ind)){BICmodel_ind=1}
      weight_vec = exp(-1/2*BIC_vec[BICmodel_ind])
      weight_vec = weight_vec/sum(weight_vec)
      theta_BIC_MA = sum(weight_vec*theta_vec[BICmodel_ind])
      r_BIC_MA = sum(weight_vec*r_vec[BICmodel_ind])
      c_BIC_MA = sum(weight_vec*c_vec[BICmodel_ind])
      if(length(BICmodel_ind)>1){
        tau_BIC_MA = colSums(weight_vec*tau_1_mat[BICmodel_ind,])
      } else if(is.vector(tau_1_mat)){
        tau_BIC_MA = tau_1_mat
      } else{
        tau_BIC_MA = tau_1_mat[BICmodel_ind,]
      }
      p_BIC_MA = mean(tau_BIC_MA>0.5)
      se_BIC_MA = sum(weight_vec*sqrt(se_vec[BICmodel_ind]^2 +(theta_vec[BICmodel_ind] -theta_BIC_MA)^2))
      ser_BIC_MA = sum(weight_vec*sqrt(ser_vec[BICmodel_ind]^2 +(r_vec[BICmodel_ind] -r_BIC_MA)^2))
      pval_BIC_MA = 2*pnorm(abs(theta_BIC_MA/se_BIC_MA),lower.tail=FALSE)
      pvalr_BIC_MA = 2*pnorm(abs(r_BIC_MA/ser_BIC_MA),lower.tail=FALSE)
    }

}
  return(list(theta_BIC_MA=theta_BIC_MA,se_BIC_MA=se_BIC_MA,pval_BIC_MA=pval_BIC_MA,
              r_BIC_MA=r_BIC_MA,p_BIC_MA=p_BIC_MA,c_BIC_MA=c_BIC_MA,
              tau_BIC_MA=tau_BIC_MA,
              theta_vec=theta_vec[BICmodel_ind]))
}



