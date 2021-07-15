#' mixIE-MA with data perturbation (mixIE-MA-DP)
#'
#' This is the main function of mixIE-MA method with data perturbation.
#'
#' @param b_exp A vector of SNP effects on the exposure variable, usually obtained from a GWAS.
#' @param b_out A vector of SNP effects on the outcome variable, usually obtained from a GWAS.
#' @param se_exp A vector of standard errors of \code{b_exp}.
#' @param se_out A vector of standard errors of \code{b_out}.
#' @param n Sample size of either one of the GWAS dataset.
#' @param flip Whether to reorient the SNPs like Egger regression?
#' @param B Number of data perturbation, default is 200.
#' @param thres_e Could ignore this.
#' @param diagnostic_plot Should the function returns diagnostic plots? Default is FALSE
#' @param ... Arguments to be passed to \code{mixIE_MA}
#'
#' @import ggplot2 stats
#' @return A list
#' \describe{
#' \item{mixIE_MA_theta}{Estimated causal effect from mixIE-MA}
#' \item{mixIE_MA_se}{Estimate standard error for \code{mixIE_MA_theta}}
#' \item{mixIE_MA_pval}{Two-sided p-value of \code{mixIE_MA_theta}}
#' \item{mixIE_MA_pi}{Estimated proportion of invalid IVs from mixIE-MA}
#' \item{mixIE_MA_tau}{A vector of (averaged) posterior probabilities of each IVs being invalid from mixIE-MA}
#' \item{mixIE_MA_DP_theta}{Estimated causal effect from mixIE-MA-DP}
#' \item{mixIE_MA_DP_se}{Estimate standard error for \code{mixIE_MA_DP_theta}}
#' \item{mixIE_MA_DP_pval}{Two-sided p-value of \code{mixIE_MA_DP_theta}}
#' \item{mixIE_MA_DP_pi}{Estimated proportion of invalid IVs from mixIE-MA-DP}
#' \item{mixIE_MA_DP_tau}{A vector of posterior probabilities of each IVs being invalid from mixIE-MA-DP}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' result = mixIE_MA_DP(b_exp=FG_T2D$beta_hat_1,
#'                      b_out=FG_T2D$beta_hat_2,
#'                      se_exp=FG_T2D$seb1,
#'                      se_out = FG_T2D$seb2,
#'                      flip=1,
#'                      n=69033)
#' result
#'
#' ## Diagnostic plot
#' set.seed(1)
#' result = mixIE_MA_DP(b_exp=FG_T2D$beta_hat_1,
#'                      b_out=FG_T2D$beta_hat_2,
#'                      se_exp=FG_T2D$seb1,
#'                      se_out = FG_T2D$seb2,
#'                      flip=1,
#'                      n=69033,diagnostic_plot=TRUE)
#'
#' ## Compare IVs' posterior probability from mixIE-MA and mixIE-MA-DP
#' result$iv_barplot
#'
#' ## Histogram of estimates from data perturbation
#' result$est_hist
#'
#' ## GWAS scatter plot (red points are valid IVs, blue points are invalid IVs)
#' result$scatter_og.plot # mixIE-MA
#' result$scatter_dp.plot # mixIE-MA-DP
#'
mixIE_MA_DP<- function(b_exp,b_out,se_exp,se_out,n,flip=1,
                       B=200,thres_e=Inf,diagnostic_plot=FALSE,...){
  if(flip==1){
    sign0 <- function(x) {
      x[x == 0] <- 1
      return(sign(x))
    }
    to_flip <- sign0(b_exp) == -1
    b_out = b_out * sign0(b_exp)
    b_exp = abs(b_exp)
  }
  if(length(b_exp)==2){
    ivw_res = mr_ivw_fe(b_exp,b_out,se_exp,se_out)
    out = list(mixIE_MA_theta=ivw_res$b, mixIE_MA_se=ivw_res$se, mixIE_MA_pval=ivw_res$pval, mixIE_MA_pi=0, mixIE_MA_tau=rep(0,2),
               mixIE_MA_DP_theta=ivw_res$b, mixIE_MA_DP_se=ivw_res$se, mixIE_MA_DP_pval=ivw_res$pval, mixIE_MA_DP_pi=0, mixIE_MA_DP_tau=rep(0,2))
  } else{
    theta_f_B  = se_f_B = rep(NaN,B)
    m = length(b_exp)
    f_og_result = mixIE_MA(b_exp=b_exp,b_out=b_out,se_exp=se_exp,se_out=se_out,n=n,flip=flip,...)
    invalid_count_B = rep(0,m)
    invalid_p_B = rep(0,m)
    for(i in 1:B){
      e_out_dp = unlist(lapply(1:m,function(x){rnorm(1,0,se_out[x])}))
      b_out_dp = b_out + e_out_dp
      f_result = mixIE_MA(b_exp=b_exp,b_out=b_out_dp,se_exp=se_exp,se_out=sqrt(2)*se_out,n=n,flip=flip,...)
      invalid_p_B = invalid_p_B + f_result$tau_BIC_MA
      invalid_ind=which(f_result$tau_BIC_MA>0.5)
      valid_ind = setdiff(1:m,invalid_ind)
      invalid_count_B[invalid_ind] = invalid_count_B[invalid_ind]+1
      if(length(valid_ind)==0){
        egger_b = mr_egger_ll(b_exp[invalid_ind],b_out_dp[invalid_ind],se_exp[invalid_ind],se_out[invalid_ind])
        b_out_dp_star = b_out_dp[invalid_ind] + egger_b$sig*e_out_dp[invalid_ind]
        egger_b_star = mr_egger_ll(b_exp[invalid_ind],b_out_dp_star,se_exp[invalid_ind],se_out[invalid_ind])
        theta_b_star = egger_b_star$b
      } else if(length(valid_ind)==1){
        ivw_theta_b = b_out_dp[valid_ind]/b_exp[valid_ind]
        ivw_se_b = sqrt(2)*se_out[valid_ind]/b_exp[valid_ind]
        egger_b = mr_egger_ll(b_exp[invalid_ind],b_out_dp[invalid_ind],se_exp[invalid_ind],sqrt(2)*se_out[invalid_ind])
        b_out_dp_star = b_out_dp[invalid_ind] + egger_b$sig*e_out_dp[invalid_ind]
        egger_b_star = mr_egger_ll(b_exp[invalid_ind],b_out_dp_star,se_exp[invalid_ind],se_out[invalid_ind])
        theta_b_star = sum(ivw_theta_b/ivw_se_b^2+egger_b_star$b/egger_b$se^2)/
          sum(1/ivw_se_b^2+1/egger_b$se^2)
      } else if(length(invalid_ind)<3){
        ivw_b = mr_ivw_fe(b_exp[valid_ind],b_out_dp[valid_ind],se_exp[valid_ind],sqrt(2)*se_out[valid_ind])
        theta_b_star = f_result$theta_BIC_MA
      } else{
        ivw_b = mr_ivw_fe(b_exp[valid_ind],b_out_dp[valid_ind],se_exp[valid_ind],sqrt(2)*se_out[valid_ind])
        egger_b = mr_egger_ll(b_exp[invalid_ind],b_out_dp[invalid_ind],se_exp[invalid_ind],sqrt(2)*se_out[invalid_ind])
        b_out_dp_star = b_out_dp[invalid_ind] + egger_b$sig*e_out_dp[invalid_ind]
        egger_b_star = mr_egger_ll(b_exp[invalid_ind],b_out_dp_star,se_exp[invalid_ind],se_out[invalid_ind])
        theta_b_star = sum(ivw_b$b/ivw_b$se^2+egger_b_star$b/egger_b$se^2)/
          sum(1/ivw_b$se^2+1/egger_b$se^2)
      }
      theta_f_B[i] = theta_b_star
    }

    fdp_theta = mean(theta_f_B[which(theta_f_B<thres_e)],na.rm=T)
    fdp_se = sd(theta_f_B[which(theta_f_B<thres_e)],na.rm=T)
    fdp_pval = 2*pnorm(abs(fdp_theta/fdp_se),lower.tail=FALSE)
    invalid_count_B = invalid_count_B/B

    out.summary = list(mixIE_MA_theta=f_og_result$theta_BIC_MA,
                       mixIE_MA_se=f_og_result$se_BIC_MA,
                       mixIE_MA_pval=f_og_result$pval_BIC_MA,
                       mixIE_MA_pi=mean(f_og_result$tau_BIC_MA>=0.5),
                       mixIE_MA_tau=f_og_result$tau_BIC_MA,
                       mixIE_MA_DP_theta=fdp_theta,
                       mixIE_MA_DP_se=fdp_se,
                       mixIE_MA_DP_pval=fdp_pval,
                       mixIE_MA_DP_pi=mean(invalid_count_B>=0.5),
                       mixIE_MA_DP_tau=invalid_count_B)
    out = out.summary
    if(diagnostic_plot==T){
      ##JUST FOR R CHECK
      theta_b <- IV <- Prob <- method <- invalid <- NULL

      plot_hist = data.frame(theta_b=theta_f_B)
      histest= ggplot(plot_hist, aes(x=theta_b)) +
        geom_histogram(bins = 100)

      iv_hist = data.frame(IV=rep(factor(1:m),2),
                           Prob=c(invalid_count_B,f_og_result$tau_BIC_MA),
                           method=rep(c('mixIE-MA-DP','mixIE-MA'),each=m))
      bariv = ggplot(data=iv_hist,aes(x=IV,y=Prob,fill=method,group=method)) +
        geom_bar(stat="identity",position = 'dodge',alpha=0.8)+
        scale_fill_brewer(palette = "Paired")+
        geom_hline(yintercept=0.5,linetype='dashed') +
        theme(axis.text=element_text(size=7))+
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE))

      plot_og.df = data.frame(b_exp=b_exp,
                              b_out=b_out,
                              invalid=factor(f_og_result$tau_BIC_MA>0.5,levels=c(T,F)))

      scatter_og.plot = ggplot(data=plot_og.df, aes(x=b_exp, y=b_out,color=invalid)) +
        geom_point(size=1) +
        theme_minimal() +
        scale_colour_manual(name="invalid",values=c( "#377EB8","#E41A1C"))+
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        geom_abline(slope=f_og_result$theta_BIC_MA,intercept = 0,color="#A6CEE3") +
        scale_size_manual(values=2)


      plot_dp.df = data.frame(b_exp=b_exp,b_out=b_out,invalid=factor(invalid_count_B>0.5,levels=c(T,F)))

      scatter_dp.plot = ggplot(plot_dp.df, aes(x=b_exp, y=b_out,color=invalid)) +
        geom_point(size=1) +
        theme_minimal() +
        scale_colour_manual(name="invalid",values=c( "#377EB8","#E41A1C"))+
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        geom_abline(slope=fdp_theta,intercept=0,color="#1F78B4")

      out = append(out.summary,list(est_hist=histest,iv_barplot=bariv,scatter_og.plot=scatter_og.plot,scatter_dp.plot=scatter_dp.plot))
    }
  }
  return(out)
}
