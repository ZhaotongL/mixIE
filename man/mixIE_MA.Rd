% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixIE_MA.R
\name{mixIE_MA}
\alias{mixIE_MA}
\title{mixIE with model averaging}
\usage{
mixIE_MA(b_exp, b_out, se_exp, se_out, n, n_model = 5, ...)
}
\arguments{
\item{b_exp}{A vector of SNP effects on the exposure variable, usually obtained from a GWAS.}

\item{b_out}{A vector of SNP effects on the outcome variable, usually obtained from a GWAS.}

\item{se_exp}{A vector of standard errors of \code{b_exp}.}

\item{se_out}{A vector of standard errors of \code{b_out}.}

\item{n}{Sample size of either one of the GWAS dataset.}

\item{n_model}{Number of top models in the final candidate list, default is 5.}

\item{...}{Arguments to be passed to \code{mixIE_multiple_start}}
}
\value{
A list
\describe{
\item{theta_BIC_MA}{Estimated causal effect from mixIE-MA}
\item{se_BIC_MA}{Estimate standard error for \code{theta_BIC_MA}}
\item{pval_BIC_MA}{Two-sided p-value of \code{theta_BIC_MA}}
\item{c_BIC_MA}{Estimated overdispersion parameter from mixIE-MA}
\item{r_BIC_MA}{Estimated pleiotropic effect from mixIE-MA}
\item{p_BIC_MA}{Estimated proportion of invalid IVs from mixIE-MA}
\item{tau_BIC_MA}{A vector of (averaged) posterior probabilities of each IVs being invalid from mixIE-MA}
\item{theta_vec}{A vector of estimated causal effects from \code{n_model} models}
}
}
\description{
This is the main function of mixIE-MA
}
