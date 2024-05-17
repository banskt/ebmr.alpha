#' @title Update by fitting GRR model
#'
#' @description Updates parameters of an EBMR fit by fitting the full GRR model,
#' including estimating hyperparameters (residual variance and prior scaling factor).
#' The method uses an SVD on Xtilde = XW^0.5, after which
#' other computations (EM algorithm to maximize hyperparameters, and posterior computations)
#' are cheap.
#'
#' @details The code is based on the ideas at https://stephens999.github.io/misc/ridge_em_svd.html
#' but extended to deal with case where n>p as well as n<=p.
#'
#' @param fit a previous ebmr fit
#'
#' @param tol a tolerance for the EM algorithm
#'
#' @param maxiter maximum number of iterations for EM algorithm
#'
#' @param compute_Sigma_diag boolean flag; set to TRUE to compute the diagonal of the posterior variance, otherwise this variance is set to 0.
#'
#' @param compute_Sigma_full boolean flag; set to TRUE to compute full Sigma matrix for testing purposes
#'
#' @param update_residual_variance boolean flag; set to FALSE to keep residual variance fixed
#'
#' @return an updated ebmr fit
ebmr.update.grr.svd = function(fit, tol = 1e-8, maxiter = 1000, compute_Sigma_diag = TRUE, compute_Sigma_full=FALSE, update_residual_variance=TRUE){

  # svd computation
  Xtilde = t(fit$wbar^0.5 * t(fit$X))
  Xtilde.svd = svd(Xtilde)

  # maximum likelihood estimation
  ytilde = drop(t(Xtilde.svd$u) %*% fit$y)
  df = length(fit$y) - length(ytilde)
  ss = sum(fit$y^2) - sum(ytilde^2)

  yt.em3 = ridge_indep_em3(ytilde, Xtilde.svd$d^2, ss, df, tol, maxiter, s2.init = fit$residual_variance, sb2.init = fit$sb2, update_s2 = update_residual_variance)

  fit$residual_variance = yt.em3$s2
  fit$sb2 = yt.em3$sb2

  if(compute_Sigma_diag){
    fit$Sigma_diag = Sigma1_diag_woodbury_svd(fit$wbar,Xtilde.svd,fit$sb2)
  } else {
    fit$Sigma_diag = rep(0, fit$p)
  }

  if(compute_Sigma_full){
    fit$Sigma_full = Sigma1_woodbury_svd(fit$wbar,Xtilde.svd,fit$sb2)
  }

  fit$mu = mu1_woodbury_svd(fit$y, fit$wbar, Xtilde.svd, fit$sb2)

  fit$h2_term = -0.5*fit$p + 0.5*sum(log(fit$sb2*fit$wbar)) - 0.5*sum(log(1+fit$sb2*Xtilde.svd$d^2))

  #fit = ebmr.update.ebnv(fit, ebnv.pm) # done to add the correct logdet_KL_term

  return(fit)
}



#' @title Function for estimating hyper-parameters for "independent" ridge regressions (ie the orthogonal case)
#'
#' @description
#' Obtains maximum likelihood estimates for (s2,sb2)
#' in the model \deqn{y_i | theta_i \sim N(theta_i, s2), theta_i \sim N(0,sb2 s2 d2_i),} 
#' where y,d2 are given and (s2,sb2) are to be estimated. However it uses a different parameterization
#' to improve convergence (see details). It also incorporates additional sum of squares (ss) and degrees of freedom (df)
#' which effectively specify a prior on s2 (necessary to deal with fitting ridge regression when n>p above).
#'
#' @details This code is based on https://stephens999.github.io/misc/ridge_em_svd.html
#' extended to take account of additional residual sum of squares (ss)
#' and degrees of freedom parameter (df), to deal with the case n>p.
#' To improve convergence (see URL above) it fits the model by using the over-parameterized model:
#' \deqn{y_i | theta_i \sim N(sb theta_i,s^2), theta_i \sim N(0,l2 d2_i).}
#' So it returns s2hat = s2 and sb2hat = sb2 *l2 / s2
#'
#' @param y numeric vector
#'
#' @param d2 numeric vector
#'
#' @param ss sum of squares for prior on s2
#'
#' @param df degrees of freedom for ss
#'
#' @param tol real scalar, specifying convergence tolerance for log-likelihood
#'
#' @param maxiter integer, maximum number of iterations
#'
#' @param s2.init initial value for s2
#'
#' @param sb2.init initial value for sb2
#'
#' @param l2.init initial value for l2
#'
#' @param update_s2 boolean; if false then s2 is fixed to s2.init
#'
ridge_indep_em3 = function(y, d2, ss = 0, df = 0, tol=1e-3, maxiter=1000, s2.init = 1, sb2.init=1, l2.init=1, update_s2 = TRUE){
  k = length(y)
  s2 = s2.init
  sb2 = sb2.init
  l2 = l2.init

  ll = sum(dnorm(y,mean=0,sd = sqrt(sb2*l2*d2 + s2),log=TRUE))
  # plus the part due to additional sum of squares
  ll = ll - (df/2) * log(2*pi*s2) - ss/(2*s2)

  loglik = ll

  for(i in 1:maxiter){

    prior_var = d2*l2 # prior variance for theta
    data_var = s2/sb2 # variance of y/sb, which has mean theta
    post_var = 1/((1/prior_var) + (1/data_var)) #posterior variance of theta
    post_mean =  post_var * (1/data_var) * (y/sqrt(sb2)) # posterior mean of theta

    sb2 = (sum(y*post_mean)/sum(post_mean^2 + post_var))^2
    l2 = mean((post_mean^2 + post_var)/d2)

    if(update_s2){
      r = y - sqrt(sb2) * post_mean # residuals
      s2 = (sum(r^2 + sb2 * post_var ) + ss) / (length(y) + df) # adjusted sum of squares and df
    }

    ll = sum(dnorm(y,mean=0,sd = sqrt(sb2*l2*d2 + s2),log=TRUE))
    ll = ll - (df/2)*log(2*pi*s2) - ss/(2*s2)
    loglik= c(loglik,ll)

    if(is.nan(ll) || (loglik[i+1]-loglik[i])<tol) break
  }
  return(list(s2=s2,sb2=sb2*l2/s2,loglik=loglik))
}

# original code from URL above, used for testing
ridge_em3 = function(y,X, s2, sb2, l2, niter=10){
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  yty = t(y) %*% y
  n = length(y)
  p = ncol(X)
  loglik = rep(0,niter)
  for(i in 1:niter){
    V = chol2inv(chol(XtX+ diag(s2/(sb2*l2),p)))

    SigmaY = l2*sb2 *(X %*% t(X)) + diag(s2,n)
    loglik[i] = mvtnorm::dmvnorm(as.vector(y),sigma = SigmaY,log=TRUE)

    Sigma1 = (s2/sb2)*V  # posterior variance of b
    mu1 = (1/sqrt(sb2))*as.vector(V %*% Xty) # posterior mean of b


    sb2 = (sum(mu1*Xty)/sum(diag(XtX %*% (mu1 %*% t(mu1) + Sigma1))))^2
    s2 = as.vector((yty + sb2*sum(diag(XtX %*% (mu1 %*% t(mu1) + Sigma1)))- 2*sqrt(sb2)*sum(Xty*mu1))/n)

    l2 = mean(mu1^2+diag(Sigma1))

  }
  return(list(s2=s2,sb2=sb2,l2=l2,loglik=loglik,postmean=mu1*sqrt(sb2)))
}


