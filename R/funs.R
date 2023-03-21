

############################################################

# takes as input a 1d vector x and outputs a matrix with nrow rows, with each row containing x
.row_rep = function(x,nrow){
  rbind(x)[rep(1,nrow),]
}


############################################################

# cosine basis functions
# note: nn=1 returns a constant function
.my_cosine = function(xx,nn){
  sqrt(2)*cos(xx*(nn-1)*pi)
}

############################################################

# simulate from a raised cosine type distribution with density (up to vertical shift) equal to functions in the
# cosine basis.
# The cosine function is raised by 1 (so that it takes minimum value zero)

.dcos = function(x,nn){(cos(nn*pi*x)+1)}
.pcos = function(x,nn){(pi*nn*x + sin(pi*nn*x))/(pi*nn)}
.qcos = function(q,nn){sapply(q,function(curr_q){uniroot(function(x){.pcos(x,nn)-curr_q},c(0,1),extendInt="upX")$root})}
.rcos = function(n,nn){.qcos(runif(n),nn)}

# .pcos(quantile(.qcos(runif(1e4),10)),10)

############################################################

# sinc function
.sinc = function(x){
  (sin(x)+(x==0))/(x + (x==0)) # by adding x==0, get the right behavior at the singularity at 0 --- .sinc(0)=1
}

# define a distribution with density proportional to f(x)=.sinc(x)^4
.dsinc4 = function(x){.sinc(x)^4/(2*pi/3)}
.psinc4 = function(p){require(pracma);1/2 + (-4*pracma::Si(2*p) + 8*pracma::Si(4*p) - ((sin(p)^2) * (4*p^2 + (8*p^2-1)*cos(2*p) + 2*p*sin(2*p)+1))/(p^3 + (p==0)))/(4*pi)}
.qsinc4 = function(q){sapply(q,function(curr_q){uniroot(function(x){.psinc4(x)-curr_q},c(-10,10),extendInt="upX")$root})}
.rsinc4 = function(n){.qsinc4(runif(n))}


############################################################
# # code to simulate data

# propensity
.g0 = function(a,w) {(2*a-1)*(0.05 + 0.9*plogis(w$W1)) + 1-a}


#' Simulate data from distributions without bandlimited counterfactual density
#'
#' @description
#' Simulates data from counterfactual density scenarios considered in Luedtke and Chung (2023).
#' @details
#' The counterfactual density under A=1 is defined via the \code{setting} input. Details on the three scenarios listed here can be found in Luedtke and Chung (2023).
#' 
#' The counterfactual density under A=0 is defined via the \code{cos_parameter} input. When \code{cos_parameter=0}, the counterfactual density under A=0 is the same as that under A=1. When \code{cos_parameter} is equal to \code{k>0} and \code{setting='nonzeroBothSides'}, data are simulated from the scenario labeled 'Alt k' in Luedtke and Chung (2023).
#' @param n sample size
#' @param setting one of \code{'spikeLeft'}, \code{'nonzeroBothSides'}, and \code{'zeroBothSides'}.
#' @param cos_parameter Nonnegative integer used to define distribution of the counterfactual density under A=0.
#' @return a list containing
#' \item{W}{data frame with observations in the rows and baseline covariates in the columns.}
#' \item{A}{binary treatment vector, where 1 indicates that an individual was treated.}
#' \item{Y}{vector of outcomes.}
#' \item{Y0}{vector of counterfactual outcomes under A=0.}
#' \item{Y1}{vector of counterfactual outcomes under A=1.}
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @examples
#' sim_data(100)
#' @export
# setting should be one of the names of the items in the true_shape lists
sim_data = function(n,setting='nonzeroBothSides',cos_parameter=1){
  # dimension of feature W
  p = 5
  # shape parameters for uniform mixture of beta distributions used to define density function on [0,1]
  true_shape1 = list(
    'spikeLeft'=c(1,5,4),
    'nonzeroBothSides'=c(1,8,4),
    'zeroBothSides'=c(2,3,4)
  )

  true_shape2 = list(
    'spikeLeft'=c(5,2,8),
    'nonzeroBothSides'=c(1,4,8),
    'zeroBothSides'=c(2,3,4)
  )

  mixture_inds = sample(1:3,n,replace=T)

  # counterfactual outcome under A=1
  Y1 = rbeta(n,true_shape1[[setting]][mixture_inds],true_shape2[[setting]][mixture_inds])

  # counterfactual outcome under A=0
  Y0 = rbeta(n,true_shape1[[setting]][mixture_inds],true_shape2[[setting]][mixture_inds])
  if(cos_parameter>0){
    Y0 = .rcos(n,cos_parameter)
  }
  # Y0 = (1-local_alt_ind)*rbeta(n,true_shape1[[setting]][mixture_inds],true_shape2[[setting]][mixture_inds]) + local_alt_ind*rbeta(n,3,1)

  W = as.data.frame(matrix(rnorm(n*p)*(0.5+Y1) + Y1,ncol=p))
  colnames(W) = paste0("W",1:p)
  # Y1 = with(W,rnorm(n,mean=W1+W2,sd=0.5*sqrt(W2^2+W3^2)))
  A = rbinom(n,1,.g0(1,W))
  Y = A*Y1 + (1-A)*Y0

  return(list(W=W,A=A,Y=Y,Y0=Y0,Y1=Y1))
}

#' Simulate data from distribution with bandlimited counterfactual density
#'
#' @description
#' Simulates data from bandlimited counterfactual density scenario considered in Luedtke and Chung (2023).
#' @details
#' The Fourier transform of the counterfactual density under A=1 has support \code{[-2,2]}.
#' @param n sample size
#' @return a list containing
#' \item{W}{data frame with observations in the rows and baseline covariates in the columns.}
#' \item{A}{binary treatment vector, where 1 indicates that an individual was treated.}
#' \item{Y}{vector of outcomes.}
#' \item{Y0}{vector of counterfactual outcomes under A=0.}
#' \item{Y1}{vector of counterfactual outcomes under A=1.}
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @examples
#' sim_data_bandlimited(100)
#' @export
sim_data_bandlimited = function(n){
  # dimension of feature W
  p = 5
  # parameters to define counterfactual densities
  true_sinc4_shifts = c(-6,-2,5)
  true_sinc4_scales = c(2,3,2)

  Y0 = rnorm(n)
  mixture_inds = sample(1:3,n,replace=T)
  Y1 = true_sinc4_scales[mixture_inds]*.rsinc4(n)+true_sinc4_shifts[mixture_inds]
  W = as.data.frame(matrix(rnorm(n*p)*(0.5+(Y1>0)) + (Y1>0),ncol=p))
  colnames(W) = paste0("W",1:p)
  # Y1 = with(W,rnorm(n,mean=W1+W2,sd=0.5*sqrt(W2^2+W3^2)))
  A = rbinom(n,1,.g0(1,W))
  Y = A*Y1 + (1-A)*Y0

  return(list(W=W,A=A,Y=Y,Y0=Y0,Y1=Y1))
}

##########################################################################

# define a function that can be used to compute the adjoint
# INPUTS
# h_y: vector of length n of evaluations of a function h in L^2 at the observed data points
# h_ygrid: evaluation of a function h in L^2 at the points in ygrid
# outcome_dens_ygrid: an n x ngrid matrix, where entry (i,j) contains the evaluation of the conditional density (x,y)\mapsto p(y|A=1,x) at covariate value X_i in the dataset and entry j in ygrid
.compute_adjoint = function(h_y,h_ygrid,outcome_dens_ygrid_val,outcome_dens_ygrid_train,gn,A_val){
  n_val = length(A_val)
  h_condmean = rowMeans(.row_rep(h_ygrid,n_val)*outcome_dens_ygrid_val)
  h_mean = rowMeans(.row_rep(h_ygrid,nrow(outcome_dens_ygrid_train))*outcome_dens_ygrid_train)
  return((A_val/gn) * (h_y - h_condmean) + h_condmean - mean(h_condmean))
}

############################################################

# if boot_mat is not NULL, then will bootstrap the data according to the values in that matrix.
# boot_mat should be a matrix of dimension n_val x num_boot, with each column denoting the indices of a bootstrap sample
.one_split = function(W_train,A_train,Y_train,W_val,A_val,Y_val,ygrid,dens_init_method,maxK_legendre,maxK_cos,gn_val=NULL,outcome_dens_val=NULL,outcome_dens_train=NULL,beta_params=NULL){
  require(ranger)
  require(orthopolynom)
  n_train = nrow(W_train)
  n_val = nrow(W_val)
  ngrid = length(ygrid)

  ############################################################
  # estimate the propensity
  if(is.null(gn_val)){
    gn_val = predict(ranger(A_train ~ ., data=W_train),W_val)$predictions
  }

  ############################################################
  # get initial fit of outcome density
  if(is.null(outcome_dens_val) | is.null(outcome_dens_train)){
    if(dens_init_method=='kde'){
      require(np)
      bws = npcdensbw(xdat=W_train[A_train==1,],ydat=data.frame(Y=Y_train[A_train==1]),bwmethod='normal-reference')
      # estimate density and then evaluate it on grid
      outcome_dens = matrix(npcdens(bws,exdat=rbind(W_val,W_train)[rep(1:n,ngrid),],eydat=data.frame(Y=rep(ygrid,each=n)))$condens,nrow=n)
    } else if(dens_init_method=='ranger'){
      h = bw.nrd0(Y_train[A_train==1])
      outcome_dens = sapply(ygrid,function(curr_y){
        kern = dnorm((Y_train-curr_y)/h)/h
        return(predict(ranger(kern[A_train==1] ~ ., data=data.frame(W_train[A_train==1,])),rbind(W_val,W_train))$predictions)
      })
    }
    # normalize outcome density estimate so that the rows have mean 1
    outcome_dens = diag(1/rowMeans(outcome_dens)) %*% outcome_dens

    outcome_dens_val = outcome_dens[1:n_val,]
    outcome_dens_train = outcome_dens[(n_val+1):(n_val+n_train),]
  }
  # compute plug-in density on this grid
  plugin_dens = colMeans(outcome_dens_train)

  ############################################################
  # For each estimator, will return the regularized EIF evaluations and values of the estimator (evaluated through maxK for each basis), to be used for confidence set construction
  out_list_regularized_eifs = list()
  out_list_regularized_maxK_ests = list()

  ############################################################

  out_list_ests = list(plugin=plugin_dens)

  ############################################################

  # compute one-step for the specified value of K using the Legendre basis
  legendre_polynomials = polynomial.functions(legendre.polynomials(maxK_legendre, normalized=TRUE ))
  if(!is.null(beta_params)){
    beta = 1/(1+((1:maxK_legendre)/beta_params[2])^(beta_params[1]))
  } else {
    beta = rep(1,maxK_legendre)
  }
  regularized_eif_legendre = Reduce('+',
                                   lapply(1:maxK_legendre,function(k){
                                     basis_fun = function(x){sqrt(2)*legendre_polynomials[[k]](2*x-1)} # multiplying by sqrt(2) to keep normalized (needed because of the change of variables)
                                     basis_fun_y = basis_fun(Y_val)
                                     basis_fun_ygrid = basis_fun(ygrid)
                                     # first row of output is the k-th term of the sum for the regularized parameter evaluation, and other rows are the regularized eif evaluations at each observation
                                     rbind(rbind(beta[k] * mean(plugin_dens*basis_fun_ygrid) * basis_fun_ygrid),
                                           (beta[k] * cbind(.compute_adjoint(basis_fun_y,basis_fun_ygrid,outcome_dens_val,outcome_dens_train,gn_val,A_val))) %*% rbind(basis_fun_ygrid))
                                   }),accumulate=T)
  # density estimators for 1:maxK_legendre basis functions
  onestep_legendre_dens = lapply(regularized_eif_legendre,function(z){colMeans(z[-1,]) + plugin_dens})
  names(onestep_legendre_dens) = paste0('legendre_',1:maxK_legendre)
  out_list_ests = append(out_list_ests,onestep_legendre_dens)
  # regularized density estimator for maxK
  onestep_legendre_regularized_dens = onestep_legendre_dens[[length(onestep_legendre_dens)]] - plugin_dens + regularized_eif_legendre[[length(onestep_legendre_dens)]][1,]
  out_list_regularized_maxK_ests = append(out_list_regularized_maxK_ests,list(legendre=onestep_legendre_regularized_dens))
  # only save regularized EIF of maxK estimator
  regularized_eif_legendre = regularized_eif_legendre[[maxK_legendre]][-1,]
  out_list_regularized_eifs = append(out_list_regularized_eifs,list(legendre=regularized_eif_legendre))

  ############################################################

  # compute one-step for the specified value of K using the cosine basis
  if(!is.null(beta_params)){
    beta = 1/(1+((1:maxK_cos)/beta_params[2])^(beta_params[1]))
  } else {
    beta = rep(1,maxK_cos)
  }
  regularized_eif_cos = Reduce('+',
                            lapply(1:maxK_cos,function(k){
                              basis_fun = function(x){.my_cosine(x,k)}
                              basis_fun_y = basis_fun(Y_val)
                              basis_fun_ygrid = basis_fun(ygrid)
                              rbind(rbind(beta[k] * mean(plugin_dens*basis_fun_ygrid) * basis_fun_ygrid),
                                    (beta[k] * cbind(.compute_adjoint(basis_fun_y,basis_fun_ygrid,outcome_dens_val,outcome_dens_train,gn_val,A_val))) %*% rbind(basis_fun_ygrid))
                              }),accumulate=T)
  # density estimators for 1:maxK_cos basis functions
  onestep_cos_dens = lapply(regularized_eif_cos,function(z){colMeans(z[-1,]) + plugin_dens})
  names(onestep_cos_dens) = paste0('cos_',1:maxK_cos)
  out_list_ests = append(out_list_ests,onestep_cos_dens)
  # regularized density estimator for maxK
  onestep_cos_regularized_dens = onestep_cos_dens[[length(onestep_cos_dens)]] - plugin_dens + regularized_eif_cos[[length(onestep_cos_dens)]][1,]
  out_list_regularized_maxK_ests = append(out_list_regularized_maxK_ests,list(cos=onestep_cos_regularized_dens))
  # only save regularized EIF of maxK estimator
  regularized_eif_cos = regularized_eif_cos[[maxK_cos]][-1,]
  out_list_regularized_eifs = append(out_list_regularized_eifs,list(cos=regularized_eif_cos))

  return(list(ests=out_list_ests,regularized_eifs=out_list_regularized_eifs,regularized_maxK_ests=out_list_regularized_maxK_ests))
}

# if supplied, the fold argument be a list of length equal to n=nrow(W)=length(A)=length(Y). Will subsume
# num_folds if supplied
# objective should be either 'estimation' or 'inference'.
#   If 'estimation', then returns a matrix with estimators in rows and evaluations of the estimates of the density at the values in ygrid in the columns
#   If 'inference', then returns a list of lists. Each element in the outer list pertains to a fold, inner list contains output of .one_split on that fold. Key quantities of interest in that list are a matrix of regularized EIF evaluations (observations in rows, ygrid values in columns) and a vector of regularized parameter estimates along the values in ygrid
.cf_est = function(W,A,Y,grid=500,dens_init_method = 'ranger',maxK_legendre=10,maxK_cos=10,num_fold=5,fold=NULL,gn=NULL,outcome_dens_list=NULL,objective='estimation',beta_params=NULL){
  n = nrow(W)
  if(length(grid)==1){ # if grid is of length 1, then it's assumed that it's an integer, and a grid on [0,1] is defined base don this integer
    ygrid = seq(0,1,length=ngrid+2)[2:(ngrid+1)]
  } else {
    ygrid = grid
  }
  if(is.null(fold)){
    fold = sample(rep(1:num_fold,ceiling(n/num_fold))[1:n])
  } else{
    num_fold = length(unique(fold))
  }
  out = lapply(sort(unique(fold)),function(fold_ind){
    if(!is.null(gn)){
      gn_val = gn[fold==fold_ind]
    } else{
      gn_val = NULL
    }
    if(!is.null(outcome_dens_list)){
      outcome_dens_val = outcome_dens_list[[as.character(fold_ind)]][fold==fold_ind,]
      outcome_dens_train = outcome_dens_list[[as.character(fold_ind)]][fold!=fold_ind,]
    } else {
      outcome_dens_val = NULL
      outcome_dens_train = NULL
    }
    .one_split(W[fold!=fold_ind,],A[fold!=fold_ind],Y[fold!=fold_ind],W[fold==fold_ind,],A[fold==fold_ind],Y[fold==fold_ind],ygrid,dens_init_method,maxK_legendre=maxK_legendre,maxK_cos=maxK_cos,gn_val=gn_val,outcome_dens_val=outcome_dens_val,outcome_dens_train=outcome_dens_train,beta_params=beta_params)
  })

  if(objective=='estimation'){
    out = Reduce('+',lapply(out,function(z){
      do.call(rbind,z$ests)
    }))/num_fold
    # normalize output so that each row has mean 1
    tmp = rownames(out)
    out = diag(1/rowMeans(out)) %*% out
    rownames(out) = tmp
  }

  return(out)
}

#' Regularized one-step estimator of counterfactual density function
#'
#' @description
#' Regularized one-step estimator of counterfactual density function under treatment A=1. Number of basis terms to include in the bias correction step is selected by cross-validation.
#' @details
#' Propensity score is estimated using the ranger package.
#' 
#' An alternative approach for estimating this function is implemented in the npcausal R package (Kennedy et al., 2021).
#' @param W data frame with observations in the rows and baseline covariates in the columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y outcome taking values between 0 and 1. If Y is not bounded in this interval, then should rescale it before calling this function so that it is.
#' @param ngrid number of grid points used in approximation procedure. Selecting larger values of ngrid should yield better performance, but slower runtime.
#' @param dens_init_method method used to obtain initial estimates of the counterfactual density. Should be either 'ranger', in which case random forest is used, or 'kde', in which case a kernel density estimator is used.
#' @param maxK_legendre Maximum number of functions in the Legendre basis to consider during cross-validation.
#' @param maxK_cos Maximum number of functions in the cosine basis to consider during cross-validation.
#' @param num_fold_final Number of folds used during cross-fitting of final estimator. Should be at least 2.
#' @return a list containing
#' \item{best_fun}{A function taking as input a vector of outcome values and outputting one-step estimates of the counterfactual density function at those values. Estimated by cross-validation to be the best of any of the estimates considered: the plug-in estimator, and one-step estimators that use the Legendre or cosine basis.}
#' \item{best_legendre_fun}{Same as the above, but only focusing on one-step estimators that use the Legendre basis.}
#' \item{best_cos_fun}{Same as the above, but only focusing on one-step estimators that use the cosine basis.}
#' \item{which_best}{Estimator selected for best_fun.}
#' \item{which_best_legendre}{Estimator selected for best_legendre_fun}
#' \item{which_best_cos}{Estimator selected for best_cos_fun}
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @references E. H. Kennedy, S. Balakrishnan, and L. Wasserman, ``Semiparametric counterfactual density estimation,'' \emph{arXiv}, 2021.
#' @examples
#' # sample size
#' n = 500
#'
#' # simulate data
#' dat = sim_data(n,setting='zeroBothSides')
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' cv_out = cv_density(W,A,Y,ngrid=500,num_fold_final=2)
#' cv_out$which_best
#' ygrid = seq(0.0025,0.9975,by=0.0025)
#' plot(ygrid,cv_out$best_fun(ygrid),xlab='y',ylab='Counterfactual Density Estimate',type='l')
#' @export

cv_density = function(W,A,Y,ngrid=500,dens_init_method = 'ranger',maxK_legendre=16,maxK_cos=16,num_fold_final=5){
  require(ranger)
  require(combinat)
  if(max(Y)>1 | min(Y)<0){
    stop("Y values should all fall in the interval [0,1].")
  }
  n = nrow(W)
  num_folds_cv = 4
  fold_cv = sample(rep(1:num_folds_cv,ceiling(n/num_folds_cv))[1:n])
  ygrid = seq(0,1,length=ngrid+2)[2:(ngrid+1)]

  print("computing nuisance functions used to evaluate cross-validated risk")
  # nuisance estimation
  # propensity
  gn_list = lapply(1:num_folds_cv,function(fold_ind){
    predict(ranger(A[fold_cv==fold_ind] ~ ., data=W[fold_cv==fold_ind,]),W)$predictions})
  # outcome density
  outcome_dens_list = lapply(1:num_folds_cv,function(fold_ind){
    W_train = W[fold_cv==fold_ind,]
    A_train = A[fold_cv==fold_ind]
    Y_train = Y[fold_cv==fold_ind]
    if(dens_init_method=='kde'){
      require(np)
      bws = npcdensbw(as.formula(paste0('Y ~ ',paste0(colnames(W),collapse="+"))),data=data.frame(W_train,Y=Y_train)[A_train==1,],bwmethod='normal-reference')
      dens_fit = npcdens(bws)
      # evaluate density on this grid
      outcome_dens_out = matrix(predict(dens_fit,newdata=data.frame(W[rep(1:n,ngrid),],Y=rep(ygrid,each=n))),nrow=n)
    } else if(dens_init_method=='ranger'){
      h = bw.nrd0(Y_train[A_train==1])
      outcome_dens_out = sapply(ygrid,function(curr_y){
        kern = dnorm((Y_train-curr_y)/h)/h
        return(predict(ranger(kern[A_train==1] ~ ., data=data.frame(W_train[A_train==1,])),W)$predictions)
      })
    }

    # normalize outcome density estimate so that the rows have mean 1
    outcome_dens_out = diag(1/rowMeans(outcome_dens_out)) %*% outcome_dens_out

    return(outcome_dens_out)
  })


  print("computing cross-validated risk")
  # permute over the 24 combinations of c(1,2,3,4)
  cv_risk = rowMeans(sapply(permn(1:4),function(fold_order){
      # indices to use for estimation
      use_for_est = fold_cv %in% c(fold_order[1],fold_order[3])

      # nuisance functions to use for estimation
      gn_est = rep(NA,n)
      gn_est[fold_cv==fold_order[1]] = gn_list[[fold_order[3]]][fold_cv==fold_order[1]]
      gn_est[fold_cv==fold_order[3]] = gn_list[[fold_order[1]]][fold_cv==fold_order[3]]
      gn_est = gn_est[use_for_est]
      outcome_dens_list_curr = list()
      outcome_dens_list_curr[[as.character(fold_order[[3]])]] = outcome_dens_list[[fold_order[3]]][use_for_est,]
      outcome_dens_list_curr[[as.character(fold_order[[1]])]] = outcome_dens_list[[fold_order[1]]][use_for_est,]
      # cross-fitted estimators
      ests = .cf_est(W[use_for_est,],A[use_for_est],Y[use_for_est],grid=ygrid,dens_init_method=dens_init_method,maxK_legendre=maxK_legendre,maxK_cos=maxK_cos,fold=fold_cv[use_for_est],gn=gn_est,outcome_dens_list=outcome_dens_list_curr)

      # compute plugin based on data in fold_order[2]
      plugin_dens = colMeans(outcome_dens_list[[fold_order[2]]][fold_cv==fold_order[2],])
      return(apply(ests,1,function(est){
        h_ygrid = est-plugin_dens
        # since we've only evaluated the densities on a grid, approximate the evaluations at the points in Y[fold_cv==fold_order[4]] by just evaluating at the nearest point on the grid
        h_y = sapply(Y[fold_cv==fold_order[4]],function(y_curr){h_ygrid[which.min(abs(y_curr-ygrid))]})
        0.5*mean((h_ygrid)^2) - mean(.compute_adjoint(h_y,h_ygrid,outcome_dens_list[[fold_order[2]]][fold_cv==fold_order[4],],outcome_dens_list[[fold_order[2]]][fold_cv==fold_order[2],],gn_list[[fold_order[2]]][fold_cv==fold_order[4]],A[fold_cv==fold_order[4]]))
      }))
  }))
  which_best = names(cv_risk)[which.min(cv_risk)]
  legendre_inds = which(substr(names(cv_risk),0,8)=='legendre')
  which_best_legendre = names(cv_risk)[legendre_inds][which.min(cv_risk[legendre_inds])]
  cos_inds = which(substr(names(cv_risk),0,3)=='cos')
  which_best_cos = names(cv_risk)[cos_inds][which.min(cv_risk[cos_inds])]

  print(paste0("the estimator selected by cross-validation is ",which_best))

  print("computing estimators based on the full dataset")
  all_ests = .cf_est(W,A,Y,grid=ygrid,dens_init_method = dens_init_method,maxK_legendre=maxK_legendre,maxK_cos=maxK_cos,num_fold=num_fold_final,fold=NULL,gn=NULL,outcome_dens_list=NULL)
  best_fun = approxfun(c(0,ygrid,1),c(0,all_ests[which_best,],0),rule=2)
  best_legendre_fun = approxfun(c(0,ygrid,1),c(0,all_ests[which_best_legendre,],0),rule=2)
  best_cos_fun = approxfun(c(0,ygrid,1),c(0,all_ests[which_best_cos,],0),rule=2)

  return(list(best_fun=best_fun,best_legendre_fun=best_legendre_fun,best_cos_fun=best_cos_fun,which_best=which_best,which_best_legendre=which_best_legendre,which_best_cos=which_best_cos))
}


#' Test for equality of counterfactual densities
#'
#' @description
#' Test for the equality of counterfactual densities under treatments A=1 and A=0.
#' @details
#' Defined based on a one-step estimator using the regularized spherical or Wald-type test defined in the reference. These tests either take a Wald-type or spherical form, and defining them requires selecting an orthonormal basis of L2(\code{[0,1]}). Here, we report results for the cosine basis and a transformation of the Legendre basis.
#' 
#' We recommend basing inference on the Wald-type test as defined with the cosine basis.
#' 
#' Propensity score is estimated using the ranger package.
#' @param W data frame with observations in the rows and baseline covariates in the columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y outcome taking values between 0 and 1. If Y is not bounded in this interval, then should rescale it before calling this function so that it is.
#' @param ngrid number of grid points used in approximation procedure. Selecting larger values of ngrid should yield better performance, but slower runtime.
#' @param dens_init_method method used to obtain initial estimates of the two counterfactual densities. Should be either 'ranger', in which case random forest is used, or 'kde', in which case a kernel density estimator is used.
#' @param maxK_legendre Number of functions in the Legendre basis to use when constructing the test. In practice we've observed numerical instability when this value is taken much larger than 50.
#' @param maxK_cos Number of functions in the cosine basis to use when constructing the test. We have not observed any numerical instability when this value is taken large, though selecting very large values will increase the runtime and memory requirements.
#' @param num_fold Number of folds used during cross-fitting. Should be at least 2.
#' @param num_boot Number of bootstrap replications used to calculate p-value.
#' @param beta_params Regularization parameters used when defining the test statistic. The k-th term of the basis expansion is weighted by \code{1/(1+(k/beta_params[2])^beta_params[1])}. Test only asymptotically valid if \code{beta_params[2]>0} and \code{beta_params[1]>1/2}.
#' @param lambda Parameter used to define the standardization operator used to define the correlation-based Wald-type test. Should take a strictly positive value that's less than 1. Values too close to zero may lead to numerical instability, and values close to 1 will lead a test that's similar to the spherical test. See the reference for details.
#' @return The p-values of the spherical and Wald-type tests based on the cosine and Legendre bases.
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @examples
#' # sample size
#' n = 250
#'
#' # simulate data from alternative
#' dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=1)
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' density_test(W,A,Y,ngrid=500,num_fold=2,num_boot=2000)
#' 
#' 
#' # simulate data from null
#' dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=0)
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' density_test(W,A,Y,ngrid=500,num_fold=2,num_boot=2000)
#' @export

density_test = function(W,A,Y,ngrid=500,dens_init_method = 'ranger',maxK_legendre=50,maxK_cos=200,num_fold=2,num_boot=2000,beta_params=c(2,5),lambda=0.5){
  require(ranger)
  if(max(Y)>1 | min(Y)<0){
    stop("Y values should all fall in the interval [0,1]. If Y is a bounded random variable, then please linearly rescale the Y values to fall in this range before calling this function.")
  }
  n = nrow(W)
  ygrid = seq(0,1,length=ngrid+2)[2:(ngrid+1)]
  fold = sample(rep(1:num_fold,ceiling(n/num_fold))[1:n])

  cf_out = lapply(0:1,function(a){
    .cf_est(W,a*A+(1-a)*(1-A),Y,grid=ygrid,dens_init_method=dens_init_method,maxK_legendre=maxK_legendre,maxK_cos=maxK_cos,num_fold=num_fold,fold=fold,gn=NULL,outcome_dens_list=NULL,objective='inference',beta_params=beta_params)
  })
  basis_names = names(cf_out[[1]][[1]]$regularized_eifs)

  Omega = lapply(basis_names,function(curr_basis_name){
    solve(Reduce('+',lapply(1:num_fold,function(fold_ind){
      # (1-lambda)*cor(cf_out[[2]][[fold_ind]]$regularized_eifs[[curr_basis_name]] - cf_out[[1]][[fold_ind]]$regularized_eifs[[curr_basis_name]])+lambda*diag(ngrid)
      # cov(cf_out[[2]][[fold_ind]]$regularized_eifs[[curr_basis_name]] - cf_out[[1]][[fold_ind]]$regularized_eifs[[curr_basis_name]])*(matrix(1-lambda,nrow=ngrid,ncol=ngrid)+lambda*diag(ngrid))
      out = cor(cf_out[[2]][[fold_ind]]$regularized_eifs[[curr_basis_name]] - cf_out[[1]][[fold_ind]]$regularized_eifs[[curr_basis_name]])
      (1-lambda)*out+lambda*diag(ngrid)
      }))/num_fold)
  })
  names(Omega) = basis_names

  test_statistics = n*sapply(basis_names,function(curr_basis_name){
    point_ests = Reduce('+',lapply(1:num_fold,function(fold_ind){
      cf_out[[2]][[fold_ind]]$regularized_maxK_ests[[curr_basis_name]] - cf_out[[1]][[fold_ind]]$regularized_maxK_ests[[curr_basis_name]]
    }))/num_fold
    return(c(
      mean(point_ests^2), # squared Hilbert norm
      rbind(point_ests)%*%Omega[[curr_basis_name]]%*%cbind(point_ests)/length(ygrid))) # regularized Wald-type
  })
  rownames(test_statistics) = c('spherical','Wald')

  p_vals = sapply(basis_names,function(curr_basis_name){
    rowMeans(test_statistics[,curr_basis_name,drop=F][,rep(1,num_boot)] <= n*sapply(1:num_boot,function(boot_ind){
      boot_draw = Reduce('+',lapply(1:num_fold,function(fold_ind){
        n_val = sum(fold==fold_ind)
        wgt = c(rmultinom(1,n_val,rep(1/n_val,n_val))-1)
        colMeans((cf_out[[2]][[fold_ind]]$regularized_eifs[[curr_basis_name]] - cf_out[[1]][[fold_ind]]$regularized_eifs[[curr_basis_name]])*cbind(wgt)[,rep(1,ngrid)])
      }))/num_fold
      return(c(
        mean(boot_draw^2), # squared Hilbert norm
        rbind(boot_draw)%*%Omega[[curr_basis_name]]%*%cbind(boot_draw)/length(ygrid))) # regularized Wald-type threshold
    }))
  })
  rownames(p_vals) = c('spherical','Wald')
  
  return(p_vals)
}

#' Maximum mean discrepancy between counterfactual distributions
#'
#' @description
#' Maximum mean discrepancy between counterfactual distributions under treatments A=1 and A=0.
#' @details
#' Maximum mean discrepancy estimator defined by computing the distance between the one-step estimators of the counterfactual kernel mean embeddings under A=1 and A=0.
#' 
#' All needed nuisances are estimated using the ranger package.
#' @param W data frame with observations in the rows and baseline covariates in the columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y outcome taking values between 0 and 1. If Y is not bounded in this interval, then should rescale it before calling this function so that it is.
#' @param ngrid number of grid points used in approximation procedure. Selecting larger values of ngrid should yield better performance, but slower runtime.
#' @param num_fold Number of folds used during cross-fitting. Should be at least 2.
#' @param num_boot Number of bootstrap replications used to calculate p-value.
#' @param bw_mult The reported MMD is based on the Gaussian kernel with bandwidth equal to the median heuristic -- namely, the median of the off-diagonal elements of the distance matrix dist(Y) -- times bw_mult.
#' @return The p-value of the test.
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @references J. Fawkes, R. Hu, R. J. Evans, and D. Sejdinovic, ``Doubly Robust Kernel Statistics for Testing Distributional Treatment Effects Even Under One Sided Overlap,'' \emph{arXiv}, 2022.
#' @references M. Muandet, M. Kanagawa, S.  Saengkyongam, and S. Marukatat, ``Counterfactual mean embeddings,'' \emph{J Mach Learn Res}, 2021.
#' @examples
#' # sample size
#' n = 250
#'
#' # simulate data from alternative
#' dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=1)
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' mmd_test(W,A,Y,ngrid=500,num_fold=2,num_boot=2000,bw_mult=1)
#' 
#' 
#' # simulate data from null
#' dat = sim_data(n,setting='nonzeroBothSides',cos_parameter=0)
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' mmd_test(W,A,Y,ngrid=500,num_fold=2,num_boot=2000,bw_mult=1)
#' @export

mmd_test = function(W,A,Y,ngrid=500,num_fold=2,num_boot=2000,bw_mult=1){
  require(ranger)
  if(max(Y)>1 | min(Y)<0){
    stop("Y values should all fall in the interval [0,1]. If Y is a bounded random variable, then please linearly rescale the Y values to fall in this range before calling this function.")
  }
  n = nrow(W)
  ygrid = seq(0,1,length=ngrid+2)[2:(ngrid+1)]
  fold = sample(rep(1:num_fold,ceiling(n/num_fold))[1:n])

  # define Gram matrix of c(Y,ygrid)
  G = as.matrix(dist(c(Y,ygrid)))
  G = dnorm(G,mean=0,sd=bw_mult*median(G[1:n,1:n][upper.tri(G[1:n,1:n])]))

  boot_out = Reduce('+',lapply(1:num_fold,function(fold_ind){
    # training and "validation" data to be used for cross-fitting
    W_train = W[fold==fold_ind,]
    A_train = A[fold==fold_ind]
    Y_train = Y[fold==fold_ind]
    val_inds = which(fold!=fold_ind)
    W_val = W[val_inds,]
    A_val = A[val_inds]
    Y_val = Y[val_inds]

    # sample sizes
    n_train = nrow(W_train)
    n_val = nrow(W_val)

    # estimate the propensity
    gn_val = predict(ranger(A_train ~ ., data=W_train),W_val)$predictions

    # estimate outcome density using training data
    outcome_dens_01 = lapply(0:1,function(a){
      h = bw.nrd0(Y_train[A_train==a])
      return(sapply(ygrid,function(curr_y){
        kern = dnorm((Y_train-curr_y)/h)/h
        return(predict(ranger(kern[A_train==a] ~ ., data=data.frame(W_train[A_train==a,])),rbind(W_val,W_train))$predictions)
      }))
    })

    # number of times each observation was sampled by bootstrap (or, in the first column, a vector of 1s used to compute the nonbootstrapped estimator)
    wgt_matrix = cbind(rep(1,n_val),rmultinom(num_boot,n_val,rep(1/n_val,n_val))-1)

    est_vecs = apply(wgt_matrix,2,function(wgt){
      # will be used to compute the quadratic form that returns the MMD^2 (and the H-squared norm of the EIF)
      est_vec = rep(0,n+ngrid)
      # define a vector that can be used to evaluate estimator of MMD^2. In particular, an estimate of MMD^2 is given by the quadratic form rbind(est_vec) %*% G $*$ cbind(est_vec)
      est_vec[1:n][val_inds] = est_vec[1:n][val_inds] + 1/n_val * wgt * (2*A_val-1)/(A_val*gn_val + (1-A_val)*(1-gn_val))
      est_vec[(n+1):(n+ngrid)] = est_vec[(n+1):(n+ngrid)] - 1/(n_val*ngrid) * rbind(wgt*(A_val/gn_val - 1)) %*% outcome_dens_01[[2]][1:n_val,]
      est_vec[(n+1):(n+ngrid)] = est_vec[(n+1):(n+ngrid)] + 1/(n_val*ngrid) * rbind(wgt*((1-A_val)/(1-gn_val) - 1)) %*% outcome_dens_01[[1]][1:n_val,]
      return(est_vec)
    })

    # the below would mean center the EIF, but it isn't needed for bootstrap (since the bootstrapped empirical distribution is centered by the empirical, and so this constant would cancel out)
    # bootstrapped_mean_eif = est_vecs
    # bootstrapped_mean_eif[(n+1):(n+ngrid),] = bootstrapped_mean_eif[(n+1):(n+ngrid),] - cbind(1/ngrid * colMeans((outcome_dens_01[[2]]-outcome_dens_01[[1]])[(n_val+1):(n_val+n_train),]))[,rep(1,ncol(bootstrapped_mean_eif))]

    # return(cbind(est_vecs[,1],bootstrapped_mean_eif[,-1]))
    return(est_vecs)
  }))/num_fold

  est_vec = boot_out[,1]
  boot_out = boot_out[,-1]

  test_statistic = c(n*rbind(est_vec)%*%G%*%cbind(est_vec))
  p_val = mean(test_statistic <= n*apply(boot_out,2,function(boot_est_vec){rbind(boot_est_vec)%*%G%*%cbind(boot_est_vec)}))
  names(p_val) = 'p-value'
  return(p_val)
}




#' One-Step Estimator of a Bandlimited Counterfactual Density
#'
#' @description
#' This function estimates a bandlimited counterfactual density function under treatment A=1.
#' @details
#' If the counterfactual density function is not bandlimited at the specified value \code{b}, then a \code{b}-bandlimiting of the density function is estimated instead. See the reference for details.
#'
#' All needed nuisances are estimated using the ranger package.
#' @param W data frame with observations in the rows and baseline covariates in the columns.
#' @param A binary treatment vector, where 1 indicates that an individual was treated.
#' @param Y real-valued outcome.
#' @param ngrid number of grid points used in approximation procedure. Selecting larger values of ngrid should yield better performance, but slower runtime.
#' @param num_fold Number of folds used during cross-fitting. Should be at least 2.
#' @param alpha 1-alpha is the desired level of confidence set. Can be a single number between 0 and 1, or a vector of such numbers.
#' @param num_boot Number of bootstrap replications used to estimate appropriate L2 radius of confidence set.
#' @param b Upper bound on the support of the Fourier transform of the counterfactual density function. If this Fourier transform has support outside of \code{[-b,b]}, then the estimand is instead the b-bandlimiting of this function.
#' @return a list containing
#' \item{onestep_fun}{A function taking as input a vector of outcome values and outputting one-step estimates of the b-bandlimited counterfactual density function at those values.}
#' \item{plugin_fun}{Same as the above, but the estimates are obtained based on a plug-in estimator. Provided for comparison purposes only --- \code{onestep_fun} is the preferred estimator.}
#' \item{L2_radius_sq}{A vector of the squared L2 radiuses of alpha-level confidence sets for the b-bandlimited counterfactual density function. These confidence sets consist of all b-bandlimited functions contained in L2 balls centered at \code{onestep_fun} with the specified radii.}
#' @references A. Luedtke and I. Chung, ``One-Step Estimation of Differentiable Hilbert-Valued Parameters,'' \emph{arXiv}, 2023.
#' @examples
#' # sample size
#' n = 500
#'
#' # simulate data
#' dat = sim_data_bandlimited(n)
#' W = dat$W
#' A = dat$A
#' Y = dat$Y
#'
#' bl_out = bandlimited_density(W,A,Y,ngrid=500,num_fold=2,num_boot=2000,alpha=seq(0.01,0.99,by=0.01),b=2)
#' ygrid = seq(-15,15,by=0.01)
#' plot(ygrid,bl_out$onestep_fun(ygrid),xlab='y',ylab='Counterfactual Density Estimate',type='l')
#' @export

bandlimited_density = function(W,A,Y,ngrid=1000,num_fold=2,alpha=0.05,num_boot=2000,b=10){
  require(ranger)
  n = nrow(W)
  ygrid = qnorm((1:ngrid)/(ngrid+1),mean=mean(Y[A==1]),sd=4*sd(Y[A==1]))
  ygrid_wgts = 1/dnorm(ygrid,mean=mean(Y[A==1]),sd=4*sd(Y[A==1]))
  fold = sample(rep(1:num_fold,ceiling(n/num_fold))[1:n])

  # define Gram matrix of c(Y,ygrid)
  G = (b/pi)*.sinc(b*as.matrix(dist(c(Y,ygrid))))

  boot_out = Reduce('+',lapply(1:num_fold,function(fold_ind){
    # training and "validation" data to be used for cross-fitting
    W_train = W[fold==fold_ind,]
    A_train = A[fold==fold_ind]
    Y_train = Y[fold==fold_ind]
    val_inds = which(fold!=fold_ind)
    W_val = W[val_inds,]
    A_val = A[val_inds]
    Y_val = Y[val_inds]

    # sample sizes
    n_train = nrow(W_train)
    n_val = nrow(W_val)

    # estimate the propensity
    gn_val = predict(ranger(A_train ~ ., data=W_train),W_val)$predictions

    # estimate outcome density using training data
    h = bw.nrd0(Y_train[A_train==1])
    outcome_dens = sapply(ygrid,function(curr_y){
      kern = dnorm((Y_train-curr_y)/h)/h
      return(predict(ranger(kern[A_train==1] ~ ., data=data.frame(W_train[A_train==1,])),rbind(W_val,W_train))$predictions)
    })

    init_dens_wgts = c(rep(0,n),c(ygrid_wgts/(n_val*ngrid) * rbind(rep(1,n_val))%*%outcome_dens[1:n_val,]))

    # number of times each observation was sampled by bootstrap (or, in the first column, a vector of 1s used to compute the nonbootstrapped estimator)
    wgt_matrix = cbind(rep(1,n_val),rmultinom(num_boot,n_val,rep(1/n_val,n_val))-1)

    est_vecs = apply(wgt_matrix,2,function(wgt){
      # will be used to compute the quadratic form that returns the squared L^2 distance from the estimated density to the true density (and the H-squared norm of the EIF)
      est_vec = rep(0,n+ngrid)
      est_vec[1:n][val_inds] = est_vec[1:n][val_inds] + 1/n_val * wgt * A_val/gn_val
      est_vec[(n+1):(n+ngrid)] = est_vec[(n+1):(n+ngrid)] - ygrid_wgts/(n_val*ngrid) * rbind(wgt*(A_val/gn_val - 1)) %*% outcome_dens[1:n_val,]
      return(est_vec)
    })

    # the below would mean center the EIF, but it isn't needed for bootstrap (since the bootstrapped empirical distribution is centered by the empirical, and so this constant would cancel out)
    # bootstrapped_mean_eif = est_vecs
    # bootstrapped_mean_eif[(n+1):(n+ngrid),] = bootstrapped_mean_eif[(n+1):(n+ngrid),] - cbind(ygrid_wgts/ngrid * colMeans(outcome_dens[(n_val+1):(n_val+n_train),]))[,rep(1,ncol(bootstrapped_mean_eif))]

    return(cbind(est_vecs[,1],init_dens_wgts,est_vecs[,-1]))
  }))/num_fold

  est_vec = boot_out[,1]
  plugin_vec = boot_out[,2]
  boot_out = boot_out[,-c(1,2)]

  est_ygrid = (G%*%cbind(est_vec))[(n+1):(n+ngrid)]
  plugin_ygrid = (G%*%cbind(plugin_vec))[(n+1):(n+ngrid)]

  onestep_fun = approxfun(ygrid,est_ygrid,rule=2)
  plugin_fun = approxfun(ygrid,plugin_ygrid,rule=2)

  # squared L^2 radius of confidence set
  alpha = sort(alpha,decreasing=T)
  L2_radius_sq = quantile(apply(boot_out,2,function(boot_est_vec){rbind(boot_est_vec)%*%G%*%cbind(boot_est_vec)}),1-alpha)
  names(L2_radius_sq) = paste0('alpha=',alpha)
  return(list(onestep_fun=onestep_fun,plugin_fun=plugin_fun,L2_radius_sq=L2_radius_sq))
}


