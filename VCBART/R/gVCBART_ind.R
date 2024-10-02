gVCBART_ind <- function(Y_train,
                        subj_id_train, ni_train, X_train,
                        W_train = matrix(0, nrow = 1, ncol = 1),
                        Z_cont_train = matrix(0, nrow = 1, ncol = 1),
                        Z_cat_train = matrix(0L, nrow = 1, ncol = 1),
                        X_test = matrix(0, nrow = 1, ncol = 1),
                        W_test = matrix(0, nrow = 1, ncol = 1),
                        Z_cont_test = matrix(0, nrow = 1, ncol = 1),
                        Z_cat_test = matrix(0, nrow = 1, ncol = 1),
                        unif_cuts = rep(TRUE, times = ncol(Z_cont_train)),
                        cutpoints_list = NULL,
                        cat_levels_list = NULL,
                        edge_mat_list = NULL,
                        graph_split = rep(FALSE, times = ncol(Z_cat_train)),
                        sparse = TRUE, 
                        M = 50,
                        mu0 = NULL, tau = NULL, nu = NULL, lambda = NULL,
                        gamma0 = NULL, xi = NULL,
                        nd = 1000, burn = 1000, thin = 1,
                        save_samples = TRUE, save_trees = TRUE,
                        verbose = TRUE, print_every = floor( (nd*thin + burn)/10))
{
  # standardize Y
  y_mean <- mean(Y_train)
  y_sd <- stats::sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd 
  
  if(!is.matrix(X_train)) stop("X_train must be a matrix!")
  if(!is.matrix(W_train)) stop("W_train must be a matrix!")
  if(!is.matrix(Z_cont_train)) stop("Z_cont_train must be a matrix!")
  if(!is.matrix(Z_cat_train)) stop("Z_cat_train must be an integer matrix!")
  
  if(all(X_train[,1] == 1)) stop("Do not include intercept as first column of X_train")
  
  # ASSUME: first column of X_train is not the intercept
  std_X_train <- matrix(NA, nrow = nrow(X_train), ncol = 1 + ncol(X_train))
  std_X_train[,1] <- 1
  x_mean <- rep(NA, times = 1+ncol(X_train))
  x_sd <- rep(NA, times = 1+ncol(X_train))
    
  x_mean[1] <- 0
  x_sd[1] <- 1
  for(j in 1:ncol(X_train)){
    x_mean[j+1] <- mean(X_train[,j])
    x_sd[j+1] <- stats::sd(X_train[,j])
    std_X_train[,j+1] <- (X_train[,j] - x_mean[j+1])/x_sd[j+1]
  }
  
  if(length(X_test) > 1){
    # There is a valid X_test
    if(ncol(X_test) != ncol(X_train)){
      print(paste("X_train has", ncol(X_train), "columns but X_test has", ncol(X_test), "columns"))
      stop("X_train and X_test must have the same number of columns")
    } else{
      # standardize X_test using the column means and sd's from the training data
      std_X_test <- matrix(NA, nrow = nrow(X_test), ncol = 1 + ncol(X_test))
      std_X_test[,1] <- 1
      for(j in 1:ncol(X_test)) std_X_test[,j+1] <- (X_test[,j] - x_mean[j+1])/x_sd[j+1]
    }
  } else{
    # No X_test provided
    std_X_test <- matrix(0, nrow = 1, ncol = 1)
  }
  
  
  std_W_train <- matrix(NA, nrow = nrow(W_train), ncol = ncol(W_train))
  w_mean <- rep(NA, times = ncol(W_train))
  w_sd <- rep(NA, times = ncol(W_train))
  
  for(j in 1:ncol(W_train)){
    w_mean[j] <- mean(W_train[,j])
    w_sd[j] <- stats::sd(W_train[,j])
    std_W_train[,j] <- (W_train[,j] - w_mean[j])/w_sd[j]
  }
  
  if(length(W_test) > 1){
    # There is a valid W_test
    if(ncol(W_test) != ncol(W_train)){
      print(paste("W_train has", ncol(W_train), "columns but W_test has", ncol(W_test), "columns"))
      stop("W_train and W_test must have the same number of columns")
    } else{
      # standardize W_test using the column means and sd's from the training data
      std_W_test <- matrix(NA, nrow = nrow(W_test), ncol = ncol(W_test))
      for(j in 1:ncol(W_test)) std_W_test[,j] <- (W_test[,j] - w_mean[j])/w_sd[j]
    }
  } else{
    # No W_test provided
    std_W_test <- matrix(0, nrow = 1, ncol = 1)
  }

  # Set some hyperparameters
  
  if(is.null(mu0)) mu0 <- rep(0, times = ncol(std_X_train))
  else if(length(mu0) != ncol(std_X_train)) stop("mu0 needs to have length 1 + ncol(X_train)")
  
  if(is.null(tau)) tau <- rep(0.5/sqrt(M), times = ncol(std_X_train))
  else if(length(tau) != ncol(std_X_train)) stop("tau needs to have length 1 + ncol(X_train)")
  
  if(is.null(nu)) nu <- 3
  if(is.null(lambda)) lambda <- stats::qchisq(0.1, df = nu)/nu
  
  if(is.null(gamma0)) gamma0 <- rep(0, times = ncol(W_train))
  if(is.null(xi)) xi <- rep(100, times = ncol(W_train))
  
  # hyperparameters for the Dirichlet prior on splitting probs
  a_u <- 0.5
  b_u <- 1
  
 
  fit <- .gvcbart_ind_fit(Y_train = std_Y_train,
                          subj_id_train = subj_id_train-1, # remember C++ is 0-indexed
                          ni_train = ni_train,
                          tX_train = t(std_X_train),
                          tW_train = t(std_W_train),
                          tZ_cont_train = t(Z_cont_train),
                          tZ_cat_train = t(Z_cat_train),
                          tX_test = t(std_X_test),
                          tW_test = t(std_W_test),
                          tZ_cont_test = t(Z_cont_test),
                          tZ_cat_test = t(Z_cat_test),
                          unif_cuts = unif_cuts,
                          cutpoints_list = cutpoints_list,
                          cat_levels_list = cat_levels_list,
                          edge_mat_list = edge_mat_list,
                          graph_split = graph_split,
                          graph_cut_type = 0, # hardcoded!
                          rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                          sparse = sparse, a_u = a_u, b_u = b_u,
                          mu0 = mu0, tau = tau,
                          gamma0 = gamma0, xi = xi,
                          lambda = lambda, nu = nu,
                          M = M, 
                          nd = nd, burn = burn, thin = thin,
                          save_samples = save_samples, save_trees = save_trees,
                          verbose = verbose, print_every = print_every)
  
  sigma_samples <- y_sd * fit$sigma
  gamma <- sweep(fit$gamma * y_sd, 2, w_sd, '/')
  gamma_mean <- fit$gamma_mean * y_sd / w_sd
  
  yhat_train_mean <- y_mean + y_sd * fit$fit_train_mean
  beta_train_mean <- rescale_beta_mean(fit$beta_train_mean, y_mean, y_sd, x_mean, x_sd)
  
  # Adjust intercept
  for (j in 1:ncol(W_train)) {
    beta_train_mean[,1] <- beta_train_mean[,1] - w_mean[j] * gamma_mean[j]
  }
  
  if(save_samples){
    yhat_train <- y_mean + y_sd * fit$fit_train
    beta_train <- rescale_beta(fit$beta_train, y_mean, y_sd, x_mean, x_sd)
    
    # Adjust intercept
    for (j in 1:ncol(W_train)) {
      beta_train[,,1] <- sweep(beta_train[,,1], 1, w_mean[j] * gamma[,j], '-')
    }
  }
  
  if(!is.null(fit$fit_test_mean)){
    yhat_test_mean <- y_mean + y_sd * fit$fit_test_mean
    beta_test_mean <- rescale_beta_mean(fit$beta_test_mean, y_mean, y_sd, x_mean, x_sd)
    for (j in 1:ncol(W_train)) {
      beta_test_mean[,1] <- beta_test_mean[,1] - w_mean[j] * gamma_mean[j];
    }
    if(save_samples){
      yhat_test <- y_mean + y_sd * fit$fit_test
      beta_test <- rescale_beta(fit$beta_test, y_mean, y_sd, x_mean, x_sd)
      for (j in 1:ncol(W_train)) {
        beta_test[,,1] <- sweep(beta_test[,,1], 1, w_mean[j] * gamma[,j], '-')
      }
    }
  }
  
  results <- list()
  results[["y_mean"]] <- y_mean
  results[["y_sd"]] <- y_sd
  results[["x_mean"]] <- x_mean
  results[["x_sd"]] <- x_sd
  results[["w_mean"]] <- w_mean
  results[["w_sd"]] <- w_sd
  results[["yhat.train.mean"]] <- yhat_train_mean
  results[["betahat.train.mean"]] <- beta_train_mean
  results[["gammahat.mean"]] <- gamma_mean
  results[["gammahat"]] <- gamma
  if(save_samples){
    results[["yhat.train"]] <- yhat_train
    results[["betahat.train"]] <- beta_train
  }
  if(!is.null(fit$fit_test_mean)){
    results[["yhat.test.mean"]] <- yhat_test_mean
    results[["betahat.test.mean"]] <- beta_test_mean
    if(save_samples){
      results[["yhat.test"]] <- yhat_test
      results[["betahat.test"]] <- beta_test
    }
  }
  results[["sigma"]] <- y_sd * fit$sigma
  results[["varcounts"]] <- fit$var_count
  if(sparse) results[["theta"]] <- fit$theta
  if(save_trees) results[["trees"]] <- fit$trees
  return(results)
}
