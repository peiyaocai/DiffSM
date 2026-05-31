packages = c("mvtnorm", "Matrix", "MASS", "clime", "maotai","genscore")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
message("All packages loaded successfully.")
{
  get_element_block_small <- function(data, trace = TRUE, weights = NULL, boot = FALSE) {
    if (is.null(dim(data))) {
      data <- matrix(data, nrow = 1)
    }
    
    n <- nrow(data)
    p <- ncol(data)
    
    if (is.null(weights)) weights <- rep(1, n)
    if (length(weights) != n) stop("weights must have length n.")
    
    Gamma_list <- vector("list", p)
    G_list     <- vector("list", p)
    deri_list  <- vector("list", p)
    
    # Precompute once
    X  <- data
    X2 <- X^2
    
    for (j in seq_len(p)) {
      if (trace) cat("Building block", j, "of", p, "\n")
      
      xj  <- X[, j]
      xj2 <- xj^2
      
      # Guard against exact zero if your formula truly needs division
      if (any(abs(xj) < 1e-12)) {
        stop(sprintf("Column %d contains values too close to zero; derivatives divide by x_j.", j))
      }
      
      # For each sample i:
      # linmat_j      = x_i * x_ij
      # linmat_deri_j = 2 * (linmat_j / x_ij) = 2 * x_i
      #
      # quamat_j      = x_i^2 * x_ij^2
      # quamat_deri_j = 4 * quamat_j / x_ij = 4 * x_i^2 * x_ij
      #
      # So we can vectorize all rows at once.
      
      lin_part <- 2 * X                                 # n x p
      qua_part <- 4 * X2 * xj                           # n x p, broadcasts xj across columns
      qua_part[, j] <- 1
      
      deri_raw <- cbind(lin_part, qua_part)             # n x (2p)
      
      if (boot) {
        # If you really need bootstrap weights, keep them real-valued.
        # Using sqrt(as.complex(...)) is usually a red flag.
        if (any(weights < 0)) {
          stop("Negative bootstrap weights encountered. Current implementation expects nonnegative weights.")
        }
        deri_weighted <- deri_raw * sqrt(weights)
      } else {
        deri_weighted <- deri_raw
      }
      
      # Gamma_j = (1/n) sum_i d_i d_i^T
      # crossprod gives t(A) %*% A, so transpose to match your old layout
      # old code used deri_raw %*% t(deri_raw), where deri_raw was (2p x n)
      # here deri_raw is (n x 2p), so use crossprod(deri_weighted)
      Gamma_j <- crossprod(deri_weighted) / n
      
      # Build G_j
      # Original:
      # lin_deri_ii[j] = 2
      # qua_deri_ii    = 4 * quamat_j / x_ij^2 = 4 * x_i^2
      # qua_deri_ii[j] = 0
      #
      # So average over i:
      G_j <- numeric(2 * p)
      G_j[j] <- 2 * mean(weights)
      
      qua_mean <- 4 * colMeans(X2 * weights)
      qua_mean[j] <- 0
      G_j[(p + 1):(2 * p)] <- qua_mean
      
      Gamma_list[[j]] <- Gamma_j
      G_list[[j]] <- G_j
      deri_list[[j]] <- deri_weighted
    }
    
    list(
      Gamma = Gamma_list,
      G = G_list,
      deri = deri_list,
      p = p,
      n = n
    )
  }
  ADMM_gen_block_small <- function(rho, lambda, alpha, threshold, maxit, elem_x, elem_y, trace = TRUE) {
    
    p <- elem_x$p
    if (is.null(p)) stop("elem_x$p is missing.")
    if (elem_y$p != p) stop("elem_x and elem_y must have the same p.")
    
    # Each block has length 2p
    s_block <- 2 * p
    
    # Diagonal indices to leave unpenalized in theta_x:
    # first p entries = linear part, next p entries = quadratic part
    tmp_mat = diag(1,p)
    diag_set <- list()
    for(j in 1:p){
      diag_set[[j]] <- which(c(tmp_mat[,j],tmp_mat[,j])!=0)
    }
    
    # Initialize lists
    theta_x <- vector("list", p)
    delta   <- vector("list", p)
    Z1      <- vector("list", p)
    Z2      <- vector("list", p)
    U1      <- vector("list", p)
    U2      <- vector("list", p)
    
    Inv1 <- vector("list", p)
    Inv2 <- vector("list", p)
    
    for (j in seq_len(p)) {
      gamma_x_j <- elem_x$Gamma[[j]]
      gamma_y_j <- elem_y$Gamma[[j]]
      
      if (!all(dim(gamma_x_j) == c(s_block, s_block))) stop("Wrong block size in elem_x$Gamma.")
      if (!all(dim(gamma_y_j) == c(s_block, s_block))) stop("Wrong block size in elem_y$Gamma.")
      
      Inv1[[j]] <- solve(gamma_x_j + gamma_y_j + diag(rho, s_block))
      Inv2[[j]] <- solve(gamma_y_j + diag(rho, s_block))
      
      theta_x[[j]] <- numeric(s_block)
      delta[[j]]   <- numeric(s_block)
      Z1[[j]]      <- numeric(s_block)
      Z2[[j]]      <- numeric(s_block)
      U1[[j]]      <- numeric(s_block)
      U2[[j]]      <- numeric(s_block)
    }
    
    soft_thresh <- function(x, t) sign(x) * pmax(0, abs(x) - t)
    
    error <- Inf
    iter <- 0
    
    while (error > threshold && iter < maxit) {
      old_theta <- unlist(theta_x, use.names = FALSE)
      old_delta <- unlist(delta,   use.names = FALSE)
      
      for (j in seq_len(p)) {
        gamma_x_j <- elem_x$Gamma[[j]]
        gamma_y_j <- elem_y$Gamma[[j]]
        g_x_j     <- elem_x$G[[j]]
        g_y_j     <- elem_y$G[[j]]
        
        # theta_x update
        theta_x[[j]] <- Inv1[[j]] %*% (
          gamma_y_j %*% delta[[j]] - g_x_j - g_y_j + rho * (Z1[[j]] - U1[[j]])
        )
        theta_x[[j]] <- as.numeric(theta_x[[j]])
        
        # delta update
        delta[[j]] <- Inv2[[j]] %*% (
          gamma_y_j %*% theta_x[[j]] + g_y_j + rho * (Z2[[j]] - U2[[j]])
        )
        delta[[j]] <- as.numeric(delta[[j]])
        
        # Z updates
        A1 <- theta_x[[j]] + U1[[j]]
        A2 <- delta[[j]] + U2[[j]]
        
        Z1[[j]] <- soft_thresh(A1, lambda / rho)
        Z1[[j]][diag_set[[j]]] <- A1[diag_set[[j]]]   # keep diagonal entries unpenalized
        
        Z2[[j]] <- soft_thresh(A2, lambda * alpha / rho)
        
        # dual updates
        U1[[j]] <- U1[[j]] + theta_x[[j]] - Z1[[j]]
        U2[[j]] <- U2[[j]] + delta[[j]] - Z2[[j]]
      }
      
      new_theta <- unlist(theta_x, use.names = FALSE)
      new_delta <- unlist(delta,   use.names = FALSE)
      
      diff_norm <- sqrt(sum((new_theta - old_theta)^2) + sum((new_delta - old_delta)^2))
      base_norm <- sqrt(sum(old_theta^2) + sum(old_delta^2))
      if (base_norm < 1e-12) base_norm <- 1
      
      error <- diff_norm / base_norm
      iter <- iter + 1
      
      if (trace) cat("iter =", iter, " error =", signif(error, 6), "\n")
    }
    
    list(
      theta_x = Z1,
      delta   = Z2,
      theta_x_vec = unlist(Z1, use.names = FALSE),
      delta_vec   = unlist(Z2, use.names = FALSE),
      iter = iter,
      error = error
    )
  }
  clime_block_small <- function(mat_list, lambda, p, trace = TRUE,dist_type = 'nc',perb = T) {
    if (!is.list(mat_list)) {
      stop("mat_list must be a list of length p.")
    }
    if (length(mat_list) != p) {
      stop("length(mat_list) must equal p.")
    }
    if(dist_type == 'nc'){
      s_block <- 2 * p
    }else if (dist_type == 1){
      s_block <- p
    }else if (dist_type == 2){
      s_block <- p
    }
    
    ret_list <- vector("list", p)
    
    for (j in seq_len(p)) {
      if (trace) print(j)
      
      submat <- mat_list[[j]]
      
      if (!is.matrix(submat)) {
        stop(sprintf("mat_list[[%d]] is not a matrix.", j))
      }
      if (!all(dim(submat) == c(s_block, s_block))) {
        stop(sprintf("mat_list[[%d]] must be a %d x %d matrix.", j, s_block, s_block))
      }
      
      clime_obj <- tryCatch(
        {
          clime(
            submat,
            lambda = lambda,
            sigma = TRUE,
            standardize = FALSE,
            perturb = perb,
            linsolver = 'primaldual'
          )
        },
        error = function(e) {
          clime(
            submat,
            lambda = lambda,
            sigma = TRUE,
            standardize = FALSE,
            perturb = TRUE,
            linsolver = 'primaldual'
          )
        }
      )
      
      ret_list[[j]] <- clime_obj$Omegalist[[1]]
    }
    
    return(ret_list)
  }
  debiased_delta_block <- function(eles_x, eles_y, sm_est, kx, ky, trace = TRUE) {
    
    p <- eles_x$p
    if (is.null(p)) stop("eles_x$p is missing.")
    if (is.null(eles_y$p)) stop("eles_y$p is missing.")
    if (eles_y$p != p) stop("eles_x and eles_y must have the same p.")
    
    if (!is.list(eles_x$Gamma) || !is.list(eles_y$Gamma)) {
      stop("eles_x$Gamma and eles_y$Gamma must be lists of block matrices.")
    }
    if (!is.list(eles_x$G) || !is.list(eles_y$G)) {
      stop("eles_x$G and eles_y$G must be lists of block vectors.")
    }
    if (!is.list(kx) || !is.list(ky)) {
      stop("kx and ky must be lists of CLIME block matrices.")
    }
    
    if (length(eles_x$Gamma) != p || length(eles_y$Gamma) != p) {
      stop("Gamma lists must have length p.")
    }
    if (length(eles_x$G) != p || length(eles_y$G) != p) {
      stop("G lists must have length p.")
    }
    if (length(kx) != p || length(ky) != p) {
      stop("kx and ky must have length p.")
    }
    
    s_block <- dim(eles_x$Gamma[[1]])[1]
    
    # allow either list output from ADMM_gen_block_small or old long vectors
    if (is.list(sm_est$theta_x) && is.list(sm_est$delta)) {
      theta_x_list <- sm_est$theta_x
      delta_list   <- sm_est$delta
    } else if (!is.null(sm_est$theta_x_vec) && !is.null(sm_est$delta_vec)) {
      theta_x_list <- lapply(seq_len(p), function(j) {
        idx <- ((j - 1) * s_block + 1):(j * s_block)
        sm_est$theta_x_vec[idx]
      })
      delta_list <- lapply(seq_len(p), function(j) {
        idx <- ((j - 1) * s_block + 1):(j * s_block)
        sm_est$delta_vec[idx]
      })
    } else {
      stop("sm_est must contain either list-form theta_x/delta or theta_x_vec/delta_vec.")
    }
    
    if (length(theta_x_list) != p || length(delta_list) != p) {
      stop("theta_x and delta must each have length p.")
    }
    
    theta_d_list <- vector("list", p)
    delta_d_list <- vector("list", p)
    
    for (j in seq_len(p)) {
      if (trace) cat("Debiasing block", j, "of", p, "\n")
      
      Gamma_x_j <- eles_x$Gamma[[j]]
      Gamma_y_j <- eles_y$Gamma[[j]]
      G_x_j     <- eles_x$G[[j]]
      G_y_j     <- eles_y$G[[j]]
      
      theta_x_j <- theta_x_list[[j]]
      delta_j   <- delta_list[[j]]
      
      Kx_j <- kx[[j]]
      Ky_j <- ky[[j]]
      
      if (!all(dim(Gamma_x_j) == c(s_block, s_block))) {
        stop(sprintf("eles_x$Gamma[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (!all(dim(Gamma_y_j) == c(s_block, s_block))) {
        stop(sprintf("eles_y$Gamma[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (!all(dim(Kx_j) == c(s_block, s_block))) {
        stop(sprintf("kx[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (!all(dim(Ky_j) == c(s_block, s_block))) {
        stop(sprintf("ky[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (length(G_x_j) != s_block) {
        stop(sprintf("eles_x$G[[%d]] must have length %d.", j, s_block))
      }
      if (length(G_y_j) != s_block) {
        stop(sprintf("eles_y$G[[%d]] must have length %d.", j, s_block))
      }
      if (length(theta_x_j) != s_block) {
        stop(sprintf("theta_x[[%d]] must have length %d.", j, s_block))
      }
      if (length(delta_j) != s_block) {
        stop(sprintf("delta[[%d]] must have length %d.", j, s_block))
      }
      
      # Upper and lower residual blocks of:
      # Gamma_aug %*% c(theta_x_j, delta_j) + G_aug
      r_x_j <- (Gamma_x_j + Gamma_y_j) %*% theta_x_j - Gamma_y_j %*% delta_j + G_x_j + G_y_j
      r_d_j <- -Gamma_y_j %*% theta_x_j + Gamma_y_j %*% delta_j - G_y_j
      
      r_x_j <- as.numeric(r_x_j)
      r_d_j <- as.numeric(r_d_j)
      
      # Debiased upper and lower blocks without building giant augmented matrices
      theta_d_j <- theta_x_j - as.numeric(Kx_j %*% r_x_j + Kx_j %*% r_d_j)
      delta_d_j <- delta_j   - as.numeric(Kx_j %*% r_x_j + (Kx_j + Ky_j) %*% r_d_j)
      
      theta_d_list[[j]] <- theta_d_j
      delta_d_list[[j]] <- delta_d_j
    }
    
    list(
      theta_d = theta_d_list,
      delta_d = delta_d_list,
      theta_d_vec = unlist(theta_d_list, use.names = FALSE),
      delta_d_vec = unlist(delta_d_list, use.names = FALSE),
      p = p
    )
  }
  get_gradients <- function(eles_x, eles_y, sm_est,kx,ky ,trace = TRUE) {
    
    p <- eles_x$p
    if (is.null(p)) stop("eles_x$p is missing.")
    if (is.null(eles_y$p)) stop("eles_y$p is missing.")
    if (eles_y$p != p) stop("eles_x and eles_y must have the same p.")
    
    if (!is.list(eles_x$Gamma) || !is.list(eles_y$Gamma)) {
      stop("eles_x$Gamma and eles_y$Gamma must be lists of block matrices.")
    }
    if (!is.list(eles_x$G) || !is.list(eles_y$G)) {
      stop("eles_x$G and eles_y$G must be lists of block vectors.")
    }
    if (length(eles_x$Gamma) != p || length(eles_y$Gamma) != p) {
      stop("Gamma lists must have length p.")
    }
    if (length(eles_x$G) != p || length(eles_y$G) != p) {
      stop("G lists must have length p.")
    }
    
    s_block <- dim(eles_x$Gamma[[1]])[1]
    
    # allow either list output from ADMM_gen_block_small or old long vectors
    if (is.list(sm_est$theta_x) && is.list(sm_est$delta)) {
      theta_x_list <- sm_est$theta_x
      delta_list   <- sm_est$delta
    } else if (!is.null(sm_est$theta_x_vec) && !is.null(sm_est$delta_vec)) {
      theta_x_list <- lapply(seq_len(p), function(j) {
        idx <- ((j - 1) * s_block + 1):(j * s_block)
        sm_est$theta_x_vec[idx]
      })
      delta_list <- lapply(seq_len(p), function(j) {
        idx <- ((j - 1) * s_block + 1):(j * s_block)
        sm_est$delta_vec[idx]
      })
    } else {
      stop("sm_est must contain either list-form theta_x/delta or theta_x_vec/delta_vec.")
    }
    
    if (length(theta_x_list) != p || length(delta_list) != p) {
      stop("theta_x and delta must each have length p.")
    }
    
    gradient_thetax_list <- vector("list", p)
    gradient_delta_list <- vector("list", p)
    
    for (j in seq_len(p)) {
      if (trace) cat("Computing block", j, "of", p, "\n")
      
      Gamma_x_j <- eles_x$Gamma[[j]]
      Gamma_y_j <- eles_y$Gamma[[j]]
      G_x_j     <- eles_x$G[[j]]
      G_y_j     <- eles_y$G[[j]]
      Kx_j <- kx[[j]]
      Ky_j <- ky[[j]]
      theta_x_j <- theta_x_list[[j]]
      delta_j   <- delta_list[[j]]
      
      if (!all(dim(Gamma_x_j) == c(s_block, s_block))) {
        stop(sprintf("eles_x$Gamma[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (!all(dim(Gamma_y_j) == c(s_block, s_block))) {
        stop(sprintf("eles_y$Gamma[[%d]] must be %d x %d.", j, s_block, s_block))
      }
      if (length(G_x_j) != s_block) {
        stop(sprintf("eles_x$G[[%d]] must have length %d.", j, s_block))
      }
      if (length(G_y_j) != s_block) {
        stop(sprintf("eles_y$G[[%d]] must have length %d.", j, s_block))
      }
      if (length(theta_x_j) != s_block) {
        stop(sprintf("theta_x[[%d]] must have length %d.", j, s_block))
      }
      if (length(delta_j) != s_block) {
        stop(sprintf("delta[[%d]] must have length %d.", j, s_block))
      }
      
      # Upper and lower residual blocks of:
      # Gamma_aug %*% c(theta_x_j, delta_j) + G_aug
      r_x_j <- (Gamma_x_j + Gamma_y_j) %*% theta_x_j - Gamma_y_j %*% delta_j + G_x_j + G_y_j
      r_d_j <- -Gamma_y_j %*% theta_x_j + Gamma_y_j %*% delta_j - G_y_j
      
      r_x_j <- as.numeric(r_x_j)
      r_d_j <- as.numeric(r_d_j)
      
      gradient_thetax_list[[j]] <- as.numeric(Kx_j %*% r_x_j + Kx_j %*% r_d_j)
      gradient_delta_list[[j]] <- as.numeric(Kx_j %*% r_x_j + (Kx_j + Ky_j) %*% r_d_j)
    }
    
    final_list <- list(theta_grad = unlist(gradient_thetax_list), 
                       delta_grad = unlist(gradient_delta_list))
    return(final_list)
  }
  get_leading_matrix <- function(x_data,y_data,kx,ky,sm_est,dist_type = 'nc',type = NULL){
    dim0 <- dim(kx[[1]])[1]
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    full_dim <- length(sm_est$delta_vec)
    leading_mat <- matrix(0,full_dim,n)
    for(i in 1:n){
      print(i)
      if(dist_type == 'nc'){
        ele_xi <- get_element_block_small(x_data[i,],trace = F)
        ele_yi <- get_element_block_small(y_data[i,],trace = F)
      }else {
        ele_xi <- get_element_pos(x_data[i,],trace = F, h_type = type, dist_type = dist_type)
        ele_yi <- get_element_pos(y_data[i,],trace = F, h_type = type, dist_type = dist_type)
      }
      grad_i <- get_gradients(eles_x = ele_xi,eles_y = ele_yi,
                              sm_est = sm_est,kx = kx, ky = ky ,trace = F)
      leading_mat[,i] <- grad_i$delta_grad
    }
    return(leading_mat)
  }
  bootstrap_max <- function(leading_matrix,crit_position,B = 100){
    n <- dim(leading_matrix)[2]
    TB <- rep(0,B)
    for(b in 1:B){
      print(b)
      leading_boot <- leading_matrix %*% rnorm(n,0,1)/sqrt(n)
      TB[b] <- max(abs(leading_boot[crit_position]))
    }
    TB
  }
  get_data_aug = function(n, beta, beta2, burnin, skip){
    
    p = dim(beta)[1]
    X = matrix(0, nrow = n, ncol = p)
    temp = rnorm(p)
    
    # burn in period
    for(t in 1:burnin){
      for(l in 1:p){
        vec_condition = temp[-l]
        vec_beta = beta[,l][-l]
        vec_beta2 = beta2[,l][-l]
        A = 2*t(vec_condition^2)%*%vec_beta + beta2[l,l]
        B = 2*t(vec_condition)%*%vec_beta2 + beta[l,l]
        mu  = -B/(2*A); sig = sqrt(-1/(2*A))
        temp[l] = rnorm(1, mean = mu, sd = sig)
      }
      if(t%%10000==0){
      }
    }
    
    #samples obtained by skipping 500 samples
    k = 1
    samples = 1
    while(samples <= n){
      for(l in 1:p){
        vec_condition = temp[-l]
        vec_beta = beta[,l][-l]
        vec_beta2 = beta2[,l][-l]
        A = 2*t(vec_condition^2)%*%vec_beta + beta2[l,l]
        B = 2*t(vec_condition)%*%vec_beta2 + beta[l,l]
        mu  = -B/(2*A); sig = sqrt(-1/(2*A))
        temp[l] = rnorm(1, mean = mu, sd = sig)
      }
      if(k%%skip == 1){
        X[samples,] <- temp
        samples<-samples+1
      }
      k = k+1
    }
    return(X)
  }
  highd_nc_I <- function(n,p,alpha,seed,m,sig,k = 5,l,l2){
    #basic settings
    rho = 1
    lambda = l*sqrt(log(p)/n)
    lambda2 = l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l = l)
    
    beta = matrix(0,p,p)
    for(j in 2:p){
      beta[j,j-1] = -0.04
      beta[j-1, j] = -0.04
      if(j>=3){
        beta[j, j-2] = -0.04
        beta[j-2, j] = -0.04
      }
    }
    diag(beta) = 8/50
    
    beta2 = matrix(0,p,p)
    for(j in 6:p){
      beta2[j-5,j] = -0.04
      beta2[j,j-5] = -0.04
    }
    diag(beta2) = -1
    
    Delta_true = matrix(0,p,p)
    for(j in 1:m){
      Delta_true[1,3*j] <- sig
      Delta_true[3*j,1] <- sig
    }
    beta_y = beta - Delta_true
    
    betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
    betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
    Delta_v = betax_vec - betay_vec
    
    
    
    set.seed(seed)
    x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
    y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
    
    
    eles_x <- get_element_block_small(x_data[1:n,])
    eles_y <- get_element_block_small(y_data[1:n,])
    eles_x1 <- get_element_block_small(x_data[(n+1):(2*n),])
    eles_y1 <- get_element_block_small(y_data[(n+1):(2*n),])
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_small(rho,lambda, alpha,1e-3,1000,eles_x,eles_y)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    t = 1
    tmp_mat = diag(1,p)
    position_set_raw = NULL
    for(j in 1:p){
      position_set_raw = c(position_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
    }
    position_set = which(position_set_raw!=0)
    position_all = 1:(2*p^2)
    inf_graph_aug = matrix(0,p,p)
    t5 <- Sys.time()
    leading_matrix <- get_leading_matrix(x_data = x_data[1:n,], y_data = y_data[1:n,], kx = kx1, ky = ky1,sm_est = sm_est,
                                         dist_type = 'nc', type = NULL)
    repeat{
      inf_graph_aug_tmp <- inf_graph_aug
      crit_position <- setdiff(position_all, position_set)
      
      boot_seq <- bootstrap_max(leading_matrix = leading_matrix,crit_position = crit_position)
      
      crit_aug <- quantile(boot_seq,0.95)
      
      
      delta_matlin = matrix(0,p,p)
      delta_matqua = matrix(0,p,p)
      for(j in 1:p){
        pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
        vec_sub = delta_d_aug[pos]
        delta_matlin[,j] = vec_sub[1:p]
        delta_matqua[,j] = vec_sub[(p+1):(2*p)]
      }
      
      rej_setlin = which(sqrt(n)*abs(delta_matlin)>=crit_aug, arr.ind = T)
      rej_setqua = which(sqrt(n)*abs(delta_matqua)>=crit_aug, arr.ind = T)
      
      inf_graph_aug[rej_setlin] = 1
      inf_graph_aug[rej_setqua] = 1
      
      inf_graph_aug = inf_graph_aug + t(inf_graph_aug)
      inf_graph_aug = (inf_graph_aug!=0)
      diag(inf_graph_aug) = 1
      
      for(j in 1:p){
        pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
        position_set_raw[pos][1:p] = inf_graph_aug[,j]
        position_set_raw[pos][(p+1):(2*p)] = inf_graph_aug[,j]
      }
      position_set = which(position_set_raw!=0)
      
      max_deg_aug = max(colSums(inf_graph_aug))-1
      indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
      t = t+1
      
      if((max_deg_aug>k)|(indicator==0) ){break}
      
    }
    t6 <- Sys.time()
    time_vec <- c(t2 - t1,t4 - t3, t6-t5)
    typeI_rej = as.numeric(max_deg_aug > k)
    result = list(inf_graph = inf_graph_aug, typeI_rej = typeI_rej, m = m, 
                  crit = crit_aug, attributes = attributes,time = time_vec)
    return(result)
  }
  generate_exponential_graph <- function(
    n,
    Theta,
    burnin = 10000,
    thin = 10,
    x_init = NULL,
    check_inputs = TRUE
  ) {
    
    Theta <- as.matrix(Theta)
    
    if (nrow(Theta) != ncol(Theta)) {
      stop("Theta must be a square matrix.")
    }
    
    p <- nrow(Theta)
    
    if (check_inputs) {
      if (max(abs(Theta - t(Theta))) > 1e-10) {
        stop("Theta must be symmetric.")
      }
      if (any(diag(Theta) <= 0)) {
        stop("All diagonal entries Theta[j,j] must be strictly positive.")
      }
      off_diag <- Theta[row(Theta) != col(Theta)]
      if (any(off_diag < 0)) {
        stop("All off-diagonal entries Theta[j,k] must be nonnegative.")
      }
    }
    
    if (is.null(x_init)) {
      x <- rexp(p, rate = diag(Theta))
    } else {
      x <- as.numeric(x_init)
      if (length(x) != p) {
        stop("x_init must have length p.")
      }
      if (any(x < 0)) {
        stop("x_init must be nonnegative.")
      }
    }
    
    total_iter <- burnin + n * thin
    out <- matrix(0, nrow = n, ncol = p)
    colnames(out) <- paste0("V", seq_len(p))
    
    save_idx <- 0
    
    for (iter in seq_len(total_iter)) {
      for (j in seq_len(p)) {
        rate_j <- Theta[j, j] + sum(Theta[j, -j] * 2*x[-j])
        
        if (rate_j <= 0) {
          stop("A conditional rate became nonpositive. Check Theta.")
        }
        
        x[j] <- rexp(1, rate = rate_j)
      }
      
      if (iter > burnin && ((iter - burnin) %% thin == 0)) {
        save_idx <- save_idx + 1
        out[save_idx, ] <- x
      }
    }
    
    list(
      X = out,
      Theta = Theta,
      burnin = burnin,
      thin = thin,
      last_state = x
    )
  }
  
  get_pos_ab <- function(n,p, Theta,eta,a = '1/2', b = '1/2'){
    domain <- make_domain("R+", p=p)
    ab_idx <- paste0('ab_',a,'_',b)
    data <-  gen(n, setting="exp", abs=FALSE, eta=eta, K=Theta, 
                 domain=domain, finite_infinity=100, burn_in=1000, thinning=1000, verbose=TRUE, 
                 remove_outofbound=TRUE)
    return(data)
  }
  
  
  get_h_functions <- function(type) {
    
    if (type == 1) {
      # h(x) = x^2
      h <- function(x) {
        x^2
      }
      
      h_prime <- function(x) {
        2 * x
      }
      hsq <- function(x){ abs(x)}
      
    } else if (type == 2) {
      # h(x) = log(x + 1)
      h <- function(x) {
        log(x + 1)
      }
      
      h_prime <- function(x) {
        1 / (x + 1)
      }
      hsq <- function(x){sqrt(log(x+1))}
      
    } else if (type == 3) {
      # h(x) = min(x, 3)^1.5
      h <- function(x) {
        pmin(x, 3)^1.5
      }
      
      h_prime <- function(x) {
        ifelse(x < 3,
               1.5 * sqrt(x),
               0)
      }
      hsq <- function(x){sqrt(pmin(x, 3)^1.5)}
    } else {
      stop("type must be 1, 2, or 3")
    }
    
    return(list(h = h, h_prime = h_prime,hsq = hsq))
  }
  get_element_pos <- function(data, trace = TRUE,h_type, dist_type) {
    if (is.null(dim(data))) {
      data <- matrix(data, nrow = 1)
    }
    
    h_all <- get_h_functions(h_type)
    h <- h_all$h
    hp <- h_all$h_prime
    hsq <- h_all$hsq
    
    n <- nrow(data)
    p <- ncol(data)
    
    
    Gamma_list <- vector("list", p)
    G_list     <- vector("list", p)
    deri_list  <- vector("list", p)
    
    # Precompute once
    X  <- data
    
    for (j in seq_len(p)) {
      if (trace) cat("Building block", j, "of", p, "\n")
      
      if( dist_type == 1){
        Xj <- X[,j]
        deri_j <- 2*X
        deri_j[,j] <- 1
        hsq_j <- hsq(Xj)
        deri_j_tilde <- deri_j*hsq_j
        
        hp_j <- hp(Xj)
        term1_j <- colMeans(hp_j*deri_j)
        term2_j <- rep(0,p)
        Gamma_j <- crossprod(deri_j_tilde) / n
        G_j <- -(term1_j + term2_j)
      }else if(dist_type == 2){
        Xj <- X[,j]
        deri_j_raw <- -sqrt(X/Xj)
        #deri_j <- cbind(1/sqrt(Xj), deri_j_raw)
        deri_j <- deri_j_raw
        hsq_j <- hsq(Xj)
        deri_j_tilde <- deri_j*hsq_j
        
        hp_j <- hp(Xj)
        term1_j <- colMeans(hp_j*deri_j)
        deri2_j_raw <- 1/2*sqrt(X/(Xj^3))
        deri2_j_raw[,j]<- 0
        #deri2_j <- cbind(-1/2*sqrt(1/Xj^3), deri2_j_raw)
        deri2_j <- deri2_j_raw
        h_j <- h(Xj)
        term2_j <- colMeans(h_j*deri2_j)
        
        Gamma_j <- crossprod(deri_j_tilde) / n
        G_j <- (term1_j + term2_j)
      }
      
      
      Gamma_list[[j]] <- Gamma_j
      G_list[[j]] <- G_j
    }
    
    list(
      Gamma = Gamma_list,
      G = G_list,
      p = p,
      n = n
    )
  }
  ADMM_gen_block_pos <- function(rho, lambda, alpha, threshold, maxit, elem_x, elem_y,dist_type,trace = TRUE) {
    
    p <- elem_x$p
    if (is.null(p)) stop("elem_x$p is missing.")
    if (elem_y$p != p) stop("elem_x and elem_y must have the same p.")
    
    if(dist_type == 1){
      s_block <- p
      tmp_mat <- diag(1,p)
      diag_set <- list()
      for(j in 1:p){
        diag_set[[j]] <- which(c(tmp_mat[,j])!=0)
      }
    }else if(dist_type == 2){
      s_block <- p
      tmp_mat <- diag(1,p)
      diag_set <- list()
      for(j in 1:p){
        diag_set[[j]] <- which(c( tmp_mat[,j])!=0)
      }
    }
    
    
    # Initialize lists
    theta_x <- vector("list", p)
    delta   <- vector("list", p)
    Z1      <- vector("list", p)
    Z2      <- vector("list", p)
    U1      <- vector("list", p)
    U2      <- vector("list", p)
    
    Inv1 <- vector("list", p)
    Inv2 <- vector("list", p)
    
    for (j in seq_len(p)) {
      gamma_x_j <- elem_x$Gamma[[j]]
      gamma_y_j <- elem_y$Gamma[[j]]
      
      if (!all(dim(gamma_x_j) == c(s_block, s_block))) stop("Wrong block size in elem_x$Gamma.")
      if (!all(dim(gamma_y_j) == c(s_block, s_block))) stop("Wrong block size in elem_y$Gamma.")
      
      Inv1[[j]] <- solve(gamma_x_j + gamma_y_j + diag(rho, s_block))
      Inv2[[j]] <- solve(gamma_y_j + diag(rho, s_block))
      
      theta_x[[j]] <- numeric(s_block)
      delta[[j]]   <- numeric(s_block)
      Z1[[j]]      <- numeric(s_block)
      Z2[[j]]      <- numeric(s_block)
      U1[[j]]      <- numeric(s_block)
      U2[[j]]      <- numeric(s_block)
    }
    
    soft_thresh <- function(x, t) sign(x) * pmax(0, abs(x) - t)
    
    error <- Inf
    iter <- 0
    
    while (error > threshold && iter < maxit) {
      old_theta <- unlist(theta_x, use.names = FALSE)
      old_delta <- unlist(delta,   use.names = FALSE)
      
      for (j in seq_len(p)) {
        gamma_x_j <- elem_x$Gamma[[j]]
        gamma_y_j <- elem_y$Gamma[[j]]
        g_x_j     <- elem_x$G[[j]]
        g_y_j     <- elem_y$G[[j]]
        
        # theta_x update
        theta_x[[j]] <- Inv1[[j]] %*% (
          gamma_y_j %*% delta[[j]] - g_x_j - g_y_j + rho * (Z1[[j]] - U1[[j]])
        )
        theta_x[[j]] <- as.numeric(theta_x[[j]])
        
        # delta update
        delta[[j]] <- Inv2[[j]] %*% (
          gamma_y_j %*% theta_x[[j]] + g_y_j + rho * (Z2[[j]] - U2[[j]])
        )
        delta[[j]] <- as.numeric(delta[[j]])
        
        # Z updates
        A1 <- theta_x[[j]] + U1[[j]]
        A2 <- delta[[j]] + U2[[j]]
        
        Z1[[j]] <- soft_thresh(A1, lambda / rho)
        Z1[[j]][diag_set[[j]]] <- A1[diag_set[[j]]]   # keep diagonal entries unpenalized
        
        Z2[[j]] <- soft_thresh(A2, lambda * alpha / rho)
        
        # dual updates
        U1[[j]] <- U1[[j]] + theta_x[[j]] - Z1[[j]]
        U2[[j]] <- U2[[j]] + delta[[j]] - Z2[[j]]
      }
      
      new_theta <- unlist(theta_x, use.names = FALSE)
      new_delta <- unlist(delta,   use.names = FALSE)
      
      diff_norm <- sqrt(sum((new_theta - old_theta)^2) + sum((new_delta - old_delta)^2))
      base_norm <- sqrt(sum(old_theta^2) + sum(old_delta^2))
      if (base_norm < 1e-12) base_norm <- 1
      
      error <- diff_norm / base_norm
      iter <- iter + 1
      
      if (trace) cat("iter =", iter, " error =", signif(error, 6), "\n")
    }
    
    list(
      theta_x = Z1,
      delta   = Z2,
      theta_x_vec = unlist(Z1, use.names = FALSE),
      delta_vec   = unlist(Z2, use.names = FALSE),
      iter = iter,
      error = error
    )
  }
  highd_pos_I <- function(n,p,alpha,seed,m,sig,type,dist_type,k = 5,l1=5,l2=0.5){
    #basic settings
    rho <- 1
    lambda1 <- l1*sqrt(log(p)/n)
    lambda2 <- l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
    
    if(dist_type == 1){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 1
        ThetaX[j-1,j] <- 1
        if(j >= 3){
          ThetaX[j,j-2] <- 1
          ThetaX[j-2,j] <- 1
        }
      }
      diag(ThetaX)<- 1
      
      Delta_true = matrix(0,p,p)
      for(j in 1:m){
        Delta_true[1,3*j] <- sig
        Delta_true[3*j,1] <- sig
      }
      ThetaY <- ThetaX  - Delta_true
      
      set.seed(seed)
      x_data <- generate_exponential_graph(2*n,ThetaX)$X
      y_data <- generate_exponential_graph(2*n,ThetaY)$X
    }else if (dist_type == 2){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 0.5
        ThetaX[j-1,j] <- 0.5
        if(j >= 3){
          ThetaX[j,j-2] <- 0.5
          ThetaX[j-2,j] <- 0.5
        }
      }
      diag(ThetaX)<- 6
      etaX<- rep(0,p)
      
      Delta_true <- matrix(0,p,p)
      for(j in 1:m){
        Delta_true[1,3*j] <- sig
        Delta_true[3*j,1] <- sig
      }
      ThetaY <- ThetaX  - Delta_true
      etaY <- etaX
      set.seed(seed)
      x_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaX, eta = etaX)
      y_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaY, eta = etaY)
    }
    
    eles_x <- get_element_pos(x_data[1:n,],h_type = type, dist_type = dist_type)
    eles_y <- get_element_pos(y_data[1:n,],h_type = type, dist_type = dist_type)
    eles_x1 <- get_element_pos(x_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    eles_y1 <- get_element_pos(y_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_pos(rho = rho, lambda = lambda1, alpha = alpha, threshold = 1e-3, maxit =1000,
                                 elem_x = eles_x, elem_y = eles_y, dist_type = dist_type)
    dmat<- matrix(0,p,p)
    for(j in 1:p){
      dmat[,j] <- sm_est$delta[[j]]
    }
    #View(dmat)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p,dist_type = dist_type)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p,dist_type = dist_type)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    t <- 1
    if(dist_type == 1){
      tmp_mat <- diag(1,p)
      position_set_raw <- NULL
      for(j in 1:p){
        position_set_raw <- c(position_set_raw, c(tmp_mat[j,]))
      }
      position_set <- which(position_set_raw!=0)
      position_all <- 1:(p^2)
    }else if (dist_type == 2){
      tmp_mat <- diag(1,p)
      position_set_raw<- NULL
      for(j in 1:p){
        position_set_raw <- c(position_set_raw, c(tmp_mat[j,]))
      }
      position_set <- which(position_set_raw!=0)
      position_all <- 1:((p)*p)
    }
    inf_graph_aug <- matrix(0,p,p)
    t5 <- Sys.time()
    leading_matrix <- get_leading_matrix(x_data = x_data[1:n,], y_data = y_data[1:n,], kx = kx1, ky = ky1,sm_est = sm_est,
                                         dist_type = dist_type, type = type)
    repeat{
      inf_graph_aug_tmp <- inf_graph_aug
      crit_position <- setdiff(position_all, position_set)
      
      boot_seq <- bootstrap_max(leading_matrix = leading_matrix,crit_position = crit_position)
      
      crit_aug <- quantile(boot_seq,0.95)
      
      if(dist_type == 1){
        delta_mat <- matrix(0,p,p)
        for(j in 1:p){
          pos <- (p*(j-1) + 1) :(p*(j-1) + p)
          delta_mat[,j] <- delta_d_aug[pos]
        }
      }else if(dist_type == 2){
        delta_mat <- matrix(0,p,p)
        for(j in 1:p){
          pos <- ((p)*(j-1) + 1) :((p)*(j-1) + (p))
          delta_mat[,j] <- delta_d_aug[pos]
        }
      }
      
      rej_set <- which(sqrt(n)*abs(delta_mat)>=crit_aug, arr.ind = T)
      
      
      inf_graph_aug[rej_set] <- 1
      
      
      inf_graph_aug <- inf_graph_aug + t(inf_graph_aug)
      inf_graph_aug <- (inf_graph_aug!=0)
      diag(inf_graph_aug) <- 1
      
      if(dist_type == 1){
        for(j in 1:p){
          pos <- (p*(j-1) + 1) :(p*(j-1) + p)
          position_set_raw[pos] <- inf_graph_aug[,j]
        }
      }else if (dist_type == 2){
        for(j in 1:p){
          pos <- ((p)*(j-1) + 1) :((p)*(j-1) + (p))
          position_set_raw[pos] <- c(inf_graph_aug[,j])
        }
      }
      
      position_set <- which(position_set_raw!=0)
      
      max_deg_aug <- max(colSums(inf_graph_aug))-1
      indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
      t = t+1
      
      if((max_deg_aug>k)|(indicator==0) ){break}
      
    }
    t6 <- Sys.time()
    time_vec <- c(t2 - t1,t4 - t3, t6-t5)
    typeI_rej = as.numeric(max_deg_aug > k)
    result = list(inf_graph = inf_graph_aug, typeI_rej = typeI_rej, m = m, 
                  crit = crit_aug, attributes = attributes,time = time_vec)
    return(result)
  }
  
  tune_lam1 <- function(x_data, y_data, lambda_seq,type, dist_type, alpha, rho){
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    
    cv_pos <- sample(1:n, replace = F)
    cv_pos1 <- cv_pos
    v <- 5
    m_tau = list()
    eles_x_list = list()
    eles_y_list = list()
    len = floor(n/v)
    
    if (dist_type == 'nc'){
      tmp_mat = diag(1,p)
      diag_set_raw = NULL
      for(j in 1:p){
        diag_set_raw = c(diag_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
      }
      diag_set = which(diag_set_raw!=0)
      s = 2*p^2
    } else if (dist_type %in% c(1,2)){
      tmp_mat = diag(1,p)
      diag_set_raw = NULL
      for(j in 1:p){
        diag_set_raw = c(diag_set_raw, c(tmp_mat[j,]))
      }
      diag_set = which(diag_set_raw!=0)
      s = 1*p^2
    }
    for(k in 1:v){
      cv_seq = cv_pos[(1 + (k - 1)*len) : (k*len)]
      cv_seq1 = cv_pos1[(1 + (k - 1)*len) : (k*len)]
      
      if(dist_type == 'nc'){
        eles_val_x <- get_element_block_small(x_data[-cv_seq,], boot = F)
        eles_val_y <- get_element_block_small(y_data[-cv_seq,], boot = F)
      }else if (dist_type %in%c(1,2)){
        eles_val_x <- get_element_pos(x_data[-cv_seq,], T, h_type = type ,dist_type = dist_type)
        eles_val_y <- get_element_pos(y_data[-cv_seq,], T, h_type = type ,dist_type = dist_type)
      }
      eles_x_list[[k]] = eles_val_x
      eles_y_list[[k]] = eles_val_y
    }
    
    for(i in 1:length(lambda_seq)){
      lambda <- lambda_seq[i]
      est_i <- matrix(0, nrow = s, ncol = v)
      for(k in 1:v){
        
        eles_val_x = eles_x_list[[k]]
        eles_val_y = eles_y_list[[k]]
        
        if(dist_type == 'nc'){
          sm_est <- ADMM_gen_block_small(rho = rho, lambda = lambda, alpha = alpha,
                                         threshold =1e-3,maxit = 1000,elem_x = eles_val_x,elem_y = eles_val_y)
        }else if(dist_type %in%c(1,2)){
          sm_est <- ADMM_gen_block_pos(rho = rho, lambda = lambda, alpha = alpha, threshold = 1e-3, maxit =1000,
                                       elem_x = eles_val_x, elem_y = eles_val_y, dist_type = dist_type)
        }
        
        est_i[,k] <- sm_est$delta_vec
      }
      
      m_tau[[i]] = est_i
    }
    
    es <- NULL
    for(i in 1:length(lambda_seq)){
      m_mat <- m_tau[[i]]
      m_lambda <- rowMeans(m_mat)
      v_lambda <- mean(diag(t(m_mat - m_lambda) %*% (m_mat - m_lambda)))
      es = c(es, v_lambda/(m_lambda%*%m_lambda))
    }
    return(cbind(lambda_seq, es))
  }
  
  highd_nc_II <- function(n,p,alpha,seed,m,sig,k = 5,l,l2){
    #basic settings
    rho = 1
    lambda = l*sqrt(log(p)/n)
    lambda2 = l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l = l)
    
    beta = matrix(0,p,p)
    for(j in 2:p){
      beta[j,j-1] = -0.04
      beta[j-1, j] = -0.04
      if(j>=3){
        beta[j, j-2] = -0.04
        beta[j-2, j] = -0.04
      }
    }
    diag(beta) = 8/50
    
    beta2 = matrix(0,p,p)
    for(j in 6:p){
      beta2[j-5,j] = -0.04
      beta2[j,j-5] = -0.04
    }
    diag(beta2) = -1
    
    Delta_true = matrix(0,p,p)
    sig_coef <- 1 + 0.5*(sig-1)
    for(s in 1:sig){
      for(j in 1:m){
        Delta_true[s, s + 2 + 3*(j-1)] <- sig_coef
        Delta_true[s + 2 + 3*(j-1), s] <- sig_coef
      }
    }
    
    
    beta_y = beta - Delta_true
    
    betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
    betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
    Delta_v = betax_vec - betay_vec
    
    
    
    set.seed(seed)
    x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
    y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
    
    
    eles_x <- get_element_block_small(x_data[1:n,])
    eles_y <- get_element_block_small(y_data[1:n,])
    eles_x1 <- get_element_block_small(x_data[(n+1):(2*n),])
    eles_y1 <- get_element_block_small(y_data[(n+1):(2*n),])
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_small(rho,lambda, alpha,1e-3,1000,eles_x,eles_y)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    t = 1
    tmp_mat = diag(1,p)
    position_set_raw = NULL
    for(j in 1:p){
      position_set_raw = c(position_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
    }
    position_set = which(position_set_raw!=0)
    position_all = 1:(2*p^2)
    inf_graph_aug = matrix(0,p,p)
    t5 <- Sys.time()
    leading_matrix <- get_leading_matrix(x_data = x_data[1:n,], y_data = y_data[1:n,], kx = kx1, ky = ky1,sm_est = sm_est,
                                         dist_type = 'nc', type = NULL)
    repeat{
      inf_graph_aug_tmp <- inf_graph_aug
      crit_position <- setdiff(position_all, position_set)
      
      boot_seq <- bootstrap_max(leading_matrix = leading_matrix,crit_position = crit_position)
      
      crit_aug <- quantile(boot_seq,0.95)
      
      
      delta_matlin = matrix(0,p,p)
      delta_matqua = matrix(0,p,p)
      for(j in 1:p){
        pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
        vec_sub = delta_d_aug[pos]
        delta_matlin[,j] = vec_sub[1:p]
        delta_matqua[,j] = vec_sub[(p+1):(2*p)]
      }
      
      rej_setlin = which(sqrt(n)*abs(delta_matlin)>=crit_aug, arr.ind = T)
      rej_setqua = which(sqrt(n)*abs(delta_matqua)>=crit_aug, arr.ind = T)
      
      inf_graph_aug[rej_setlin] = 1
      inf_graph_aug[rej_setqua] = 1
      
      inf_graph_aug = inf_graph_aug + t(inf_graph_aug)
      inf_graph_aug = (inf_graph_aug!=0)
      diag(inf_graph_aug) = 1
      
      for(j in 1:p){
        pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
        position_set_raw[pos][1:p] = inf_graph_aug[,j]
        position_set_raw[pos][(p+1):(2*p)] = inf_graph_aug[,j]
      }
      position_set = which(position_set_raw!=0)
      
      max_deg_aug = max(colSums(inf_graph_aug))-1
      indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
      t = t+1
      
      if((max_deg_aug>k)|(indicator==0) ){break}
      
    }
    t6 <- Sys.time()
    time_vec <- c(t2 - t1,t4 - t3, t6-t5)
    typeII = as.numeric(max_deg_aug <= k)
    result = list(inf_graph = inf_graph_aug, typeII = typeII, m = m, 
                  crit = crit_aug, attributes = attributes,time = time_vec)
    return(result)
  }
  highd_pos_II <- function(n,p,alpha,seed,m,sig,type,dist_type,k = 5,l1=5,l2=0.5){
    #basic settings
    rho <- 1
    lambda1 <- l1*sqrt(log(p)/n)
    lambda2 <- l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
    
    if(dist_type == 1){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 1
        ThetaX[j-1,j] <- 1
        if(j >= 3){
          ThetaX[j,j-2] <- 1
          ThetaX[j-2,j] <- 1
        }
      }
      diag(ThetaX)<- 1
      
      Delta_true = matrix(0,p,p)
      sig_coef <- -(abs(sig)+1)
      for(s in 1:sig){
        for(j in 1:m){
          Delta_true[s, s + 2 + 3*(j-1)] <- sig_coef
          Delta_true[s + 2 + 3*(j-1), s] <- sig_coef
        }
      }
      ThetaY <- ThetaX  - Delta_true
      
      set.seed(seed)
      x_data <- generate_exponential_graph(2*n,ThetaX)$X
      y_data <- generate_exponential_graph(2*n,ThetaY)$X
    }else if (dist_type == 2){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 0.5
        ThetaX[j-1,j] <- 0.5
        if(j >= 3){
          ThetaX[j,j-2] <- 0.5
          ThetaX[j-2,j] <- 0.5
        }
      }
      diag(ThetaX)<- 6
      etaX<- rep(0,p)
      
      Delta_true <- matrix(0,p,p)
      for(j in 1:m){
        Delta_true[1,3*j] <- sig
        Delta_true[3*j,1] <- sig
      }
      ThetaY <- ThetaX  - Delta_true
      etaY <- etaX
      set.seed(seed)
      x_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaX, eta = etaX)
      y_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaY, eta = etaY)
    }
    
    eles_x <- get_element_pos(x_data[1:n,],h_type = type, dist_type = dist_type)
    eles_y <- get_element_pos(y_data[1:n,],h_type = type, dist_type = dist_type)
    eles_x1 <- get_element_pos(x_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    eles_y1 <- get_element_pos(y_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_pos(rho = rho, lambda = lambda1, alpha = alpha, threshold = 1e-3, maxit =1000,
                                 elem_x = eles_x, elem_y = eles_y, dist_type = dist_type)
    dmat<- matrix(0,p,p)
    for(j in 1:p){
      dmat[,j] <- sm_est$delta[[j]]
    }
    #View(dmat)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p,dist_type = dist_type)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p,dist_type = dist_type)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    t <- 1
    if(dist_type == 1){
      tmp_mat <- diag(1,p)
      position_set_raw <- NULL
      for(j in 1:p){
        position_set_raw <- c(position_set_raw, c(tmp_mat[j,]))
      }
      position_set <- which(position_set_raw!=0)
      position_all <- 1:(p^2)
    }else if (dist_type == 2){
      tmp_mat <- diag(1,p)
      position_set_raw<- NULL
      for(j in 1:p){
        position_set_raw <- c(position_set_raw, c(tmp_mat[j,]))
      }
      position_set <- which(position_set_raw!=0)
      position_all <- 1:((p)*p)
    }
    inf_graph_aug <- matrix(0,p,p)
    t5 <- Sys.time()
    leading_matrix <- get_leading_matrix(x_data = x_data[1:n,], y_data = y_data[1:n,], kx = kx1, ky = ky1,sm_est = sm_est,
                                         dist_type = dist_type, type = type)
    repeat{
      inf_graph_aug_tmp <- inf_graph_aug
      crit_position <- setdiff(position_all, position_set)
      
      boot_seq <- bootstrap_max(leading_matrix = leading_matrix,crit_position = crit_position)
      
      crit_aug <- quantile(boot_seq,0.95)
      
      if(dist_type == 1){
        delta_mat <- matrix(0,p,p)
        for(j in 1:p){
          pos <- (p*(j-1) + 1) :(p*(j-1) + p)
          delta_mat[,j] <- delta_d_aug[pos]
        }
      }else if(dist_type == 2){
        delta_mat <- matrix(0,p,p)
        for(j in 1:p){
          pos <- ((p)*(j-1) + 1) :((p)*(j-1) + (p))
          delta_mat[,j] <- delta_d_aug[pos]
        }
      }
      
      rej_set <- which(sqrt(n)*abs(delta_mat)>=crit_aug, arr.ind = T)
      
      
      inf_graph_aug[rej_set] <- 1
      
      
      inf_graph_aug <- inf_graph_aug + t(inf_graph_aug)
      inf_graph_aug <- (inf_graph_aug!=0)
      diag(inf_graph_aug) <- 1
      
      if(dist_type == 1){
        for(j in 1:p){
          pos <- (p*(j-1) + 1) :(p*(j-1) + p)
          position_set_raw[pos] <- inf_graph_aug[,j]
        }
      }else if (dist_type == 2){
        for(j in 1:p){
          pos <- ((p)*(j-1) + 1) :((p)*(j-1) + (p))
          position_set_raw[pos] <- c(inf_graph_aug[,j])
        }
      }
      
      position_set <- which(position_set_raw!=0)
      
      max_deg_aug <- max(colSums(inf_graph_aug))-1
      indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
      t = t+1
      
      if((max_deg_aug>k)|(indicator==0) ){break}
      
    }
    t6 <- Sys.time()
    time_vec <- c(t2 - t1,t4 - t3, t6-t5)
    typeII = as.numeric(max_deg_aug <= k)
    result = list(inf_graph = inf_graph_aug, typeII = typeII, m = m, 
                  crit = crit_aug, attributes = attributes,time = time_vec)
    return(result)
  }
  
  highd_nc_single <- function(n,p,alpha,seed,m,sig,l,l2){
    #basic settings
    rho = 1
    lambda = l*sqrt(log(p)/n)
    lambda2 = l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l = l)
    
    beta = matrix(0,p,p)
    for(j in 2:p){
      beta[j,j-1] = -0.04
      beta[j-1, j] = -0.04
      if(j>=3){
        beta[j, j-2] = -0.04
        beta[j-2, j] = -0.04
      }
    }
    diag(beta) = 8/50
    
    beta2 = matrix(0,p,p)
    for(j in 6:p){
      beta2[j-5,j] = -0.04
      beta2[j,j-5] = -0.04
    }
    diag(beta2) = -1
    
    Delta_true = matrix(0,p,p)
    for(j in 1:m){
      Delta_true[1,3*j] <- sig
      Delta_true[3*j,1] <- sig
    }
    beta_y = beta - Delta_true
    
    betax_vec = c(diag(beta), diag(beta2), beta2[lower.tri(beta2)], beta[lower.tri(beta)])
    betay_vec = c(diag(beta_y), diag(beta2), beta2[lower.tri(beta2)], beta_y[lower.tri(beta_y)])
    Delta_v = betax_vec - betay_vec
    
    
    
    set.seed(seed)
    x_data = get_data_aug(2*n, beta, beta2, 1000, 10)
    y_data = get_data_aug(2*n, beta_y, beta2, 1000, 10)
    
    
    eles_x <- get_element_block_small(x_data[1:n,])
    eles_y <- get_element_block_small(y_data[1:n,])
    eles_x1 <- get_element_block_small(x_data[(n+1):(2*n),])
    eles_y1 <- get_element_block_small(y_data[(n+1):(2*n),])
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_small(rho,lambda, alpha,1e-3,1000,eles_x,eles_y)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p,perb = F)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p,perb = F)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    
    
    var_e1 <- get_var(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 1, dist_type = 'nc',type = NULL)
    var_e2 <- get_var(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 2, dist_type = 'nc',type = NULL)
    result1 <- get_cov(delta_d = delta_d_aug, var_e = var_e1, Delta_true = Delta_true, dist_type = 'nc',n = n,p = p)
    result2 <- get_cov(delta_d = delta_d_aug, var_e = var_e2, Delta_true = Delta_true, dist_type = 'nc',n = n,p = p)
    return(list(r1 = result1, r2 = result2))
  }
  get_var <- function(x_data, y_data, kx,ky,sm_est,var_type,dist_type, type){
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    if(var_type == 1){
      leading_matrix <- get_leading_matrix(x_data = x_data, y_data = y_data, kx = kx, ky = ky,sm_est = sm_est,
                                           dist_type = dist_type, type = type)
      var_e <- rowMeans(leading_matrix^2)
    }else if (var_type == 2){
      s_block <- length(sm_est$delta[[1]])
      var_e <- rep(0,s_block*p)
      for(i in 1:n){
        print(i)
        if(dist_type == 'nc'){
          ele_xi <- get_element_block_small(x_data[i,],trace = F)
          ele_yi <- get_element_block_small(y_data[i,],trace = F)
        }else {
          ele_xi <- get_element_pos(x_data[i,],trace = F, h_type = type, dist_type = dist_type)
          ele_yi <- get_element_pos(y_data[i,],trace = F, h_type = type, dist_type = dist_type)
        }
        for(j in 1:p){
          thetaxj <- sm_est$theta_x[[j]]
          deltaj <- sm_est$delta[[j]]
          kxj <- kx[[j]]
          kyj <- ky[[j]]
          gammaxj <- ele_xi$Gamma[[j]]
          gxj <- ele_xi$G[[j]]
          gammayj <- ele_yi$Gamma[[j]]
          gyj <- ele_yi$G[[j]]
          
          z1j <- kxj%*%(gammaxj%*%thetaxj + gxj)
          z2j <- kyj%*%(gammayj%*%(thetaxj - deltaj) + gyj)
          var_e[(s_block*(j-1) + 1): (s_block*(j-1) + s_block )] <- var_e[(s_block*(j-1) + 1): (s_block*(j-1) + s_block )] +
            z1j^2 + z2j^2
        }
      }
      var_e <- var_e/n
    }
    return(var_e)
  }
  get_cov <- function(delta_d, var_e, Delta_true,dist_type,n,p){
    tmp_delta <- NULL
    if(dist_type == 'nc'){
      for(j in 1:p){
        tmp_delta <- c(tmp_delta,rep(0,p), Delta_true[,j])
      }
    }else{
      for(j in 1:p){
        tmp_delta <- c(tmp_delta, Delta_true[,j])
      }
    }
    Dif_Supp <- (tmp_delta!=0)
    
    clb <- delta_d + sqrt(var_e)*qnorm(0.025)/sqrt(n)
    cub <- delta_d + sqrt(var_e)*qnorm(0.975)/sqrt(n)
    
    supp <- which(tmp_delta!=0)
    supp_c <- which(tmp_delta == 0)
    typeI <- sum((tmp_delta[supp_c]<clb[supp_c]) | (tmp_delta[supp_c] > cub[supp_c]))/length(supp_c)
    typeII <- sum((0>=clb[supp])*(0<=cub[supp]))/length(supp)
    res <- list(supp = supp, supp_c = supp_c, typeI = typeI, typeII = typeII, ci = list(clb, cub))
    return(res)
  }
  highd_pos_single <- function(n,p,alpha,seed,m,sig,type,dist_type,l1=5,l2=0.05){
    #basic settings
    rho <- 1
    lambda1 <- l1*sqrt(log(p)/n)
    lambda2 <- l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
    
    if(dist_type == 1){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 1
        ThetaX[j-1,j] <- 1
        if(j >= 3){
          ThetaX[j,j-2] <- 1
          ThetaX[j-2,j] <- 1
        }
      }
      diag(ThetaX)<- 1
      
      Delta_true = matrix(0,p,p)
      for(j in 1:m){
        Delta_true[1,3*j] <- sig
        Delta_true[3*j,1] <- sig
      }
      ThetaY <- ThetaX  - Delta_true
      
      set.seed(seed)
      x_data <- generate_exponential_graph(2*n,ThetaX)$X
      y_data <- generate_exponential_graph(2*n,ThetaY)$X
    }else if (dist_type == 2){
      ThetaX <- matrix(0,p,p)
      for(j in 2:p){
        ThetaX[j,j-1] <- 0.5
        ThetaX[j-1,j] <- 0.5
        if(j >= 3){
          ThetaX[j,j-2] <- 0.5
          ThetaX[j-2,j] <- 0.5
        }
      }
      diag(ThetaX)<- 6
      etaX<- rep(0,p)
      
      Delta_true <- matrix(0,p,p)
      for(j in 1:m){
        Delta_true[1,3*j] <- sig
        Delta_true[3*j,1] <- sig
      }
      ThetaY <- ThetaX  - Delta_true
      etaY <- etaX
      set.seed(seed)
      x_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaX, eta = etaX)
      y_data <- get_pos_ab(n = 2*n, p =p,Theta = ThetaY, eta = etaY)
    }
    
    eles_x <- get_element_pos(x_data[1:n,],h_type = type, dist_type = dist_type)
    eles_y <- get_element_pos(y_data[1:n,],h_type = type, dist_type = dist_type)
    eles_x1 <- get_element_pos(x_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    eles_y1 <- get_element_pos(y_data[(n+1):(2*n),],h_type = type, dist_type = dist_type)
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_pos(rho = rho, lambda = lambda1, alpha = alpha, threshold = 1e-3, maxit =1000,
                                 elem_x = eles_x, elem_y = eles_y, dist_type = dist_type)
    dmat<- matrix(0,p,p)
    for(j in 1:p){
      dmat[,j] <- sm_est$delta[[j]]
    }
    #View(dmat)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p,dist_type = dist_type,perb = F)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p,dist_type = dist_type,perb = F)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    
    
    var_e1 <- get_var(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 1, dist_type = dist_type,type = type)
    var_e2 <- get_var(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 2, dist_type = dist_type,type = type)
    var_e1 <- matrix(var_e1,nrow = p, ncol = p)
    var_e1 <- pmax(var_e1, t(var_e1))
    var_e1 <- as.vector(var_e1)
    var_e2 <- matrix(var_e2,nrow = p, ncol = p)
    var_e2 <- pmax(var_e2, t(var_e2))
    var_e2 <- as.vector(var_e2)
    result1 <- get_cov(delta_d = delta_d_aug, var_e = var_e1, Delta_true = Delta_true, dist_type = dist_type,n = n,p = p)
    result2 <- get_cov(delta_d = delta_d_aug, var_e = var_e2, Delta_true = Delta_true, dist_type = dist_type,n = n,p = p)
    return(list(r1 = result1, r2 = result2))
  }
  
  get_ab_data <- function(n,p,Theta, a,b){
    eta <- rep(0,p)
    domain <- make_domain("R+", p=p)
    if ((a == '1/2')&(b == '1/2') ){
      data <-  gen(n, setting="exp", abs=FALSE, eta=eta, K=Theta, 
                   domain=domain, finite_infinity=100, burn_in=1000, thinning=100, verbose=TRUE, 
                   remove_outofbound=TRUE)
    }else{
      ab_idx <- paste0('ab_',a,'_',b)
      data <- gen(n, setting=ab_idx, abs=FALSE, eta=eta, K=Theta, 
                  domain=domain, finite_infinity=100, burn_in=1000, thinning=100, verbose=TRUE, 
                  remove_outofbound=TRUE)
    }
  }
  get_elements_ab <- function(data, hfunc, hp1, hp2,a,b){
    a_val <- as.numeric(eval(parse(text = a)))
    b_val <- as.numeric(eval(parse(text = b)))
    ab_idx <- paste0('ab_',a,'_',b)
    Gamma <- list()
    G <- list()
    n <- dim(data)[1]
    p <- dim(data)[2]
    domain <- make_domain("R+", p=p)
    dm <- 1 + (1-1/(1+4*exp(1)*max(6*log(p)/n, sqrt(6*log(p)/n))))
    h_hp <- get_h_hp(mode=hfunc, para=hp1, para2=hp2)
    h_hp_dx <- h_of_dist(h_hp, data,domain)
    if(n == 1){
      h_hp_dx$hdx <- matrix(h_hp_dx$hdx, nrow = 1)
      h_hp_dx$hpdx <- matrix( h_hp_dx$hpdx, nrow = 1)
    }
    elts <- get_elts_ab(h_hp_dx$hdx, h_hp_dx$hpdx, data, a=a_val, b=b_val, setting=ab_idx,
                        centered=TRUE, scale=NULL, diag=dm)
    for(j in 1:p){
      Gamma[[j]] <- elts$Gamma_K[1:p, (p*(j-1)+1 ):(p*(j-1) + p)]
      G[[j]] <- -elts$g_K[(p*(j-1)+1 ):(p*(j-1) + p)]
    }
    return(list(Gamma = Gamma, G = G, n = n, p = p))
  }
  get_leading_matrix_ab <- function(x_data,y_data,kx,ky,sm_est,hfunc, hp1, hp2){
    dim0 <- dim(kx[[1]])[1]
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    full_dim <- length(sm_est$delta_vec)
    leading_mat <- matrix(0,full_dim,n)
    for(i in 1:n){
      print(i)
      ele_xi <- get_elements_ab(matrix(x_data[i, ],nrow = 1), hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
      ele_yi <- get_elements_ab(matrix(y_data[i, ],nrow = 1), hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
      grad_i <- get_gradients(eles_x = ele_xi,eles_y = ele_yi,
                              sm_est = sm_est,kx = kx, ky = ky ,trace = F)
      leading_mat[,i] <- grad_i$delta_grad
    }
    return(leading_mat)
  }
  get_var_ab <- function(x_data, y_data, kx,ky,sm_est,var_type, hfunc, hp1,hp2){
    n <- dim(x_data)[1]
    p <- dim(x_data)[2]
    if(var_type == 1){
      leading_matrix <- get_leading_matrix_ab(x_data = x_data, y_data = y_data, kx =kx, ky = ky,
                                              sm_est = sm_est, hfunc = hfunc, hp1 = hp1, hp2 = hp2)
      var_e <- rowMeans(leading_matrix^2)
    }else if (var_type == 2){
      s_block <- length(sm_est$delta[[1]])
      var_e <- rep(0,s_block*p)
      for(i in 1:n){
        print(i)
        ele_xi <- get_elements_ab(matrix(x_data[i, ],nrow = 1), hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
        ele_yi <- get_elements_ab(matrix(y_data[i, ],nrow = 1), hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
        for(j in 1:p){
          thetaxj <- sm_est$theta_x[[j]]
          deltaj <- sm_est$delta[[j]]
          kxj <- kx[[j]]
          kyj <- ky[[j]]
          gammaxj <- ele_xi$Gamma[[j]]
          gxj <- ele_xi$G[[j]]
          gammayj <- ele_yi$Gamma[[j]]
          gyj <- ele_yi$G[[j]]
          
          z1j <- kxj%*%(gammaxj%*%thetaxj + gxj)
          z2j <- kyj%*%(gammayj%*%(thetaxj - deltaj) + gyj)
          var_e[(s_block*(j-1) + 1): (s_block*(j-1) + s_block )] <- var_e[(s_block*(j-1) + 1): (s_block*(j-1) + s_block )] +
            z1j^2 + z2j^2
        }
      }
      var_e <- var_e/n
    }
    return(var_e)
  }
  highd_ab_I <- function(n,p,alpha,seed,m,sig,hfunc,hp1,hp2,a,b,k = 5,l1,l2){
    #basic settings
    rho <- 1
    lambda1 <- l1*sqrt(log(p)/n)
    lambda2 <- l2*sqrt(log(p)/n)
    attributes = list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
    
    Kx <- cov_cons(mode="sub", p=p, seed=999, spars=0.2, eig= 2, subgraphs=5)
    Delta_true <- matrix(0,p,p)
    for(j in 1:m){
      Delta_true[1,3*j] <- sig
      Delta_true[3*j,1] <- sig
    }
    Ky <- Kx  - Delta_true
    diag(Ky) <- diag(Ky) + abs(min(eigen(Ky)$values)) + 1
    diag(Kx) <- diag(Ky)
    set.seed(seed)
    x_data <- get_ab_data(n = 2*n, p =p,Theta = Kx, a = a, b = b)
    y_data <- get_ab_data(n = 2*n, p =p,Theta = Ky, a = a, b = b)
    
    eles_x <- get_elements_ab(x_data[1:n,], hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
    eles_y <- get_elements_ab(y_data[1:n,], hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
    eles_x1 <- get_elements_ab(x_data[(n+1):(2*n),], hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
    eles_y1 <- get_elements_ab(y_data[(n+1):(2*n),], hfunc = hfunc, hp1 = hp1, hp2 = hp2, a = a, b = b)
    
    t1 <- Sys.time()
    sm_est <- ADMM_gen_block_pos(rho = rho, lambda = lambda1, alpha = alpha, threshold = 1e-3, maxit =1000,
                                 elem_x = eles_x, elem_y = eles_y, dist_type = 2)
    dmat<- matrix(0,p,p)
    for(j in 1:p){
      dmat[,j] <- sm_est$delta[[j]]
    }
    #View(dmat)
    t2 <- Sys.time()
    
    t3 <- Sys.time()
    kx1 <- clime_block_small(eles_x1$Gamma,lambda2,p,dist_type = 2)
    ky1 <- clime_block_small(eles_y1$Gamma,lambda2,p,dist_type = 2)
    t4 <- Sys.time()
    
    
    debiase_obj <- debiased_delta_block(eles_x = eles_x,eles_y = eles_y,
                                        sm_est = sm_est, kx = kx1, ky = ky1)
    delta_d_aug <- debiase_obj$delta_d_vec
    
    var_e1 <- get_var_ab(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 1, hfunc = hfunc, hp1 = hp1, hp2 = hp2)
    var_e2 <- get_var_ab(x_data[1:n,], y_data[1:n,], kx = kx1, ky = ky1, sm_est = sm_est,var_type = 2, hfunc = hfunc, hp1 = hp1, hp2 = hp2)
    result1 <- get_cov(delta_d = delta_d_aug, var_e = var_e1, Delta_true = Delta_true, dist_type = 2,n = n,p = p)
    result2 <- get_cov(delta_d = delta_d_aug, var_e = var_e2, Delta_true = Delta_true, dist_type = 2,n = n,p = p)
    
    t <- 1
    tmp_mat <- diag(1,p)
    position_set_raw<- NULL
    for(j in 1:p){
      position_set_raw <- c(position_set_raw, c(tmp_mat[j,]))
    }
    position_set <- which(position_set_raw!=0)
    position_all <- 1:((p)*p)
    inf_graph_aug <- matrix(0,p,p)
    t5 <- Sys.time()
    leading_matrix <- get_leading_matrix_ab(x_data = x_data[1:n,], y_data = y_data[1:n,], kx =kx1, ky = ky1,
                                            sm_est = sm_est, hfunc = hfunc, hp1 = hp1, hp2 = hp2)
    var_e1 <- rowMeans(leading_matrix^2)
    result1 <- get_cov(delta_d = delta_d_aug, var_e = var_e1, Delta_true = Delta_true, dist_type = 2,n = n,p = p)
    repeat{
      inf_graph_aug_tmp <- inf_graph_aug
      crit_position <- setdiff(position_all, position_set)
      
      boot_seq <- bootstrap_max(leading_matrix = leading_matrix,crit_position = crit_position)
      
      crit_aug <- quantile(boot_seq,0.95)
      delta_mat <- matrix(0,p,p)
      for(j in 1:p){
        pos <- (p*(j-1) + 1) :(p*(j-1) + p)
        delta_mat[,j] <- delta_d_aug[pos]
      }
      
      
      rej_set <- which(sqrt(n)*abs(delta_mat)>=crit_aug, arr.ind = T)
      
      
      inf_graph_aug[rej_set] <- 1
      
      
      inf_graph_aug <- inf_graph_aug + t(inf_graph_aug)
      inf_graph_aug <- (inf_graph_aug!=0)
      diag(inf_graph_aug) <- 1
      for(j in 1:p){
        pos <- (p*(j-1) + 1) :(p*(j-1) + p)
        position_set_raw[pos] <- inf_graph_aug[,j]
      }
      
      position_set <- which(position_set_raw!=0)
      
      max_deg_aug <- max(colSums(inf_graph_aug))-1
      indicator = max(abs(inf_graph_aug - inf_graph_aug_tmp))
      t = t+1
      
      if((max_deg_aug>k)|(indicator==0) ){break}
      
    }
    t6 <- Sys.time()
    time_vec <- c(t2 - t1,t4 - t3, t6-t5)
    typeI_rej = as.numeric(max_deg_aug > k)
    result = list(inf_graph = inf_graph_aug, typeI_rej = typeI_rej, m = m, 
                  crit = crit_aug, attributes = attributes,time = time_vec,inf = list(result1, result2))
    return(result)
  }
  
}
