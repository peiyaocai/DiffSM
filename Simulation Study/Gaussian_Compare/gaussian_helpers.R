#library(genscore)
library(mvtnorm)
library(Matrix)
library(MASS)
library(clime)
library(maotai)
get_elements = function(x_data, y_data){
  nobs = dim(x_data)[1]
  gamma_x = t(x_data)%*%x_data/nobs
  gamma_y = t(y_data)%*%y_data/nobs
  return(list(gamma_x = gamma_x, gamma_y = gamma_y))
}
ADMM_sym1 = function(rho, lambda, alpha, threshold, maxit, elem){
  
  
  #initialize
  gamma_x = elem$gamma_x
  gamma_y = elem$gamma_y
  theta_x = matrix(0,p,p)
  delta = matrix(0,p,p)
  Z1 = matrix(0,p,p)
  Z2 = matrix(0,p,p)
  U1 = matrix(0,p,p)
  U2 = matrix(0,p,p)
  error = 1 + threshold
  i = 0
  while((error > threshold)&(i<=maxit)){
    
    theta_0 = rbind(theta_x, delta)
    
    C1 = diag(2,p) + (gamma_y%*%delta)/2 + (delta%*%gamma_y)/2 + rho*(Z1 - U1)
    A1 = (gamma_x + gamma_y)/2
    B1 = A1 + diag(1,p)
    theta_x = sylvester(A1, B1, C1)
    C2 = diag(-1,p) + (gamma_y%*%theta_x)/2 + (theta_x%*%gamma_y)/2 + rho*(Z2 - U2)
    A2 = gamma_y/2
    B2 = A2 + diag(1,p)
    delta = sylvester(A2, B2, C2)
    
    A1 = theta_x + U1
    A2 = delta + U2
    
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    diag(Z1) = diag(A1)
    Z2 = sign(A2)*pmax(0, abs(A2)-lambda*alpha/rho)
    
    U1 = U1 + theta_x - Z1
    U2 = U2 + delta - Z2
    
    theta_1 = rbind(theta_x, delta)
    dif = theta_1 - theta_0
    error = norm(dif, type = 'F')/norm(theta_0, type = 'F')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(list(theta_x = Z1, delta = Z2))
}
matdif = function(A,B){
  pos_vec = c()
  if(is.null(A)){
    return(B)
  }else
  {for(i in 1:(dim(A)[1])){
    for(j in 1:(dim(B)[1])){
      if (max(abs(A[i,] - B[j,])) ==0){
        pos_vec = c(pos_vec,j)
      }
    }
  }
    pos_vec = unique(pos_vec)
    if(is.null(pos_vec)){return(B)}else{
      return(B[-pos_vec,])}}
}
gaussian_I <- function(n,p,alpha,seed,m,sig,k = 5,l1,l2){
  #basic settings
  rho <- 1
  lambda1 <- l1*sqrt(log(p)/n)
  lambda2 <- l2*sqrt(log(p)/n)
  attributes <- list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
  
  set.seed(9999)
  A <- matrix(0,p,p)
  A[upper.tri(A)] <- rbinom(length(A[upper.tri(A)]),1,0.2)
  A <- A + t(A)
  diag(A) <- 1
  
  Delta_true <- matrix(0,p,p)
  for(j in 1:m){
    Delta_true[1,3*j] <- sig
    Delta_true[j*3,1] <- sig
  }
  B <- A - Delta_true
  e1 <- min(eigen(A)$value)
  e2 <- min(eigen(B)$value)
  
  ThetaX <- A + diag(1,p)*(0.5 + abs(min(e1,e2)))
  ThetaY <- B + diag(1,p)*(0.5 + abs(min(e1,e2)))
  Delta_true <- ThetaX - ThetaY
  Delta_v <- Delta_true[lower.tri(Delta_true)]
  Dif_Supp <- (Delta_true!=0)
  
  set.seed(seed = seed)
  #data generation
  x_data_all = rmvnorm(2*n, mean = rep(0,p), sigma = solve(ThetaX))
  y_data_all = rmvnorm(2*n, mean = rep(0,p), sigma = solve(ThetaY))
  
  x_data <- x_data_all[1:n,]
  y_data <- y_data_all[1:n,]
  if(!dir.exists('Data')){dir.create('Data')}
  xdir <- paste0('Data/','XNull_n',n,'p',p,'sig',sig,'seed',seed, '.csv')
  ydir <- paste0('Data/','YNull_n',n,'p',p,'sig',sig,'seed',seed, '.csv')
  write.csv(t(x_data), file = xdir, row.names = FALSE)
  write.csv(t(y_data), file = ydir, row.names = FALSE)
  x_data1 <- x_data_all[(n+1):(2*n),]
  y_data1 <- y_data_all[(n+1):(2*n),]
  eles <- get_elements(x_data, y_data)
  eles1 <- get_elements(x_data1, y_data1)
  
  sm_est = ADMM_sym1(rho, lambda1, alpha, 1e-3, 1000, eles)
  theta_hat = rbind(sm_est$theta_x, sm_est$delta)
  
  # debiasing
  Gamma_hat <- matrix(0, 2*p, 2*p)
  Gamma_hat1 <- matrix(0, 2*p, 2*p)
  
  Gamma_hat[1:p, 1:p] <- eles$gamma_x + eles$gamma_y
  Gamma_hat[1:p, (p+1):(2*p)] <-  -eles$gamma_y
  Gamma_hat[(p+1):(2*p), 1:p] <- -eles$gamma_y
  Gamma_hat[(p+1):(2*p), (p+1):(2*p)] <- eles$gamma_y
  
  Gamma_hat1[1:p, 1:p] <- eles1$gamma_x + eles1$gamma_y
  Gamma_hat1[1:p, (p+1):(2*p)] <-  -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), 1:p] <- -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), (p+1):(2*p)] <- eles1$gamma_y
  
  Ghat <- rbind(diag(-2,p), diag(1,p))
  
  clime_objx <- clime(eles1$gamma_x, lambda = l2*sqrt(log(p)/n), sigma = T, standardize = F,perturb = T)
  clime_objy <- clime(eles1$gamma_y, lambda = l2*sqrt(log(p)/n), sigma = T, standardize = F,perturb = T)
  kx <- clime_objx$Omegalist[[1]]
  ky <- clime_objy$Omegalist[[1]]
  Mhat <- rbind(cbind(kx,kx), cbind(kx, kx+ky))
  
  theta_d = theta_hat - Mhat%*%(Gamma_hat%*%theta_hat + Ghat)
  
  theta_true = rbind(ThetaX, Delta_true)
  
  #calculate the variance of graph nodes - edge(2,3)
  var_e = matrix(0,p,p)
  gamma_x = eles$gamma_x
  gamma_y = eles$gamma_y
  
  for(j in 1:p){
    for(k in 1:p){
      
      
      Mx = kx[j,]
      My = ky[j,]
      thetax_k = sm_est$theta_x[,k]
      delta_k =  sm_est$delta[,k]
      
      var_jk = Mx%*%gamma_x%*%Mx*thetax_k[k] + kx[j,k]^2 + My%*%gamma_y%*%My*(thetax_k[k] - delta_k[k]) + ky[j,k]^2
      
      
      var_e[j,k] = var_jk
      
    }
  }
  #variance symmetrization
  var_e <- pmax(var_e, t(var_e))
  
  
  delta_d = theta_d[(p+1):(p*2),]
  clb = delta_d + sqrt(var_e)*qnorm(0.025)/sqrt(n)
  cub = delta_d + sqrt(var_e)*qnorm(0.975)/sqrt(n)
  
  supp = which(Dif_Supp!=0, arr.ind = T)
  supp_c = which(Dif_Supp==0, arr.ind = T)
  #supp_c = supp_c[which(supp_c[,1]!=supp_c[,2]),]
  
  typeI = sum((Delta_true[supp_c]<clb[supp_c]) | (Delta_true[supp_c] > cub[supp_c]))/dim(supp_c)[1]
  typeII = sum((0>=clb[supp])*(0<=cub[supp]))/dim(supp)[1]
  
  t <- 1
  edges_set <- NULL
  edge_all <- which(diag(-1,p)>=0, arr.ind = T)#positions of off-diagonal terms
  edges_set_list <- list()
  inf_graph <- matrix(0,p,p)
  leading_list <- list()
  for(i in 1:n){
    xi <- x_data[i,]
    yi <- y_data[i,]
    gamma_xi <- xi%*%t(xi)
    gamma_yi <- yi%*%t(yi)
    
    gamma_i <- rbind(cbind(gamma_xi + gamma_yi, -gamma_yi), cbind(-gamma_yi, gamma_yi))
    leading_i <- -Mhat%*%(gamma_i%*%theta_hat + Ghat)/sqrt(n)
    leading_sub_i <- leading_i[(p+1):(p*2),]
    diag(leading_sub_i) <- 0
    leading_list[[i]] <- leading_sub_i
  }
  
  repeat{
    inf_graph_tmp = inf_graph
    edge_crit = matdif(edges_set, edge_all)
    TB = NULL
    for(c in 1:100){ #loop for Monte Carlo
      boot_weights <- rnorm(n)
      boot_res <- Reduce(`+`, Map(`*`, leading_list, boot_weights))
      print(c)
      TB = c(TB, max(abs(boot_res[edge_crit])))
    }
    crit = quantile(TB, 0.95)
    
    delta_d = theta_d[(p+1):(2*p),]
    delta_d[edges_set] = 0
    rej_set <- which(sqrt(n)*abs(delta_d) >= crit, arr.ind = T)
    edges_set <- rbind(edges_set, rej_set)
    inf_graph <- matrix(0,p,p)
    inf_graph[edges_set] <- 1
    inf_graph <- inf_graph + t(inf_graph)
    inf_graph <- (inf_graph !=0)
    diag(inf_graph) <- 0
    max_deg <- max(colSums(inf_graph))
    indicator <- max(abs(inf_graph_tmp - inf_graph))
    t <- t+1
    if((indicator==0)|(max_deg>k)){break}
  }
  
  typeI_rej <- as.numeric(max_deg > k)
  
  results <- list(inf_graph = inf_graph, typeI_rej = typeI_rej, m = m, crit = crit, attributes = attributes,
                  typeI = typeI, typeII = typeII, ci = list(clb,cub),est = sm_est$delta)
  
  results
}

gaussian_II <- function(n,p,alpha,seed,m,sig,k = 5,l1,l2){
  #basic settings
  rho <- 1
  lambda1 <- l1*sqrt(log(p)/n)
  lambda2 <- l2*sqrt(log(p)/n)
  attributes <- list(n = n, p = p, alpha = alpha, seed = seed, m = m,l1 = l1)
  
  set.seed(9999)
  A <- matrix(0,p,p)
  A[upper.tri(A)] <- rbinom(length(A[upper.tri(A)]),1,0.2)
  A <- A + t(A)
  diag(A) <- 1
  
  Delta_true = matrix(0,p,p)
  sig_coef <- 0.5*abs(sig) + 0.5
  for(s in 1:sig){
    for(j in 1:m){
      Delta_true[s, s + 2 + 3*(j-1)] <- sig_coef
      Delta_true[s + 2 + 3*(j-1), s] <- sig_coef
    }
  }
  B <- A - Delta_true
  e1 <- min(eigen(A)$value)
  e2 <- min(eigen(B)$value)
  
  ThetaX <- A + diag(1,p)*(0.5 + abs(min(e1,e2)))
  ThetaY <- B + diag(1,p)*(0.5 + abs(min(e1,e2)))
  Delta_true <- ThetaX - ThetaY
  Delta_v <- Delta_true[lower.tri(Delta_true)]
  Dif_Supp <- (Delta_true!=0)
  
  set.seed(seed = seed)
  #data generation
  x_data_all = rmvnorm(2*n, mean = rep(0,p), sigma = solve(ThetaX))
  y_data_all = rmvnorm(2*n, mean = rep(0,p), sigma = solve(ThetaY))
  
  x_data <- x_data_all[1:n,]
  y_data <- y_data_all[1:n,]
  if(!dir.exists('Data')){dir.create('Data')}
  xdir <- paste0('Data/','XAlt_n',n,'p',p,'sig',sig,'seed',seed, '.csv')
  ydir <- paste0('Data/','YAlt_n',n,'p',p,'sig',sig,'seed',seed, '.csv')
  write.csv(t(x_data), file = xdir, row.names = FALSE)
  write.csv(t(y_data), file = ydir, row.names = FALSE)
  x_data1 <- x_data_all[(n+1):(2*n),]
  y_data1 <- y_data_all[(n+1):(2*n),]
  eles <- get_elements(x_data, y_data)
  eles1 <- get_elements(x_data1, y_data1)
  
  sm_est = ADMM_sym1(rho, lambda1, alpha, 1e-3, 1000, eles)
  theta_hat = rbind(sm_est$theta_x, sm_est$delta)
  
  # debiasing
  Gamma_hat <- matrix(0, 2*p, 2*p)
  Gamma_hat1 <- matrix(0, 2*p, 2*p)
  
  Gamma_hat[1:p, 1:p] <- eles$gamma_x + eles$gamma_y
  Gamma_hat[1:p, (p+1):(2*p)] <-  -eles$gamma_y
  Gamma_hat[(p+1):(2*p), 1:p] <- -eles$gamma_y
  Gamma_hat[(p+1):(2*p), (p+1):(2*p)] <- eles$gamma_y
  
  Gamma_hat1[1:p, 1:p] <- eles1$gamma_x + eles1$gamma_y
  Gamma_hat1[1:p, (p+1):(2*p)] <-  -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), 1:p] <- -eles1$gamma_y
  Gamma_hat1[(p+1):(2*p), (p+1):(2*p)] <- eles1$gamma_y
  
  Ghat <- rbind(diag(-2,p), diag(1,p))
  
  clime_objx <- clime(eles1$gamma_x, lambda = l2*sqrt(log(p)/n), sigma = T, standardize = F,perturb = T)
  clime_objy <- clime(eles1$gamma_y, lambda = l2*sqrt(log(p)/n), sigma = T, standardize = F,perturb = T)
  kx <- clime_objx$Omegalist[[1]]
  ky <- clime_objy$Omegalist[[1]]
  Mhat <- rbind(cbind(kx,kx), cbind(kx, kx+ky))
  
  
  theta_d = theta_hat - Mhat%*%(Gamma_hat%*%theta_hat + Ghat)
  
  theta_true = rbind(ThetaX, Delta_true)
  
  t <- 1
  edges_set <- NULL
  edge_all <- which(diag(-1,p)>=0, arr.ind = T)#positions of off-diagonal terms
  edges_set_list <- list()
  inf_graph <- matrix(0,p,p)
  leading_list <- list()
  for(i in 1:n){
    xi <- x_data[i,]
    yi <- y_data[i,]
    gamma_xi <- xi%*%t(xi)
    gamma_yi <- yi%*%t(yi)
    
    gamma_i <- rbind(cbind(gamma_xi + gamma_yi, -gamma_yi), cbind(-gamma_yi, gamma_yi))
    leading_i <- -Mhat%*%(gamma_i%*%theta_hat + Ghat)/sqrt(n)
    leading_sub_i <- leading_i[(p+1):(p*2),]
    diag(leading_sub_i) <- 0
    leading_list[[i]] <- leading_sub_i
  }
  
  repeat{
    inf_graph_tmp = inf_graph
    edge_crit = matdif(edges_set, edge_all)
    TB = NULL
    for(c in 1:100){ #loop for Monte Carlo
      boot_weights <- rnorm(n)
      boot_res <- Reduce(`+`, Map(`*`, leading_list, boot_weights))
      print(c)
      TB = c(TB, max(abs(boot_res[edge_crit])))
    }
    crit = quantile(TB, 0.95)
    
    delta_d = theta_d[(p+1):(2*p),]
    delta_d[edges_set] = 0
    rej_set <- which(sqrt(n)*abs(delta_d) >= crit, arr.ind = T)
    edges_set <- rbind(edges_set, rej_set)
    inf_graph <- matrix(0,p,p)
    inf_graph[edges_set] <- 1
    inf_graph <- inf_graph + t(inf_graph)
    inf_graph <- (inf_graph !=0)
    diag(inf_graph) <- 0
    max_deg <- max(colSums(inf_graph))
    indicator <- max(abs(inf_graph_tmp - inf_graph))
    #rej_set_raw = which(sqrt(n)*abs(delta_d) >= crit, arr.ind = T)
    
    #has_reverse <- function(row, mat) {
    #  j <- row[1]
    #  k <- row[2]
    #  any(mat[,1] == k & mat[,2] == j)
    #}
    
    #filtered_rows <- apply(rej_set_raw, 1, function(row) has_reverse(row, rej_set_raw))
    #rej_set = rej_set_raw[filtered_rows,]
    
    #edges_set = rbind(edges_set, rej_set)
    #edges_set_list[[t]] = edges_set
    #inf_graph = matrix(0,p,p)
    #inf_graph[edges_set] = 1
    #diag(inf_graph) = 0
    #max_deg = max(colSums(inf_graph))
    #indicator = max(abs(inf_graph_tmp - inf_graph))
    t <- t+1
    if((indicator==0)|(max_deg>k)){break}
  }
  
  typeII <- as.numeric(max_deg <= k)
  
  results <- list(inf_graph = inf_graph, typeII = typeII, m = m, crit = crit, attributes = attributes)
  
  results
}