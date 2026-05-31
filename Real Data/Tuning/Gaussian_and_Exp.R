#This script reproduces the parameter tuning results for Gaussian and non-positive exponential
#graphical models.

packages <- c("mvtnorm", "Matrix", "MASS", "clime","maotai")
# Loop to check and install missing packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Verification
message("All packages loaded successfully.")


#functions used
{
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
get_element = function(data){
  # this function is used to get gamma and g for single normal conditional graph
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  
  n = dim(data)[1]
  p = dim(data)[2]
  s = 2*p + choose(p,2)
  
  pos_mat = matrix(1, p,p)
  pos_mat[upper.tri(pos_mat)] = 0
  diag(pos_mat) = 0
  tmat_lower = which(pos_mat !=0, arr.ind = T)
  tmat = tmat_lower[,c(2,1)]
  
  tmp_deri = function(x, j){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2, x^2_jx^2_k), j is the corresponding coordinate
    p = length(x)
    deri_upper = rep(0,p)
    deri_upper[j] = 1
    
    deri_mid = rep(0,p)
    deri_mid[j] = 2*x[j]
    
    
    tx_left = x[tmat[,1]]^2
    tx_right = x[tmat[,2]]^2
    
    deri_left = tx_left
    deri_left[which(tmat[,1]==j)] = 2*x[j]
    deri_left[which(tmat[,1]!=j)] = 0
    deri_left = deri_left*tx_right
    
    deri_right = tx_right
    deri_right[which(tmat[,2] == j)] = 2*x[j]
    deri_right[which(tmat[,2] != j)] = 0
    deri_right = deri_right*tx_left
    
    deri_lower = deri_left + deri_right
    
    deri_full = c(deri_upper, deri_mid, 2*deri_lower)
    return(deri_full)
  }
  
  tmp_deri_ii = function(x, deri){
    # a tmp function to get the derivatives of the vector(x1,..xp, x1^2.. xp^2, x^2_jx^2_k), j is the corresponding coordinate
    # deri is the vector we got in tmp_deri
    p = length(x)
    tmp_sub = deri[1:p]
    j = which(tmp_sub!=0)
    deri_ii = deri/x[j]
    deri_ii[1:p] = rep(0,p)
    return(deri_ii)
  }
  
  # now calculate the empirical Gamma and G for the graph
  
  Gamma = matrix(0, s, s)
  G = rep(0, s)
  
  for(i in 1:n){
    x_i = data[i,]
    Gamma_i = matrix(0, s, s)
    G_i = rep(0, s)
    for(j in 1:p){
      deri_i = tmp_deri(x_i,j)
      deri_ii = tmp_deri_ii(x_i, deri_i)
      
      Gamma_i = Gamma_i + deri_i %*% t(deri_i)
      G_i = G_i + deri_ii
      print(paste('i = ', i, '; j = ', j, sep = ''))
    }
    Gamma = Gamma + Gamma_i
    G = G + G_i
  }
  
  return(list(Gamma = Gamma/n, G = G/n))
}
tr = function(mat){
  a = sum(diag(mat))
  return(a)
}
ADMM_gen_block = function(rho, lambda, alpha, threshold, maxit, elem_x, elem_y){
  
  tmp_mat = diag(1,p)
  diag_set_raw = NULL
  for(j in 1:p){
    diag_set_raw = c(diag_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
  }
  diag_set = which(diag_set_raw!=0)
  
  
  gamma_x = elem_x$Gamma
  gamma_y = elem_y$Gamma
  g_x = elem_x$G
  g_y = elem_y$G
  s = dim(gamma_x)[1]
  
  Inv_mat1 = solve(gamma_x + gamma_y + diag(rho, s))
  Inv_mat2 = solve(gamma_y + diag(rho,s))
  
  theta_x = rep(0, s)
  delta = rep(0, s)
  Z1 = rep(0, s)
  Z2 = rep(0, s)
  U1 = rep(0, s)
  U2 = rep(0, s)
  
  error = 1 + threshold
  i = 0
  while((error > threshold)&(i<=maxit)){
    
    theta_0 = c(theta_x, delta)
    theta_x = Inv_mat1%*%(gamma_y%*%delta - g_x - g_y + rho*(Z1 - U1))
    delta = Inv_mat2%*%(gamma_y%*%theta_x + g_y + rho*(Z2 - U2))
    
    A1 = theta_x + U1
    A2 = delta + U2
    
    Z1 = sign(A1)*pmax(0, abs(A1)-lambda/rho)
    Z1[diag_set] = A1[diag_set]
    Z2 = sign(A2)*pmax(0, abs(A2)-lambda*alpha/rho)
    
    U1 = U1 + (theta_x - Z1)
    U2 = U2 + (delta - Z2)
    
    theta_1 = c(theta_x, delta)
    dif = theta_1 - theta_0
    error = norm(dif, type = '2')/norm(theta_0, type = '2')
    i = i+1
    print(paste(lambda, '+',alpha, '+',i, '+', error, sep = ' '))
  }
  return(list(theta_x = Z1, delta = Z2))
}
get_element_block = function(data, trace = T,weights,boot){
  if(is.null(dim(data))){data = matrix(data, nrow = 1, ncol = length(data))}
  n = dim(data)[1]
  p = dim(data)[2]
  
  
  
  Gamma = matrix(0, 2*(p^2),2*(p^2))
  G = rep(0, 2*p^2)
  deri_list=list()
  for(j in 1:p){
    if(trace ==T){
      print(j)}
    pos = ((2*p*(j-1))+1) : ((2*p*(j-1)) + 2*p)
    
    deri_raw = NULL
    for(i in 1:n){
      datai = data[i,]
      dataij = datai[j]
      
      linmat = (datai)%*%t(datai)
      linmat_j = linmat[,j]
      linmat_deri_j = 2*(linmat_j/dataij)
      
      quamat = (datai^2)%*%t(datai^2)
      quamat_j = quamat[,j]
      quamat_deri_j = 2*2*quamat_j/(dataij)
      quamat_deri_j[j] = 1
      
      if(boot==T){
        deri_j = c(linmat_deri_j,quamat_deri_j)*sqrt(as.complex(weights[i]))
      }else{
        deri_j = c(linmat_deri_j,quamat_deri_j)
      }
      
      deri_raw = cbind(deri_raw, deri_j)
      
      lin_deri_ii = rep(0,p)
      lin_deri_ii[j] = 2
      qua_deri_ii = 2*2*quamat_j/(dataij^2)
      qua_deri_ii[j] = 0
      
      deri_j_ii = c(lin_deri_ii, qua_deri_ii)*weights[i]
      G[pos] = G[pos]+deri_j_ii/n
    }
    deri_list[[j]] = deri_raw
    Gamma[pos,pos] = deri_raw%*%t(deri_raw)/n
  }
  return(list(Gamma = Gamma,G = G, deri= deri_list))
}
load('../RealData.Rdata')
load('../breast40.Rdata')
}

#Gaussian Graphical Model: Lambda1, tuning result: 0.4
{
  set.seed(10)
  alpha = 0.5
  rho = 1
  dat_X <- exp_breast_sex_dif[,2:(expset_breast_men_dim+1)]
  dat_Y <- exp_breast_sex_dif[,(2+expset_breast_men_dim):(1+expset_breast_men_dim+expset_breast_women_dim)]
  dat_X = log(dat_X+1)
  dat_Y = log(dat_Y+1)
  dat_X_final = t((dat_X - rowMeans(dat_X)))
  dat_Y_final = t((dat_Y - rowMeans(dat_Y)))
  dat_X_final[which(dat_X_final==0)] = 1e-5
  dat_Y_final[which(dat_Y_final==0)] = 1e-5
  
  x_data = dat_X_final[1:100,]
  x_data1 = dat_X_final[-(1:100),]
  y_data = dat_Y_final[1:100,]
  y_data1 = dat_Y_final[-(1:100),]
  
  n = dim(x_data)[1]
  p = dim(x_data)[2]
  
  
  # es tuning
  lambda_seq = seq(from = 0.2, to = 1.8, by = 0.2)
  eles = get_elements(x_data, y_data)
  n = dim(x_data)[1]
  p = dim(x_data)[2]
  
  
  cv_pos = sample(1:n, replace = F)
  v = 5
  m_tau = list()
  for(i in 1:length(lambda_seq)){
    lambda = lambda_seq[i]
    len = floor(n/v)
    est_i = matrix(0, nrow = p*(p-1)/2, ncol = v)
    for(k in 1:v){
      cv_seq = cv_pos[(1 + (k - 1)*len) : (k*len)]
      x_val = x_data[-cv_seq,]
      y_val = y_data[-cv_seq,]
      eles_val = get_elements(x_val, y_val)
      
      sm_est = ADMM_sym1(rho, lambda, alpha, 1e-3, 1000, eles_val)
      
      
      theta = sm_est$theta_x
      delta = sm_est$delta
      
      est_hat = delta[lower.tri(delta)]
      est_i[,k] = est_hat
    }
    m_tau[[i]] = est_i
  }
  
  es = NULL
  for(i in 1:length(lambda_seq)){
    m_mat = m_tau[[i]]
    m_lambda = rowMeans(m_mat)
    v_lambda = mean(diag(t(m_mat - m_lambda) %*% (m_mat - m_lambda)))
    es = c(es, v_lambda/(m_lambda%*%m_lambda))
  }
  plot(lambda_seq, es)
  lambda = lambda_seq[which.min(es)]
  print(paste0('The optimal lambda is, ', lambda))
}

#Gaussian Graphical Model: Lambda 2, tuning result: 0.2
{
  set.seed(10)
  rho = 1
  dat_X <- exp_breast_sex_dif[,2:(expset_breast_men_dim+1)]
  dat_Y <- exp_breast_sex_dif[,(2+expset_breast_men_dim):(1+expset_breast_men_dim+expset_breast_women_dim)]
  dat_X = log(dat_X+1)
  dat_Y = log(dat_Y+1)
  dat_X_final = t((dat_X - rowMeans(dat_X)))
  dat_Y_final = t((dat_Y - rowMeans(dat_Y)))
  dat_X_final[which(dat_X_final==0)] = 1e-5
  dat_Y_final[which(dat_Y_final==0)] = 1e-5
  
  x_data = dat_X_final[1:100,]
  x_data1 = dat_X_final[-(1:100),]
  y_data = dat_Y_final[1:100,]
  y_data1 = dat_Y_final[-(1:100),]
  
  n = dim(x_data)[1]
  p = dim(x_data)[2]
  
  m_seq = seq(from = 0.2, to = 0.9,by = 0.1)
  
  clime_x = clime(x_data1, lambda = m_seq,standardize = F)
  clime_y = clime(y_data1, lambda = m_seq,standardize = F)
  
  cv_x = cv.clime(clime_x, loss = 'tracel2', fold = 5)
  cv_y = cv.clime(clime_y, loss = 'tracel2', fold = 5)
  lambda = (cv_x$lambdaopt + cv_y$lambdaopt)/2
  print(paste0('The optimal lambda is, ', lambda))
}


#Exponential Graphical Model: Lambda1, tuning result: 1.6
{
  source('../highd_helpers.R')
  rho <- 1
  dat_X <- exp_breast_sex_dif[,2:(expset_breast_men_dim+1)]
  dat_Y <- exp_breast_sex_dif[,(2+expset_breast_men_dim):(1+expset_breast_men_dim+expset_breast_women_dim)]
  dat_X <- log(dat_X+1)+0.25
  dat_Y <- log(dat_Y+1)+0.25
  dat_X_final <- t((dat_X))
  dat_Y_final <- t((dat_Y))
  
  
  x_data <- dat_X_final[1:100,]
  x_data1 <- dat_X_final[-(1:100),]
  y_data <- dat_Y_final[1:100,]
  y_data1 <- dat_Y_final[-(1:100),]
  
  n <- dim(x_data)[1]
  p <- dim(x_data)[2]
  alpha <- 0.5
  type <- 1
  dist_type <- 1
  

  set.seed(10)
  lam_seq <- seq(0.2,1.8, by = 0.2)
  l1_res <- tune_lam1(x_data, y_data, lam_seq, type = type, dist_type = dist_type,
                      alpha = 0.5,rho = 1)
  plot(l1_res)
  lambda_fin <- l1_res[,1][which.min(l1_res[,2])]
  print(paste0('The optimal lambda is, ', lambda_fin))
}


#Exponential Graphical Model: Lambda2, tuning result: 0.2
{
  m_seq <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7, 0.8,0.9)
  eles_x <- get_element_pos(x_data1,h_type = type, dist_type = dist_type)
  eles_y <- get_element_pos(y_data1,h_type = type, dist_type = dist_type)
  
  lambda_x = c()
  lambda_y = c()
  
  x_loss = NULL
  y_loss = NULL
  set.seed(10)
  for(j in 1:p){
    x_sub = eles_x$deri[[j]]
    y_sub = eles_y$deri[[j]]
    
    clime_x = clime(x_sub, lambda = m_seq, standardize = F)
    clime_y = clime(y_sub, lambda = m_seq, standardize = F)
    
    cv_x = cv.clime(clime_x, loss = 'tracel2', fold = 5)
    cv_y = cv.clime(clime_y, loss = 'tracel2', fold = 5)
    
    lambda_x = c(lambda_x, cv_x$lambdaopt)
    lambda_y = c(lambda_y, cv_y$lambdaopt)
    print(j)
    print(cv_x$lambdaopt)
    print(cv_y$lambdaopt)
  }
  
  lambda_fin <- round(mean(c(lambda_x,lambda_y)),1)
  print(paste0('The optimal lambda is, ', lambda_fin))
}