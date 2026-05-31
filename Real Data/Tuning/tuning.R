args = commandArgs(TRUE)
seed = as.numeric(args[[1]])
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

lambda_seq = seq(from = 0.1, to = 1.8, by = 0.1)
set.seed(seed)
cv_pos = sample(1:n, replace = F)
cv_pos1 = sample(1:dim(y_data)[1], replace = F)
v = 5
m_tau = list()
eles_x_list = list()
eles_y_list = list()
len = floor(n/v)

tmp_mat = diag(1,p)
diag_set_raw = NULL
for(j in 1:p){
  diag_set_raw = c(diag_set_raw, c(tmp_mat[j,], tmp_mat[j,]))
}
diag_set = which(diag_set_raw!=0)
s = 2*p^2
for(k in 1:v){
  cv_seq = cv_pos[(1 + (k - 1)*len) : (k*len)]
  cv_seq1 = cv_pos1[(1 + (k - 1)*len) : (k*len)]
  
  eles_val_x = get_element_block(x_data[-cv_seq,],F,rep(1, dim(x_data[-cv_seq,])[1]),F)
  eles_val_y = get_element_block(y_data[-cv_seq1,],F,rep(1, dim(y_data[-cv_seq1,])[1]),F)
  eles_x_list[[k]] = eles_val_x
  eles_y_list[[k]] = eles_val_y
}

for(i in 1:length(lambda_seq)){
  lambda = lambda_seq[i]
  est_i = matrix(0, nrow = s, ncol = v)
  for(k in 1:v){
    
    eles_val_x = eles_x_list[[k]]
    eles_val_y = eles_y_list[[k]]
    
    sm_est = ADMM_gen_block(rho, lambda, alpha, 1e-2, 1000, eles_val_x, eles_val_y)
    
    
    delta = sm_est$delta
    est_i[,k] = delta
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

res1 = cbind(lambda_seq,es)

fin = list(res1 = res1)

if(!dir.exists('results'){dir.create('results')})
saveRDS(fin, paste0('results/','NC_Tuning', seed, '.rds'))

