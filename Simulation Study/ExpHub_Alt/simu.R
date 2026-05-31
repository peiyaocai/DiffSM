args <- commandArgs(TRUE)

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
sig <- as.numeric(args[[3]])
s <- as.numeric(args[[4]])
type <- as.numeric(args[[5]])

source('../highd_helpers.R')
highd_pos_II <- function(n,p,alpha,seed,m,sig,type,dist_type,k = 5){
  #basic settings
  l1<-5
  l2<-0.5
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
    sig_coef <- -(1 +0.5*(abs(sig)-1))
    #sig_coef <- -(abs(sig)+1)
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
if(!dir.exists('../Results')){dir.create('../Results')}
base_res <- paste0('../Results/ExpAlt_n',n,'p',p,'sig',sig)
if(!dir.exists(base_res)){dir.create(base_res)}
for( seed in ((s-1)*20 + 1): ((s-1)*20 + 20)){
res_dir <- paste0(base_res, '/seed', seed, 'type',type,'.rds')
res <- highd_pos_II(n = n,p = p,alpha = 0.5,seed = seed,
                   m = 6,sig = sig, type = type, dist_type = 1, k=5)
saveRDS(res,res_dir)
}
