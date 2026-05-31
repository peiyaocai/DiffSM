library(tidyr)
library(ggplot2)
#Reproducing Table 1
{
#Normal Conditional Part
n<-c(400,500,600)
p<-c(50,60,70)
sig<-1
m<-5
nc_rej <- rep(0,length(n))
for(i in 1:length(n)){
  ni <- n[i]
  pi <- p[i]
base_res <- paste0('Results/Null_n',ni,'p',pi,'sig',sig,'m',m)
rej_vec <- rep(0,300)
for(seed in 1:300){
  res_dir <- paste0(base_res, '/seed', seed, '.rds')
  tmp <- readRDS(res_dir)
  rej_vec[seed] <- tmp$typeI_rej
}
nc_rej[i] <-mean(rej_vec)
}
#Reproducing Table 1: Exponential Graphical Model Part
n<-c(400,500,600)
p<-c(20,30,40)
sig<--2
type <- 1
exp_rej <- rep(0,3)
for(i in 1:3){
  ni <- n[i]
  pi <- p[i]
base_res <- paste0('Results/ExpNull_n',ni,'p',pi,'sig',sig)
rej_vec <- rep(0,300)
deg_vec <- rep(0,300)
for(seed in 1:300){
  res_dir <- paste0(base_res, '/seed', seed,'type',type ,'.rds')
  tmp <- readRDS(res_dir)
  rej_vec[seed] <- tmp$typeI_rej
}
exp_rej[i]<-mean(rej_vec)
}
exp_rej

#Output Table1
print(rbind(nc_rej,exp_rej))
}




#Reproducing Figure 1(a)
{
n<-c(400,500,600)
p<-50
sig<-c(1,2,3,4,5)
m <- 6
power_all <- rep(0,length(sig))
nc_power <- NULL
for(nn in n){
for(i in 1:length(sig)){
  si <- sig[i]
  base_res <- paste0('Results/Alt_n',nn,'p',p,'sig',si,'m',m)
  t2_vec <- rep(0,300)
  for(seed in (1:300)){
    res_dir <- paste0(base_res, '/seed', seed ,'.rds')
    tmp <- readRDS(res_dir)
    t2_vec[seed] <- tmp$typeII
  }
  power_all[i]<-1-mean(t2_vec)
}
nc_power <- rbind(nc_power,power_all)
}
nc_power
df <- data.frame(
  sig = sig,
  `n = 400` = nc_power[1,],
  `n = 500` = nc_power[2,],
  `n = 600` = nc_power[3,],
  check.names = FALSE
)

df_long <- pivot_longer(
  df,
  cols = c(`n = 400`, `n = 500`, `n = 600`),
  names_to = "Sample_Size",
  values_to = "Power"
)


ggplot(df_long, aes(x = sig, y = Power, color = Sample_Size, linetype = Sample_Size)) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 4) +
  scale_color_manual(
    name = "Sample size",
    breaks = c("n = 400", "n = 500", "n = 600"),
    values = c(
      "n = 400" = "red",
      "n = 500" = "darkgreen",
      "n = 600" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Sample size",
    breaks = c("n = 400", "n = 500", "n = 600"),
    values = c(
      "n = 400" = "solid",
      "n = 500" = "dotted",
      "n = 600" = "dashed"
    )
  ) +
  labs(
    x = "Signal Strength",
    y = "Power"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.position = "bottom",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 3
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 5,
        size = 5
      )
    )
  )
}

#Reproducing Figure 1(b)
{
n<-c(400,500,600)
p<-30
sig<-c(2,3,4,5,6)
m <- 6
power_all <- rep(0,length(sig))
exp_power <- NULL
for(nn in n){
for(i in 1:length(sig)){
  si <- sig[i]
  base_res <- paste0('Results/ExpAlt_n',nn,'p',p,'sig',si)
  t2_vec <- rep(0,300)
  for(seed in 1:300){
    res_dir <- paste0(base_res, '/seed', seed ,'type1','.rds')
    tmp <- readRDS(res_dir)
    t2_vec[seed] <- tmp$typeII
  }
  power_all[i]<-1-mean(t2_vec)
}
exp_power <- rbind(exp_power, power_all)
}
exp_power
df <- data.frame(
  sig = sig,
  `n = 400` = exp_power[1,],
  `n = 500` = exp_power[2,],
  `n = 600` = exp_power[3,],
  check.names = FALSE
)

df_long <- pivot_longer(
  df,
  cols = c(`n = 400`, `n = 500`, `n = 600`),
  names_to = "Sample_Size",
  values_to = "Power"
)

ggplot(df_long, aes(x = sig, y = Power, color = Sample_Size, linetype = Sample_Size)) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 4) +
  scale_color_manual(
    name = "Sample size",
    breaks = c("n = 400", "n = 500", "n = 600"),
    values = c(
      "n = 400" = "red",
      "n = 500" = "darkgreen",
      "n = 600" = "blue"
    )
  ) +
  scale_linetype_manual(
    name = "Sample size",
    breaks = c("n = 400", "n = 500", "n = 600"),
    values = c(
      "n = 400" = "solid",
      "n = 500" = "dotted",
      "n = 600" = "dashed"
    )
  ) +
  labs(
    x = "Signal Strength",
    y = "Power"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.position = "bottom",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 3
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 5,
        size = 5
      )
    )
  )

}


#Reproducing Table 1 in online supplementary material
{
#Normal Conditional Part
n<- c(400,500,600)
p <- c(50,60,70)
sig <- 1
m <- 5
t1_vec <- rep(0,length(n))
t2_vec <- rep(0,length(n))
cis_vec <- rep(0,length(n))
cisc_vec <- rep(0,length(n))
for(i in 1:length(n)){
  ni <- n[i]
  pi <- p[i]
  tmp_mat <- diag(1,pi)
  tmp_vec <- c()
  for(j in 1:pi){
    tmp_vec <- c(tmp_vec, tmp_mat[j,])
  }
  nuisance <- which(tmp_vec!=0)
  base_res <- paste0('Results/Single_n',ni,'p',pi,'sig',sig,'m',m)
  t1<- rep(0,300)
  t2 <- rep(0,300)
  cis <- rep(0,300)
  cisc <- rep(0,300)
  for(seed in 1:300){
    res_dir <- paste0(base_res, '/seed', seed ,'.rds')
    tmp <- readRDS(res_dir)
    
    supp <- tmp$r2$supp
    supp_c <- setdiff(tmp$r2$supp_c,nuisance)
    
    clb <- tmp$r2$ci[[1]]
    cub <- tmp$r2$ci[[2]]
    
    ci_len <- cub - clb
    cis[seed] <- sum(ci_len[supp])/length(supp)
    cisc[seed] <- sum(ci_len[supp_c])/length(supp_c)
    
    t1[seed] <- tmp$r2$typeI
    t2[seed] <- tmp$r2$typeII
  }
  t1_vec[i] <- mean(t1)
  t2_vec[i] <- mean(t2)
  cis_vec[i] <- mean(cis)
  cisc_vec[i] <- mean(cisc)
}
nc_res <- round(cbind(t1_vec,1-t2_vec,cis_vec,cisc_vec),3)
nc_res
#Reproducing Table 1 in online supplementary material: Exponential Graphical Model Part
n<- c(400,500,600)
p <- c(20,30,40)
sig <- -2
t1_vec <- rep(0,length(n))
t2_vec <- rep(0,length(n))
cis_vec <- rep(0,length(n))
cisc_vec <- rep(0,length(n))
type <- 1
for(i in 1:length(n)){
  ni <- n[i]
  pi <- p[i]
  tmp_mat <- as.vector(diag(1,pi))
  nuisance <- which(tmp_mat!=0)
  base_res <- paste0('Results/ExpSig_n',ni,'p',pi,'sig',sig)
  t1<- rep(0,300)
  t2 <- rep(0,300)
  cis <- rep(0,300)
  cisc <- rep(0,300)
  for(seed in 1:300){
    res_dir <- paste0(base_res, '/seed', seed ,'type',type,'.rds')
    tmp <- readRDS(res_dir)
    
    supp <- tmp$r2$supp
    supp_c <- setdiff(tmp$r2$supp_c,nuisance)
    
    clb <- tmp$r2$ci[[1]]
    cub <- tmp$r2$ci[[2]]
    
    ci_len <- cub - clb
    cis[seed] <- sum(ci_len[supp])/length(supp)
    cisc[seed] <- sum(ci_len[supp_c])/length(supp_c)
    
    #res <- get_cov_sym(clb, cub, Delta_true)
    t1[seed] <- tmp$r2$typeI
    t2[seed] <- tmp$r2$typeII
  }
  t1_vec[i] <- mean(t1)
  t2_vec[i] <- mean(t2)
  cis_vec[i] <- mean(cis)
  cisc_vec[i] <- mean(cisc)
}
exp_res <- round(cbind(t1_vec,
1-t2_vec,
cis_vec,
cisc_vec),3)
exp_res

#Output Table 1 in online supplementary material
print(rbind(nc_res,exp_res))}




#Reproducing Tables 3--5 in online supplementary material
{
n <- c(400,500,600)
p <- c(50,60,70)
sig <- 2
m<-5
rej_all_sm <- rep(0,length(n))
rej_all_jl <-  rep(0,length(n))
degs_all <- matrix(0,300,length(n))
t1_jl <- rep(0,length(n))
t2_jl <- rep(0,length(n))
t1_sm <- rep(0,length(n))
t2_sm <- rep(0,length(n))
ci_s_sm <- rep(0,length(n))
ci_sc_sm <- rep(0,length(n))
ci_s_jl <- rep(0,length(n))
ci_sc_jl <- rep(0,length(n))
err_sm <- rep(0,length(n))
err_jl <- rep(0,length(n))

for(i in 1:length(n)){
  ni <- n[i]
  pi <- p[i]
  Delta_true <- matrix(0,pi,pi)
  for(j in 1:m){
    Delta_true[1,3*j] <- sig
    Delta_true[j*3,1] <- sig
  }
  Dif_Supp <- (Delta_true!=0)
  supp <- which(Dif_Supp!=0, arr.ind = T)
  supp_c <- which(Dif_Supp==0, arr.ind = T)
  supp_c <- supp_c[which(supp_c[,1]!=supp_c[,2]),]
  
  base_res <- paste0('Gaussian/R_res/Null_n',ni,'p',pi,'sig',sig,'m5')
  base_res_jl <- paste0('Gaussian/Res')
  deg_seq <- rep(0,300)
  rej_sm <- rep(0,300)
  rej_jl <- rep(0,300)
  t1_seq <- rep(0,300)
  t2_seq <- rep(0,300)
  t1_seq_sm <- rep(0,300)
  t2_seq_sm <- rep(0,300)
  
  s_sm <- rep(0,300)
  sc_sm <- rep(0,300)
  s_jl <- rep(0,300)
  sc_jl <- rep(0,300)
  
  err_sq_sm <- rep(0,300)
  err_sq_jl <- rep(0,300)
  for(seed in 1:300){
    print(seed)
    res_dir <- paste0(base_res, '/seed', seed ,'.rds')
    tmp <- readRDS(res_dir)
    est_sm <- tmp$est
    deg_seq[seed] <- max(colSums(tmp$inf_graph))
    rej_sm[seed] <- (max(colSums(tmp$inf_graph))>m)
    
    jl_dir <- paste0(base_res_jl, '/inf_graph_n', ni,'p',pi,'sig',sig,'seed',seed,'.csv')
    inf_jl <- read.csv(jl_dir)
    diag(inf_jl) <- 0
    max_deg <- max(colSums(inf_jl))
    rej_jl[seed] <- (max_deg > m)
    t1_seq_sm[seed] <- tmp$typeI
    t2_seq_sm[seed] <- tmp$typeII
    
    
    clb <- read.csv(paste0(base_res_jl, '/clb_n', ni,'p',pi,'sig',sig,'seed',seed,'.csv'))
    cub <- read.csv(paste0(base_res_jl, '/cub_n', ni,'p',pi,'sig',sig,'seed',seed,'.csv'))
    est_jl <- as.matrix(read.csv(paste0(base_res_jl, '/est_n', ni,'p',pi,'sig',sig,'seed',seed,'.csv')))
    
    err_sq_sm[seed] <- norm(est_sm - Delta_true, type = 'F')/pi
    err_sq_jl[seed] <- norm(est_jl - Delta_true,type = 'F')/pi
    typeI <- sum((Delta_true[supp_c]<clb[supp_c]) | (Delta_true[supp_c] > cub[supp_c]))/dim(supp_c)[1]
    typeII <- sum((0>=clb[supp])*(0<=cub[supp]))/dim(supp)[1]
    t1_seq[seed] <- typeI
    t2_seq[seed] <- typeII
    
    s_sm[seed] <- mean((tmp$ci[[2]] - tmp$ci[[1]])[supp])
    sc_sm[seed] <- mean((tmp$ci[[2]] - tmp$ci[[1]])[supp_c])
    s_jl[seed] <- mean((cub - clb)[supp])
    sc_jl[seed] <- mean((cub-clb)[supp_c])
  }
  degs_all[,i] <- deg_seq
  rej_all_sm[i] <- mean(rej_sm)
  rej_all_jl[i] <- mean(rej_jl)
  t1_jl[i] <- mean(t1_seq,na.rm = T)
  t2_jl[i] <- mean(t2_seq)
  t1_sm[i] <- mean(t1_seq_sm)
  t2_sm[i] <- mean(t2_seq_sm)
  
  ci_s_sm[i] <- mean(s_sm)
  ci_sc_sm[i] <- mean(sc_sm)
  ci_s_jl[i] <- mean(s_jl)
  ci_sc_jl[i] <- mean(sc_jl)
  
  err_sm[i] <- mean(err_sq_sm)
  err_jl[i] <- mean(err_sq_jl)
}
t1_jl
1 - t2_jl
t1_sm
1 - t2_sm
ci_s_sm
ci_sc_sm
ci_s_jl
ci_sc_jl
#Output Table 3 in online supplementary material
table_3 <- cbind(c(t1_sm[1],(1-t2_sm)[1],ci_s_sm[1], ci_sc_sm[1]),
                 c(t1_jl[1],(1-t2_jl)[1],ci_s_jl[1], ci_sc_jl[1]),
                 c(t1_sm[2],(1-t2_sm)[2],ci_s_sm[2], ci_sc_sm[2]),
                 c(t1_jl[2],(1-t2_jl)[2],ci_s_jl[2], ci_sc_jl[2]),
                 c(t1_sm[3],(1-t2_sm)[3],ci_s_sm[3], ci_sc_sm[3]),
                 c(t1_jl[3],(1-t2_jl)[3],ci_s_jl[3], ci_sc_jl[3]))
table_3

#Output Table 4 in online supplementary material
rbind(rej_all_sm, rej_all_jl)



#Ouput Table 5 in online supplementary material
rbind(err_sm, err_jl)
}