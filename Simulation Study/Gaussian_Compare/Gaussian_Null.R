args <- commandArgs(TRUE)

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
sig <- as.numeric(args[[3]])
s <- as.numeric(args[[4]])
m <- as.numeric(args[[5]])
l1 <- 5
l2 <- 0.5
alpha <- 0.5
rho <- 1


source('gaussian_helpers.R')
if(!dir.exists('R_res')){dir.create('R_res')}
base_res <- paste0('R_res/Null_n',n,'p',p,'sig',sig,'m',m)
if(!dir.exists(base_res)){dir.create(base_res)}
for(seed in ((s-1)*1+1):((s-1)*1+1 )){
res <- gaussian_I(n =n, p = p, alpha = alpha, seed = seed, sig = sig,l1 = l1, l2 = l2, m = m)
res_dir <- paste0(base_res, '/seed', seed, '.rds')
saveRDS(res,res_dir)
}
