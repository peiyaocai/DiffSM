args <- commandArgs(TRUE)

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
sig <- as.numeric(args[[3]])
seed <- as.numeric(args[[4]])
m <- as.numeric(args[[5]])



source('../highd_helpers.R')
if(!dir.exists('../Results')){dir.create('../Results')}
base_res <- paste0('../Results/Null_n',n,'p',p,'sig',sig,'m',m)
if(!dir.exists(base_res)){dir.create(base_res)}
res_dir <- paste0(base_res, '/seed', seed, '.rds')
res <- highd_nc_I(n = n,p = p,alpha = 0.5,seed = seed,m = m,sig = sig,k = 5,l = 5,l2 = 0.5)
saveRDS(res,res_dir)

