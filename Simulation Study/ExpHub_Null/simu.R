args <- commandArgs(TRUE)

n <- as.numeric(args[[1]])
p <- as.numeric(args[[2]])
sig <- as.numeric(args[[3]])
s <- as.numeric(args[[4]])
type <- as.numeric(args[[5]])


source('../highd_helpers.R')
if(!dir.exists('../Results')){dir.create('../Results')}
base_res <- paste0('../Results/ExpNull_n',n,'p',p,'sig',sig)
if(!dir.exists(base_res)){dir.create(base_res)}
for( seed in ((s-1)*20 + 1): ((s-1)*20 + 20)){
res_dir <- paste0(base_res, '/seed', seed, 'type',type,'.rds')
res <- highd_pos_I(n = n,p = p,alpha = 0.5,seed = seed,
                   m = 5,sig = sig, type = type, dist_type = 1, k=5)
saveRDS(res,res_dir)
}

