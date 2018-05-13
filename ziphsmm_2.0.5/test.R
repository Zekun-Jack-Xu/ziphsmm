library(ziphsmm)
glogit <- function(p){
  k <- length(p) - 1
  if(k==0) {x <- log(p) - log(1-p)}else{
    x <- rep(NA, k)
    for(j in 1:k) x[j] <- log(p[j+1]) - log(p[1])}
  return(x)
}


ginvlogit <- function(x){
  k <- length(x) + 1
  p <- rep(NA,k)
  den <- 1+sum(exp(x))
  p[1] <- 1 / den
  for(j in 2:k) p[j] <- exp(x[j-1]) / den
  return(p)
} 




####################
###############################
##########################################################
require("Rcpp")
setwd("~/Desktop/myRpackage")
compileAttributes("ziphsmm_2.0.2")
setwd("~/Desktop/myRpackage/ziphsmm_2.0.2")
devtools::document()
devtools::load_all()

devtools::use_build_ignore("cran-comments.md")
devtools::use_build_ignore("test.R")
devtools::use_build_ignore("testcont.R")
devtools::check()
devtools::build()

for(i in 1:20){
   print(pmf_expbase(10,1,i))
   print(pmf_expbase(10,1.5,i))
}
for(i in 1:20){
  print(cdf_expbase(10,1,i))
  print(cdf_expbase(10,2,i))
}
slow <- sapply(1:30,function(h) cdf_expbase(10,1,h))
fast <- sapply(1:30,function(h) cdf_expbase(10,2,h))
plot(slow, type='l')
lines(fast,type='l',col='green')
u <- sapply(1:50,function(k) random_expbase(5,1, 30))

########################
M <- 3
n <- 2000
prior <- c(0.5,0.3,0.2)
dtrate <- c(6,5,4)
dtparm <- matrix(c(0.2,0.1,0.2),nrow=3)
zeroparm <- c(0,-0.2)
emitparm <- matrix(c(4,0.3,5,0.2,6,-0.1),3,2,byrow=TRUE)
tpmparm <- c(1,0.2,0.5,-0.2,0,0.2)

emit_x <- matrix(c(rep(1,1000),rep(0,1000)),nrow=2000,ncol=1)
dt_x <- emit_x
tpm_x <- emit_x
zeroinfl_x <- emit_x
trunc <- c(18,15,10)

re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                   trunc, M, n, dt_x,tpm_x, emit_x, zeroinfl_x)
y <- re$series

rrr <- hsmmfit_exp(y,M,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
            dt_x,zeroinfl_x,emit_x,tpm_x,method="BFGS",control=list(trace=1))
retrieve_hsmm_aft(rrr$working_parm,trunc,M,1,1,1,1)

decode <- hsmmviterbi_exp(y,M, trunc,dtrate,dtparm,
                          prior,zeroparm,emitparm,tpmparm,
                          dt_x, zeroinfl_x, emit_x, tpm_x)
sum(decode!=re$state)

#only cov in aft
tpmparm <- c(1,0.5,0)
zeroparm <- 0
emitparm <- matrix(c(4,5,6),nrow=3)
emit_x <- NULL
tpm_x <- NULL
zeroinfl_x <- NULL
re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                   trunc, M, n, dt_x,tpm_x, emit_x, zeroinfl_x)
y <- re$series
hsmm_cov_loglik_aft(y,prior,dtrate,dtparm,tpmparm,zeroparm,emitparm,
                    trunc,M,0,covdt, matrix(0,2000,1),matrix(1,2000,1), 
                    matrix(1,2000,1))
allparm <- c(log(dtrate),dtparm,c(0,0),
             zeroparm, c(t(emitparm)),tpmparm)
retrieve_hsmm_aft(allparm,trunc,M,0,1,0,0)
hsmm_cov_nllk_aft(y,allparm,trunc,M,0,1,0,0,
                  covdt,matrix(0,2000,1),matrix(1,2000,1), 
                  matrix(1,2000,1))
rrr <- hsmmfit_exp(y,M,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
                   dt_x,zeroinfl_x,emit_x,method="BFGS",control=list(trace=1))


#only cov in tpm
tpm_x <- matrix(c(rep(1,1000),rep(0,1000)),nrow=2000,ncol=1)
tpmparm <- c(1,0.2,0.5,-0.2,0,0.2)
zeroparm <- 0
dtparm <- matrix(c(0,0,0),nrow=3)
emitparm <- matrix(c(4,5,6),nrow=3)
emit_x <- NULL
dt_x <- NULL
zeroinfl_x <- NULL
re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                   trunc, M, n, dt_x,tpm_x, emit_x, zeroinfl_x)
y <- re$series
hsmm_cov_loglik_aft(y,prior,dtrate,dtparm,tpmparm,zeroparm,emitparm,
                    trunc,M,1,matrix(0,2000,1),covtpm,matrix(1,2000,1), 
                    matrix(1,2000,1))
allparm <- c(log(dtrate),c(0,0), #remove dtparm since empty
             zeroparm, c(t(emitparm)),tpmparm)
retrieve_hsmm_aft(allparm,trunc,M,1,0,0,0)
hsmm_cov_nllk_aft(y,allparm,trunc,M,1,0,0,0,
                  matrix(0,2000,1),covtpm,matrix(1,2000,1), 
                  matrix(1,2000,1))
rrr <- hsmmfit_exp(y,M,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
                   dt_x,zeroinfl_x,emit_x,method="BFGS",control=list(trace=1))

#4 states cov in tpm
prior <- c(0.25,0.25,0.25,0.25)
dtrate <- c(5,4,3,2)
trunc <- rep(10,4)
tpm_x <- matrix(c(rep(1,1000),rep(0,1000)),nrow=2000,ncol=1)
tpmparm <- rep(c(0,0.1),8)
zeroparm <- 0
dtparm <- matrix(c(0,0,0,0),nrow=4)
emitparm <- matrix(c(4,5,6,7),nrow=4)
emit_x <- NULL
dt_x <- NULL
zeroinfl_x <- NULL
re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                   trunc, 4, 2000, dt_x,tpm_x, emit_x, zeroinfl_x)
y <- re$series
hsmm_cov_loglik_aft(y,prior,dtrate,dtparm,tpmparm,zeroparm,emitparm,
                    trunc,4,1,matrix(0,2000,1),covtpm,matrix(1,2000,1), 
                    matrix(1,2000,1))
allparm <- c(log(dtrate),c(0,0,0), #remove dtparm since empty
             zeroparm, c(t(emitparm)),tpmparm)
retrieve_hsmm_aft(allparm,trunc,M,2,0,0,0)
hsmm_cov_nllk_aft(y,allparm,trunc,M,2,0,0,0,
                  matrix(0,2000,1),covtpm,matrix(1,2000,1), 
                  matrix(1,2000,1))
rrr <- hsmmfit_exp(y,4,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
                   dt_x,zeroinfl_x,emit_x,tpm_x,method="BFGS",control=list(trace=1))

######
#only cov in aft
#if 2 state, give 0 to tpmparm
prior <- c(0.5,0.5)
tpmparm <- 0
dtrate <- c(6,4)
dtparm <- matrix(c(0.2,0.1),nrow=2)
zeroparm <- 0
emitparm <- matrix(c(4,6),nrow=2)
emit_x <- NULL
tpm_x <- NULL
zeroinfl_x <- NULL
dt_x <- matrix(c(rep(1,1000),rep(0,1000)),nrow=2000,ncol=1)
trunc <- c(15,13)
re <- hsmmsim2_exp(prior,dtrate,dtparm,zeroparm,emitparm,tpmparm,
                   trunc, 2, n, dt_x,tpm_x, emit_x, zeroinfl_x)

y <- re$series
hsmm_cov_loglik_aft(y,prior,dtrate,dtparm,tpmparm,zeroparm,emitparm,
                    trunc,2,0,covdt, matrix(0,2000,1),matrix(1,2000,1), 
                    matrix(1,2000,1))

allparm <- c(log(dtrate),dtparm,0,
             zeroparm, c(t(emitparm))) #2 state , no tpm in all parm
retrieve_hsmm_aft(allparm,trunc,2,0,1,0,0)
hsmm_cov_nllk_aft(y,allparm,trunc,2,0,1,0,0,
                  covdt,matrix(0,2000,1),matrix(1,2000,1), 
                  matrix(1,2000,1))
rrr <- hsmmfit_exp(y,2,trunc,dtrate,dtparm,prior,zeroparm,emitparm,tpmparm,
                   dt_x,zeroinfl_x,emit_x,method="BFGS",control=list(trace=1))

#####################################################################

