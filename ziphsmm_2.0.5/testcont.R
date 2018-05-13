#devtools::use_package("glmnet") #already run
require("Rcpp")
setwd("~/Desktop/myRpackage")
compileAttributes("ziphsmm_2.0.5")
setwd("~/Desktop/myRpackage/ziphsmm_2.0.5")
devtools::document()
devtools::load_all()
#devtools::use_package("pracma")
###############################################
devtools::use_build_ignore("cran-comments.md")
devtools::use_build_ignore("test.R")
devtools::use_build_ignore("testcont.R")
devtools::check()
devtools::build()
devtools::build_win() #check windows
#upload to http://win-builder.r-project.org/ for additional checks
devtools::release("~/Desktop/myRpackage/ziphsmm_2.0.5")


convr <- function(vec1, vec2, m){
  result <- 0
  for(i in 1:(m-1)) result <- result + vec1[m-i] * vec2[i]
  return(result)
}


smp <- function(parmmat, M, time, ngrid, lower){
  
  pp <- array(0,dim=c(M,M,ngrid))
  pp[,,1] <- diag(1,M) #initialize pp at 0
  
  grids <- seq(lower, time, length=ngrid)
  h <- median(diff(grids))
  
  densities <- array(0, dim=c(M,M,ngrid))
  cdfs <- array(0, dim=c(M,M,ngrid))
  
  for(m in 1:ngrid){
    nextindex <- 1
    for(i in 1:M) {
      for(j in 1:M){
        if(i!=j){
          densities[i,j,m] <- dgamma(grids[m],
                                     shape=parmmat[nextindex,1],
                                     scale=parmmat[nextindex,2])
          cdfs[i,j,m] <- pgamma(grids[m],
                                shape=parmmat[nextindex,1],
                                scale=parmmat[nextindex,2])
          nextindex <- nextindex + 1
        }
      }
    }
  }
  
  
  
  for(m in 2:ngrid){
    for(i in 1:M){
      for(j in 1:M){
        if(i!=j){
          
          pp[i,j,m] <- 0
          
          for(k in 1:M){
            if(k!=i){
              vec1 <- pp[k,j,1:(m-1)]
              #must normalize!!!
              vec2 <- rep(NA, m-1)
              for(n in 1:(m-1)){
                vec2[n] <- densities[i,k,n]*cdfs[i,k,n]/sum(cdfs[i,,n])
              }
              
              pp[i,j,m] <- pp[i,j,m] + h*convr(vec1,vec2,m)
            }
          }
        }
      }
    }
    for(i in 1:M) pp[i,i,m] <- max(0,1 - sum(pp[i,,m]))
    
  }
  return(pp[,,ngrid])
  
}

parmmat <- matrix(c(2,1,2,1,2,1,2,1,2,1,2,1),6,2,byrow=T)
system.time(smp(parmmat,3,time=1,50,0.05))
system.time(smprcpp(parmmat,3,time=1,50,0.05))

smp(parmmat,3,time=4,50,0.05)
smprcpp(parmmat,3,time=4,50,0.05)

parmmat <- matrix(c(2,0.5,1,2,4,2,2,1,3,0.5,2,2),6,2,byrow=T)
smp(parmmat,3,time=1,50,0.05)
smprcpp(parmmat,3,time=1,50,0.05)
smp(parmmat,3,time=10,100,0.05)
smprcpp(parmmat,3,time=1,50,0.05)

smp(matrix(c(3,1,4,2),2,2,byrow=T),2,10,50,0.05)
smprcpp(matrix(c(3,1,4,2),2,2,byrow=T),2,10,50,0.05)



#########################################################################
#consensus
#########################################################################
devtools::use_build_ignore("cran-comments.md")
devtools::use_build_ignore("test.R")
devtools::use_build_ignore("testcont.R")
devtools::check()
devtools::build()
#without covariates
#devtools::use_package("glmnet") #already run
rm(list=ls())
require("Rcpp")
setwd("~/Desktop/myRpackage")
compileAttributes("ziphsmm_2.0.5")
setwd("~/Desktop/myRpackage/ziphsmm_2.0.5")
devtools::document()
devtools::load_all()

#import from Matrix bdiag!!!

set.seed(930518)
nsubj <- 10
ns <- 5040
ylist <- vector(mode="list",length=nsubj)
timelist <- vector(mode="list",length=nsubj)

prior1 <- c(0.5,0.2 ,0.3 )
omega1 <- matrix(c(-0.3,0.2,0.1,
                   0.1,-0.2,0.1,
                   0.15,0.2,-0.35),3,3,byrow=TRUE)
prior2 <- c(0.3,0.3 ,0.4 )
omega2 <- matrix(c(-0.5,0.25,0.25,
                   0.2,-0.4,0.2,
                   0.15,0.3,-0.45),3,3,byrow=TRUE)
emit <- c(50,200,600)
zero <- c(0.2,0,0)

for(n in 1:nsubj){
  
  timeindex <- rep(1,ns)
  for(i in 2:ns) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
  timelist[[n]] <- timeindex
  
  if(n<=5){
    result <- hmmsim.cont(ns, 3, prior1, omega1, emit, zero, timeindex)
    ylist[[n]] <- result$series
  }else{
    result <- hmmsim.cont(ns, 3, prior2, omega2, emit, zero, timeindex)
    ylist[[n]] <- result$series
  }
}

prior_init <- c(0.5,0.2,0.3)
emit_init <- c(50, 225, 650)
zero_init <- 0.2
tpm_init <- matrix(c(-0.3,0.2,0.1,0.1,-0.2,0.1,0.15,0.2,-0.35),3,3,byrow=TRUE)


####two group
M <- 3
priorclust <- NULL
priorclust <- c(1,1,1,1,1,2,2,2,2,2)

tpmclust <- c(1,1,1,1,1,2,2,2,2,2)
zeroclust <- rep(1,10)
emitclust <- rep(1,10)
group <- vector(mode="list",length=2)
group[[1]] <- 1:5; group[[2]] <- 6:10
###


####three group
M <- 3
priorclust <- NULL
tpmclust <- c(1,1,1,2,2,2,3,3,3,3)
zeroclust <- rep(1,10)
emitclust <- rep(1,10)
group <- vector(mode="list",length=3)
group[[1]] <- 1:3; group[[2]] <- 4:6; group[[3]] <- 7:10
###

###only common
M <- 3
tpmclust <- rep(1,10)
tpmclust <- NULL

priorclust <- rep(1,10)
zeroclust <- rep(1,10)
emitclust <- rep(1,10)
group <- vector(mode="list",length=2)
group[[1]] <- 1:5; group[[2]] <- 6:10

##error
M <- 3
priorclust <- NULL
tpmclust <- NULL
tpmclust <- c(rep(1,5),rep(2,5))
zeroclust <- NULL
emitclust <- NULL

group <- vector(mode="list",length=2)
group[[1]] <- 1:5; group[[2]] <- 6:10


#######################################
time <- proc.time()
result <- dist_learn(ylist, timelist, prior_init, tpm_init, 
                     emit_init, zero_init,NULL, rho=100,priorclust,tpmclust,
                     emitclust,zeroclust,group,
                     maxit=10, tol=1e-4, method="CG", print=TRUE)
proc.time() - time

time <- proc.time()
result <- dist_learn(ylist, timelist, prior_init, tpm_init, 
                     emit_init, zero_init,NULL, rho=1000,priorclust,tpmclust,
                     emitclust,zeroclust,group,
                     maxit=10, tol=1e-4, method="CG", print=TRUE)
proc.time() - time

time <- proc.time()
result <- dist_learn(ylist, timelist, prior_init, tpm_init, 
                     emit_init, zero_init,NULL, rho=1,priorclust,tpmclust,
                     emitclust,zeroclust,group,
                     maxit=50, tol=1e-4, ncores=2, method="CG", print=TRUE)
proc.time() - time

################################################
#with covariates, hypothesis testing on a nonsignificant effect
##############################################
rm(list=ls())
require("Rcpp")
setwd("~/Desktop/myRpackage")
compileAttributes("ziphsmm_2.0.5")
setwd("~/Desktop/myRpackage/ziphsmm_2.0.5")
devtools::document()
devtools::load_all()
set.seed(12933)
nsubj <- 20
ns <- 4000
ylist <- vector(mode="list",length=nsubj)
xlist <- vector(mode="list",length=nsubj)
timelist <- vector(mode="list",length=nsubj)

priorparm1 <- 0
priorparm2 <- 1
tpmparm1 <- c(-2,-2)
tpmparm2 <- c(0,0)
zeroparm <- c(-2,0)
emitparm <- c(4,0, 6,0)
zeroindex <- c(1,0)

for(n in 1:nsubj){
  
  xlist[[n]] <- matrix(rep(c(0,1,0,1),rep(1000,4)),nrow=4000,ncol=1)
  
  timeindex <- rep(1,4000)
  for(i in 2:4000) timeindex[i] <- timeindex[i-1] + sample(1:4,1)
  timelist[[n]] <- timeindex
  
  if(n<=10){
    workparm <- c(priorparm1,tpmparm1,zeroparm,emitparm)
  }else{
    workparm <- c(priorparm2,tpmparm2,zeroparm,emitparm)
  }
  
  result <- hmmsim2.cont(workparm,2,4000,zeroindex,emit_x=xlist[[n]],
                         zeroinfl_x=xlist[[n]],timeindex=timeindex)
  ylist[[n]] <- result$series
}

prior_init=c(0.5,0.5)
tpm_init=matrix(c(-0.1,0.1,0.1,-0.1),2,2,byrow=TRUE)
zero_init=0.2
emit_init=c(50,400)

####
M <- 2
priorclust <- NULL
tpmclust <- c(rep(1,10),rep(2,10))
zeroclust <- rep(1,20)
emitclust <- rep(1,20)
slopeclust <- rep(1,20)

priorclust <- c(rep(1,10),rep(2,10))
zeroclust <- NULL
emitclust <- NULL

group <- vector(mode="list",length=2)
group[[1]] <- 1:10; group[[2]] <- 11:20
###

#####4 groups
M <- 2
tpmclust <- rep(1,20)
tpmclust <- NULL
priorclust <- rep(c(1,2,3,4),rep(5,4))
zeroclust <- rep(c(1,2,3,4),rep(5,4))
emitclust <- rep(c(1,2,3,4),rep(5,4))
slopeclust <- rep(1,10)
group <- vector(mode="list",length=4)
group[[1]] <- 1:5; group[[2]] <- 6:10; group[[3]] <- 11:15; group[[4]] <- 16:20

#####common
M <- 2
tpmclust <- rep(1,20)
tpmclust <- NULL
priorclust <- rep(1,20)
zeroclust <- rep(1,20)
emitclust <- rep(1,20)
slopeclust <- rep(1,20)
group <- vector(mode="list",length=1)
group[[1]] <- 1:20


###ERROR case
priorclust <- NULL
tpmclust <- NULL
emitclust <- NULL
zeroclust <- NULL
slopeclust <- NULL
group <- vector(mode="list",length=1)
group[[1]] <- 1:20
###

time <- proc.time()
result <- dist_learn2(ylist, xlist, timelist, prior_init, tpm_init, 
                      emit_init, zero_init, NULL, rho=0.001, priorclust,tpmclust,
                      emitclust,zeroclust,slopeclust,group,
                      maxit=30, tol=1e-4, method="CG",print=TRUE)
proc.time() - time




