#' needed R packages
library("EMMREML")
library("Matrix")
args <- commandArgs(T)

#' get system time for seed and then generate random index
randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}
#' generate train idx and test idx
cvSampleIndex <- function( sampleNum, cross = 5, seed = 1,randomSeed = FALSE ) {
  if(randomSeed == TRUE){
    seed <- randomSeed()
  }
  cv <- cross
  resList <- list()
  
  # leave-one-out
  if( cv == sampleNum ){
    vec <- 1:sampleNum
    for( i in 1:sampleNum ){
      resList[[i]] <- list( trainIdx = vec[-i], testIdx = i, cvIdx = i)
    }
  }else {
    #random samples 
    set.seed(seed)
    index <- sample(1:sampleNum, sampleNum, replace = FALSE )
    step = floor( sampleNum/cv )
    
    start <- NULL
    end <- NULL
    train_sampleNums <- rep(0, cv)
    for( i in c(1:cv) ) {
      start <- step*(i-1) + 1
      end <- start + step - 1
      if( i == cv ) 
        end <- sampleNum
      
      testIdx <- index[start:end]
      trainIdx <- index[-c(start:end)]
      resList[[i]] <- list( trainIdx = trainIdx, testIdx = testIdx, cvIdx = i)
    }
  }
  names(resList) <- paste0("cv",1:cross)
  resList
}
#' main
pheno <- args[1]
G_MATRIX <- args[2]
idd_count <- as.integer(args[3])
output <- args[4]
feature <- unlist(strsplit(pheno,split="[/.]"))
feature <- feature[length(feature)-1]
G <- as.matrix(read.table(file = G_MATRIX))
y <- as.matrix(read.table(file = pheno))
y <- scale(y,center=T,scale=T)
nRep <- 100
corre <- matrix(nrow=nRep+1,ncol=1,NA)#' modified
heritablity <- matrix(nrow=nRep,ncol=1,NA)
var <- matrix(nrow=nRep,ncol=2,NA)
for(i in 1:nRep){
  cvSampleList <- cvSampleIndex(length(y), 10, randomSeed == TRUE)
  cvIdx <- 1
  trainIdx <- cvSampleList[[cvIdx]]$trainIdx
  testIdx <- cvSampleList[[cvIdx]]$testIdx
  trainPheno = y[trainIdx]
  testPheno = y[testIdx]
  funout <- emmreml(y=trainPheno, X=matrix(rep(1, idd_count)[trainIdx], ncol=1), Z=diag(idd_count)[trainIdx,], K=G, varbetahat=TRUE, varuhat=TRUE)
  corre[i,1] <- cor(testPheno-as.vector(funout$betahat), funout$uhat[testIdx])
  var[i,] <- c(funout$Ve,funout$Vu)
  heritablity[i,1] <- funout$Vu/(funout$Ve+funout$Vu)
  rm(funout)
}
corre[nRep+1,1] <- mean(corre[1:nRep,1])
write.table(corre, file = paste0(output,"/",feature,"_corre-gblup.txt"), append = FALSE, row.names = FALSE, col.names = c("correlation"), quote = FALSE)
write.table(var, file = paste0(output,"/",feature,"_var-gblup.txt"), append = FALSE, row.names = FALSE, col.names = c("vare","varu"), quote = FALSE)
write.table(heritablity, file = paste0(output,"/",feature,"_heritablity-gblup.txt"), append = FALSE, row.names = FALSE, col.names = c("heritablity"), quote = FALSE)

