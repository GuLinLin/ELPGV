#' needed R packages
library("BGLR")
args<-commandArgs(T)

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
SNP <- args[2]
output <- args[3]
temp <- args[4]

feature <- unlist(strsplit(pheno,split="[/.]"))
feature <- feature[length(feature)-1]
Markers <- as.matrix(read.table(file = SNP))
y <- as.matrix(read.table(file = pheno))
y <- scale(y,center=T,scale=T)
nIter=5000;
burnIn=1000;
thin=3;
saveAt=paste0(temp,"/",feature,"/byes_");
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=Markers,model='BayesA')) #' model = c(BayesA, BayesB, BayesC, RKHS) etc.
nQTL<-30
nRep <- 100
corre <- matrix(nrow=nRep+1,ncol=1,NA)#' modified

for(i in 1:nRep){
  cvSampleList <- cvSampleIndex(length(y),10,randomSeed == TRUE)
  cvIdx <- 1
  trainIdx <- cvSampleList[[cvIdx]]$trainIdx
  testIdx <- cvSampleList[[cvIdx]]$testIdx
  yNA <- y
  yNA[testIdx] <- NA #' set missing value for prediction
  fit <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  corre[i,1] <- cor(y[testIdx], fit$yHat[testIdx])
  rm(fit)
}
corre[nRep+1,1] <- mean(corre[1:nRep,1])
write.table(corre, file = paste0(output,"/",feature,"_corre.txt"), append = FALSE, row.names = FALSE, col.names = c("correlation"), quote = FALSE)

