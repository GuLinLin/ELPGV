# ELPGV [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/GuLinLin-JMU/ELPGV/issues) [![](https://img.shields.io/badge/Release-v1.1.0-important.svg)](https://github.com/GuLinLin-JMU/ELPGV/blob/master/ELPGV_0.1.0.tar.gz) [![](https://img.shields.io/badge/license-GPL3.0-blue.svg)](https://github.com/GuLinLin-JMU/ELPGV/blob/master/LICENSE)<br>

## *[E](https://github.com/GuLinLin-JMU/ELPGV)nsemble [L](https://github.com/GuLinLin-JMU/ELPGV)earning method for [P](https://github.com/GuLinLin-JMU/ELPGV)rediction for [G](https://github.com/GuLinLin-JMU/ELPGV)enetic [V](https://github.com/GuLinLin-JMU/ELPGV)alue*<br>

![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
![](https://encrypted-tbn3.gstatic.com/images?q=tbn:ANd9GcS3RzhXKSfXpWhWhvClckwi1Llj1j3HvjKpjvU8CQv4cje23TwS "windows logo")

## Brief introduction <br>
**Background:** Whole genome variants offer sufficient information for genetic prediction of human disease risk, and prediction of animal and plant breeding values. Many sophisticated statistical algorithms have been developed for enhancing the prediction ability; however, each method has its own advantages and disadvantages, so far, no one method can beat others.<br>
**Results:** We herein propose an Ensemble Learning method for Prediction for Genetic Value (ELPGV), which assembles predictions from several basic models such as GBLUP, BayesA, BayesB and BayesCπ, to more accurate predictions. We validated ELPGV with a variety of well-known datasets and a serious of simulated datasets. All revealed that ELPGV was able to significantly enhance the prediction ability than any basic models.<br>
**Conclusions:** We here introduced a new method ELPGV, to predict more accurate genetic value from several different basic models in the context of GS. It performed as well as or better than several commonly used methods in a broad range of tasks. The insights gained through the development, testing, and application of ELPGV in the present study can benefit genomics application of machine learning in the future.<br>

## Version and download <br>
* [Linux Version 0.1.0](https://github.com/GuLinLin-JMU/ELPGV/blob/master/ELPGV_0.1.0.tar.gz) -First version released on Jan, 6th, 2021<br>
* [Windows Version 0.1.0](https://github.com/GuLinLin-JMU/ELPGV/blob/master/ELPGV_0.1.0.zip) -First version released on Jan, 6th, 2021<br>
## Running build-in data (runing in linux)
```R
library("ELPGV")
data(exdata)
train_PredMat = mkg$train
test_PredMat = mkg$test[,-1]
TBV = mkg$test[,1] # True breeding value of testing group
ELPGV_pred = ELPGV(rep_times = 100, interation_times=20, weight_min=0, weight_max=1, 
                   weight_dimension=ncol(test_PredMat), rate_min=-0.01, rate_max=0.01, 
                   paticle_number=20, CR = 1.0, train_PredMat=as.matrix(train_PredMat),
                   test_PredMat=as.matrix(test_PredMat), R = 0.5,IW = 1, AF1 = 2, 
                   AF2 = 2, type ="pcc", select="auto", index=NULL)
PredMat = cbind(TBV,ELPGV_pred)
colnames(PredMat) = c("obs", "BayesA", "BayesB", "BayesCpai", "GBLUP", "ELPGV")
cor(PredMat)
```
## Quick running your data (runing in windows)
```R
library("ELPGV")
#loading training predictions
train_PredMat = read.table(system.file("exampleFile/exampleTrainingFile.txt", package = "ELPGV"), header = T)
#loading testing predictions
test_PredMat = read.table(system.file("exampleFile/exampleTestingFile.txt", package = "ELPGV"), header = T)
ELPGV_pred = ELPGV(rep_times = 100, interation_times=20, weight_min=0, weight_max=1, 
                   weight_dimension=ncol(test_PredMat), rate_min=-0.01, rate_max=0.01, 
                   paticle_number=20, CR = 1.0, train_PredMat=as.matrix(train_PredMat),
                   test_PredMat=as.matrix(test_PredMat), R = 0.5,IW = 1, AF1 = 2, 
                   AF2 = 2, type ="pcc", select="auto", index=NULL)
colnames(ELPGV_pred) = c("BayesA", "BayesB", "BayesCpai", "GBLUP", "ELPGV")
head(ELPGV_pred)
```
## Genomic Best Linear Unbiased Prediction <br>
```bash
$ Rscript Genomic_Best_Linear_Unbiased_Prediction.R phenotype.txt G_Matrix.txt 599 /DRNGS_out_path/
```
## An Rscript for (Bayesian) High-Dimensional Regression <br>
```bash
$ Rscript Bayesian.R phenotype.txt SNP_Matrix.txt /out_path/ /tmp_path/
```

## How to access help <br>
If users have any bugs or issues or any suggestions are available, feel free to contact:<br>
:yum: **Linlin Gu:** linlin-gu@outlook.com <br>
:yum: **Prof. Ming Fang:** fangming618@126.com <br>

## Citation <br>
