# Set working directory
setwd("~/GoogleDrive/LabFolderJohnny")

# Source the functions with the rankSumGibbsSampler and signRankGibbsSampler
source('signRankSampler.R')
source('rankSumSampler.R')

library(effsize)
library(BayesFactor)
library(sn)
library(parallel)
library(logspline)
library(foreach)
library(doMC)


no_cores <- detectCores() - 1

# Set to FALSE if you want to run the full simulation study
testFase <- TRUE

nIter <- 1e4
nBurnin <- 5e3
nReps <- nRuns <- 1e3

allMu <- c(0, 0.5, 1.5)
allN <- c(10, 20, 50)

allScenarios <- c("normal", "cauchy", "cauchyNormal", "uniform", "uniformNormal", "skew",
                  "skewNormal", "logistic", "logisticNormal", "tTwo", "tThree", "tTwoNormal", "tThreeNormal", "mixture")


if (testFase) {
  nIter <- 3e2
  nBurnin <- 1e2
  nReps <- nRuns <- 10
  
  allMu <- c(0,1)
  allN <- c(10,50)
  allScenarios <- c("logistic")
}

nMu <- length(allMu)
nN <- length(allN)

myColnames <- c("N", "mu","obsMeanX", "obsMeanY", "obsMedianY", "obsMedianY", "obsMedianDiff",
                "obsVarX", "obsVarY", "myBF", "moreyBF", "W", "pval_W", "pval_T",
                "cohenD", paste(100*c(0.025, 0.25, 0.5, 0.75, 0.975),"% Q", sep =  ""))
nColsResults <- length(myColnames)
resultsSignRank <- resultsRankSum <- matrix(ncol = nColsResults, nrow = 0, dimnames = list(NULL,  myColnames))

registerDoMC(core=no_cores)

allResults <- list()



analyzeSamples <- function(thisN, thisMu, thisSig = 1, paired = TRUE, nIter = 1e3, nBurnin = 1e3, myCase = "mixture",
                           decorrelateStep = TRUE) {
  
  switch(myCase,
         mixture = {
           myProbs <- c(0.3, 0.7)
           mixMus <- c(0,5)
           popMean <- sum(myProbs * mixMus)
           components <- sample(1:2, prob=myProbs, size=thisN, replace=TRUE)
           x <- rnorm(thisN, mixMus[components], thisSig)
           y <- rnorm(thisN, popMean - thisMu, thisSig)
           },
         cauchy = {
           x <- rcauchy(thisN, scale = 10, location = thisMu)
           y <- rcauchy(thisN, location = 0)
         },
         skew = {
           x <- rsn(thisN, thisMu, alpha = 20)
           y <- rsn(thisN, 0, alpha = 20)
         },
         uniform = {
           x <- runif(thisN, thisMu - 2, thisMu + 2)
           y <- runif(thisN, -2, 2)
         }, 
         logistic = {      
           x <- rlogis(thisN, thisMu, 1)
           y <- rlogis(thisN, 0, 1)
         },
         categorical = {      
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rnorm(thisN, 0, thisSig)
           x <- as.numeric(cut(x, breaks = c(-Inf,seq(-3,3, by = 1), Inf)))
           y <- as.numeric(cut(y, breaks = c(-Inf,seq(-3,3, by = 1), Inf)))
           },
         cauchyNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rcauchy(thisN, location = 0)
         },
         logisticNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rlogis(thisN, location = 0)
         },
         uniformNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- runif(thisN, -2,2)
         },
         skewNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rsn(thisN, alpha= 20)
         },
         tTwo = {
           x <- rt(thisN, df = 2, ncp = thisMu)
           y <- rt(thisN, df = 2)
         },
         tThree = {
           x <- rt(thisN, df = 3, ncp = thisMu)
           y <- rt(thisN, df = 3)
         },
         tTwoNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rt(thisN, df = 2)
         },
         tThreeNormal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rt(thisN, df = 3)
         },
         normal = {
           x <- rnorm(thisN, thisMu, thisSig)
           y <- rnorm(thisN, 0, thisSig)
         })

    thisRow <- try( expr = {

    if (paired) {
      
      mySamples <- unlist(lapply(list(1,2,3), function(a) signRankGibbsSampler(x, y, nSamples = nIter+nBurnin, progBar = FALSE, decorrelate = decorrelateStep)$deltaSamples[nBurnin:(nIter+nBurnin)]))

    } else {
      
      mySamples <- unlist(lapply(list(1,2,3), function(a) rankSumGibbsSampler(x, y, nSamples = nIter+nBurnin, progBar = FALSE, decorrelate = decorrelateStep)$deltaSamples[nBurnin:(nIter+nBurnin)]))
      
    }

    rsQuants <- quantile(mySamples, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    bf01 <- giveBayesFactorZeroOne(mySamples)
    moreyBf <- 1/exp(ttestBF(x, y, rscale = 1/sqrt(2), paired = paired)@bayesFactor$bf)
    
    freqRes <- wilcox.test(x, y, paired = paired)
    testStat <- unname(freqRes$statistic)
    pVal <- unname(freqRes$p.value)
    pValT <- t.test(x, y, paired = paired)$p.value
    cohenD <- cohen.d(x, y, paired = paired)$estimate
    
    obsMeanX <- mean(x) 
    obsMeanY <- mean(y)
    obsMedianX <- median(x)
    obsMedianY <- median(y)
    obsMedianDiff <- median(x-y)
    obsVarX <- var(x)
    obsVarY <- var(y)
    
    return(unname(c(obsMeanX, obsMeanY, obsMedianX, obsMedianY, obsMedianDiff, obsVarX,
                    obsVarY, bf01, moreyBf, testStat, pVal, pValT, cohenD, rsQuants)))
  } , silent = FALSE)
  
  if(!is.numeric(thisRow)) {
    thisRow <- rep(0, (nColsResults-2))
  }
  
  return(thisRow)
}

#########################
######## Run It! ########
#########################
for (thisScenario in allScenarios) {
  
  myFilename <- paste(thisScenario,".Rdata", sep = "")
  resultsSignRank <- resultsRankSum <- matrix(ncol = nColsResults, nrow = 0, dimnames = list(NULL,  myColnames))
  
  for (i in 1:nMu) {
    
    for (j in 1:nN) {
  
      thisMu <- allMu[i]
      thisN <- allN[j]
      thisSig <- 1
      
      myRankSumResult <- foreach(k = 1:nRuns, .combine = 'rbind') %dopar% {
        analyzeSamples(thisN = thisN, thisMu = thisMu, paired = FALSE, myCase = thisScenario, nIter = nIter, nBurnin = nBurnin)
      }
    
      mySignRankResult <- foreach(k = 1:nRuns, .combine = 'rbind') %dopar% {
        analyzeSamples(thisN = thisN, thisMu = thisMu, paired = TRUE, myCase = thisScenario, nIter = nIter, nBurnin = nBurnin)
      }
      
      resultsRankSum <- rbind(resultsRankSum, cbind(rep(thisN, nRuns), rep(thisMu, nRuns), (myRankSumResult)))
      resultsSignRank <- rbind(resultsSignRank, cbind(rep(thisN, nRuns), rep(thisMu, nRuns), (mySignRankResult)))
      
      rownames(resultsRankSum) <- rownames(resultsSignRank) <- NULL
      print(c(thisScenario, thisMu, thisN))
    }
    
  }
  
  resultList <- list(resultsRankSum = resultsRankSum, resultsSignRank = resultsSignRank)
  save(resultList, file = myFilename)
  allResults[[thisScenario]] <- resultList
}

save(allResults, file = "latentTtestSimuResults.Rdata")

