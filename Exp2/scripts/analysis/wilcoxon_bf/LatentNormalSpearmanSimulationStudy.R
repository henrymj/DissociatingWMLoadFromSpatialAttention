# Set working directory
setwd("~/GoogleDrive/LabFolderJohnny")

# Source the functions with the rankSumGibbsSampler and signRankGibbsSampler
source('spearmanSampler.R')

library(effsize)
library(BayesFactor)
library(sn)
library(logspline)
library(parallel)
library(foreach)
library(doMC)
library(copula)
library(hypergeo)

pearsonBayesFactor10 <- function(n, r, kappa=1) {
  logHyperTerm <- log(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=((n+2/kappa)/2), z=r^2)) 
  logResult <- log(2^(1-2/kappa))+0.5*log(pi)-lbeta(1/kappa, 1/kappa) +
    lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+logHyperTerm
  realResult <- exp(Re(logResult))
  return(realResult)
}

no_cores <- detectCores() - 1

# Set to FALSE if you want to run the full simulation study
testFase <- TRUE

nIter <- 1e4
nBurnin <- 5e3
nReps <- nRuns <- 1e3

allTau <- c(0, 0.2, 0.4, 0.6, 0.8)
allN <- c(10, 20, 50)

allScenarios <- c("normal", "clayton", "gumbel", "frank", "t")

if (testFase) {
  nIter <- 50
  nBurnin <- 10
  nReps <- nRuns <- 5
  
  allTau <- c(0, 0.2)
  allN <- c(10)
  allScenarios <- c("normal", "clayton", "gumbel", "frank", "t")
}

nTau <- length(allTau)
nN <- length(allN)
nQuants <- 41
allQuants <- seq(0, 1, length.out = nQuants)

myColnames <- c("N", "tau", "obsPearson", "obsKendall", "obsSpearman","myBF","paraBF",
                "pVal", "pValPara", paste(100*allQuants,"% Q", sep =  ""))
nColsResults <- length(myColnames)
resultsKendall <- resultsSpearman <- matrix(ncol = nColsResults, nrow = 0, dimnames = list(NULL,  myColnames))

registerDoMC(core=no_cores)

allResults <- list()

analyzeSamples <- function(thisN, thisTau, thisSig = 1, useKendall = TRUE, nIter = 1e3, nBurnin = 1e3, myCase = "arch", myCopula,
                           nColsResults, allQuants) {
  
  if (myCase == "arch") {
    
    thisTheta <- ifelse(useKendall, iTau(archmCopula(myCopula, dim = 2), thisTau), iRho(archmCopula(myCopula, dim = 2), thisTau))
    thisTheta <- ifelse(is.na(thisTheta), 0, thisTheta)
    thisCopula <- archmCopula(myCopula, param = thisTheta, dim = 2)
    mydat <- rCopula(n =  thisN, copula = thisCopula)
    
  } else if (myCopula == "normal") {
    
    thisTheta <- ifelse(useKendall, iTau(normalCopula(dim = 2), thisTau), iRho(normalCopula(dim = 2), thisTau))
    mydat <- rCopula(n =  thisN, copula = normalCopula(dim = 2, param = thisTheta))
    
  } else if (myCopula == "t") {
    
    thisTheta <- ifelse(useKendall, iTau(normalCopula(dim = 2), thisTau), iRho(normalCopula(dim = 2), thisTau))
    mydat <- rCopula(n =  thisN, copula = tCopula(dim = 2, param = thisTheta))
    
  }
  
  x <- mydat[, 1]
  y <- mydat[, 2]
  
  thisRow <- try( expr = {
    
    mySamples <- unlist(lapply(list(1,2,3), function(a) spearmanGibbsSampler(x, y, nSamples = nIter+nBurnin, 
                                                                             progBar = FALSE)$rhoSamples[nBurnin:(nIter+nBurnin)]))

    if (useKendall) {
      mySamples <- pearsonToKendall(mySamples)
    } else {
      mySamples <- pearsonToSpearman(mySamples)
    }
    
    rsQuants <- quantile(mySamples, probs = allQuants)
    myLogspline <- logspline(mySamples)
    bf01 <- dlogspline(0, fit = myLogspline) / 0.5
    paraBf <- pearsonBayesFactor10(r = cor(x,y), n = thisN)
    
    freqRes <- cor.test(x, y, method = ifelse(useKendall, "kendall", "spearman"))
    testStat <- unname(freqRes$estimate)
    pVal <- unname(freqRes$p.value)
    pValPara <-  cor.test(x, y)$p.value
    
    obsPearson <- cor(x,y)
    obsKendall <- cor(x,y, method="kendall")
    obsSpearman <- cor(x,y, method="spearman")
    
    return((c(obsPearson, obsKendall, obsSpearman, 
              bf01, paraBf, pVal, pValPara, rsQuants)))
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
  
  myFilename <- paste(thisScenario,"Simulation.Rdata", sep = "")
  resultsKendall <- resultsSpearman <- matrix(ncol = nColsResults, nrow = 0, dimnames = list(NULL,  myColnames))
  
  thisCase <- ifelse((thisScenario == "normal" | thisScenario == "t"), "bla", "arch")
  
  for (i in 1:nTau) {
    
    for (j in 1:nN) {
      
      thisTau <- allTau[i]
      thisN <- allN[j]
      thisSig <- 1
      
      mySpearResult <- foreach(k = 1:nRuns, .combine = 'rbind') %dopar% {
        analyzeSamples(thisN = thisN, thisTau = thisTau, useKendall = FALSE, myCopula = thisScenario, 
                       myCase = thisCase, nIter = nIter, nBurnin = nBurnin, nColsResults = nColsResults,
                       allQuants = allQuants)
      }
      
      
      resultsSpearman <- rbind(resultsSpearman, cbind(rep(thisN, nRuns), rep(thisTau, nRuns), (mySpearResult)))

      rownames(resultsSpearman) <- NULL
      
      print(c(thisScenario, thisTau, thisN))
    }
    
  }
  
  resultList <- list(resultsSpearman = resultsSpearman)
  save(resultList, file = myFilename)
  allResults[[thisScenario]] <- resultList
}

save(allResults, file = "latentCorrelationResults.Rdata")

