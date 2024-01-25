# R script for reproducing the analyses in 
#       van Doorn, J.B., Ly, A., Marsman, M., & Wagenmakers, E.-J.  (2020). 
#           Bayesian Rank-Based Hypothesis Testing for the Rank Sum Test, the Signed Rank Test, 
#           and Spearman's rho. Journal of Applied Statistics, 47, 2984-3006.
#
# Online supplementary material at https://osf.io/gny35/
#
# First source the files with the relevant functions
source('rankBasedCommonFunctions.R')
source('rankSumSampler.R') # Wilcoxon rank sum function
source('signRankSampler.R') # Wilcoxon signed-rank function
source('spearmanSampler.R')# Spearman's rho function

# Load the example data
mathDat <- read.csv("DataExampleStudentMath.csv") # for rank sum  & Spearman's rho
progDat <- read.csv("DataExampleProgabide.csv") # For signed-rank

# Example 1: Wilcoxon Rank Sum
# Do students who flunk a math course report drinking more alcohol?
x <- mathDat$Dalc[mathDat$Math.T3 < 10] # Alcohol intake of students who failed
y <- mathDat$Dalc[mathDat$Math.T3 >= 10] # Alcohol intake of students who passed

# Default Cauchy prior width is set to 1/sqrt(2) for the sampler 
rankSumSamples <- rankSumGibbsSampler(xVals = x,
                                      yVals = y, 
                                      nSamples = 5e3)

# Posterior distribution
hist(rankSumSamples$deltaSamples)

# Give the posterior samples for delta to the function below to compute BF01
computeBayesFactorOneZero(rankSumSamples$deltaSamples, whichTest = "Wilcoxon", 
                          priorParameter = 1 / sqrt(2))

# One-sided test, positive values for delta only
computeBayesFactorOneZero(rankSumSamples$deltaSamples, 
                          whichTest = "Wilcoxon", 
                          oneSided = "right",
                          priorParameter = 1 / sqrt(2))

# One-sided test, negative values for delta only
computeBayesFactorOneZero(rankSumSamples$deltaSamples, 
                          whichTest = "Wilcoxon", 
                          oneSided = "left",
                          priorParameter = 1 / sqrt(2))

# Example 2: Wilcoxon Signed Rank
# Does progabide reduce the frequency of epileptic seizures? 
x <- progDat$p1 # Alcohol intake of students who failed
y <- progDat$p4 # Alcohol intake of students who passed

# Default Cauchy prior width is set to 1/sqrt(2) for the sampler 
signedRankSamples <- signRankGibbsSampler(xVals = x,
                                          yVals = y, 
                                          nSamples = 5e3)
# OR
signedRankSamples <- signRankGibbsSampler(xVals = x - y,
                                          testValue = 0,
                                          nSamples = 5e3)


# Posterior distribution
hist(signedRankSamples$deltaSamples, freq = FALSE)

# Give the posterior samples for delta to the function below to compute BF01
computeBayesFactorOneZero(signedRankSamples$deltaSamples, 
                          whichTest = "Wilcoxon",
                          priorParameter = 1 / sqrt(2))

# 

# Example 3: Spearman's rho
# Is performance on a math exam associated with the quality of family relations?
x <- mathDat$Math.T3 # Score for math exam 3
y <- mathDat$famrel # Quality of family relations

# Default beta prior width is set to a = b = 1 for the sampler 
rhoSamples <- spearmanGibbsSampler(xVals = x,
                                   yVals = y, 
                                   nSamples = 5e3)

# Posterior distribution
hist(rhoSamples$rhoSamples, freq = FALSE)

# Give the posterior samples for rho to the function below to compute BF01
computeBayesFactorOneZero(rhoSamples$rhoSamples, 
                          whichTest = "Spearman",
                          priorParameter = 1)