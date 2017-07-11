##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP") # make this more flexible

required.packages <- c("dplyr", "lmtest", "sandwich", "stargazer", "ggplot2", "xtable", "pwr", "wesanderson", "extrafont", "car", "knitr")
packages.missing <- required.packages[!required.packages %in% installed.packages()[,"Package"]]

if(length(packages.missing) > 0) {install.packages(required.packages, repo="https://cran.cnr.berkeley.edu/")}
lapply(required.packages, library, character.only = TRUE)

######################
## Define functions ##
######################

## MainFx tests primary hypotheses from linear model ##

MainFx <- function(model, clustvars, hypotheses) {

    require(lmtest)
    require(car)

    if (missing(clustvars)) model$vcov <- vcov(model)
    else model$vcov <- vcovHC(model, type = "HC1", cluster = clustvars)

    lapply(hypotheses, linearHypothesis(model, H))

    for (hyp in hypotheses) {

        linearHypothesis(model, hyp, vcov. = model$vcov)

    }

    return(tests)

}

## PermFx returns approximations of the exact p-value ##

PermFx <- function(model, data, subset, weights, na.action, clustvars, hypotheses, iterations) {

}

## CovFx returns cov. adjusted ATE estimates (subsume by MainFx) ##

## HetFx returns het treatment effects (subsume by MainFx) ##

## FDR returns minimum q-values ##

FDR <- function(pvals) {

    return pvals

}

###############
## Load data ##
###############

cars_data <- read_dta("../data/auto.dta")

## Append all survey versions ##

################
## Clean data ##
################

################
## Estimation ##
################

## Balance checks

## Primary hypotheses tests

# Model1 <- lm(depvar ~ treatvars)
# MainFx(Model1, clustvars = surveyid, hypotheses = c("Ind = 0", "Col = 0", "Ind - Col = 0"))

## Covariate adjustment

# Model1 <- lm(depvar ~ treatvars + controlvars)
# MainFx(Model1, clustvars = surveyid, hypotheses = c("Ind = 0", "Col = 0", "Ind - Col = 0"))

## Randomization inference

# Model1 <- lm(depvar ~ treatvars)
# PermFx(Model1, clustvars = surveyid, hypotheses = c("Ind = 0", "Col = 0", "Ind - Col = 0"))

## Het effects

# Model1 <- lm(depvar ~ treatvars*hetvar)
# MainFx(Model1, clustvars = surveyid, hypotheses = c(""))
