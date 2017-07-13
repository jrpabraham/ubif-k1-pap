##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP") # make this more flexible

required.packages <- c("dplyr", "lmtest", "sandwich", "stargazer", "ggplot2", "xtable", "pwr", "wesanderson", "extrafont", "multiwayvcov", "multcomp", "knitr")
packages.missing <- required.packages[!required.packages %in% installed.packages()[,"Package"]]

if(length(packages.missing) > 0) {install.packages(required.packages, repo="https://cran.cnr.berkeley.edu/")}
lapply(required.packages, library, character.only = TRUE)

######################
## Define functions ##
######################

## RegTest tests primary hypotheses from linear model ##

RegTest <- function(equation, clustvars, hypotheses) {

    model <- lm(equation)

    if (missing(clustvars)) model$vcov <- vcov(model)
    else model$vcov <- cluster.vcov(model, cluster = clustvars)

    model$test <- summary(glht(model, linfct = hypotheses, vcov = model$vcov))$test

    numhyp <- length(hypotheses)

    EST <- matrix(nrow = numhyp, ncol = 3)

    for (i in 1:numhyp) {

        EST[i, 1] <- model$test$coefficients[i]
        EST[i, 2] <- model$test$sigma[i]
        EST[i, 3] <- model$test$pvalues[i]

    }

    colnames(EST) <- c("Estimate", "SE", "P")

    return(EST)

}

## PermTest returns approximations of the exact p-value ##

PermTest <- function(equation, clustvars, hypotheses) {

    model <- lm(equation)

    if (missing(clustvars)) model$vcov <- vcov(model)
    else model$vcov <- cluster.vcov(model, cluster = clustvars)

    model$test <- summary(glht(model, linfct = hypotheses, vcov = model$vcov))$test

    numhyp <- length(hypotheses)

    EST <- matrix(nrow = numhyp, ncol = 3)

    for (i in 1:numhyp) {

        EST[i, 1] <- model$test$coefficients[i]
        EST[i, 2] <- model$test$sigma[i]
        EST[i, 3] <- model$test$pvalues[i]

    }

    colnames(EST) <- c("Estimate", "SE", "P")

    return(EST)

}

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
