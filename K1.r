##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP") # make this more flexible
set.seed(47269801)

required.packages <- c("dplyr", "lmtest", "sandwich", "stargazer", "ggplot2", "xtable", "pwr", "wesanderson", "extrafont", "multiwayvcov", "multcomp", "knitr")
packages.missing <- required.packages[!required.packages %in% installed.packages()[,"Package"]]

if(length(packages.missing) > 0) {install.packages(required.packages, repo="https://cran.cnr.berkeley.edu/")}
lapply(required.packages, library, character.only = TRUE)

######################
## Define functions ##
######################

## RegTest tests primary hypotheses from linear model ##

RegTest <- function(equation, clustvars, hypotheses, data) {

    model <- lm(equation, data = data)

    if (missing(clustvars)) model$vcov <- vcov(model)
    else model$vcov <- cluster.vcov(model, cluster = clustvars)

    model$test <- summary(glht(model, linfct = hypotheses, vcov = model$vcov))$test

    numhyp <- length(hypotheses)

    EST <- matrix(nrow = numhyp, ncol = 4)

    for (i in 1:numhyp) {

        EST[i, 1] <- model$test$coefficients[i]
        EST[i, 2] <- model$test$tstat[i]
        EST[i, 3] <- model$test$sigma[i]
        EST[i, 4] <- model$test$pvalues[i]

    }

    colnames(EST) <- c("Estimate", "Tstat", "SE", "P")

    return(EST)

}

## PermTest returns approximations of the exact p-value ##

PermTest <- function(equation, treatvars, clustvars, hypotheses, iterations, data) {

    stopifnot(length(hypotheses) <= 1)

    obsEST <- RegTest(equation, clustvars, hypotheses, data)
    obsStat <- obsEST[1, 2]

    simEST <- matrix(ncol = 4)

    for (i in 1:iterations) {

        simTreat <- data[, treatvars, drop = FALSE]
        simTreat <- simTreat[sample(nrow(simTreat)),]

        simData <- cbind(simTreat, data[, !(names(data) %in% treatvars), drop = FALSE])
        colnames(simData) <- c(treatvars, colnames(simData)[(1+length(treatvars)):ncol(simData)])

        simEST <- rbind(simEST, RegTest(equation, clustvars, hypotheses, data = simData))

    }

    simSTAT <- simEST[2:nrow(simEST), 2]
    countSTAT <- matrix(abs(simSTAT) >= abs(obsStat), ncol = 1)

    ExactP <- matrix(1, nrow = 1, ncol = nrow(countSTAT)) %*% countSTAT
    ExactP <- ExactP / iterations

    EST <- cbind(obsEST, ExactP)

    colnames(EST) <- c("Estimate", "Tstat", "SE", "P", "ExactP")

    return(EST)

}

## FDR returns minimum q-values ##

FDR <- function(pvals) {

    return pvals

}

################
## Clean data ##
################

#Create locals for simulation
  OBS <- 1000
  TE  <- .8
  HET <- .4

#Generate treatment
  Treat <- sample(0:1,OBS,rep = TRUE,prob = c(.5,.5)) %>%
  factor(levels = c(0,1), labels = c("Control","Treatment"))

#Generate gender
 Gen <- sample(0:1,OBS,rep = TRUE,prob = c(.5,.5))  %>%
  factor(levels = c(0,1), labels = c("Male","Female"))

#Generate factor variable measuring highest level of education
 Edu <- sample(1:3,OBS,rep = TRUE,prob = c(.5,.3,.2)) %>%
  factor(levels = c(1,2,3), labels = c("Primary school","High school","University & above"))

#Generate income
 LnInc <- rnorm(OBS, mean = 5, sd = 1)
 Inc <- exp(LnInc)

#Generate y with notreatment effect
  y_nottreat <- rnorm(OBS, 0, 1)

#Generate outcome with noisy treatment effect of ___
  y_Teffect <- rnorm(OBS, TE, 1)
  y_Teffect[Treat == "Control"] <- 0
  y_treated = y_nottreat + y_Teffect

#Generate outcome with noisy treatment effect of ___ and noisy het of ___ gender
  y_GenTeffect <- rnorm(OBS, HET, 1)
  y_GenTeffect[Treat == "Control"] <- 0
  y_GenTeffect[Gen == "Male"] <- 0
  y_HetTreated <- y_treated + y_GenTeffect

#Generate id
  ID <- matrix(1:OBS, ncol = 1)

#Create, save dataframe
  TestData <- data.frame(ID, Treat, Gen, Edu, Inc, y_nottreat, y_treated, y_HetTreated)

################
## Estimation ##
################

depvars <- c("y_nottreat", "y_treated", "y_HetTreated")
hypotheses <- c("TreatTreatment = 0", "TreatTreatment = 1", "TreatTreatment = 0.8")

RegTest(y_nottreat ~ Treat, clustvars = TestData$ID, hypotheses = c("TreatTreatment = 0"), data = TestData)
PermTest(y_nottreat ~ Treat, treatvars = c("Treat"), clustvars = TestData$ID, hypotheses = c("TreatTreatment = 0"), iterations = 100, data = TestData)
