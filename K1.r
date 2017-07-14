##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP") # make this more flexible
set.seed(47269801)

required.packages <- c("dplyr", "multiwayvcov", "multcomp", "knitr")
packages.missing <- required.packages[!required.packages %in% installed.packages()[,"Package"]]

if(length(packages.missing) > 0) {install.packages(required.packages, repo="https://cran.cnr.berkeley.edu/")}
lapply(required.packages, library, character.only = TRUE)

######################
## Define functions ##
######################

## RegTest conducts asymptotic test from linear model ##

RegTest <- function(equation, clustvars, hypotheses, data) {

    model <- lm(equation, data = data, na.action = na.omit)

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

## PermTest returns MC approximations of the exact p-value ##

PermTest <- function(equation, treatvars, clustvars, hypotheses, iterations, data) {

    stopifnot(length(hypotheses) <= 1)

    obsEST <- RegTest(equation, clustvars, hypotheses, data)
    obsStat <- obsEST[1, 2]

    simEST <- matrix(ncol = 4)

    for (i in 1:iterations) {

        simTreat <- data[, treatvars, drop = FALSE]
        simTreat <- simTreat[sample(nrow(simTreat)),]

        simData <- cbind(simTreat, data[, !(names(data) %in% treatvars), drop = FALSE])
        colnames(simData)[1:2] <- treatvars

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

    return(pvals)

}

################
## Clean data ##
################

## Create locals for simulation ##

  OBS <- 510

## Generate treatment ##

  Treat <- sample(0:2,OBS, rep = TRUE, prob = c(.33, .33, 0.33)) %>%
  factor(levels = c(0, 1, 2), labels = c("Poverty", "Ind.", "Col."))

  Pov <- (Treat == "Poverty") * 1
  Ind <- (Treat == "Ind.") * 1
  Col <- (Treat == "Col.") * 1

## Generate gender ##

 Gen <- sample(0:1,OBS,rep = TRUE,prob = c(.5,.5))  %>%
  factor(levels = c(0,1), labels = c("Male","Female"))

## Generate factor variable measuring highest level of education ##

 Edu <- sample(1:3,OBS,rep = TRUE,prob = c(.5,.3,.2)) %>%
  factor(levels = c(1,2,3), labels = c("Primary school","High school","University & above"))

## Generate income ##

 LnInc <- rnorm(OBS, mean = 5, sd = 1)
 Inc <- exp(LnInc)

## Generate y with notreatment effect ##

  yNull <- rnorm(OBS, 0, 1)

## Generate outcome with effects
  yInd <- (0.8 * Ind) + rnorm(OBS, 0, 1)
  yCol <- (0.4 * Col) + rnorm(OBS, 0, 1)

## Generate id ##

  ID <- matrix(1:OBS, ncol = 1)

## Create, save dataframe ##

  TestData <- data.frame(ID, Treat, Pov, Ind, Col, Gen, Edu, Inc, yNull, yInd, yCol)

################
## Estimation ##
################

## Plain OLS ##

hypotheses <- c("Ind = 0", "Col = 1", "Ind - Col = 0")
equations <- c("yNull ~ Ind + Col", "yInd ~ Ind + Col", "yCol ~ Ind + Col")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (eqn in equations) {

        # RES <- rbind(RES, RegTest(eqn, clustvars = TestData$ID, hypotheses = c(h), data = TestData))
        RES <- rbind(RES, PermTest(eqn, treatvars = c("Treat"), clustvars = TestData$ID, hypotheses = c(h), iterations = 100, data = TestData))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- equations
    colnames(RES)[6] <- "Min. Q"

    print("--------------------------------------------------------------------", quote = FALSE)
    print(paste("H_0:", h), quote = FALSE)
    print(RES, quote = FALSE)
    print("--------------------------------------------------------------------", quote = FALSE)

}

## Covariate adjustment ##

hypotheses <- c("Ind = 0", "Col = 1", "Ind - Col = 0")
equations <- c("yNull ~ Ind + Col + Gen + LnInc", "yInd ~ Ind + Col + Gen + LnInc", "yCol ~ Ind + Col + Gen + LnInc")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (eqn in equations) {

        # RES <- rbind(RES, RegTest(eqn, clustvars = TestData$ID, hypotheses = c(h), data = TestData))
        RES <- rbind(RES, PermTest(eqn, treatvars = c("Treat"), clustvars = TestData$ID, hypotheses = c(h), iterations = 100, data = TestData))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- equations
    colnames(RES)[6] <- "Min. Q"

    print("--------------------------------------------------------------------", quote = FALSE)
    print(paste("H_0:", h), quote = FALSE)
    print(RES, quote = FALSE)
    print("--------------------------------------------------------------------", quote = FALSE)

}

## Het effects ##
