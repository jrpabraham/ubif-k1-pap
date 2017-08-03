##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP")
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
        colnames(simData)[1:length(treatvars)] <- treatvars

        simEST <- rbind(simEST, RegTest(equation, clustvars, hypotheses, data = simData))

    }

    simSTAT <- simEST[2:nrow(simEST), 2]
    countSTAT <- matrix(abs(simSTAT) >= abs(obsStat), ncol = 1)

    ExactP <- colSums(countSTAT) / iterations

    EST <- cbind(obsEST, ExactP)

    colnames(EST) <- c("Estimate", "Tstat", "SE", "P", "ExactP")

    return(EST)

}

## FDR returns minimum q-values ##

FDR <- function(pvals, step) {

    if (sum(is.na(pvals) == FALSE) <= 1) {return(pvals)}
    if (missing(step)) {step <- 0.001}

    allpvals <- cbind(as.matrix(pvals), matrix(1:nrow(as.matrix(pvals)), ncol = 1))

    pvals <- na.omit(allpvals)
    nump <- nrow(pvals)

    pvals <- pvals[order(pvals[, 1]), ]
    rank <- matrix(1:nump, ncol = 1)
    pvals <- cbind(pvals, rank, matrix(0, nrow = nump, ncol = 1))

    qval <- 1

    while (qval > 0) {

        qfirst <- qval / (1 + qval)
        fdrtemp <- (qfirst * rank) / nump

        subrank <- which(fdrtemp >= as.matrix(pvals[, 1]))

        if (length(subrank) < 1) {
            numreject <- 0
        } else numreject <- max(subrank)

        qsec <- qfirst * (nump / (nump - numreject))
        fdrtemp <- (qsec * rank) / nump

        subrank <- which(fdrtemp >= as.matrix(pvals[, 1]))

        if (length(subrank) < 1) {
            numreject <- 0
        } else numreject <- max(subrank)

        pvals[which(pvals[, 3] <= numreject), 4] <- qval

        qval <- qval - step

    }

    pvals <- pvals[order(pvals[, 2]), ]

    qvals <- matrix(nrow = nrow(allpvals), ncol = 1)
    qvals[match(pvals[, 2], allpvals[, 2]), 1] <- pvals[, 4]

    return(as.matrix(qvals))

}

################
## Clean data ##
################

data <- read.delim(file = "K1__Field_Survey_v34.csv")
data <- as.data.frame(data[2:nrow(data), ])

## Participant ID ##

data$survey.id <- as.numeric(as.character(data$numb1))
data$survey.id[data$survey.id == NA] <- as.numeric(as.character(data$numb2))
data <- data[complete.cases(data$survey.id), ]

## Treatment assignment ##

data$treatment[data$condition == "poor"] <- 0
data$treatment[data$condition == "individual"] <- 1
data$treatment[data$condition == "community"] <- 2

data$poor <- ifelse(data$condition == "poor", 1, 0)
data$ind <- ifelse(data$condition == "individual", 1, 0)
data$com <- ifelse(data$condition == "community", 1, 0)

## Self-efficacy ##

data$sel.score <- as.numeric(as.character(data$sel.con)) + as.numeric(as.character(data$sel.pers)) + as.numeric(as.character(data$sel.com)) + as.numeric(as.character(data$sel.prob)) + as.numeric(as.character(data$sel.bett))
data$sel.score.z <- (data$sel.score - mean(data$sel.score)) / sd(data$sel.score)

## Judgement ##

data$jud.score <- as.numeric(as.character(jud.judg)) + as.numeric(as.character(jud.emb)) + as.numeric(as.character(jud.ups)) + as.numeric(as.character(jud.fam)) + as.numeric(as.character(jud.com))
data$jud.score.z <- (data$jud.score - mean(data$jud.score)) / sd(data$jud.score)

## Sociodemographics ##

data$soc.age <- as.numeric(as.character(data$soc.age))
data$soc.fem <- ifelse(data$soc.gen == "2", 1, 0)
data$soc.chr <- ifelse(data$soc.rel == "1", 1, 0)

################
## Estimation ##
################

attach(data)

## Plain OLS ##

hypotheses <- c("ind = 0", "com = 1", "ind - com = 0")
equations <- c("sel.score ~ ind + com", "jud.score ~ ind + com")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (eqn in equations) {

        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = survey.id, hypotheses = c(h), iterations = 1000, data = data))

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

hypotheses <- c("ind = 0", "com = 1", "ind - com = 0")
equations <- c("sel.score ~ ind + com + soc.age", "jud.score ~ ind + com + soc.age")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (eqn in equations) {

        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = survey.id, hypotheses = c(h), iterations = 1000, data = data))

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
