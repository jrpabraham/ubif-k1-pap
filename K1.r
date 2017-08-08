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

data <- read.delim(file = "K1__Field_Survey_v34.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE, na.strings = "")
data <- as.data.frame(data[2:nrow(data), ])
attach(data)

## Participant ID ##

data$survey.id <- as.numeric(data$numb1)
data$survey.id[is.na(data$survey.id)] <- as.numeric(data$numb2)
data <- data[complete.cases(data$survey.id), ]

## Treatment assignment ##

data$treatment[data$condition == "poor"] <- 0
data$treatment[data$condition == "individual"] <- 1
data$treatment[data$condition == "community"] <- 2

factor(data$treatment)

data$poor <- ifelse(data$condition == "poor", 1, 0)
data$ind <- ifelse(data$condition == "individual", 1, 0)
data$com <- ifelse(data$condition == "community", 1, 0)

data$msg1 <- recode(as.numeric(data$ORG_MESSAGE), `2` = 0, `3` = 2)
data$msg2 <- recode(as.numeric(data$ORG_MESSAGE_2), `2` = 0)
data$msg3 <- recode(as.numeric(data$ORG_MESSAGE_3), `3` = 2)

## Self-efficacy ##

selvars <- c(data$sel.con, data$sel.pers, data$sel.com, data$sel.prob, data$sel.bett)

for (var in selvars) {

    var[var == "-99"] <- NA

}

data$sel.score <- as.numeric(data$sel.con) + as.numeric(data$sel.pers) + as.numeric(data$sel.com) + as.numeric(data$sel.prob) + as.numeric(data$sel.bett)
data$sel.score.z <- (data$sel.score - mean(data$sel.score)) / sd(data$sel.score)

## Judgement ##

judvars <- c(data$jud.fam, data$jud.com, data$jud.judg, data$jud.emb, data$jud.ups)

for (var in judvars) {

    var[var == "-99"] <- NA

}

data$jud.score <- as.numeric(data$jud.fam) + as.numeric(data$jud.com) + (6 - as.numeric(data$jud.judg)) + (6 - as.numeric(data$jud.emb)) + (6 - as.numeric(data$jud.ups))
data$jud.score.z <- (data$jud.score - mean(data$jud.score)) / sd(data$jud.score)

## Affect ##

affvars <- c(data$aff.pos, data$aff.ash, data$aff.pow, data$aff.fina)

for (var in judvars) {

    var[var == "-99"] <- NA

}

data$jud.score <- as.numeric(data$aff.pos) + as.numeric(data$aff.pow) + (7 - as.numeric(data$aff.ash)) + (7 - as.numeric(data$aff.fina))
data$jud.score.z <- (data$jud.score - mean(data$jud.score)) / sd(data$jud.score)

## Video selection ##

data$vid.imp1 <- ifelse(data$vid.dec1 == "3" | data$vid.dec1 == "5", 1, 0)
data$vid.imp2 <- ifelse(data$vid.dec2 == "3" | data$vid.dec2 == "5", 1, 0)
data$vid.num <- data$vid.imp1 + data$vid.imp2

## Intertemporal choice ##

data$sav.save = ifelse(data$sav.dec == "2" | data$sav.dec == "3", 1, 0)
data$sav.amt[data$sav.dec == "1"] = 0
data$sav.amt[data$sav.dec == "2"] = 100
data$sav.amt[data$sav.dec == "3"] = 200

# variable for inconsistent choice, interval estimates of discounting parameter
# recode refusals?

## Query theory ##

quevars <- c(data$que.rat1, data$que.rat2, data$que.rat3, data$que.rat4, data$que.rat5)

for (var in judvars) {

    var[var == "-99"] <- NA

}

quedf <- cbind(data$que.rat1, data$que.rat2, data$que.rat3, data$que.rat4, data$que.rat5)
data$que.nonm <- apply(quedf, 1, function(x) length(x[is.na(x)]))

# data$que.mrp <- apply(quedf, 1
# data$que.mri <- apply(quedf, 1
#
# data$que.smrd <- (2 * (data$que.mrp - data$que.mri)) / data$que.nonm

## Frame evaluation ##

# should we analyze categorical with multilogit, dummies, continuous?

data$msg.dec <- recode(data$msg.dec, "2" = 1, "1" = 0)

data$eva.poor[data$msg1 == 0] <- as.numeric(data$eva.msg1[data$msg1 == 0])
data$eva.poor[data$msg2 == 0] <- as.numeric(data$eva.msg2[data$msg2 == 0])
data$eva.poor[data$msg3 == 0] <- as.numeric(data$eva.msg3[data$msg3 == 0])

data$eva.ind[data$msg1 == 1] <- as.numeric(data$eva.msg1[data$msg1 == 1])
data$eva.ind[data$msg2 == 1] <- as.numeric(data$eva.msg2[data$msg2 == 1])
data$eva.ind[data$msg3 == 1] <- as.numeric(data$eva.msg3[data$msg3 == 1])

data$eva.com[data$msg1 == 2] <- as.numeric(data$eva.msg1[data$msg1 == 2])
data$eva.com[data$msg2 == 2] <- as.numeric(data$eva.msg2[data$msg2 == 2])
data$eva.com[data$msg3 == 2] <- as.numeric(data$eva.msg3[data$msg3 == 2])

data$eva.vid.poor <- as.numeric(data$eva.rank.vid_8)
data$eva.vid.ind <- as.numeric(data$eva.rank.vid_9)
data$eva.vid.com <- as.numeric(data$eva.rank.vid_10)

data$eva.conf <- as.numeric(data$eva.conf)

data$eva.emp.poor <- as.numeric(data$eva.rank.emp_5)
data$eva.emp.ind <- as.numeric(data$eva.rank.emp_6)
data$eva.emp.com <- as.numeric(data$eva.rank.emp_7)

## Ladder scales ##

data$ses.lad.now <- as.numeric(data$ses.lad.now)
data$ses.lad.y2 <- as.numeric(data$ses.lad.y2)

## Sociodemographics ##

data$soc.age <- as.numeric(as.character(data$soc.age))
data$soc.pri[is.na(data$soc.edu) == FALSE] <- ifelse(as.numeric(data$soc.edu[is.na(data$soc.edu) == FALSE]) >= 4, 1, 0)
data$soc.fem[is.na(data$soc.gen) == FALSE] <- ifelse(data$soc.gen[is.na(data$soc.gen) == FALSE] == "2", 1, 0)
data$soc.chr[is.na(data$soc.rel) == FALSE] <- ifelse(data$soc.rel[is.na(data$soc.rel) == FALSE] == "1", 1, 0)
data$ses.unemp <- ifelse(data$ses.emp == "1", 1, 0)
data$soc.inc <- as.numeric(data$soc.inc)
data$soc.con <- as.numeric(data$soc.con)
data$soc.sav <- as.numeric(data$soc.sav) - 1

## Survey validity ##

data$end.hear <- as.numeric(data$end.hear) - 1

################
## Estimation ##
################

## Plain OLS ##

hypotheses <- c("ind = 0", "com = 1", "ind - com = 0")
depvars <- c("vid.num", "sav.save", "msg.dec")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (depvar in depvars) {

        eqn <- paste(depvar, "~ ind  + com", sep = " ")
        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = data$survey.id, hypotheses = c(h), iterations = 100, data = data))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- depvars
    colnames(RES)[6] <- "Min. Q"

    print("----------------------------------------------------------------", quote = FALSE)
    print(paste("H_0:", h), quote = FALSE)
    print(RES, quote = FALSE)

}

## Covariate adjustment ##

covariates <- c("soc.fem", "soc.age", "soc.pri", "soc.chr", "soc.sav", "ses.unemp", "soc.inc")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (depvar in depvars) {

        eqn <- paste(depvar, "~ ind  + com +", covariates, sep = " ")
        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = data$survey.id, hypotheses = c(h), iterations = 100, data = data))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- depvars
    colnames(RES)[6] <- "Min. Q"

    print("----------------------------------------------------------------", quote = FALSE)
    print(paste("H_0:", h), quote = FALSE)
    print(RES, quote = FALSE)

}

## Het effects ##

hetvars <- c("soc.fem", "soc.age", "soc.pri", "soc.chr", "soc.sav", "ses.unemp")

for (hetvar in hetvars) {

    hypotheses <- c(paste("ind:", hetvar, " = 0"), paste(hetvar, ":com", " = 0"), paste("ind:", hetvar, " - ", hetvar, ":com", " = 0"))

    for (h in hypotheses) {

        RES <- matrix(nrow = 1, ncol = 5)

        for (depvar in depvars) {

            eqn <- paste(depvar, " ~ ind*", hetvar, " + com*", hetvar)
            RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = data$survey.id, hypotheses = c(h), iterations = 100, data = data))

        }

        RES <- RES[2:nrow(RES), 1:ncol(RES)]
        RES <- cbind(RES, FDR(RES[, 4]))

        rownames(RES) <- depvars
        colnames(RES)[6] <- "Min. Q"

        print("----------------------------------------------------------------", quote = FALSE)
        print(paste("H_0:", h), quote = FALSE)
        print(RES, quote = FALSE)

    }

}
