##############################
## Install missing packages ##
##############################

setwd("/Users/Justin/Google Drive/UBIF/UBIF_Deliverables/UBIF_PAP/K1_PAP")
set.seed(47269801)

required.packages <- c("dplyr", "multiwayvcov", "multcomp", "reshape2", "knitr")
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

k1_df <- read.delim(file = "K1__Field_Survey_v34.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE, na.strings = "")
k1_df <- as.data.frame(k1_df[2:nrow(k1_df), ])
attach(k1_df)

## Participant ID ##

k1_df$survey.id <- as.numeric(as.character(k1_df$numb1))
k1_df$survey.id[is.na(k1_df$survey.id)] <- as.numeric(as.character(k1_df$numb2))
data <- data[complete.cases(k1_df$survey.id), ]

# sum(duplicated(k1_df$survey.id))

## Treatment assignment ##

k1_df$treatment[k1_df$condition == "poor"] <- 0
k1_df$treatment[k1_df$condition == "individual"] <- 1
k1_df$treatment[k1_df$condition == "community"] <- 2

factor(k1_df$treatment)

k1_df$poor <- ifelse(k1_df$condition == "poor", 1, 0)
k1_df$ind <- ifelse(k1_df$condition == "individual", 1, 0)
k1_df$com <- ifelse(k1_df$condition == "community", 1, 0)

k1_df$msg1 <- recode(as.numeric(k1_df$ORG_MESSAGE), `2` = 0, `3` = 2)
k1_df$msg2 <- recode(as.numeric(k1_df$ORG_MESSAGE_2), `2` = 0)
k1_df$msg3 <- recode(as.numeric(k1_df$ORG_MESSAGE_3), `3` = 2)

## Self-efficacy ##

selvars <- c(k1_df$sel.con, k1_df$sel.pers, k1_df$sel.com, k1_df$sel.prob, k1_df$sel.bett)

for (var in selvars) {

    var[var == "-99"] <- NA

}

k1_df$sel.score <- as.numeric(k1_df$sel.con) + as.numeric(k1_df$sel.pers) + as.numeric(k1_df$sel.com) + as.numeric(k1_df$sel.prob) + as.numeric(k1_df$sel.bett)
k1_df$sel.score.z <- (k1_df$sel.score - mean(k1_df$sel.score)) / sd(k1_df$sel.score)

## Judgement ##

judvars <- c(k1_df$jud.fam, k1_df$jud.com, k1_df$jud.judg, k1_df$jud.emb, k1_df$jud.ups)

for (var in judvars) {

    var[var == "-99"] <- NA

}

k1_df$jud.score <- as.numeric(k1_df$jud.fam) + as.numeric(k1_df$jud.com) + (6 - as.numeric(k1_df$jud.judg)) + (6 - as.numeric(k1_df$jud.emb)) + (6 - as.numeric(k1_df$jud.ups))
k1_df$jud.score.z <- (k1_df$jud.score - mean(k1_df$jud.score)) / sd(k1_df$jud.score)

## Affect ##

affvars <- c(k1_df$aff.pos, k1_df$aff.ash, k1_df$aff.pow, k1_df$aff.fina)

for (var in judvars) {

    var[var == "-99"] <- NA

}

k1_df$aff.score <- as.numeric(k1_df$aff.pos) + as.numeric(k1_df$aff.pow) + (7 - as.numeric(k1_df$aff.ash)) + (7 - as.numeric(k1_df$aff.fina))
k1_df$aff.score.z <- (k1_df$aff.score - mean(k1_df$aff.score)) / sd(k1_df$aff.score)

## Video selection ##

k1_df$vid.imp1 <- ifelse(k1_df$vid.dec1 == "3" | k1_df$vid.dec1 == "5", 1, 0)
k1_df$vid.imp2 <- ifelse(k1_df$vid.dec2 == "3" | k1_df$vid.dec2 == "5", 1, 0)
k1_df$vid.num <- k1_df$vid.imp1 + k1_df$vid.imp2

## Intertemporal choice ##

k1_df$sav.save = ifelse(k1_df$sav.dec == "2" | k1_df$sav.dec == "3", 1, 0)
k1_df$sav.amt[k1_df$sav.dec == "1"] = 0
k1_df$sav.amt[k1_df$sav.dec == "2"] = 100
k1_df$sav.amt[k1_df$sav.dec == "3"] = 200

# variable for inconsistent choice, interval estimates of discounting parameter
# recode refusals?

## Query theory (savings) ##

que_df <- k1_df[names(k1_df) %in% c("survey.id", "que.rat1", "que.rat2", "que.rat3", "que.rat4", "que.rat5")]

k1_df$que.nonm <- apply(que_df[, 1:5], 1, function(x) length(x[is.na(x) == FALSE]))

que_df <- melt(que_df, id = c("survey.id"))
que_df$variable <- as.numeric(que_df$variable)
que_df$value <- as.numeric(que_df$value)
que_df <- dcast(que_df[is.na(que_df$value) == FALSE, ], survey.id ~ value, median, value.var = "variable")
names(que_df) <- c("survey.id", "que.mri", "que.mrp")
k1_df <- merge(k1_df, que_df, all.x = TRUE)

k1_df$que.smrd <- (2 * (k1_df$que.mrp - k1_df$que.mri)) / k1_df$que.nonm
k1_df$que.smrd[is.na(k1_df$que.mrp)] <- 1
k1_df$que.smrd[is.na(k1_df$que.mri)] <- -1

# dealing with missing by filling in bounds for now, obs with values over this have duplicates

## Frame evaluation ##

k1_df$msg.dec <- recode(k1_df$msg.dec, "2" = 1, "1" = 0)

k1_df$eva.poor[k1_df$msg1 == 0] <- as.numeric(k1_df$eva.msg1[k1_df$msg1 == 0])
k1_df$eva.poor[k1_df$msg2 == 0] <- as.numeric(k1_df$eva.msg2[k1_df$msg2 == 0])
k1_df$eva.poor[k1_df$msg3 == 0] <- as.numeric(k1_df$eva.msg3[k1_df$msg3 == 0])

k1_df$eva.ind[k1_df$msg1 == 1] <- as.numeric(k1_df$eva.msg1[k1_df$msg1 == 1])
k1_df$eva.ind[k1_df$msg2 == 1] <- as.numeric(k1_df$eva.msg2[k1_df$msg2 == 1])
k1_df$eva.ind[k1_df$msg3 == 1] <- as.numeric(k1_df$eva.msg3[k1_df$msg3 == 1])

k1_df$eva.com[k1_df$msg1 == 2] <- as.numeric(k1_df$eva.msg1[k1_df$msg1 == 2])
k1_df$eva.com[k1_df$msg2 == 2] <- as.numeric(k1_df$eva.msg2[k1_df$msg2 == 2])
k1_df$eva.com[k1_df$msg3 == 2] <- as.numeric(k1_df$eva.msg3[k1_df$msg3 == 2])

k1_df$eva.vid.poor <- as.numeric(k1_df$eva.rank.vid_8)
k1_df$eva.vid.ind <- as.numeric(k1_df$eva.rank.vid_9)
k1_df$eva.vid.com <- as.numeric(k1_df$eva.rank.vid_10)

k1_df$eva.conf <- as.numeric(k1_df$eva.conf)

k1_df$eva.emp.poor <- as.numeric(k1_df$eva.rank.emp_5)
k1_df$eva.emp.ind <- as.numeric(k1_df$eva.rank.emp_6)
k1_df$eva.emp.com <- as.numeric(k1_df$eva.rank.emp_7)

## Ladder scales ##

k1_df$ses.lad.now <- as.numeric(k1_df$ses.lad.now)
k1_df$ses.lad.y2 <- as.numeric(k1_df$ses.lad.y2)

## Sociodemographics ##

k1_df$soc.age <- as.numeric(as.character(k1_df$soc.age))
k1_df$soc.pri[is.na(k1_df$soc.edu) == FALSE] <- ifelse(as.numeric(k1_df$soc.edu[is.na(k1_df$soc.edu) == FALSE]) >= 4, 1, 0)
k1_df$soc.fem[is.na(k1_df$soc.gen) == FALSE] <- ifelse(k1_df$soc.gen[is.na(k1_df$soc.gen) == FALSE] == "2", 1, 0)
k1_df$soc.chr[is.na(k1_df$soc.rel) == FALSE] <- ifelse(k1_df$soc.rel[is.na(k1_df$soc.rel) == FALSE] == "1", 1, 0)
k1_df$ses.unemp <- ifelse(k1_df$ses.emp == "1", 1, 0)
k1_df$soc.inc <- as.numeric(k1_df$soc.inc)
k1_df$soc.con <- as.numeric(k1_df$soc.con)
k1_df$soc.sav <- as.numeric(k1_df$soc.sav) - 1

## Survey validity ##

k1_df$end.hear <- as.numeric(k1_df$end.hear) - 1

################
## Estimation ##
################

## Plain OLS ##

hypotheses <- c("ind = 0", "com = 1", "ind - com = 0")
depvars <- c("vid.num", "sav.save", "msg.dec")
mechvars <- c("sel.score", "jud.score", "aff.score", "que.smrd", "ses.lad.now", "ses.lad.y2")

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (depvar in depvars) {

        eqn <- paste(depvar, "~ ind  + com", sep = " ")
        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = k1_df$survey.id, hypotheses = c(h), iterations = 100, data = k1_df))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- depvars
    colnames(RES)[6] <- "Min. Q"

    print("----------------------------------------------------------------", quote = FALSE)
    print(paste("H_0:", h), quote = FALSE)
    print(RES, quote = FALSE)

}

for (h in hypotheses) {

    RES <- matrix(nrow = 1, ncol = 5)

    for (depvar in mechvars) {

        eqn <- paste(depvar, "~ ind  + com", sep = " ")
        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = k1_df$survey.id, hypotheses = c(h), iterations = 100, data = k1_df))

    }

    RES <- RES[2:nrow(RES), 1:ncol(RES)]
    RES <- cbind(RES, FDR(RES[, 4]))

    rownames(RES) <- mechvars
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
        RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = k1_df$survey.id, hypotheses = c(h), iterations = 100, data = k1_df))

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

    hypotheses <- c(paste("ind:", hetvar, " = 0", sep = ""), paste(hetvar, ":com", " = 0", sep = ""), paste("ind:", hetvar, " - ", hetvar, ":com", " = 0", sep = ""))

    for (h in hypotheses) {

        RES <- matrix(nrow = 1, ncol = 5)

        for (depvar in depvars) {

            eqn <- paste(depvar, " ~ ind*", hetvar, " + com*", hetvar)
            RES <- rbind(RES, PermTest(eqn, treatvars = c("treatment", "poor", "ind", "com"), clustvars = k1_df$survey.id, hypotheses = c(h), iterations = 100, data = k1_df))

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
