# Cox regression modelling and survival analysis
# 


# Reading in data ---------------------------------------------------------


##At Home
setwd("~/NGS/PDAC.CTC")

pdac <- data.table::fread("pdac.study.20170310.txt")
#pdac <- data.table::fread("pdac.study.20170308.txt")
#151 patients

spaceless <- function(x) {colnames(x) <- gsub(" ", "_", colnames(x));x}
pdac <- spaceless(pdac)
#pdac$tenplusCTC <- cut(pdac$CTC_Count, breaks = c(-Inf, 9.9, Inf), labels = c("less than 10 CTCs", "10+ CTCs"), right = F)
#pdac$threeplusCTC <- cut(pdac$CTC_Count, breaks = c(-Inf, 2.9, Inf), labels = c(0, 1), right = F)
#pdac$CTC.strat <- cut(pdac$CTC_Count, breaks = c(-Inf, 0.9, 2.1, 5.1, 9.1, Inf), labels = c("0 CTCs", "1-2 CTCs", "3-5 CTCs", "6-9 CTCs", "10+ CTCs"), right = F)
#pdac$optCA199 <- cut(pdac$CA19_9_At_Draw, breaks = c(-Inf, 175, Inf), labels = c(0, 2), right = F)
#pdac$OS_mos <- pdac$OS_Days / 30.4167
#table(pdac$OS_mos)
#pdac$Recurrence <- as.numeric(pdac$Recurrence)
#pdac$Stage <- pdac$Preop_Stage
#pdac$Stage[pdac$Stage == 0] <- NA
#table(pdac$Stage)
#Sex as numbers
#pdac$Sex_num <- pdac$Sex
#pdac$Sex_num[pdac$Sex_num == "M"] <- 1
#pdac$Sex_num[pdac$Sex_num == "F"] <- 0
#table(pdac$Sex_num)
#Just the study patients
pdac.study <- filter(pdac, study_patients == 1)
#127 patients
table(pdac.study$Preop_Stage)
#Just the PDAC patients
pdac.only <- filter(pdac.study, PDAC == 1)
#100 patients
#Stage for pdac pts only
pdac.only$Stage <- pdac.only$Preop_Stage
pdac.only$Stage[pdac.only$Stage == 0] <- NA
table(pdac.only$Stage)



# Cox regression univariate exploratory -----------------------------------

# * * * * All PDAC Patients (100) ----------------------------------------------------
colnames(pdac)
covariates <- c("Age","Sex","Preop_Stage", "Stage_preop_for_recur", "Downstaged", "Metastatic", "Mesenteric.vv.inv", "CTC_Count", "CEA_CTC_present", "CTCs_Present", "CA19_9_At_Draw", "Rad_Tumor_Size", "Tumor_Size", "Tumor_Grade", "Poor_Grade", "Margins", "PNI", "LVI", "Node_Status", "oneplusCTC", "log_CTC", "log_CA199", "log_tumor_size", "tenplusCTC", "threeplusCTC")

univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS_Days, Mortality) ~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = pdac.only)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         beta<-signif(x$coef[1], digits=3);#coeficient beta
                         HR <-signif(x$coef[2], digits=3);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#                           beta HR (95% CI for HR) wald.test  p.value
Age                      0.0338   1.03 (1.01-1.06)      6.35   0.0118
Sex                      -0.119 0.888 (0.538-1.47)      0.22    0.642
Preop_Stage               0.478    1.61 (1.2-2.17)      10.1  0.00148
Stage_preop_for_recur     0.519   1.68 (1.27-2.22)      13.4 0.000252
Downstaged              -0.0207 0.979 (0.508-1.89)         0    0.951
Mesenteric_or_mets        0.427  1.53 (0.913-2.58)      2.61    0.106
Mesenteric_involvement   -0.445 0.641 (0.358-1.15)      2.24    0.134
Tumor_location_head_1    -0.391 0.676 (0.303-1.51)      0.91     0.34
Metastatic                 1.05   2.86 (1.67-4.88)      14.7 0.000125
CTC_Count                0.0303   1.03 (1.01-1.05)      7.79  0.00524
CEA_CTC_present          0.0887   1.09 (0.662-1.8)      0.12    0.729
CTCs_Present              0.982   2.67 (1.27-5.62)      6.69   0.0097
CA19_9_At_Draw         1.78e-05            1 (1-1)      14.4 0.000151
Rad_Tumor_Size          0.00364     1 (0.861-1.17)         0    0.963
Tumor_Size               0.0243  1.02 (0.869-1.21)      0.08    0.773
Tumor_Grade               0.582  1.79 (0.946-3.38)      3.21   0.0734
Poor_Grade               -0.839 0.432 (0.19-0.983)         4   0.0455
Margins                   0.395  1.48 (0.626-3.53)       0.8     0.37
PNI                        0.99   2.69 (0.36-20.1)      0.93    0.334
LVI                       0.587   1.8 (0.757-4.27)      1.77    0.184
Node_Status                1.06   2.9 (0.973-8.62)      3.65    0.056
oneplusCTC                0.982   2.67 (1.27-5.62)      6.69   0.0097
log_CTC                   0.527   1.69 (1.28-2.25)      13.5 0.000242
log_CA199                0.0838  1.09 (0.973-1.22)      2.19    0.139
log_tumor_size          -0.0187 0.981 (0.466-2.07)         0    0.961



# Fine and Gray -----------------------------------------------------------

##############################################################################
# Adds-on functions to crr() function in 'cmprsk' package by Gray, RJ. 
# Written by Luca Scrucca
#
# Reference: 
# Scrucca L, Santucci A, Aversa F (2009) Regression Modeling of Competing 
#   Risk Using R: An In Depth Guide for Clinicians. Submitted to Bone Marrow
#   Transplantation
##############################################################################

# This ensure that the package is loaded
if(!require(cmprsk))
{ stop("the package 'cmprsk' is required, please install it. \nSee help(install.packages).") }

factor2ind <- function(x, baseline)
{
  # Given a factor variable x, create an indicator matrix of dimension 
  # length(x) x (nlevels(x)-1) dropping the column corresponding to the 
  # baseline level (by default the first level is used as baseline).
  # Example:
  # > x = gl(4, 2, labels = LETTERS[1:4])
  # > factor2ind(x)
  # > factor2ind(x, "C")
  xname <- deparse(substitute(x))
  n <- length(x)
  x <- as.factor(x)
  if(!missing(baseline)) x <- relevel(x, baseline)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x)-1)] <- 1
  X[is.na(x),] <- NA
  dimnames(X) <- list(names(x), paste(xname, levels(x), sep = ":"))
  return(X[,-1,drop=FALSE])
}

modsel.crr <- function (object, ..., d = log(object$n)) 
{
  if(class(object) != "crr") 
    stop("object is not of class 'crr'")
  objects <- list(object, ...)
  nmodels <- length(objects)
  modnames <- paste("Model ", format(1:nmodels), ": ", 
                    lapply(objects, function(x) x$call), 
                    sep = "", collapse = "\n")
  # add null model
  mod0 <- object
  mod0$loglik <- mod0$loglik.null
  mod0$coef <- mod0$call$cov1 <- mod0$call$cov2 <- NULL
  objects <- c(list(mod0), objects)
  nmodels <- nmodels + 1
  #
  modnames <- c("Model 0: Null model", modnames)
  ns <- sapply(objects, function(x) x$n) 
  dfs <- sapply(objects, function(x) length(x$coef)) 
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same dataset")
  out <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  loglik <- sapply(objects, function(x) x$loglik)
  crit <- sapply(objects, function(x) -2*x$loglik + d*length(x$coef))
  out[,1] <- ns
  out[,2] <- loglik
  out[,3] <- dfs
  out[,4] <- crit
  out[,5] <- crit - min(crit)
  if(d==log(object$n)) critname <- "BIC"
  else if(d == 2) critname <- "AIC"
  else critname <- "Criterion"
  colnames(out) <- c("Num.obs", "logLik", "Df.fit", critname, paste(critname, "diff"))
  rownames(out) <- 0:(nmodels-1)
  title <- "Model selection table\n"
  topnote <- modnames
  structure(as.data.frame(out), heading = c(title, topnote), 
            class = c("anova", "data.frame"))
}

mod.os1 <- crr(pdac.only$Ftime.os, pdac.only$Status.os, ctc2[,1])
modsel.crr(mod.os1, mod.os2, mod.os3, mod.os4, mod.os5, mod.os6, mod.os7, mod.os8, mod.os9, mod.os10, mod.os11, mod.os12)




attach(resected)
covs <- cbind.data.frame(Age, factor2ind(Sex, "F"), factor2ind(Stage, "1"), factor2ind(Downstaged, "0"), CTC_Count, log_CTC, log_CA199, CA19_9_At_Draw, log_tumor_size, Rad_Tumor_Size, factor2ind(CTCs_Present, "No"), factor2ind(threeplusCTC, "0"))
crr_models <- lapply(covs, function(x){crr(resected$RFS_Days, resected$Recurrence_or_DOC, x)})
crr_results <- lapply(crr_models, 
                      function(x){
                        x <- summary.crr(x, conf.int = 0.90)
                        p.value <- signif(x$coef[,"p-value"], digits = 20)
                        beta <- signif(x$coef[1], digits = 3)
                        HR <- signif(x$coef[2], digits = 3)
                        HR.confint.lower <- signif(x$conf.int[,"5%"], 3)
                        HR.confint.upper <- signif(x$conf.int[,"95%"], 3)
                        HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                        number <- x$n
                        res <- c(beta, HR, number, p.value)
                        names(res) <- c("beta", "HR (95% CI for HR)", "num.cases", "p.value")
                        return(res)
                      })
res <- t(as.data.frame(crr_results, check.names = F))
as.data.frame(res)

#                      beta HR (95% CI for HR) num.cases p.value
Age                0.0307  1.03 (0.998-1.06)        40    0.12
Sex:M            -0.00641 0.994 (0.526-1.88)        40    0.99
Stage:2            -0.486 0.615 (0.321-1.18)        40    0.22
Stage:3              1.21   3.34 (1.79-6.23)        40  0.0014
Downstaged:1         1.21   3.34 (1.79-6.23)        40  0.0014
CTC_Count            0.27   1.31 (1.16-1.48)        40 0.00025
log_CTC             0.974   2.65 (1.64-4.28)        40 0.00086
log_CA199          0.0813  1.08 (0.927-1.27)        37     0.4
CA19_9_At_Draw   0.000374            1 (1-1)        37    0.17
log_tumor_size      0.523  1.69 (0.703-4.05)        40    0.33
Rad_Tumor_Size     0.0665   1.07 (0.953-1.2)        40    0.34
CTCs_Present:Yes     1.24    3.47 (1.67-7.2)        40  0.0051
threeplusCTC:1      0.991   2.69 (1.55-4.69)        40  0.0032
