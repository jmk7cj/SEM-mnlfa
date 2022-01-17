#--------------------------------------------------------------------------------------------#
# Author: Joseph Kush (jkush1@jhu.edu) 
# 
# Title: Utilizing Moderated Nonlinear Factor Analysis Models for Integrative Data Analysis
#
# Date: 1/31/2022
#
# Purpose: Master .R file to conduct integrative data analyses using moderated 
#          non-linear factor analysis

#          Step 0: Load packages, set working directory, import data, etc.

#          Step 1: Create Mplus input files for CFAs, estimate models, examine output

#          Step 2: Create Mplus input files for MNLFA model building, estimate models, 
#                  examine output

#          Step 3: Conduct LRT between each item-model and the baseline model

#          Step 4: Remove remaining non-significant parameters, estimate next-to-last and 
#                  final MNLFA model

#          Step 5: Merge estimated factor scores to be used in subsequent analyses
#--------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------#
# Step 0: Load packages, set working directory, import data, etc.
#--------------------------------------------------------------------------------------------#
# Remove working environment, close any open connections
rm(list = ls()); closeAllConnections()

# Load necessary packages
library("parallel")
library("MplusAutomation")
library("MASS")

# Set working directory to folder
setwd("/Users/myfolder") # set own path
myfolder <- getwd()

# Import and view data
data <- read.csv("data.csv")
head(data)

# Demographics
table(data$sex)
table(data$race)
table(data$study_id)

# Reproduce Table 1 (top half)
table_1a <- cbind(table(data$study_id, data$race), 
                 table(data$study_id, data$sex), 
                 table(data$study_id))
colnames(table_1a)[ncol(table_1a)] <- c("total")

table_1a

# Reproduce Table 1 (bottom half)
table_1b <- NULL 
for(i in 1:length(table(data$study_id))) {
  table_1b <- cbind(table_1b, 
  cbind(colMeans(subset(data, 
                        subset = data$study_id == levels(as.factor(data$study_id))[i])[, 5:13], 
                 na.rm=T)))
}
table_1b <- cbind(table_1b, cbind(colMeans(data[, 5:13], na.rm=T)))
colnames(table_1b) <- c(levels(as.factor(data$study_id)), "total")

table_1b

# Store Table 1 as .csv
write.csv(table_1a, "table_1a.csv")
write.csv(table_1b, "table_1b.csv")
#--------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------#
# Step 1: Create Mplus input files for CFAs, estimate models, examine output
#--------------------------------------------------------------------------------------------#
# First, prepare datafile for Mplus (all variables numeric)
data_cfa <- data
data_cfa[is.na(data_cfa)] <- -999 
data_cfa$id <- 1:nrow(data_cfa)
data_cfa$sex <- as.numeric(as.factor(data_cfa$sex))-1
data_cfa$race <- as.numeric(as.factor(data_cfa$race))-1

# Create five binary study dummy indicator variables 
data_cfa$study_id <- as.numeric(as.factor(data_cfa$study_id))
data_cfa$study_1 <- ifelse(data_cfa$study_id == 1, 1, 0) 
data_cfa$study_2 <- ifelse(data_cfa$study_id == 2, 1, 0) 
data_cfa$study_3 <- ifelse(data_cfa$study_id == 3, 1, 0) 
data_cfa$study_4 <- ifelse(data_cfa$study_id == 4, 1, 0) 
data_cfa$study_5 <- ifelse(data_cfa$study_id == 5, 1, 0) 

data_cfa <- data_cfa[, c(1,2,15:19,3:14)]

# View Mplus data
head(data_cfa)

# Create new 'cfa' sub-folder within original working directory (delete if already exists)
unlink(paste(myfolder,"/cfa", sep=""), recursive=T)

# Create new folder
dir.create(paste(myfolder,"/cfa", sep=""))

# Export Mplus data to new folder
setwd(paste(myfolder,"/cfa", sep="")); getwd()
write.table(data_cfa, "data_cfa.dat", row.names=F, col.names=F, quote=F) 

# Determine number of processors to run in parallel
my_processors <- detectCores() - 1; my_processors

# Create Mplus input files, in which CFAs are fit separately for each study
for(i in min(data_cfa$study_id):max(data_cfa$study_id)) {
  input <- paste(
    "title: CFA for study",i,"
    
    data:
    file = data_cfa.dat;
    
    variable: 
    names = id study_id study_1-study_5 sex race x1-x9 hs;
    usevariables = x1-x9;
    categorical = x1-x9;
    useobservations = study_id == ",i,";
    missing = all (-999);
    
    analysis:
    estimator = wlsmv;
    processors = ",my_processors,";
    
    model: 
    Factor BY x1-x9;

    output:
    standardized;
    stdyx;
    ", sep="")
  write.table(input, paste("cfa_study",i,".inp", sep=""), quote=F, row.names=F, col.names=F)
}

# No variability in x3 for Study_ID = 1, so estimate CFA without x3 (but with x1, x2, x4-x9)
for(i in 1) {
  input <- paste(
    "title: CFA for study",i,"
    
    data:
    file = data_cfa.dat;
    
    variable: 
    names = id study_id study_1-study_5 sex race x1-x9 hs;
    usevariables = x1 x2 x4-x9;
    categorical = x1 x2 x4-x9;
    useobservations = study_id == ",i,";
    missing = all (-999);
    
    analysis:
    estimator = wlsmv;
    processors = ",my_processors,";
    
    model: 
    Factor BY x1 x2 x4-x9;
    
    output:
    standardized;
    stdyx;
    ", sep="")
  write.table(input, paste("cfa_study",i,".inp", sep=""), quote=F, row.names=F, col.names=F)
}

# Estimate models
runModels(replaceOutfile="never")

# Mplus input file for a final CFA for all observations (pooled across study)
input <- paste(
  "title: CFA for all observations (pooled across study)
  
  data:
  file = data_cfa.dat;
  
  variable: 
  names = id study_id study_1-study_5 sex race x1-x9 hs;
  usevariables = x1-x9;
  categorical = x1-x9;
  missing = all (-999);
  
  analysis:
  estimator = wlsmv;
  processors =",my_processors,";
  
  model: 
  Factor BY x1-x9;
  
  output:
  standardized;
  stdyx;
  ", sep="")
write.table(input, "cfa_allobs.inp", quote=F, row.names=F, col.names=F)

# Estimate final pooled model
runModels("cfa_allobs.inp", replaceOutfile="never")

# Examine output from CFAs to determine dimensionality
out_cfa_study1 <- readModels("cfa_study1.out")
out_cfa_study2 <- readModels("cfa_study2.out")
out_cfa_study3 <- readModels("cfa_study3.out")
out_cfa_study4 <- readModels("cfa_study4.out")
out_cfa_study5 <- readModels("cfa_study5.out")
out_cfa_allobs <- readModels("cfa_allobs.out")

# Fit statistics across the models
out_cfa_study1$summaries
out_cfa_study1$summaries[c("Parameters", "ChiSqM_DF", "RMSEA_Estimate", 
                           "CFI", "TLI", "SRMR")]
out_cfa_study2$summaries[c("Parameters", "ChiSqM_DF", "RMSEA_Estimate", 
                           "CFI", "TLI", "SRMR")]
# ...etc

# final pooled model with all observations
out_cfa_allobs$summaries[c("Parameters", "ChiSqM_DF", "RMSEA_Estimate", 
                           "CFI", "TLI", "SRMR")]

# standardized estimates
out_cfa_allobs$parameters$stdyx.standardized[1:9, ]
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 2: Create Mplus input files for MNLFA model building, estimate models, examine output
#--------------------------------------------------------------------------------------------#
# First, prepare datafile for Mplus
data_mnlfa <- data_cfa

# Create new 'mnlfa' sub-folder within original working directory (delete if already exists) 
unlink(paste(myfolder,"/mnlfa", sep=""), recursive=T)

# Create new folder
dir.create(paste(myfolder,"/mnlfa", sep=""))

# Export Mplus data to new folder
setwd(paste(myfolder,"/mnlfa", sep="")); getwd()
write.table(data_mnlfa, "data_mnlfa.dat", row.names=F, col.names=F, quote=F) 

# Baseline model allows covariates to moderate the factor mean & factor variance
baseline <- paste(
  "title: Moderation of factor mean and variance
  
  data:
  file = data_mnlfa.dat;
  
  variable: 
  names = id study_id study_1-study_5 sex race x1-x9 hs; 
  usevariables = study_2-study_5 sex race x1-x9; !study_1 is reference group
  categorical = x1-x9; 
  missing = all (-999);
  constraint = study_2-study_5 sex race;
  
  analysis:
  estimator = mlr;
  link = logit;
  processors = ", my_processors, ";
  !estimator = wlsmv; !cannot be used with certain model constraints
  
  model: 
  Factor BY x1-x9;!measurement model
  
  !allow covariates to moderate factor mean (linear function) 
  Factor ON study_2-study_5 sex race;
  [Factor@0]; !constrain factor mean to zero to identify model
  
  !factor variance implicitly set to one to identify model
  Factor(factor_variance); !label for factor variance
  
  model constraint:
  new (f_study_2 f_study_3 
  f_study_4 f_study_5 f_sex f_race); 
  
  !allow covariates to moderate factor variance 
  !(log-linear function to avoid negative values)
  factor_variance = EXP(f_study_2*study_2 + f_study_3*study_3 + 
  f_study_4*study_4 + f_study_5*study_5 + 
  f_sex*sex + f_race*race);
  
  output:
  sampstat; 
  svalues;
  tech1;
  ", sep="")
write.table(baseline, "baseline.inp", quote=F, row.names=F, col.names=F)


# Item-models allows covariates to moderate the factor mean & factor variance, as
# well as the item intercept and item loading (done sequentially for each item)
for(i in 1:9) {
  input <- paste(
    "title: Moderation of factor mean and variance, as well as 
    item intercept and factor loading for item x",i,"
    
    data:
    file = data_mnlfa.dat;
    
    variable: 
    names = id study_id study_1-study_5 sex race x1-x9 hs; 
    
    usevariables = study_2-study_5 sex race x1-x9; !study_1 is reference group
    
    categorical = x1-x9; 
    missing = all (-999);
    constraint = study_2-study_5 sex race;
    
    analysis:
    estimator = mlr;
    link = logit;
    processors = ", my_processors, ";
    !estimator = wlsmv; !cannot be used with certain model constraints
    
    model: 
    Factor BY x1-x9;!measurement model
    
    !allow covariates to moderate factor mean (linear function) 
    Factor ON study_2-study_5 sex race;
    [Factor@0]; !constrain factor mean to zero to identify model
    
    !factor variance implicitly set to one to identify model
    Factor(factor_variance); !label for factor variance
    
    !allow covariates to moderate item i intercept
    x",i," ON study_2-study_5 sex race;
    
    Factor BY x",i," (x",i,"_loading); !label for item i loading
    
    
    model constraint:
    new (f_study_2 f_study_3 
    f_study_4 f_study_5
    f_sex f_race); 
    
    new (x_int x_study_2 x_study_3
    x_study_4 x_study_5
    x_sex x_race);
  
    !allow covariates to moderate factor variance 
    !(log-linear function to avoid negative values)
    factor_variance = EXP(f_study_2*study_2 + f_study_3*study_3 + 
    f_study_4*study_4 + f_study_5*study_5 + 
    f_sex*sex + f_race*race);
    
    !allow covariates to moderate item i loading
    x",i,"_loading = x_int + 
    x_study_2*study_2 + x_study_3*study_3 + 
    x_study_4*study_4 + x_study_5*study_5 + 
    x_sex*sex + x_race*race;
    
    output:
    sampstat; 
    svalues;
    tech1;
    ", sep="")
  write.table(input, paste("x",i,"_model.inp", sep=""), quote=F, row.names=F, col.names=F)
}

# Estimate all of the models (may take a few minutes)
runModels(replaceOutfile="never") 
#--------------------------------------------------------------------------------------------#





#--------------------------------------------------------------------------------------------#
# Step 3: Conduct LRT between each item-model and the baseline model
#--------------------------------------------------------------------------------------------#
modelResults <- readModels(what=c("summaries", "parameters"))

modelLRT <- do.call("rbind", sapply(modelResults,"[", "summaries"))
modelLRT <- modelLRT[, c("Title", "Parameters", "LL", "LLCorrectionFactor")]
modelLRT$pval <- NA

for(i in 2:10) {
  L0 <- modelLRT[1,3]
  c0 <- modelLRT[1,4]
  p0 <- modelLRT[1,2]
  
  L1 <- modelLRT[i,3]
  c1 <- modelLRT[i,4]
  p1 <- modelLRT[i,2]
  
  cd <- ((p0*c0) - (p1*c1)) / (p0-p1)
  lrt <- -2 * (L0 - L1) / cd
  modelLRT[i,"pval"] <- dchisq(x=lrt, df=(p1-p0))
}
modelLRT[, 2:5]
# Results indicate each item-model fits the data significantly 
# better than the baseline model (sig. p-values)


# As a result of each item model producing a significant LRT, 
# item-models 1-9 will keep covariate moderation of an intercept 
# or loading, only if the covariate effect is significant
modelParms <- do.call("rbind", sapply(modelResults,"[", "parameters"))
baselineParms <- data.frame(do.call(cbind,modelParms[[1]]))
for(i in 2:10) {
  assign(paste0("x",i-1,"Parms", sep=""), data.frame(do.call(cbind,modelParms[[i]])))
}

# For the baseline model (factor moderation only):
# Rows 10-12 give moderation of factor mean 
# Rows 24-26 give moderation of factor variance
baselineParms
baselineParms[c(10:15, 27:32), ]

# For each item-model:
# Rows 10-15 give moderation of factor mean 
# Rows 33-38 give moderation of factor variance
# Rows 16-21 give moderation of item intercept (for item i)
# Rows 40-45 give moderation of factor loading (for item i)
x1Parms
x1Parms[c(16:21, 40:45), ] 
#for item 1:
#intercept moderators: none
#loading moderators: study4 & sex

x2Parms[c(16:21, 40:45), ]
#for item 2:
#intercept moderators: study2, study3, study4, & study5
#loading moderators: study3

x3Parms[c(16:21, 40:45), ]
#for item 3:
#intercept moderators: study3, study4, & study5
#loading moderators: study5

x4Parms[c(16:21, 40:45), ]
#for item 4:
#intercept moderators: study3 & study5
#loading moderators: none

x5Parms[c(16:21, 40:45), ]
#for item 5:
#intercept moderators: study3 & sex
#loading moderators:  study3

x6Parms[c(16:21, 40:45), ]
#for item 6:
#intercept moderators: study2 & sex
#loading moderators: study3, study4, study5, & race

x7Parms[c(16:21, 40:45), ]
#for item 7:
#intercept moderators: study3, study4, study5, & sex
#loading moderators: study3, study4, & race 

x8Parms[c(16:21, 40:45), ]
#for item 8:
#intercept moderators: study3, study4, study5, & sex
#loading moderators: study3, study4, & study5

x9Parms[c(16:21, 40:45), ]
#for item 9:
#intercept moderators: study3, study5, & race
#loading moderators: study3, study4, study5, & sex


# Reproduce Table 2
table_2 <- modelLRT[, 2:4]
table_2 <- table_2[rep(1:nrow(table_2), each = 7), ]
table_2[c(seq(1,70, 7)+1, 
          seq(1,70, 7)+2, 
          seq(1,70, 7)+3,
          seq(1,70, 7)+4,
          seq(1,70, 7)+5,
          seq(1,70, 7)+6), ] <- c("") 
table_2[rep(c("est", "p"), times=4)] <- c("")


# Add factor parameter estimates for each model to Table 2
table_2[c(2:7), c(4:5)] <- modelResults$baseline.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(2:7), c(6:7)] <- modelResults$baseline.out$parameters$unstandardized[c(27:32), c("est", "pval")]

table_2[c(9:14), c(4:5)] <- modelResults$x1_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(9:14), c(6:7)] <- modelResults$x1_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(16:21), c(4:5)] <- modelResults$x2_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(16:21), c(6:7)] <- modelResults$x2_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(23:28), c(4:5)] <- modelResults$x3_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(23:28), c(6:7)] <- modelResults$x3_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(30:35), c(4:5)] <- modelResults$x4_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(30:35), c(6:7)] <- modelResults$x4_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(37:42), c(4:5)] <- modelResults$x5_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(37:42), c(6:7)] <- modelResults$x5_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(44:49), c(4:5)] <- modelResults$x6_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(44:49), c(6:7)] <- modelResults$x6_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(51:56), c(4:5)] <- modelResults$x7_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(51:56), c(6:7)] <- modelResults$x7_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(58:63), c(4:5)] <- modelResults$x8_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(58:63), c(6:7)] <- modelResults$x8_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]

table_2[c(65:70), c(4:5)] <- modelResults$x9_model.out$parameters$unstandardized[c(10:15), c("est", "pval")]
table_2[c(65:70), c(6:7)] <- modelResults$x9_model.out$parameters$unstandardized[c(33:38), c("est", "pval")]


# Add item parameter estimates for item models to Table 2
table_2[c(9:14), c(8:9)] <- modelResults$x1_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(9:14), c(10:11)] <- modelResults$x1_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(16:21), c(8:9)] <- modelResults$x2_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(16:21), c(10:11)] <- modelResults$x2_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(23:28), c(8:9)] <- modelResults$x3_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(23:28), c(10:11)] <- modelResults$x3_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(30:35), c(8:9)] <- modelResults$x4_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(30:35), c(10:11)] <- modelResults$x4_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(37:42), c(8:9)] <- modelResults$x5_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(37:42), c(10:11)] <- modelResults$x5_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(44:49), c(8:9)] <- modelResults$x6_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(44:49), c(10:11)] <- modelResults$x6_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(51:56), c(8:9)] <- modelResults$x7_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(51:56), c(10:11)] <- modelResults$x7_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(58:63), c(8:9)] <- modelResults$x8_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(58:63), c(10:11)] <- modelResults$x8_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2[c(65:70), c(8:9)] <- modelResults$x9_model.out$parameters$unstandardized[c(16:21), c("est", "pval")]
table_2[c(65:70), c(10:11)] <- modelResults$x9_model.out$parameters$unstandardized[c(40:45), c("est", "pval")]

table_2

# Store Table 2 as .csv
write.csv(table_2, paste(myfolder,"/table_2.csv", sep=""))
#--------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------#
# Step 4: Remove remaining non-significant parameters, estimate next-to-last 
#         and final MNLFA model
#--------------------------------------------------------------------------------------------#
# Based on information above, construct the next-to-last model,
# keeping only significant moderators of item intercepts and 
# loadings (but always keeping moderation of factor mean and 
# variance, regardless of significance)
next_to_last_model <- paste(
  "title: next-to-last MNLFA model
  
  data:
  file = data_mnlfa.dat;
  
  variable: 
  names = id study_id study_1-study_5
  sex race x1-x9 hs;
  
  usevariables = study_2-study_5 sex race x1-x9; !study_1 is reference group
  
  categorical = x1-x9;
  missing = all (-999);
  constraint = study_2-study_5 sex race;
  
  analysis:
  estimator = mlr;
  link = logit;
  processors = ", my_processors, ";
  !estimator = wlsmv; !cannot be used with certain model constraints
  
  model:
  factor BY x1*1 (x1_loading);
  factor BY x2*1 (x2_loading);
  factor BY x3*1 (x3_loading);
  factor BY x4*1 (x4_loading);
  factor BY x5*1 (x5_loading);
  factor BY x6*1 (x6_loading);
  factor BY x7*1 (x7_loading);
  factor BY x8*1 (x8_loading);
  factor BY x9*1 (x9_loading);
  
  !allow covariates to moderate factor mean (linear function) 
  Factor ON study_2-study_5 sex race;
  [Factor@0]; !constrain factor mean to zero to identify model
  
  !factor variance implicitly set to one to identify model
  Factor(factor_variance); !label for factor variance
  
  ! Moderation of item intercepts (previously determined)
  ! no moderation of item x1 intercept 
  x2 ON study_2-study_5;
  x3 ON study_3-study_5; 
  x4 ON study_3 study_5;
  x5 ON study_3 sex;
  x6 ON study_2 sex;
  x7 ON study_3-study_5 sex;
  x8 ON study_3-study_5 sex;
  x9 ON study_3 study_5 race;
  
  model constraint:
  NEW(f_study_2 f_study_3
  f_study_4 f_study_5
  f_sex f_race);
  
  !intercepts for loading moderation equation 
  !(no sig. moderators for x4)
  NEW(int1 int2 int3 int5 int6 int7 int8 int9);  
  
  NEW(x1_study_4 x1_sex);
  NEW(x2_study_3);
  NEW(x3_study_5);
  ! no loading moderation of x4
  NEW(x5_study_3);
  NEW(x6_study_3 x6_study_4 x6_study_5 x6_race);
  NEW(x7_study_3 x7_study_4 x7_race);
  NEW(x8_study_3 x8_study_4 x8_study_5);
  NEW(x9_study_3 x9_study_4 x9_study_5 x9_sex);
  
  !allow covariates to moderate factor variance
  !(log-linear function to avoid negative values)
  factor_variance = EXP(f_study_2*study_2 + f_study_3*study_3 + 
  f_study_4*study_4 + f_study_5*study_5 + 
  f_sex*sex + f_race*race);
  
  !allow covariates to moderate item loadings
  x1_loading = int1 + x1_study_4*study_4 + x1_sex*sex;
  
  x2_loading = int2 + x2_study_3*study_3;
  
  x3_loading = int3 + x3_study_5*study_5;
  
  ! no loading moderation of x4
  
  x5_loading = int5 + x5_study_3*study_3;
  
  x6_loading = int6 + x6_study_3*study_3 + x6_study_4*study_4 + 
  x6_study_5*study_5 + x6_race*race;
  
  x7_loading = int7 + x7_study_3*study_3 + x7_study_4*study_4 + 
  x7_race*race;
  
  x8_loading = int8 + x8_study_3*study_3 + x8_study_4*study_4 + 
  x8_study_5*study_5;
  
  x9_loading = int9 + x9_study_3*study_3 + x9_study_4*study_4 + 
  x9_study_5*study_5 + x9_sex*sex; 
  
  output:
  sampstat;
  svalues;
  tech1;
  
  ", sep="")
write.table(next_to_last_model, "next_to_last_model.inp", quote=F, row.names=F, col.names=F)

# Estimate next-to-last model 
# Note: this model may take a long time to run (upwards of 1hr+)
runModels("next_to_last_model.inp", replaceOutfile="never") 

# Read in the output of the next-to-last model 
out_next_to_last_model <- readModels("next_to_last_model.out") 


# Remove any non-significant moderation of item parameters
# (leave moderation of factor parameters regardless of significance)
out_next_to_last_model$parameters$unstandardized


# All covariate moderation of the factor mean stays (regardless of significance)
out_next_to_last_model$parameters$unstandardized[10:15, ]

# All covariate moderation of the factor variance stays (regardless of significance)
out_next_to_last_model$parameters$unstandardized[51:56, ]

# Only significant covariate moderation of item intercepts stay
out_next_to_last_model$parameters$unstandardized[16:39, ]
# drop: 
# x2*study_2, x2*study_5
# x3*study_4, x3*study_5
# x6*sex
# x7*study_4
# x8*study_3
# x9*study_3

# Only significant covariate moderation of item loadings stay
out_next_to_last_model$parameters$unstandardized[65:83, ]
# drop: 
# x3*study_5
# x6*study_3, x6*study_4
# x7*study_4, x7*race
# x8*study_3, x8_study_4, x8*study_5
# x9*study_3, x9*study_4, x9*study_5, x9*sex


# NOTE: This last pruning effort will result in the final MNLFA model


# Build final MNLFA model 
final_MNLFA_model <- paste(
  "title: final MNLFA model
  
  data:
  file = data_mnlfa.dat;
  
  variable: 
  names = id study_id study_1-study_5
  sex race x1-x9 hs;
  
  idvariable = id;

  usevariables = study_2-study_5 sex race x1-x9; !study_1 is reference group
  
  categorical = x1-x9;
  missing = all (-999);
  constraint = study_2-study_5 sex race;
  
  analysis:
  estimator = mlr;
  link = logit;
  processors = ", my_processors, ";
  !estimator = wlsmv; !cannot be used with certain model constraints
  
  model:
  factor BY x1*1 (x1_loading);
  factor BY x2*1 (x2_loading);
  factor BY x3*1 (x3_loading);
  factor BY x4*1 (x4_loading);
  factor BY x5*1 (x5_loading);
  factor BY x6*1 (x6_loading);
  factor BY x7*1 (x7_loading);
  factor BY x8*1 (x8_loading);
  factor BY x9*1 (x9_loading);
  
  !allow covariates to moderate factor mean (linear function) 
  Factor ON study_2-study_5 sex race;
  [Factor@0]; !constrain factor mean to zero to identify model
  
  !factor variance implicitly set to one to identify model
  Factor(factor_variance); !label for factor variance
  
  ! Moderation of item intercepts (previously determined)
  ! no moderation of item x1 intercept 
  x2 ON study_3 study_4;
  x3 ON study_3; 
  x4 ON study_3 study_5;
  x5 ON study_3 sex;
  x6 ON study_2;
  x7 ON study_3 study_5 sex;
  x8 ON study_4 study_5 sex;
  x9 ON study_5 race;
  
  model constraint:
  NEW(f_study_2 f_study_3 
  f_study_4 f_study_5
  f_sex f_race);
  
  !intercepts for loading moderation equation 
  !(no sig. moderators for x4)
  NEW(int1 int2 int5 int6 int7);  
  
  NEW(x1_study_4 x1_sex);
  NEW(x2_study_3);
  ! no loading moderation of x3
  ! no loading moderation of x4
  NEW(x5_study_3);
  NEW(x6_study_5 x6_race);
  NEW(x7_study_3);
  ! no loading moderation of x8
  ! no loading moderation of x9
  
  !allow covariates to moderate factor variance
  !(log-linear function to avoid negative values)
  factor_variance = EXP(f_study_2*study_2 + f_study_3*study_3 + 
  f_study_4*study_4 + f_study_5*study_5 + 
  f_sex*sex + f_race*race);
  
  !allow covariates to moderate item loadings
  x1_loading = int1 + x1_study_4*study_4 + x1_sex*sex;
  x2_loading = int2 + x2_study_3*study_3;
  ! no loading moderation of x3
  ! no loading moderation of x4
  x5_loading = int5 + x5_study_3*study_3;
  x6_loading = int6 + x6_study_5*study_5 + x6_race*race;
  x7_loading = int7 + x7_study_3*study_3;
  ! no loading moderation of x8
  ! no loading moderation of x9
  
  output:
  sampstat;
  svalues;
  tech1;
  
  savedata:
  save = fscores; !save estimated factor scores
  file = est_factor_scores.csv;
  
  ", sep="")
write.table(final_MNLFA_model, "final_MNLFA_model.inp", quote=F, row.names=F, col.names=F)

# Estimate final MNLFA model (much faster than penultimate model)
runModels("final_MNLFA_model.inp", replaceOutfile="never") 
          
# Read in the output of the final MNLFA model
out_final_MNLFA_model <- readModels("final_MNLFA_model.out") 
out_final_MNLFA_model$parameters$unstandardized


# Reproduce Table 3
table_3 <- out_final_MNLFA_model$parameters$unstandardized[c(10:15, 43:48), c("est", "se", "pval")]

table_3

# Store Table 3 as .csv
write.csv(table_3, paste(myfolder,"/table_3.csv", sep=""))


# Reproduce Table 4
table_4 <- data.frame(matrix(NA, nrow=37, ncol=6))

# Fill in item intercepts for Table 4
table_4[c(1, 5, 9, 12, 16, 20, 25, 30, 35), c(1:3)] <- 
  out_final_MNLFA_model$parameters$unstandardized[c(33:41), c("est", "se", "pval")]

# Fill in item loading for Table 4
table_4[c(1, 5, 9, 12, 16, 20, 25, 30, 35), c(4:6)] <- 
  out_final_MNLFA_model$parameters$unstandardized[c(49:50, 3:4, 51:53, 8:9), c("est", "se", "pval")]

# Fill in item intercept moderation for certain covariates
table_4[c(6:7, 10, 13:14, 17:18, 21, 26:28, 31:33, 36:37), c(1:3)] <- 
  out_final_MNLFA_model$parameters$unstandardized[c(16:31), c("est", "se", "pval")]

# Fill in item loading moderation for certain covariates
table_4[c(2:3, 6, 17, 22:23, 26), c(4:6)] <- 
  out_final_MNLFA_model$parameters$unstandardized[c(54:60), c("est", "se", "pval")]

table_4[is.na(table_4)] <- ""

table_4

# Store Table 4 as .csv
write.csv(table_4, paste(myfolder,"/table_4.csv", sep=""))
#--------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------#
# Step 5: Merge estimated factor scores to be used in subsequent analyses
#--------------------------------------------------------------------------------------------#
est_factor_scores <- read.table("est_factor_scores.csv")
colnames(est_factor_scores) <- tolower(out_final_MNLFA_model$savedata_info$fileVarNames)
head(est_factor_scores)

# Keep just id and estimated factor score
est_factor_scores <- est_factor_scores[,c("id", "factor")]

# Merge estimated factor scores in with original data
merged_data <- merge(data_cfa, est_factor_scores, by=c("id"))
merged_data[merged_data == -999] <- NA

# View merged data with factor scores
head(merged_data)

# Estimate the effect of aggressive behavior factor scores on 
# probability of high school graduation (binary outcome) using 
# a logistic regression
logit_model <- glm(hs ~ study_2 + study_3 + study_4 + study_5 +
                     sex + race + factor, family = binomial(link = "logit"), merged_data)
summary(logit_model)

# Odds ratio estimate of factor score
exp(cbind("Odds ratio" = coef(logit_model), confint.default(logit_model, level = 0.95)))

# Reproduce Table 5
table_5 <- round(cbind(exp(cbind("Odds ratio" = coef(logit_model), 
                                 confint.default(logit_model, level = 0.95))), 
                       summary(logit_model)$coef[,4]), digits = 3)

# Store Table 5 as .csv
write.csv(table_5, paste(myfolder,"/table_5.csv", sep=""))
#--------------------------------------------------------------------------------------------#
