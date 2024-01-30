
# File name: lengthClustTest.R
# Author: Vanessa Tobias <vanessa_tobias@fws.gov>
# Date: 2021-08-31

# This script builds on lengthClust.R.
# The purpose is to create training and testing datasets to check that
#  the method performs consistently with data that are not used
#  to train the classifier. We don't have "true" ages for any of the fish
#  to test the classifier, but if the method gives similar results with 
#  different training datasets, that supports the trustworthiness of the 
#  classifications that it makes.

# This script was developed to support the SSA for the Bay-Delta DPS of 
# Longfin Smelt. Documentation of methods, results, and interpretation can be
# found in Appendix C of that report.

#### 0. Setup ####

#### Load libraries ####
library(lme4)
library(mgcv)
library(RColorBrewer)
library(mclust)
library(beepr)
library(dplyr)

#### 1. Functions ####
set.seed(124)

classByKD <- function(dat){
  # dat = a data frame. Do the subsetting outside this function
  #      needs to have columns for 
  #       - Length (fork length of the fish)
  
  tmp <- dat
  
  #### Kernel Density ####
  # Catch month and year combinations with zero records
  if (length(na.omit(tmp$Length)) %in% c(0, 1)) {
    groupCount <- NA
    centers <- data.frame("middle" = NA,
                          "age" = NA)
    month <- NA
    # ageClust[which(ageClust$Year == i & ageClust$month == j), 
    #          c("groupCount", "centers", "cutoffs")] <- c(NA, NA, NA)
    
    Age0MaxLength <- NA
    Age1MaxLength <- NA
    Age2MaxLength <- NA
    
    maxMin <- NA
    tpx <- NA
    
  } else{
    
    # Calculate the kernel density
    den <- density(na.omit(tmp$Length), bw = 6.5)
    # Approximate the slope of the kernel density plot
    der1 <- den$y-lead(den$y)
    # find local maxima & minima (turning points) using first derivative
    #  -- Maxes = Center of FL for groups
    #  -- Mins = Cutoff points between groups
    tpx <- den$x[which(sign(der1) != sign(lead(der1)))] 
    tpy <- den$y[which(sign(der1) != sign(lead(der1)))]
    # ID maxima and minima by sign. first pos = maximum starts
    tpSign <- sign(tpy-lead(tpy))
    # find the NAs (last in vector), assign them the opposite sign of the previous element
    #  -- if there's only 1 it has to be a maximum. This fixes the error with length 0 replacement
    tpSign[which(is.na(tpSign))] <- if(length(tpSign)==1) 1 else 0-tpSign[which(is.na(tpSign))-1]
    maxMin <- rep("Max", length(tpSign))
    maxMin[which(tpSign == -1)] <- "Min"
    
    # Use the SFBS cutoffs as rough guides to classify the centers of the high density areas
    #  -- the centers should be far enough away from the SFBS cutoffs to make this clean (I hope)
    centers <- data.frame("middle" = tpx[which(maxMin == "Max")],
                          "age" = 0)
    centers$age[which(centers$middle > lfAgeLen$Age0MaxLength[j])] <- 1
    centers$age[which(centers$middle > lfAgeLen$Age1MaxLength[j])] <- 2
    
    # Classify measurements into age groups based on distance to centers
    tmp$DensAge <- NA
    tmp$DensAge[which(!is.na(tmp$Length))] <- if(length(centers$middle) == 1) {
      rep(centers$age[1], length(tmp$Length[which(!is.na(tmp$Length))]))
    } else {
      unlist(lapply(tmp$Length, FUN = function(x) centers$age[which.min(abs(x - centers$middle))]))
    }
    
    
    ##### Classify measurements by cutoffs ####
    if(length(tpx[which(maxMin == "Min")]) ==0) {
      cutoffs <- data.frame("cutoffs" = NA,
                            "age" = NA)
      tmp$CutAge <- tmp$DensAge
      
      Age0MaxLength <- NA
      Age1MaxLength <- NA
      Age2MaxLength <- NA
    } else {
      cutoffs <- data.frame("cutoffs" = tpx[which(maxMin == "Min")],
                            "age" = 0)
      # use SBFS cutoffs as a "prior"/best guess
      #   calculate the closest SFBS thresholds to label the cutoffs
      
      tmp$CutAge <- NA
      tmp$CutAge[which(!is.na(tmp$Length))] <- 0
      for(k in 1:length(cutoffs$cutoffs)){
        tmp$CutAge[which(tmp$Length > cutoffs$cutoffs[k])] <- centers$age[k+1]
        # save cutoffs to ageClust dataframe in correct columns.
        cutoffs$age[k] <- names(which.min(abs(cutoffs$cutoffs[k] - priors[j, c("Age0MaxLength", "Age1MaxLength", "Age2MaxLength")])))
      }
      Age0MaxLength <- NA
      Age1MaxLength <- NA
      Age2MaxLength <- NA
      
      Age0MaxLength <- if(length(which((cutoffs$age == "Age0MaxLength")))>0) max(cutoffs$cutoffs[which(cutoffs$age == "Age0MaxLength")]) else NA
      Age1MaxLength <- if(length(which((cutoffs$age == "Age1MaxLength")))>0) max(cutoffs$cutoffs[which(cutoffs$age == "Age1MaxLength")]) else NA
      Age2MaxLength <- if(length(which((cutoffs$age == "Age2MaxLength")))>0) max(cutoffs$cutoffs[which(cutoffs$age == "Age2MaxLength")]) else NA
    }
  }
  return(list("dat" = tmp,
              "clusters" = data.frame("Year" = if(length(tmp$Year) == 0) NA else unique(tmp$Year),
                                      "month" = if(length(tmp$month) == 0) NA else unique(tmp$month),
                                      "groupCount" = sum(maxMin == "Max"), 
                                      "centers" = as.character(list(if(length(tpx[which(maxMin == "Max")]) == 0) {NA} else tpx[which(maxMin == "Max")])), 
                                      "cutoffs" = as.character(list(if(length(tpx[which(maxMin == "Min")]) == 0) {NA} else tpx[which(maxMin == "Min")])),
                                      "Age0MaxLength" = Age0MaxLength,
                                      "Age1MaxLength" = Age1MaxLength,
                                      "Age2MaxLength" = Age2MaxLength)))
  
}




classBymclust <- function(KDlist){
  #KDList comes from running a dataset through classByKD
  
  tmp <- KDlist$dat
  clusters <- KDlist$clusters
  #### Use mclust to find clusters ####
  
  # run mclust
  # Makes density calcs using a Gaussian finite mixture model from Mclust
  # Still making the centers using Mclust and getting unrealistic results.
  # Make the density calcs and apply logic from KD method?
  if(length(na.omit(tmp$Length)) == 0){
    groupCount.mclust <- 0
    # tmp$mclustAge <- NA
    clusters$groupCount.mclust <-NA
    clusters$centers.mclust <- NA
  } else {
    clust <- densityMclust(na.omit(tmp$Length)) 
    groupCount.mclust <- clust$G
    
    # label clusters as ages, using info from SFBS as a guide
    centers.mclust <- data.frame("middle" = clust$parameters$mean,
                                 "age" = 0,
                                 "group" = 1:clust$G)
    centers.mclust$age[which(centers.mclust$middle > lfAgeLen$Age0MaxLength[j])] <- 1
    centers.mclust$age[which(centers.mclust$middle > lfAgeLen$Age1MaxLength[j])] <- 2
    
    # assign ages to observations
    clustClass <- data.frame(classification = clust$classification)
    tmp$mclustAge <- merge(clustClass, centers.mclust, 
                           by.x = "classification", by.y = "group",
                           all.x = TRUE)$age
    
    clusters$groupCount.mclust <- as.numeric(groupCount.mclust) 
    clusters$centers.mclust <- as.character(list(clust$parameters$mean))
  }
  
  
  return(list("dat" = tmp,
              "clusters" = clusters))
}      


#### 2. Read in the data ####
allDatWY <- read.csv(".\\Data_Original\\allDatLengthsWY_20210830.csv")
allDatWY <- allDatWY[which(allDatWY$Year > 1979),]
allDatWY <- allDatWY[which(allDatWY$Year < 2019),]
allDatWY$Cohort.f <- as.factor(allDatWY$Cohort.f)

# Bay Study length classification table
# use SFBS as a starting place, but make a few changes to account for
# information in other studies. (E.g., SFBS consisent floor at 40mm is problematic
# in the classification code; some potential age-3 fish are also problematic)
# Bay Study length classes
lfAgeLen <- read.csv(".\\Data_Original\\BayStudyAgeLength.csv",
                     stringsAsFactors = FALSE)
priors <- lfAgeLen
priors$Age2MaxLength <- lfAgeLen$Age1MaxLength+45
# added value above is somewhat arbitrary. Idea is to separate double age-1 cutoffs.



# K-fold cross-validation ####
# - train model on data not in fold t
# - test the model using data in fold t
#  --> names of results follow the fold that was left out
#### Split the data ####
k <- 5
allDatWY$testGroup <- sample(1:k, length(allDatWY$Year), replace = TRUE)


#### Run the classifier on training datasets ####


allDatWY$DensAge <- NA
allDatWY$CutAge <- NA
allDatWY$mclustAge <- NA

#### Create dataset to hold clusters ####

Year = rep(unique(allDatWY$Year), each = 12) #allDat2$Year
month = rep(1:12, length(unique(allDatWY$Year)))
gamDat <- data.frame("cohort" = c(Year, Year - 1, Year - 2),
                     "calYear" = rep(Year, 3),
                     "month" = c(month, month + 12, month + 24),
                     "calMonth" = rep(month, 3))

gamDat$cohort.f <- factor(gamDat$cohort, levels = 1980:2020)


for(t in 1:k){
assign(paste0("ageClustTmp"), data.frame(Year = rep(unique(allDatWY$Year), each = 12), 
                                         month = rep(1:12, length(unique(allDatWY$Year))), 
                                         groupCount = NA,
                                         centers = NA,
                                         cutoffs = NA,
                                         Age0MaxLength = NA,
                                         Age1MaxLength = NA,
                                         Age2MaxLength = NA))

for(y in unique(allDatWY$Year)){
  for(m in 1:12){
    j <- m #workaround for the function. fix this later.
    
    assign("test", 
           classByKD(dat = allDatWY[which(allDatWY$testGroup == t & 
                                            allDatWY$Year == y &
                                            allDatWY$month == m),]))
  
  ageClustTmp[which(ageClustTmp$Year == y & ageClustTmp$month == m), 
            c("groupCount", "centers", "cutoffs",
              "Age0MaxLength", "Age1MaxLength", "Age2MaxLength")] <- test$clusters[,3:8]
  }
}
  assign(paste0("ageClust", t, "Long"), 
         with(ageClustTmp,{
           data.frame("cohort" = c(Year, Year - 1, Year - 2),
                      "cohort.f" = factor(c(Year, Year - 1, Year - 2), levels = 1978:2020),
                      "calYear" = rep(Year, 3),
                      "month" = c(month, month + 12, month + 24),
                      "calMonth" = rep(month, 3),
                      "MaxLength" = c(Age0MaxLength,
                                      Age1MaxLength,
                                      Age2MaxLength),
                      "MaxLengthAdj" = NA)
         })
)
  
  
#### Adjust the classes by regression ####
assign(paste("AgeRegGI", t, sep = "_") , 
       gam(MaxLength ~ s(month, k = 6, m = 2, bs = "tp") +
             s(month, by = cohort.f, k = 6, m = 1, bs = "tp") +
             s(cohort.f, bs = "re", k = 6), # not sure about k here
           data = eval(parse(text = paste0("ageClust", t, "Long"))),
           method = "REML",
           drop.unused.levels = TRUE))
  
  
  ageClustTmp$Age0MaxAdj <- predict(eval(parse(text = paste("AgeRegGI", t, sep = "_"))), 
                                 newdata = data.frame(
                                   "month" = ageClustTmp$month, 
                                   "cohort.f" = factor(ageClustTmp$Year,
                                                       levels = 1980:2020)))
  ageClustTmp$Age1MaxAdj <- predict(eval(parse(text = paste("AgeRegGI", t, sep = "_"))), 
                                 newdata = data.frame(
                                   "month" = ageClustTmp$month + 12, 
                                   "cohort.f" = factor(ageClustTmp$Year - 1,
                                                       levels = 1980:2020)))
  
  ageClustTmp$Age2MaxAdj <- predict(eval(parse(text = paste("AgeRegGI", t, sep = "_"))), 
                                 newdata = data.frame(
                                   "month" = ageClustTmp$month + 24, 
                                   "cohort.f" = factor(ageClustTmp$Year - 2,
                                                       levels = 1980:2020)))
  
  
assign(paste0("ageClust", t), ageClustTmp)
  
assign(paste0("MaxLengthAdj", t), 
       with(eval(parse(text = paste0("ageClust", t, "Long"))),
            as.vector(predict(eval(parse(text = paste("AgeRegGI", t, sep = "_"))), 
                              newdata = gamDat))))
rm(ageClustTmp)
}


#### Classify the testing data ####
# - classify each group against the remaining groups
# start with group 1 as the "testing" data
# look at how 1 would be classified by each of 2:5 

# Dataset to classify lengths from
testDat1 <- allDatWY[which(allDatWY$testGroup == 1),]

for(t in 2:5){
  ClusterAge <- rep(NA, length(testDat1$Age))
  ClusterMonth <- rep(NA, length(testDat1$Age))
  ClusterCohort <- rep(NA, length(testDat1$Age))
  for(y in 1980:2016){ #reliable years for data
    for(j in 1:12){
      # Dataset that has max lengths to use for classifying the test dataset
      assign("tmp", eval(parse(text = paste0("ageClust", t))))
      tmp <- tmp[which(tmp$Year == y & tmp$month == j), ]
      xs <- which(!is.na(tmp[1, c("Age0MaxAdj", "Age1MaxAdj", "Age2MaxAdj")]))
      ys <- tmp[1, c("Age0MaxAdj", "Age1MaxAdj", "Age2MaxAdj")[xs]]
      lx = length(xs) # 0 if all of the x's are NA
      
      breaks <- as.vector(unlist(c(0, 
                                   tmp[1, c("Age0MaxAdj", "Age1MaxAdj", "Age2MaxAdj")], 
                                   350)))
      naList <- which(is.na(breaks))
      if(length(naList) > 0){
        for(d in max(naList):min(naList)){
          breaks[d] <- breaks[d +1]-1
        }
      }
      
      labs <- 0:3
      ClusterAge[which(testDat1$Year == y & testDat1$month == j)] <- as.numeric(
        as.character(cut(x = testDat1$Length[which(testDat1$Year == y & testDat1$month == j)],
                         breaks = breaks,
                         labels = labs)))
      
      ClusterMonth[which(testDat1$Year == y & testDat1$month == j)] <- 
        testDat1$month[which(testDat1$Year == y & testDat1$month == j)] + 
        12*ClusterAge[which(testDat1$Year == y & testDat1$month == j)]
      
      ClusterCohort[which(testDat1$Year == y & testDat1$month == j)] <- t - 
        ClusterAge[which(testDat1$Year == y & testDat1$month == j)]
    }
  }
  # Bundle the age assignments into a data.frame
  # naming logic: classifying data from testDat1 into clusters based on fold t
  assign(paste0("testDat1Clust_", t), data.frame(ClusterAge,
                                                 ClusterMonth,
                                                 ClusterCohort))
}
rm(ClusterAge, ClusterMonth, ClusterCohort)

#### Compare age classes by KD method and SFBS cutoffs ####
# Ages by KD method
summary(as.factor(testDat1Clust_2$ClusterAge))
# Bay Study ages
summary(as.factor(allDatWY[which(allDatWY$testGroup == 1), "Age"]))

#### Make a confusion matrix against SFBS age classes ####
# - Cluster age = vertical labels down the left side (includes age 3)
# - SFBS age = horizontal labels across the top (only goes up to 2)
table(testDat1Clust_2$ClusterAge, allDatWY[which(allDatWY$testGroup == 1), "Age"])
table(testDat1Clust_3$ClusterAge, allDatWY[which(allDatWY$testGroup == 1), "Age"])
table(testDat1Clust_4$ClusterAge, allDatWY[which(allDatWY$testGroup == 1), "Age"])
table(testDat1Clust_5$ClusterAge, allDatWY[which(allDatWY$testGroup == 1), "Age"])

calcTotAcc <- function(x, y){
  tab <- table(x, y)
  dd <- sum(diag(tab)) # sum of diag = matching classifications
  tot <- sum(tab)      # total number of values summarized by the table
  acc <- dd/tot
  return(acc)
}

#### Summarize classification ability ####
# - assume in table(x, y) y = "truth", x = "prediction" (labels across the top are "truth")
# Overall Accuracy = count(same class)/total
# all have the same length = length(na.omit(testDat1Clust_2$ClusterAge)) = 72667
matrix(c(
  sum(diag(table(testDat1Clust_2$ClusterAge, testDat1Clust_2$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_2$ClusterAge, testDat1Clust_3$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_2$ClusterAge, testDat1Clust_4$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_2$ClusterAge, testDat1Clust_5$ClusterAge)))/72667,
  
  sum(diag(table(testDat1Clust_3$ClusterAge, testDat1Clust_2$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_3$ClusterAge, testDat1Clust_3$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_3$ClusterAge, testDat1Clust_4$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_3$ClusterAge, testDat1Clust_5$ClusterAge)))/72667,
  
  sum(diag(table(testDat1Clust_4$ClusterAge, testDat1Clust_2$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_4$ClusterAge, testDat1Clust_3$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_4$ClusterAge, testDat1Clust_4$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_4$ClusterAge, testDat1Clust_5$ClusterAge)))/72667,
  
  sum(diag(table(testDat1Clust_5$ClusterAge, testDat1Clust_2$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_5$ClusterAge, testDat1Clust_3$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_5$ClusterAge, testDat1Clust_4$ClusterAge)))/72667,
  sum(diag(table(testDat1Clust_5$ClusterAge, testDat1Clust_5$ClusterAge)))/72667
),
ncol = 4, byrow = TRUE)
# diag values should be 1 and they are

# Misclassification rate = 1-Accuracy

# Class-wise metrics are appropriate for multinomial classification

