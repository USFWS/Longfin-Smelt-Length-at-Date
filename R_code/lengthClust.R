# File name: lengthClust.R
# Author: Vanessa Tobias <vanessa_tobias@fws.gov>
# Date: 2021-08-31

# This script uses data on the forklengths of Longfin Smelt that were captured
# in Interagency Ecological Program (IEP) surveys to develop length-at-date
# cutoffs for specific year classes of Longfin Smelt. It builds on existing
# static length criteria to create year class-specific length criteria.

# This script was developed to support the SSA for the Bay-Delta DPS of 
# Longfin Smelt. Documentation of methods, results, and interpretation can be
# found in Appendix C of that report.


#### 1. SETUP ####

#### 1.1 Functions ####

#### 1.2 Load Libraries ####
library(lme4)
library(mgcv)
library(RColorBrewer)
library(mclust)
library(beepr)
library(tidyverse)
library(ggridges)

#### 1.3 Prep Data ####

allDatWY <- read.csv(".\\Data_Original\\allDatLengthsWY_20210830.csv", 
                     stringsAsFactors = FALSE)

allDatWY$AgeMonth.f <- as.factor(allDatWY$AgeMonth)

ageClust <- data.frame(Year = rep(unique(allDatWY$Year), each = 12), 
                       month = rep(1:12, length(unique(allDatWY$Year))), 
                       # break1 = NA,
                       # break2 = NA,
                       # break3 = NA,
                       # break4 = NA,
                       groupCount = NA,
                       centers = NA,
                       cutoffs = NA,
                       Age0MaxLength = NA,
                       Age1MaxLength = NA,
                       Age2MaxLength = NA)
ageClust <- ageClust[which(ageClust$Year > 1979), ]
ageClust <- ageClust[order(ageClust$Year, ageClust$month),]


#### 1.4 Set prior guesses for length classes ####

# use SFBS as a starting place, but make a few changes to account for
# information in other studies. (E.g., SFBS consisent floor at 40mm is problematic
# in the classification code; some potential age-3 fish are also problematic)
# Bay Study length classes

lfAgeLen <- read.csv(".\\Data_Original\\BayStudyAgeLength.csv",
                     stringsAsFactors = FALSE)

priors <- lfAgeLen
priors$Age2MaxLength <- lfAgeLen$Age1MaxLength+45
# added value above is somewhat arbitrary. Idea is to separate double age-1 cutoffs.


#### 2. CLASSIFY LENGTHS ####

#### 2.1 Do the classifications ####

# !!! BW can influence the x value for the peaks in the KD method

# This version was copied from mclust.R on 2/3/2021.
#  AND UPDATED ON 2/10/2021 to run mclust on the KD instead of raw lengths.
# It does classification using the Kernel Density (KD) method that I wrote 
#  and the Mixing method using the mclust package.
# It writes a column of predictions to the original dataset.
# Note that DensAge and CutAge are the raw versions of the KD method. You
#  still need to run the regression code to smooth out the classifications.
# The smoothed cutoffs are used to generate clustAge, which is the version that
#  should be used to correct the ages.

# https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
# -- Univariate density estimation section

allDatWY$DensAge <- NA
allDatWY$CutAge <- NA
allDatWY$mclustAge <- NA

for(i in unique(allDatWY$Year[which(allDatWY$Year > 1979 &
                                    allDatWY$Program != "SLS")])){ 
  for(j in 1:12){
    # Subset the data
    tmp <- allDatWY[which(allDatWY$month == j & allDatWY$Year ==i),]
    
    #### Kernel Density ####
    # Catch month and year combinations with zero records
    if (length(na.omit(tmp$Length)) %in% c(0, 1)) {
      
      ageClust[which(ageClust$Year == i & ageClust$month == j), 
               c("groupCount", "centers", "cutoffs")] <- c(NA, NA, NA)
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
          ageClust[which(ageClust$Year == i & ageClust$month == j),  cutoffs$age[k]] <- cutoffs$cutoffs[k]
        }
      }
      
      
      #### Use mclust to find clusters ####
      
      # run mclust
      # Makes density calcs using a Gaussian finite mixture model from Mclust
      # Still making the centers using Mclust and getting unrealistic results.
      # Make the density calcs and apply logic from KD method?
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

      
      #### Save the results ####
      allDatWY$DensAge[which(allDatWY$month == j &
                               allDatWY$Year ==i &
                               allDatWY$Program != "SLS")] <- tmp$DensAge
      allDatWY$CutAge[which(allDatWY$month == j &
                              allDatWY$Year ==i &
                              allDatWY$Program != "SLS")] <- tmp$CutAge
      allDatWY$mclustAge[which(allDatWY$month == j &
                                 allDatWY$Year ==i &
                                 allDatWY$Program != "SLS")] <- tmp$mclustAge
      
      ageClust[which(ageClust$Year == i & ageClust$month == j), 
               c("groupCount", "centers", "cutoffs",
                 "groupCount.mclust", "centers.mclust")] <- c(sum(maxMin == "Max"),
                                                              as.character(list(tpx[which(maxMin == "Max")])),
                                                              as.character(list(if(length(tpx[which(maxMin == "Min")]) == 0) {NA} else tpx[which(maxMin == "Min")])),
                                                              groupCount.mclust,
                                                              as.character(list(clust$parameters$mean)))
    }
  }
  rm(tmp)
}
ageClust$groupCount <- as.numeric(ageClust$groupCount)

ageClustLong <- data.frame("cohort" = c(ageClust$Year, ageClust$Year - 1, ageClust$Year - 2),
                           "calYear" = rep(ageClust$Year, 3),
                           "month" = c(ageClust$month, ageClust$month + 12, ageClust$month + 24),
                           "calMonth" = rep(ageClust$month, 3),
                           "MaxLength" = c(ageClust$Age0MaxLength,
                                               ageClust$Age1MaxLength,
                                               ageClust$Age2MaxLength))
ageClustLong$cohort.f <- as.factor(ageClustLong$cohort)
# Limit the analysis to before 2015 because of issues with effort
# ageClustLong <- ageClustLong[which(ageClustLong$calYear < 2016),]


#### 2.2 Adjust KD cutoffs by regression ####

# Discrete months
plot(ageClust$month, ageClust$Age0MaxLength,
     xlim = c(1, 12),
     ylim = c(0, 250))
points(ageClust$month, ageClust$Age1MaxLength,
       col = "red")
points(ageClust$month, ageClust$Age2MaxLength,
       col = "blue")

# Running months
plot(ageClust$month, ageClust$Age0MaxLength,
     xlim = c(1, 36),
     ylim = c(0, 250))
points(ageClust$month+12, ageClust$Age1MaxLength,
       col = "red")
points(ageClust$month+24, ageClust$Age2MaxLength,
       col = "blue")

plot(ageClustLong$month, ageClustLong$MaxLength,pch = 16, col = "grey")

plot(ageClustLong$month, ageClustLong$MaxLength,
     log = "y")

plot(ageClust$month, ageClust$Age0MaxLength,
     xlim = c(1, 36),
     ylim = c(20, 250),
     log = "y")
points(ageClust$month+12, ageClust$Age1MaxLength,
       col = "red")
points(ageClust$month+24, ageClust$Age2MaxLength,
       col = "blue")

# Some of the cutoff points are missing or they're unreasonable, given the
#   other points next to them. Use regression to smooth things out.

# Pedersen et al. 2019 - Model GS or GI fits our needs here
# GS = shared global trend + group level similar smoothness (shared penalty)
#    similar to a GLMM with varying slopes
AgeRegGS <- gam(MaxLength ~ s(month, k = 6) +
                  s(month, cohort.f, bs = "fs", k = 6, m = 2),
              data = ageClustLong,
              method = "REML")
summary(AgeRegGS)
gam.check(AgeRegGS)
plot(AgeRegGS)

# GI = shared global trend + group level different smoothness (individual penalty)
AgeRegGI <- gam(MaxLength ~ s(month, k = 6, m = 2, bs = "tp") +
                  s(month, by = cohort.f, k = 6, m = 1, bs = "tp") +
                  s(cohort.f, bs = "re", k = 6), # not sure about k here
                data = ageClustLong,
                method = "REML")
summary(AgeRegGI)
gam.check(AgeRegGI)
plot(AgeRegGI)
anova.gam(AgeRegGS, AgeRegGI, test = "Chisq")


ageClustLong$MaxLengthAdj <- as.vector(predict(AgeRegGI, 
                                     newdata = ageClustLong))

ageClustWide <- data.frame(
  "calYear" = ageClustLong$calYear,
  "calMonth" = ageClustLong$calMonth,
  "age0MaxAdj" = ageClustLong$MaxLengthAdj[which(ageClustLong$month < 13)])
ageClustWide <- merge(ageClustWide, 
                      ageClustLong[which(ageClustLong$month > 12 & ageClustLong$month < 25),
                                    c("calYear", "calMonth", "MaxLengthAdj")])
names(ageClustWide)[4] <- "age1MaxAdj"
ageClustWide <- merge(ageClustWide, 
                      ageClustLong[which(ageClustLong$month > 24),
                                   c("calYear", "calMonth", "MaxLengthAdj")])
names(ageClustWide)[5] <- "age2MaxAdj"

ageClustWide <- ageClustWide[!duplicated(ageClustWide), ]
ageClustWide <- ageClustWide[order(ageClustWide$calYear, ageClustWide$calMonth),]

plot(0, 0, type = "n",
     xlim = c(1, 36),
     ylim = c(0, 250),
     xaxs = "i",
     yaxs = "i",
     xlab = "Month",
     ylab = "Fork Length (mm)")
lines(priors$Month + 12*rep(0:2, each = length(priors$Month)), 
      c(priors$Age0MaxLength, priors$Age1MaxLength, priors$Age2MaxLength))

for(i in unique(ageClustWide$calYear)){
  
  tmp <- ageClustWide[which(ageClustWide$calYear == i),]
  
  with(tmp,
       lines(c(calMonth, calMonth + 12, calMonth + 24),
             c(age0MaxAdj, age1MaxAdj, age2MaxAdj),
             col = rgb(0, 0, 0, 1/4)))
}


pred.ageClust <- ageClust 
pred.ageClust$Year.f <- as.factor(pred.ageClust$Year)
pred.ageClust <- pred.ageClust[order(pred.ageClust$Year,
                                     pred.ageClust$month),]

pred.ageClust$Age0MaxAdj <- predict(AgeRegGI, 
                               newdata = data.frame("month" = pred.ageClust$month, 
                                                    "cohort.f" = pred.ageClust$Year.f))
pred.ageClust$Age1MaxAdj <- predict(AgeRegGI, 
                               newdata = data.frame(
                                 "month" = pred.ageClust$month + 12, 
                                 "cohort.f" = factor(pred.ageClust$Year - 1,
                                                      levels = levels(pred.ageClust$Year.f))))

pred.ageClust$Age2MaxAdj <- predict(AgeRegGI, 
                                    newdata = data.frame(
                                      "month" = pred.ageClust$month + 24, 
                                      "cohort.f" = factor(pred.ageClust$Year - 2,
                                                          levels = levels(pred.ageClust$Year.f))))

pred.ageClust <- merge(pred.ageClust, 
      ageClustLong[which(ageClustLong$month < 13), c("calYear", "calMonth", "MaxLengthAdj")],
      by.x = c("Year", "month"),
      by.y = c("calYear", "calMonth"),
      all.x = TRUE)
names(pred.ageClust)[15] <- "Age0MaxAdj_v2"
pred.ageClust <- merge(pred.ageClust, 
                       ageClustLong[which(ageClustLong$month > 12 &
                                            ageClustLong$month < 25), 
                                    c("calYear", "calMonth", "MaxLengthAdj")],
                       by.x = c("Year", "month"),
                       by.y = c("calYear", "calMonth"),
                       all.x = TRUE)
names(pred.ageClust)[16] <- "Age1MaxAdj_v2"
pred.ageClust <- merge(pred.ageClust, 
                       ageClustLong[which(ageClustLong$month > 24), 
                                    c("calYear", "calMonth", "MaxLengthAdj")],
                       by.x = c("Year", "month"),
                       by.y = c("calYear", "calMonth"),
                       all.x = TRUE)
names(pred.ageClust)[17] <- "Age2MaxAdj_v2"

pred.ageClust <- pred.ageClust[order(c(pred.ageClust$month,
                                       pred.ageClust$Year)),]


# Plot Cutoffs for KD method and SFBS together ####
# -- ageClustLong as the basis for the lines --> cohorts
plot(0, 0, type = "n",
     xlim = c(1, 36),
     ylim = c(0, 250),
     xaxs = "i",
     yaxs = "i",
     xlab = "Month",
     ylab = "Fork Length (mm)")
for(i in unique(ageClustLong$cohort)){
  with(ageClustLong[which(ageClustLong$cohort == i),],
       lines(month,
             MaxLengthAdj,
             col = rgb(0, 0, 0, 1/4)))
}
lines(priors$Month + 12*rep(0:2, each = length(priors$Month)),
      c(priors$Age0MaxLength, priors$Age1MaxLength, priors$Age2MaxLength))

data.frame(ageClustLong$calYear,
           ageClustLong$calMonth)



# 2.3 Use regression-adjusted cutoffs ####
allDatWY$ClusterAge <- NA
allDatWY$ClusterMonth <- NA
allDatWY$ClusterCohort <- NA
allDatWY$mclustMonth <- NA
allDatWY$mclustCohort <- NA
for(t in unique(allDatWY$Year[which(allDatWY$Year > 1980)])){
  #[which(allDatWY$Year > 1981 & allDatWY$Year < 2016)]
  for(j in 1:12){
    tmp <- pred.ageClust[which(pred.ageClust$Year == t & pred.ageClust$month == j), ]
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
      #classify based on distribution cutoffs
    allDatWY$ClusterAge[which(allDatWY$Year == t & allDatWY$month == j)] <- as.numeric(
      as.character(cut(x = allDatWY$Length[which(allDatWY$Year == t & allDatWY$month == j)],
                                  breaks = breaks,
                                  labels = labs)))
    
    allDatWY$ClusterMonth[which(allDatWY$Year == t & allDatWY$month == j)] <- 
      allDatWY$month[which(allDatWY$Year == t & allDatWY$month == j)] + 
      12*allDatWY$ClusterAge[which(allDatWY$Year == t & allDatWY$month == j)]
    
    allDatWY$ClusterCohort[which(allDatWY$Year == t & allDatWY$month == j)] <- t - 
      allDatWY$ClusterAge[which(allDatWY$Year == t & allDatWY$month == j)]
    
    # classify based on mclust ages
    # code already makes age classifications. Just need to calculate months & cohorts
    
    allDatWY$mclustMonth[which(allDatWY$Year == t & allDatWY$month == j)] <- 
      allDatWY$month[which(allDatWY$Year == t & allDatWY$month == j)] + 
      12*allDatWY$mclustAge[which(allDatWY$Year == t & allDatWY$month == j)]
    
    allDatWY$mclustCohort[which(allDatWY$Year == t & allDatWY$month == j)] <- t - 
      allDatWY$mclustAge[which(allDatWY$Year == t & allDatWY$month == j)]
    }
    }

allDatWY$ClusterCohort.f <- factor(allDatWY$ClusterCohort)

#### 2.3 Write a copy of the data with the classifications ####
rm(tmp, ys, centers, centers.mclust, clust, clustClass, cutoffs, den, 
   breaks, der1, groupCount.mclust, i, j, k, labs, lx, maxMin, 
   t, tpSign, tpx, tpy, xs)

write.csv(allDatWY, 
          ".\\Data_Derived\\allDatLengthsWYClusters_20210830.csv",
          row.names = FALSE)
write.csv(allDatWY[, c("Year", "AgeMonth", "Length", "Age", "Cohort",
                       "month", "Date", "Program",
                       "ClusterAge", "ClusterMonth", "ClusterCohort")], 
          ".\\Data_Derived\\Length_Clusters_20210830.csv",
          row.names = FALSE)

# Read a the copy back in if you need it...
# allDatWY <- read.csv(".\\Data_Derived\\allDatLengthsWYClusters_20210830.csv",
#                      #row.names = FALSE,
#                      stringsAsFactors = FALSE)
# Fix the factors...
allDatWY$AgeMonth.f <- as.factor(allDatWY$AgeMonth.f)
allDatWY$Cohort.f <- as.factor(allDatWY$Cohort.f)
allDatWY$WYearType.0.f <- factor(allDatWY$WYearType.0.f,
                                 levels = c("C", "D", "BN", "AN", "W"))


#### 2.4 Check for differences among age classification methods ####
# make a vector where the max age is 2+
allDatWY$ClusterAge2p <- allDatWY$ClusterAge
allDatWY$ClusterAge2p[which(allDatWY$ClusterAge2p > 2)] <- 2

# Table 2 from TN02/Appendix C ####
tableComp <- with(allDatWY[which(allDatWY$Year > 1981 & allDatWY$Year < 2016),], {
  data.frame("Comparison" = c("SFBS v. Kernel Density",
                              "Mixing v. Kernel Density",
                              "SFBS v. Mixing"),
             "sameCount" = c(sum(Age == ClusterAge2p),
                        sum(mclustAge == ClusterAge2p),
                        sum(Age == mclustAge)))
})
tableComp$samePct <- tableComp$sameCount/length(allDatWY$Age[which(allDatWY$Year > 1981 & allDatWY$Year < 2016)])

with(allDatWY[which(allDatWY$Year > 1981 & allDatWY$Year < 2016),], {
  print(table(Age, ClusterAge2p))
})




##### 3. Plot results ####


plot(allDatWY$ClusterMonth,
     allDatWY$Length,
     pch = 16, col = as.factor(allDatWY$Program))
legend("topleft", pch = 16, 
       col = 1:length(levels(as.factor(allDatWY$Program))),
       legend = levels(as.factor(allDatWY$Program)), 
       cex = 0.5)

for(i in unique(allDatWY$Program)){
  png(paste0("./Figures/cluster_class_", i, ".png"),
      height = 5, width = 6, units = "in", res = 300)
  with(allDatWY[which(allDatWY$Program == i),], {
    plot(ClusterMonth,
         Length,
         main = i,
         xlim = c(0, 55),
         ylim = c(0, 300))
  })
  abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)
  dev.off()
}

# -- One plot per calendar month over years ####
for(j in 1:12){
  
  print(paste("j = ", j))
  
  # Make the blank plot frame
  par(cex = 1.25)
  plot(0, 0, type = "n",
       xlim = c(1980, 2022),
       ylim = c(0, 225),
       xlab = "Calendar Year",
       ylab = "Fork Length",
       main = paste(month.name[j], "Cutoffs"))
  
# Plot SFBS cutoffs (priors) as horizontal lines
  abline(h = priors$Age0MaxLength[j], lty = 2)
  abline(h = priors$Age1MaxLength[j], lty = 2)  
  abline(h = priors$Age2MaxLength[j], lty = 2)
  
  # Label the regions
  text(x = 2022, y = (lfAgeLen$Age1MaxLength[j] + 225)/2, "Age-2+", srt = 90)
  text(x = 2022, y = (lfAgeLen$Age1MaxLength[j] + lfAgeLen$Age0MaxLength[j])/2, "Age-1", srt = 90)
  text(x = 2022, y = (lfAgeLen$Age0MaxLength[j])/2, "Age-0", srt = 90)  
  
  # Plot the raw length data
  with(allDatWY[which(allDatWY$month == j),], {
    points(Year, Length,
         pch = ".", col = rgb(0, 0, 0, 1/5))
  }
       )
  
  # Subset the calculated clusters to the given month
  tmp <- ageClust[which(ageClust$month == j),]
  
  # Plot markers for the cluster edges and centers
  for(i in 1:length(tmp$Year[which(tmp$Year > 1979)])){
    print(i)
    
    # Skip years with only 1 group for now:
    if(is.na(eval(parse(text = tmp$cutoffs[i])))) next 
    
    # Plot centers of groups
    points(rep(tmp$Year[i], length(eval(parse(text = tmp$centers[i])))),
           eval(parse(text = tmp$centers[i])),
           pch = "-", col = "grey")

    # Plot cutoffs between groups
    points(rep(tmp$Year[i], length(eval(parse(text = tmp$cutoffs[i])))),
           eval(parse(text = tmp$cutoffs[i])))
  }
  legend("topright", horiz = TRUE,
         pch = c("-", "°", "."), 
         pt.cex = 1.5,
         cex = 0.8,
         col = c("grey", "black", "black"),
         legend = c("Centers", "Cutoffs", "Data Pts"),
         bty = "n")
}

# -- One plot per year, over calendar months ####
pdf("./Figures/length_clusters_xmonth.pdf")
for(j in unique(ageClust$Year[which(ageClust$Year > 1979)])){
  
  print(paste("j = ", j))
  
  # Make the blank plot frame
  par(cex = 1.25)
  plot(0, 0, type = "n",
       xlim = c(1, 13),
       ylim = c(0, 225),
       xlab = "Month of Calendar Year",
       ylab = "Fork Length",
       main = j)
  
  # Label the regions
  text(x = 13, y = (lfAgeLen$Age1MaxLength[12] + 225)/2, "Age-2+", srt = 90)
  text(x = 13, y = (lfAgeLen$Age1MaxLength[12] + lfAgeLen$Age0MaxLength[12])/2, "Age-1", srt = 90)
  text(x = 13, y = (lfAgeLen$Age0MaxLength[12])/2, "Age-0", srt = 90)  
  

  # Subset the calculated clusters to the given year
  tmp <- ageClust[which(ageClust$Year == j),]
  
  # Plot the raw length data  
  with(allDatWY[which(allDatWY$Year == j),], {
    points(jitter(month), Length,
           pch = ".", 
           cex = 2,
           col = as.factor(ClusterAge)) #rgb(0, 0, 0, 1/5))
  }
  )  
  
  # Plot markers for the cluster edges and centers
  for(i in 1:12){
    print(i)
    
    # Plot centers of groups
    points(rep(tmp$month[i], length(eval(parse(text = tmp$centers[i])))),
           eval(parse(text = tmp$centers[i])),
           pch = 16, col = "grey")
    
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age0MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age1MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age2MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    
    
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age0MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age1MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age2MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
  }
}
dev.off()

# -- One plot per year with months on the x-axis ####
pdf("./Figures/length_clusters_xmonth.pdf")
for(j in unique(ageClust$Year[which(ageClust$Year > 1979)])){
  
  print(paste("j = ", j))
  
  # Make the blank plot frame
  par(cex = 1.25)
  plot(0, 0, type = "n",
       xlim = c(1, 13),
       ylim = c(0, 225),
       xlab = "Month of Calendar Year",
       ylab = "Fork Length",
       main = j)
  
  # Label the regions
  text(x = 13, y = (lfAgeLen$Age1MaxLength[12] + 225)/2, "Age-2+", srt = 90)
  text(x = 13, y = (lfAgeLen$Age1MaxLength[12] + lfAgeLen$Age0MaxLength[12])/2, "Age-1", srt = 90)
  text(x = 13, y = (lfAgeLen$Age0MaxLength[12])/2, "Age-0", srt = 90)  
  
  
  # Subset the calculated clusters to the given year
  tmp <- ageClust[which(ageClust$Year == j),]
  
  # Plot the raw length data  
  with(allDatWY[which(allDatWY$Year == j),], {
    points(jitter(month), Length,
           pch = ".", 
           cex = 2,
           col = as.factor(ClusterAge)) #rgb(0, 0, 0, 1/5))
  }
  )  
  
  # Plot markers for the cluster edges and centers
  for(i in 1:12){
    print(i)
    
    # Plot centers of groups
    points(rep(tmp$month[i], length(eval(parse(text = tmp$centers[i])))),
           eval(parse(text = tmp$centers[i])),
           pch = 16, col = "grey")
    
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age0MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age1MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    points(pred.ageClust$month[which(pred.ageClust$Year == j)], 
           pred.ageClust$Age2MaxLength[which(pred.ageClust$Year == j)],
           col = "blue")
    
    
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age0MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age1MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
    lines(pred.ageClust$month[which(pred.ageClust$Year == j)], 
          pred.ageClust$Age2MaxAdj[which(pred.ageClust$Year == j)],
          col = "blue", lwd = 2, lty = 2)
  }
}
dev.off()


pdf("./Figures/lengthDistsClusterCohorts.pdf")
par(mfrow = c(3, 1))
for(i in 1980:2015){
  p <- ggplot(allDatWY[which(allDatWY$ClusterCohort == i),], 
              aes(x = Length, 
                  y = fct_rev(as.factor(ClusterMonth)))) +  #reverse the order of years = more intuitive display(?)
    geom_density_ridges(bandwidth = 6.5, alpha = 0.5) +
    scale_x_continuous(limits = c(40, 150)) +
    scale_y_discrete(breaks = fct_rev(factor(levels(as.factor(allDatWY$ClusterMonth)))), #ugly but effective
                     limits = levels(fct_rev(as.factor(allDatWY$ClusterMonth)))) +
    theme_classic() +
    ggtitle(paste("Year Class:", i))
  print(p)
}
dev.off()

# -- Example Cohort for TN ####
# Figure 5 in TN2/Appendix C ####
png("./Figures/lengthDistsClusterCohorts1982.png",
    width = 4, height = 6, units = "in", res = 300)
par(mfrow = c(3, 1))
for(i in 1982){
  p <- ggplot(allDatWY[which(allDatWY$ClusterCohort == i),], 
              aes(x = Length, 
                  y = fct_rev(as.factor(ClusterMonth)))) +  #reverse the order of years = more intuitive display(?)
    geom_density_ridges(bandwidth = 6.5, alpha = 0.5) +
    scale_x_continuous(limits = c(40, 150)) +
    scale_y_discrete(breaks = fct_rev(factor(
      levels(as.factor(allDatWY$ClusterMonth)))), #ugly but effective
                     limits = levels(
                       fct_rev(as.factor(
                         allDatWY$ClusterMonth[which(allDatWY$ClusterMonth < 37)])))) +
    theme_classic() +
    ggtitle(paste("Year Class:", i)) +
    xlab("Fork Length (mm)") +
    ylab("Month of Age Class")
  print(p)
}
dev.off()


#### 4. Visualize growth patterns ####

# Fork length by age in months for each cohort
pdf("./Figures/yLength_xClusterMonth_byCCohort.pdf")
for(i in unique(allDatWY$ClusterCohort)){
  with(allDatWY[which(allDatWY$ClusterCohort == i),],
         plot(ClusterMonth, Length,
              xlim = c(0, 40),
              ylim = c(0, 200),
              pch = 16,
              xlab = "Month",
              ylab = "Fork Length",
              main = paste("Year Class", i, sep = " ")))
  abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)
}
dev.off()

# Plot of forklength by month for all data
# - Normal y axis
plot(allDatWY$ClusterMonth, 
     allDatWY$Length)
abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)
# - Logged y axis
plot(allDatWY$ClusterMonth, 
     allDatWY$Length,
     log = "y")
abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)

# Check patterns for fish classified as > 36 months of age
plot(allDatWY$ClusterMonth[which(allDatWY$ClusterMonth < 37)], 
     allDatWY$Length[which(allDatWY$ClusterMonth < 37)],
     xlim = c(1, 36),
     ylim = c(0, 200))
points(allDatWY$ClusterMonth[which(allDatWY$ClusterMonth > 36)] - 12, 
       allDatWY$Length[which(allDatWY$ClusterMonth > 36)], col = "red")
abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)


# Fork length v. age in months assigned by the cluster procedure
plot(allDatWY$mclustMonth, allDatWY$Length)
abline(v = c(12.5, 24.5, 36.5), col = "grey", lty = 2)
with(na.omit(allDatWY[, c("mclustMonth", "Length")]), 
     lines(smooth.spline(mclustMonth, Length),
           col = "red"))


# # Investigate cluster of very small lengths in months 15, 16, and 17
# #  -- should those be months 3, 4, and 5?
# plot(allDatWY$ClusterMonth, allDatWY$Length,
#      xlim = c(13, 19),
#      ylim = c(0, 40))
# with(allDatWY[which(allDatWY$ClusterMonth %in% 13:19 &
#                       allDatWY$Length < 40),],
#      text(ClusterMonth, Length, 
#           labels = paste(Year, ClusterMonth, sep = "-"),
#           pos = 4))
# # all of the problem points are in CY2016
# # not enough data points yet - Bay Study was off the water most of that year. 
# # Maybe just stop the analysis at 2015 like for the distribution paper.
# plot(allDatWY$ClusterMonth, allDatWY$Length,
#      xlim = c(13, 19),
#      ylim = c(0, 40))
# with(allDatWY[which(allDatWY$ClusterMonth %in% 13:19 &
#                       allDatWY$Length < 40),],
#      text(ClusterMonth, Length, 
#           labels = paste(ClusterCohort, ClusterMonth, sep = "-"),
#           pos = 4))
# # Cohort 2015, Age 1




