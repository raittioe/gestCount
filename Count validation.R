# Validating functions modified for count outcomes following example in page 19 at Tompsett, D., Vansteelandt, S., Dukes, O., & De Stavola, B. (2022). gesttools: General Purpose G-Estimation in R. Observational Studies 8(1), 1-28. https://doi.org/10.1353/obs.2022.0003.
 
library(gesttools)
library(geeM)
library(tibble)
library(tidyr)
library(rsample)

# load functions from Count functions.R

# Generate data with a Poisson distributed outcome

datas <- dataexamples2(n = 1000, seed = 123, Censoring = T)
data <- datas$datagestmult
summary(data$Y) #count outcome

# Format data

data <- FormatData(
  data = data, idvar = "id", timevar = "time", An = "A",Cn="C",
  varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1
)
idvar <- "id"
timevar <- "time"
Yn <- "Y"
An <- "A"
Cn <- "C"
outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
censoringmodel <- c("C~L+U+as.factor(time)")
EfmVar <- NA

#Analyze data, true effect estimates are 1 and 1, or exp(1) and exp(1)


gestMultipleCount(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
             censoringmodel , type = 3, EfmVar,
             cutoff = 2
)


gestbootCount(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
                       censoringmodel , type = 3, EfmVar,
                       bn=50,seed=1000,cutoff =2)




