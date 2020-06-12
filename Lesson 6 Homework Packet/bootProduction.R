# bootstrapping confidence intervals for most things in regression 
# isn't all that hard when using the boot package
#
# it requires writing auxilliary functions for extracting 
# the information you need
#
# bootstrapping a confidence interval for the population mean response (value of y)
# for a fixed value of the predictor variables (x) isn't hard either.  However,
# getting prediction intervals for individual responses is more difficult and we'll 
# leave that for a more advance class.


# same data as in the presentation
production = data.frame(
  NumItems = c(175,189,344, 88,114,338,271,173,284,277,337, 58,146,277,123,227, 63,337,146, 68),
  Time     = c(195,215,243,162,185,231,234,166,253,196,220,168,207,225,169,215,147,230,208,172)
)

# load the boot package, install it if you haven't already
require(boot)

# throughout, we'll use what is sometimes called case resampling
# if each row of your data frame corresponds to an individual observation
# of (x,y), then we randomly choose among the rows (cases) with replacement.
# Case resampling is very general.
# Another kind of resampling is called residual resampling
# which is both more complicated and requires more restrictive
# assumptions.

# Bootstrap the correlation coefficient
bootCorr <- function(df, i){
  # return the correlation coefficient
  cor.test( df$NumItems[i], df$Time[i] )$estimate
}

# find bootstrap distribution
bootCorr.obj <- boot( production, bootCorr, R = 1000)
# compute CI's, this gives several different confidence intervals
# but we want only the BCa interval
boot.ci(bootCorr.obj)

# compare to classic CI
cor.test( production$NumItems, production$Time )$conf.int

#######################################################

# bootstrap regression coefficients
bootCoef <- function( df, i ){
  coef(lm( Time ~ NumItems, data = df[i,]))
}

# find bootstrap distribution
bootCoef.obj <- boot( production, bootCoef, R = 1000)
# compute CI's, index argument is necessary to determine which
# coefficient to bootstrap
boot.ci(bootCeof.obj,index=1) # intercept
boot.ci(bootCoef.obj,index=2) # slope

# compare to classic CI's
lin.mod <- lm( Time ~ NumItems, data = production)
confint(lin.mod)

#################################################################

# bootstrap population mean response
# you can have multiple predictor values here, like NumItems = c(200, 225, 230)
# but then you'll need to add index = 1, 2, or 3 to extract the 
# confidence intervals in boot.ci below
new = data.frame( NumItems = 200)

# bootstrap mean response
bootMeanResp <- function( df, i, newdata ){
  lin.mod <- lm( Time ~ NumItems, data = df[i,])
  predict( lin.mod, newdata )
}

# get bootstrap distribution
# the last argument "newdata = new" gets passed to the bootMeanResp() function
# for making the predictions
bootResp.obj <- boot( production, bootMeanResp, R = 1000, newdata = new)
# get confidence interval
boot.ci( bootResp.obj ) # if newdata has more than one entry, will need to add index = .. here

# compare to classic interval
predict(lin.mod,newdata,interval="confidence")
