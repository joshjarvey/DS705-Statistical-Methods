anovaResample <- function(x,group,B=1000,method=2,var.equal=FALSE)
#
# code is not intended to be efficient and uses slow for loop for transparency
#  
# if equal.var==TRUE use standard anova test statistic
# if equal.var==FALSE use Welch corrected test stat
# method 3: permutation test
#   test the hypothesis that all the samples are from the same population
#   and exchangable.  each resample samples from the pooled samples without replacement
# method 4: simple bootstrap
#   similar to 1 but resamples occur with replacment (no longer a permuation test)
# method 1:  resample pooled residuals with replacment
# method 2:  resample unpooled residuals with replacement
  
  {  
  out <- oneway.test(x~group,var.equal=var.equal)
  f.anova.obs <- unname(out$statistic)  
  p.anova.obs <- unname(out$p.value)
  
  group_names <- levels(group)
  k <- nlevels(group)
  
  if (method==1|method==2){
    # compute residuals 
    xres <- x
    for (i in 1:k){
      ind <- group==group_names[i]
      xres[ind] <- x[ind] - mean(x[ind])
    }
  }
  
  f.anova.boot <- numeric(B)
  r.sq.boot <- numeric(B)
  for (i in 1:B){
    if (method==3){
      xshuffle <- sample(x,replace=FALSE)
    } else if (method==4){
      xshuffle <- sample(x,replace=TRUE)
    } else if (method==1){
      xshuffle <- sample(xres,replace=TRUE)
    } else if (method==2){
      xshuffle <- x
      for (j in 1:k){
        ind <- group==group_names[j]
        xshuffle[ind] <- sample(xres[ind],replace=TRUE)
      }
    }
    f.anova.boot[i] <- unname(oneway.test(xshuffle~group,var.equal=var.equal)$statistic)
    r.sq.boot[i] <- summary(lm(xshuffle~group))$r.squared
  }
  # correct any division by zero NaNs
  f.anova.boot[is.na(f.anova.boot)]=10*f.anova.obs
  p.anova.boot <- sum(f.anova.boot >=f.anova.obs )/B
  out <- list(aov,f.anova.obs,p.anova.obs,p.anova.boot)
  names(out) <- c('aov.obs','f.obs','p.obs','p.boot')
  if (var.equal){
    print('Assuming equal variances - using standard F')
  } else {
    print('Assuming unequal variances - using corrected F')
  }
  print(paste('observed F: ',f.anova.obs))
  print(paste('observed p-value: ',p.anova.obs))
  print(paste('resampled p-value: ',p.anova.boot))
  return(out)
}


