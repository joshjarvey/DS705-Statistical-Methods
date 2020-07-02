anovaResampleFast <- function(x,group,B=1000,method=2,var.equal=FALSE)
  #
  # code is not intended to be efficient and uses slow for loop for transparency
  #  
  # if equal.var==TRUE use standard anova test statistic
  # if equal.var==FALSE use Welch corrected test stat
  # method 3: permutation test
  #   test the hypothesis that all the samples are from the same population
  #   and exchangable.  each resample samples from the pooled samples without replacement
  # method 4: simple bootstrap
  #   similar to 3 but resamples occur with replacment (no longer a permuation test)
  # method 1:  resample pooled residuals with replacment
  # method 2:  resample unpooled residuals with replacement

{  
  out <- oneway.test(x~group,var.equal=var.equal)
  f.anova.obs <- unname(out$statistic)  
  p.anova.obs <- unname(out$p.value)
  
  group <- factor(group)
  
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
  N <- length(x)
  n.samp <- as.numeric(unname(tapply(x,group,length)))
  if (method==3){
    xshuffle <- t(replicate(B,sample(x,replace=FALSE)))
  } else if (method==4){
    xshuffle <- matrix( sample(x,size=N*B,replace=TRUE), B, N)
  } else if (method==1){
    xshuffle <- matrix( sample(xres,size=N*B,replace=TRUE), B, N)
  } else if (method==2){
    xshuffle <- matrix(0,B,N)
    for (j in (1:k)){
      ind <- group==group_names[j]
      xshuffle[,ind] <- matrix(sample(xres[ind],size=n.samp[j]*B,replace=TRUE),B,n.samp[j])      
    }
  }
  
  mean.samp <- matrix(0,B,k)
  var.samp <- matrix(0,B,k)
  for (j in (1:k)){
    ind <- group==group_names[j]
    mean.samp[,j]=apply(xshuffle[,ind],1,mean)
    var.samp[,j]=apply(xshuffle[,ind],1,var)
  }
  if (var.equal){
    grandmean.samp <- apply(xshuffle,1,mean)
    ssg.samp <- numeric(B)
    sse.samp <- numeric(B)
    for (j in (1:k)){
      ind <- group==group_names[j]
      ssg.samp <- ssg.samp + n.samp[j]*(mean.samp[,j]-grandmean.samp)^2
      sse.samp <- sse.samp + (n.samp[j]-1)*var.samp[,j]
    }
    r.sq.boot <- ssg.samp/(sse.samp+ssg.samp)
    f.anova.boot <- (ssg.samp/(k-1))/(sse.samp/(N-k))
  } else{
    # apply Welch correction for heterogeneity of variance
    w.samp <- t( n.samp / (t(var.samp)) )
    wrscl.samp <- w.samp / apply(w.samp,1,sum)
    grandmean.samp <- apply(wrscl.samp*mean.samp,1,sum)
    ssg.samp <- numeric(B)*0
    mse.samp <- numeric(B)*0
    for (j in (1:k)){
      ind <- group==group_names[j]
      ssg.samp <- ssg.samp + w.samp[,j]*(mean.samp[,j]-grandmean.samp)^2
      mse.samp <- mse.samp + (1-wrscl.samp[,j])^2/(n.samp[j]-1) 
    }
    df.samp <- (k^2-1)/(3*mse.samp)
    mse.samp <- 1 + 2*(k-2)/(k^2-1)*mse.samp
    sse.samp <- mse.samp*df.samp
    r.sq.boot <- (ssg.samp)/(ssg.samp+sse.samp)
    f.anova.boot <- (ssg.samp/(k-1))/mse.samp
  }
  
  # correct any division by zero NaNs
  f.anova.boot[is.na(f.anova.boot)]=10*f.anova.obs
  p.anova.boot <- sum(f.anova.boot >=f.anova.obs )/B
  out <- list(f.anova.obs,p.anova.obs,p.anova.boot)
  names(out) <- c('f.obs','p.obs','p.boot')
  if (var.equal){
    print('Assuming equal variances - using standard F')
  } else {
    print('Assuming unequal variances - using Welch corrected F')
  }
  print(paste('observed F: ',f.anova.obs))
  print(paste('observed p-value: ',p.anova.obs))
  print(paste('resampled p-value: ',p.anova.boot))
  return(out)
}


