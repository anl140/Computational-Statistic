# Andrew Le
# A11565521

# Problem 1 
# Install this package in order to make use of the contingency table
install.packages("gmodels") 
library("gmodels", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
# Grab data to use 
admissions <- fread('http://media.pearsoncmg.com/aw/aw_deveaux_stats_2/activstats/text/Ch03_Magnet_schools.txt')
# Get all of the counts
tab = table(admissions) 

# Function to permute the data
permIndepChisqTest <- function(tab,B){
  boot <- numeric(B)
  # Grab the row
  row <- margin.table(tab,1)
  # Grab the col
  col <- margin.table(tab,2)
  # Calculate the observed chi squared value 
  d_obs <- chisq.test(tab)$stat
  
    for(b in 1:B){
      r2d <- r2dtable(1,row,col)
      ntab <- as.data.frame(r2d)
      # Store the value into boot
      boot[b] <- chisq.test(ntab)$stat
    }
  # Generate the p-val by looking at the chisquared value from the boot 
  p_val <- (sum(boot >= d_obs) +1 )/ (B+1)
  return(p_val)
}
# Test the times
a <- proc.time()
permIndepChisqTest(tab,B=2000)
proc.time() - a
# user  system elpased
# 1.129 0.186  1.339
# p-val 0.000499

# Problem 2  
bootCI <- function(x,conf,B){
 B = B
 out = numeric(B)
 alpha <- 1 - conf
 # calculate the mean value of the original data
 x_obs <- mean(x)
 # get the length
 n <- length(x)
 # get the standard deviation
 sd <- sd(x)
 
    for (b in 1:B){
      # Make the bootstrap observations
      x_boot = sample(x,length(x),replace=TRUE)
      # Calc the boot strap mean value
      x_b <- mean(x_boot)
      # Compute the standard deviation
      sb <- sqrt(var(x_boot))
      # Compute t-statistic
      tstat <- (x_b - x_obs)/(sb/sqrt(n)) 
      out[b] = tstat
    }
lval <- quantile(out, alpha/2)
rval <- quantile(out,(1-(alpha/2)))
#Return the confidence interval
return(c(x_obs + (lval* sd)/sqrt(n), x_obs + (rval*sd)/sqrt(n)))  
}

b <- proc.time()
bootCI(USArrests$Assault,conf=0.99,B=1e4)
proc.time()-b
# user system elapsed
# 0.403 0.012  0.427
#CI for this particular data set of US Arrests is (136.9408, 201.825) at 99% CI

# Problem 3
install.packages("boot")
library(boot)
# Value of the mean and variance
meanfun <- function(x,f){
  return(c(mean(x*f),var(x*f)))
}

bootCIpackage <- function(x,conf=0.99,B=1e4){
  # Calculate the bootstrap vals
  b_val <- boot(x,meanfun, R=B, stype='f')
  # Calculate the boot CI
  CI <- boot.ci(b_val)
  print(CI)
}

c <- proc.time()
bootCIpackage(USArrests$Assault,conf=0.99,B=1e4)
c - proc.time()
# user  system   elapsed
# -0.591 -0.048  -0.677
# CI (162.8,181.8) at the level of 99% confidence , studentized

b <- proc.time()
bootCI(USArrests$Assault,conf=0.99,B=1e6)
proc.time()-b 
# user system elapsed
# 39.443 0.471 40.189
# CI (137.949,201.587) at the level of 99% confidence

# Here we can see that our own function took much longer than the boot package's function for computing
# the confidence interval. 
  