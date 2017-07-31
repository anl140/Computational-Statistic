# Andrew Le
# A11565521
# HW 3

#Problem 1

# This code was taken from the solutions given to us in order to prevent any possible errors in my code.

#The function bootCI produces a bootstrap Studentized confidence interval for the mean
bootCI <- function(x, conf, B){
  #calculate the mean and standard deviation of the original sample
  mean.x <- mean(x)
  sd.x <- sd(x)
  #create a vector to store the value of t test statistic
  t <- numeric(B)
  n <- length(x)
  #use a for loop to bootstrap for B times
  for (b in 1 : B){
    #generate the new sample from the original one
    boot <- sample(x, n, replace = TRUE)
    mean.boot <- mean(boot)
    sd.boot <- sd(boot)
    #calculate the t test statistic
    t[b] <- (mean.boot - mean.x)/(sd.boot/sqrt(n))
  }
  #calculate the confidence interval
  ci <- mean.x + (sd.x/sqrt(n))*quantile(t,c((1-conf)/2,1-(1-conf)/2))
  return(ci)
}

x <- rnorm(1000,mean=0,sd=1)
conf <- 0.99
B <- 1e4
bootCI(x,0.99,B=1e4)
# 0.5% = -0.06696779 99.5% = 0.09732120

CI1 <- rep(0,1000)
CI2 <- rep(0,1000)
count <- 0
D <- seq(from=0,to=100,by=10)

for (i in D){
  M <- 1000
  out <- numeric(M)
  for (m in 1:M){
    x <- rnorm(n=100,mean=0,sd=1)
    out <- t.test(x,conf.level=0.99)
    #Check and see if the confidence interval is equal to one and store it in CI1
    if(out[m] == 1){
      CI1 <- out[m]
    }
    #Check and see if the confidence interval is equal to zero and store it in CI2
    else if(out[m] == 0){
      CI2 <- out[m]
    }
  }
  
  pl1 <- sum(CI1)/1000
  pl2 <- sum(CI2)/1000

}

# I believe the issue I had in this problem was that I wasn't able to properly plot all of the graphs as desired. 
# I understand that we had to loop this twice, one through a sequence that incremented gradually, and this had
# to be repeated for M= 1000 times. Part A grabbed values from the standard normal while Part B used exponential values.
# We are supposed to compute the confidence interval for both part and have two seperate confidence intervals where
# when the Conf Int equalled one, it was stored into the CI1 , and if it was 0 it was stored in CI2. I was unable to 
# get the graph part correctly.

# Part B
CI1 <- rep(0,1000)
CI2 <- rep(0,1000)
count <- 0
D <- seq(from=0,to=100,by=10)

for (i in D){
  M <- 1000
  out <- numeric(M)
  for (m in 1:M){
    x <- rexp(1000)
    out <- t.test(x,conf.level=0.99)
    #Check and see if the confidence interval is equal to one and store it in CI1
    if(out[m] == 1){
      CI1 <- out[m]
    }
    #Check and see if the confidence interval is equal to zero and store it in CI2
    else if(out[m] == 0){
      CI2 <- out[m]
    }
  }
  
  pl1 <- sum(CI1)/1000
  pl2 <- sum(CI2)/1000
  
}

plot(pl1,pl2)


#Problem 2
# I intend to perform some Monte Carlo tests, in order to observe what happens to as n  becomes large. 
# I assume from the description that we should be able to see approximately a normal curve. I suspect
# that this comes from the idea of having a multitude of samples to observe this from. I expect these 
# results to hold true because running Monte Carlo tests from the underlying distribution will be of 
# more use than random sample data via bootstrap.

theta <- 0
B <- 1e4
Median <- numeric(B)
N <- seq(from=1, to=900, by=100) #n-value incremements 
# plot the histograms of the values into a 3x3 in order to see the effects as N increases
par(mfrow=c(3,3)) 
for (i in N){
  for( b in 1:B){
    x_b <- rnorm(1000,mean=0,sd=1) # sample data from the normal distribution
    Median[b] <- median(x_b,na.rm=FALSE) # collect the median values
  }
  hist(Median,freq=FALSE) # plot out the histogram
  var <- dnorm(Median,mean=0,sd=1, log=FALSE)
  # perform the test to check it out
  test <-ks.test(Median,"pnorm",mean=0,sd=sqrt(var))
}


# Problem 3
#Part A
# Please download and load the quake_times file if you don't have it already.
load("~/Downloads/quake_times.rda")
test <- ks.test(quake_times, "punif", min=0, max=1)
# One-sample Kolmogorov-Smirnov test
# D = 0.037652, p-value = 0.01761

#Part B
#Get the length of the data, and use it as the parameter for the loop
n <- length(quake_times)
# Arrange the observations in order
sqt <- sort(quake_times)
w_obs <- 1:n
for (i in 1:n){
  # Record the values into this vector
  w_obs[i] <- abs(sqt[i] - i/n)/sqrt(sqt[i]*(1 - sqt[i]))
}
# Determine the largest value from that vector
w_val <- max(w_obs) 
w_val 
# test stat from data is = 0.08286188 

# Anderson-Darling
B <- 5000
w_vec <- numeric(B)
w_max <- numeric(B)

for (b in 1:B){
  # Simulate values from uniform distribution
  x_b <- runif(n=n, min=0, max=1)
  # Sort the boot oberservations
  x_b <- sort(x_b) 
  # Compute test statistics from sample
  for (i in 1:n){
    w_vec[i] <- abs(x_b[i] - i/n)/sqrt(x_b[i] * (1 - x_b[i]))
  }
  # Get Max value of W over i
  w_max[b] <- max(w_vec)
}

pval <- (sum(w_max <= w_val)+1)/(B+1) 
pval