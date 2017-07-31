# Andrew Le
# A11565521

# Problem 1 
# Part A

# This is a function that calculates the subsampling bootstrap mean confidence interval
subsample.mean.CI <- function(x, conf = 0.95, sub = length(x)/10 , B = 999){
  # These functions will calculate the mean, standard deviation
  xbar <- mean(x)
  sd <- sd(x)
  # Create a vector that will store the values that will be used later
  t <- numeric(B)
  # Specify the size when we sample from our data
  sub <- length(x)/10
  # Repeat the process B times, in this case 999 times
  for (b in 1:B){
    # Sample without replacement
    boot <- sample(x,sub, replace = FALSE)
    # Calculate the bootstrap mean and standard deviation 
    mean.b <- mean(boot)
    sd.boot <- sd(boot)
    # Compute the t-ration and store it
    t[b] <- (mean.b - xbar)/(sd.boot/sqrt(sub))
  }
  
  # Produce the subsample mean CI
  lCI <- xbar + quantile(t,(1-conf)/2) * (sd/sqrt(length(x)))
  rCI <- xbar + quantile(t,1-(1-conf)/2) * (sd/sqrt(length(x)))
  return(c(lCI,rCI))
}

# Quick test for the confidence interval
x <- rnorm(100,2,5)
subsample.mean.CI(x,conf=0.95, sub , B=999)

# Part B

# This is a function that will produce the bootstrap confidence interval that will be used in part b
bootCI <- function(x, conf, B){
  # Compute mean, standard deviation
  mean.x <- mean(x)
  sd.x <- sd(x)
  t <- numeric(B)
  n <- length(x)
  for (b in 1 : B){
    # sample with replacement
    boot <- sample(x, n, replace = TRUE)
    # Compute bootstrap mean and standard deviation
    mean.boot <- mean(boot)
    sd.boot <- sd(boot)
    # Calculate t-ration and store it 
    t[b] <- (mean.boot - mean.x)/(sd.boot/sqrt(n))
  }
  # Return the confidence interval
  ci <- mean.x + (sd.x/sqrt(n))*quantile(t,c((1-conf)/2,1-(1-conf)/2))
  return(ci)
}

x <- rnorm(100,2,5)
bootCI(x,conf=0.95,B=999)

# Sample size
n = 10000
# Repeat the process 1000 times
M = 1000
# Subsample size
u <- c(n/50,n/20,n/10,n/5,n/2)
# Create vectors that will store values that we will use later on 
t_res <- numeric(5)
avg_len <- numeric(5)
b <- numeric(5)
abl <- numeric(5)
  # Loop through the function
  for(i in 1:5){
    # These will help us store values later
    L1 <- 0
    L2 <- 0
    # Repeat the process M times
    for(j in 1:M){
      for(n in 1:length(n)){
      # Take on rnorm values
      x <- rnorm(u[i])
      # Compute the studet confidence interval
      t_ci <- t.test(x, conf.level = .95)$conf.int
      # Does a check to see if covers the population mean and length, then records it
      if( 0 > t_ci[1] && 0 < t_ci[2]){
        t_res[i] <- t_res[i] + 1
      }
      L1 <- L1 + t_ci[2] - t_ci[1]
      subsampCI <- subsample.mean.CI(x, conf = 0.95, sub = length(x)/10, B=1000)
      if( 0 > subsampCI[1] && 0 < subsampCI[2]){
        b[i] <- b[i] +1
      }
      L2 <- L2 + subsampCI[2] - subsampCI[1]
    }
    # Coverage averaged over the repeats and length
    avg_len[i] = L1/M
    abl[i] = L2/M
  }
}
# Makes a plot that will show the the coverage
t_p <- t_res*(1/M)
boot_p <- b*(1/M)
plot(x = u, y = t_p, col = 'orange', cex = 1.5, pch = 20,
     main = "Coverage Over Average",
     ylab = "Proportion of CI's containing the population mean",
     ylim = c(min(c(t_p,boot_p)),max(c(t_p,boot_p))))
points(x = u, y = boot_p, pch = 20, col = 'pink')
legend("center", pch = 20:15, col =c('orange','pink'), legend = c("Subsample Mean CI","Bootstrap CI"))
# length plot
plot(x = u, y = avg_len, col = 'maroon', cex = 1.5, pch = 20,
     main = "Length Over Average",
     ylab = "Average length of CI's",
     ylim = c(min(c(avg_len,abl)),max(c(avg_len,abl))))
lines(smooth.spline(x=u,y=avg_len, spar=0.05),col = 'black')
points(x = u, y = abl, col = 'navy', pch = 15, cex = 1)
lines(smooth.spline(x=u,y=abl, spar=0.02),col = 'black')
legend("topright", pch = 16:15, col =c('navy','maroon'), legend = c("Subsample Mean CI","Bootstrap CI"))

# I had difficulty running the part where it should have been multiple n's (1000...10000), so I did it for the biggest case
# and I observed that with a sample this big the Subsample Mean CI and Bootstrap CI were spot on with each other 
# for a large sample, I also observed with a smaller sample, the points were very mismatched with each other. In
# conclusion, but a larger sample, we were able to see much better results!

# Problem 2 

permSST <- function(y,g,B=999){
  
  #Create a null vector
  dat.vect = c()
  g <- c(1,2,3,4,5)
  #Loop through the matrix
  for(i in 1:dim(y)[2])
  {
    dat.vect = c(dat.vect,y[,i])
  }
  # Calculate the mean of our data vector 
  mean = mean(dat.vect)
  sst = 0
  sstb = numeric(B)
  for(i in 1:dim(y)[2])
  {
    #Calculate the sum of squares treatment
    sst = sst + dim(y)[1]*dim(y)[2]*(mean(y[,i])-mean)^2
  }
  # Repeat this process 1000 times
  for(i in 1:B)
  {
    # Sample from our vector values
    db = sample(dat.vect,replace = T)
    # Put our values into a dataframe 
    dat.boot = as.data.frame(matrix(db,nrow = dim(y)[1],ncol = dim(y)[2]))
    
    for(j in 1:length(g))
    {
      #Calculate the sum of squares treatement boot value
      sstb[i] = sstb[i] + dim(y)[1]*dim(y)[2]*(mean(dat.boot[,j]) - mean)^2
    }
  }
  # compute the p-value
  pval = (length(which(sstb >= sst)+1))/(B+1)
  return(pval)
}

y = c(0.2, -0.5, 1.1, 1.2, 0.1,-0.7,0.1,1.5,1.9,2.0,0.3,0.6,-0.8,1.7,-1.4,3.3,-0.1,1.6,1.9,1.1,-0.1,1.8,1,5,-2.4)
g = c(1, 2,3,4,5)
y <- as.data.frame(matrix(y,nrow=5))
permSST(y,g,B=999)

# Part B
# Produces a window to view the plots
quartz()
# Places all 9 plots in the quartz window
par(mfrow = c(3,3))
# The amount of repeats
M <- 200
# Store the values later on
l <- numeric(M)
# Store the one-way anova p vals
y = c(0.2, -0.5, 1.1, 1.2, 0.1,-0.7,0.1,1.5,1.9,2.0,0.3,0.6,-0.8,1.7,-1.4,3.3,-0.1,1.6,1.9,1.1,-0.1,1.8,1,5,-2.4)
# There are five carefully chosen values for Tau, that should put us further away from the alternative as we increase
Tau <- c(1,5,20,50,100)
for (j in c(2,5,10)){
  for (m in c(10,30,100)){
    #Repeat this M times
    for (n in 1:length(M)){
      # Generate data to use 
      mat <- matrix(rnorm(j*Tau,m),j,m)
      for (t in 1:j){
        mat[t,] <- sample(matrix(1,j),replace=TRUE)
      }
  
      l[n] <- permSST(mat,g,B=999)
    }
    plot(l/M, main = paste("Avg of p-values over M",j))
    
  }
}

# I was able to produce all the graphs, but I wasn't able to adjust the graphs properly as desired.

