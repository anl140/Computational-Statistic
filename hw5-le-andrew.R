# Andrew Le
# A11565521

# Problem 1
#Bootstrap between-groups of sum of squares test. This will take on a data matrix (dat) and compute all of the means
# and it shall calculate a bootstrap p-value

boot.oneway.test <- function(dat, B=999){
  
  #Create a null vector
  dat.vect = c()
  #Loop through the matrix
  for(i in 1:dim(dat)[2])
  {
    dat.vect = c(dat.vect,dat[,i])
  }
  # Calculate the mean of our data vector 
  mean = mean(dat.vect)
  sst = 0
  sstb = numeric(B)
  for(i in 1:dim(dat)[2])
  {
    #Calculate the sum of squares treatment
    sst = sst + dim(dat)[1]*dim(dat)[2]*(mean(dat[,i])-mean)^2
  }
  #Repeat this process 1000 times
  for(i in 1:B)
  {
    #Sample from our vector values
    db = sample(dat.vect,replace = T)
    #Put our values into a dataframe 
    dat.boot = as.data.frame(matrix(db,nrow = dim(dat)[1],ncol = dim(dat)[2]))
    for(j in 1:dim(dat)[2])
    {
      #Calculate the sum of squares treatement boot value
      sstb[i] = sstb[i] + dim(dat)[1]*dim(dat)[2]*(mean(dat.boot[,j]) - mean)^2
    }
  }
  # compute the p-value
  pval = (sum(sstb >= sst)+1)/(B+1)
  return(pval)
}
load("~/Downloads/smokers.rda")
smokers = smokers
boot.oneway.test(smokers,B=999)
# When the test is ran on the smokers data, we get out a p-value of 0.024

# Problem 2

dat = matrix(rnorm(100),10,10)
friedman.test(dat)$p.value 
#The p.value is 0.28 from a stimulated matrix thats 10 x 10

#Function to compute the G value under the null
fval <- function(I,J){
  #Loop through the vector of I
  for(i in 1:length(i)){
    #Loop through the vector of J
    for(j in 1:length(j)){
    #Generate a matrix with random normal values that i.i.d and take on I,J as the row and cols values
    dat <-matrix(rnorm(I[i]*J[j],mean=0,sd=1),I[i],J[j])
    #Rank the data values
    rvec <- rank(dat)
    #Compute the sum
    r.j <- (sum(rvec))
    #Compute the G test statistic of Friedman test
    G = (12/(I[i]*J[j]*(J[j]+1))) * (sum(r.j) - (3*I[i])*(J[j]+1))
    }
    return(G)
  } 
}

I = c(5,10,30,50,100,500)
J = c(50)
fval(I,J)


#Part B  

# Values to see what happens as I gets asympototically larger
I <- c(5,10,20,50)
J <- c(2,5,10)
M = 1e4
SampleDistributionFunction = numeric(M)
#Iterate through the vector
for(i in 1:length(I)){
  #Iterate through the J vector
  for (j in 1:length(J)){
    #Repeat the process M times
    for(m in 1:M){
    #Generate a matrix with random normal values that i.i.d and take on I,J as the row and cols values
    mx <-matrix(rnorm(I[i]*J[j],mean=0,sd=1),I[i],J[j]) 
    #Rank the matrix
    r.mx <- rank(mx)
    #Format the data to be in a matrix form in order to take the friedman test statistic
    test <- matrix(r.mx,I[i],J[j])
    SampleDistributionFunction[m] <- friedman.test(test)$p.value
    }
    #Plot the sample distribution function
    plot(ecdf(SampleDistributionFunction))
    #plot the line y=x 
    abline(0,1)
  }
  
}
#From looking at the graph we can see that the as we increase I, it gets closer to 1

# Problem 3

#Part A
# This is a test for a permutation of a friedman test 
perm.friedman.test <- function(dat,B=999){
  #Grab the data
  perm.d = dat
  #Transpose the data
  dat = t(dat)
  #Calculate the friedman test p_value
  f.stat <- friedman.test(dat)$p.value
  f.boot = numeric(B)
  
  #Allows the loop to run 999 times
  for(b in 1:B){
    #Sample from the data and permutes the colums
    permb.d = apply(dat,1,function(x){sample(dat,length(dat),replace=FALSE)})
    #Transpose the data matrix
    permb.d = t(permb.d)
    #Compute the friedman test p_val vector
    f.boot[b] <- friedman.test(permb.d)$p.value
  }
  #Calculate the p_value
  p.val = 1- (sum(f.boot >= f.stat))/(B+1)
  return(p.val)
}

dat = matrix(rnorm(25),5,5)
perm.friedman.test(dat,B=999)
# When this test is ran, I get out a p_value of 0.205 from my random normal generated data matrix

# Part B
#This is the same code from number 2, part b
I <- c(5,10,20,50)
J <- c(2,5,10)
M = 1e4
boot = numeric(M)
#Loop through the vector of I
for(i in 1:length(I)){
  #Loop through the vector of J
  for (j in 1:length(J)){
    #Repeat the process M times
    for(m in 1:M){
      #Generate a matrix with random normal values that i.i.d and take on I,J as the row and cols values
      mx <- matrix(rnorm(I[i]*J[j],mean=0,sd=1),I[i],J[j]) 
      #Rank the matrix values
      r.mx <- rank(mx)
      #Put the ranked vector back in a matrix
      dat <- matrix(r.mx,I[i],J[j])
      #Compute the test value
      boot[m] <- perm.friedman.test(dat)$p.value
    }
    plot(ecdf(boot))
    abline(0,1)
  }
  
}

#My Number 2, part b works fine but I'm not able to solve my issue with having a $ operator is invalid for atomic 
#vectors when trying to compute my test statistic. The code is supposed to plot out the sample distribution function
#like in number 2, part b. It should be making a plot to show if for the permutation friedman test value will tend
#to 1 if we increase our I values.
  