# Andrew Le
# A11565521

# Problem 1 
# This function function takes in vector and the according side of the alternative hypothesis
perm.ks.test <- function(x, y, alternative = c("two.sided","less","greater"), B=999) {
#B value for Monte Carlo simulations later
B = 999
D.b = numeric(B)
#Compute ks test value from the original data that we are given, (this is using the cloudseeding data)
D = ks.test(x,y)$stat

  #Go through a loop to make multiple test statistics using the ks test
  for (b in 1:B){
   #Sample through Z
   z = sample(c(x,y),(length(x)+length(y)),replace=FALSE)
   #This is the new x vector of vals
   x.s <- z[1:length(x)]
   #This is the new y vector of vals
   y.s <- z[length(x)+1:(length(x)+length(y))]
   #Compute the test stat, 999 times from the B
   D.b[b] <- ks.test(x.s,y.s)$stat
  }
    #Go to the appropriate alternative test and calculate the p-value from there
    if(alternative == "two.sided"){
    p.val = (sum(abs(D.b) >= abs(D)) +1)/(B+1)
    }
    
    if(alternative == "less"){
      p.val = (sum(D.b >= D)+ 1)/ (B+1)
    }
    
    else if(alternative == "greater"){
      p.val = (sum(D.b <= D)+1)/ (B+1)
    }
    
  return(p.val)
}

#the cloudseeding file was in my downloads 
load("~/Downloads/cloudseeding.rda")
dat = cloudseeding
str(dat)
x = dat$seeded
y = dat$unseeded

perm.ks.test(x,y,alternative="two.sided") #p-value is 0.021
perm.ks.test(x,y,alternative="less") #p-value is 0.013
perm.ks.test(x,y,alternative="greater") #p-value is 0.993

# Problem 2

pattern.ks.stat <- function(p){
  #Take on the value of the length of the pattern p that is being used
  z = 1:length(p)
  #Condition to check if 0 is in x
  x = z[p == 0]
  #Condition to check if 1 is in y
  y = z[p == 1]
  
  #loop through the length of p
  for (i in 1:length(p)){
    #if i is equal to 0, store it in x
    if(i == 0){
      x = c(x,z[i])
    }
    #if i isn't equal to 0, then store it in y
    if(i != 0){
      y = c(y,z[i])
    }
  }
 
  x = x[-1]
  y = y[-1]
  #Perform and compute all of the alternative hypothesis sides
  a <- ks.test(x,y,alternative="less")$stat
  b <- ks.test(x,y,alternative="two.sided")$stat
  c <- ks.test(x,y,alternative="greater")$stat
  ks.vals <- c(a,b,c)
  #Return a vector with all of the values 
  return(ks.vals)
}

p <- c(0,0,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,0)
pattern.ks.stat(p)
# From this test, we can observe some very small p-vals from all different alternatives of the test. I believe the 
# purpose of this problem was allow us to see that the ks test follows very closely to a ranking system as we 
# can see from this problem


# Problem 3
# This is a seperate function from the qqtest in order to compute the quantiles for a given x and y that we have
qq.s <-function(x,y,alternative = c("two.sided","less","greater")){
  
  qx = quantile(x,probs=(1:9)/10)
  qy = quantile(y,probs=(1:9)/10) 
  #Find the particular hypothesis that we want to test by going through the conditions
  if(alternative=="two.sided"){
    D = max(abs(qx-qy))
  }
  
  if(alternative == "greater"){
    D = min(qx-qy)
  }
  
  else if(alternative =="less"){
    D = max(qx-qy)
  }
  
  return(D)
} 

qqtest <- function(x,y, alternative = c("two.sided","less","greater"), B=999){
#This is the same exact thing like in number 1, except we must gather our test statistic from using
#the quantile values to gather our test statistic.
D.b = numeric(B)
qq.s(x,y,alternative = "two.sided")
for (b in 1:B){
  z = sample(c(x,y),(length(x)+length(y)),replace=FALSE)
  x.s <- z[1:length(x)]
  y.s <- z[length(x)+1:(length(x)+length(y))]
  D.b[b] <- qq.s(x.s,y.s)
}
  if(alternative == "two.sided"){
    p.val = (sum(abs(D.b) >= abs(D)) +1)/(B+1)
  }

  if(alternative == "less"){
   p.val = (sum(D.b >= D)+ 1)/ (B+1)
  }

  else if(alternative == "greater"){
    p.val = (sum(D.b <= D)+1)/ (B+1)
  }
return(p.val)
}
#This code isn't debugged properly, but it is supposed to take the x,y values and desired alternative and compute 
#the p-value from there. The concept is still all the same but instead this time, the D value is taken from the 
#quantiles, while the first problem use the ks.test statistic. Every methodology is exactly the same from the 
#previous problem above.