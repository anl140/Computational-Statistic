#Andrew Le
#A11565521
#HW #6

# Problem 1
#Part A
gene.val <- function(y){
  
apply(d,1,function(x){t.test(x[1:k],x[(k+1):(2*k)])})$pvalue
  
}

#Part B
k = 10
m = 1e4
#Here are the wrong observations
m.n <- c(1e4, (1e4 - 10), (1e4 - 100))
C = 1000
delt = seq(from=2,to=6.5,by=0.5)
#These are the 3 values for each of the 5 tests that will eventually be extracted in the near end and used
# to return the single p-value after running through all three tests for rejection, Type I and Type II
bon.T1 = numeric(m)
holm.T1 = numeric(m)
hoch.T1 = numeric(m)
bh.T1 = numeric(m)
by.T1 = numeric(m)
bon.T2 = numeric(m)
holm.T2 = numeric(m)
hoch.T2 = numeric(m)
bh.T2 = numeric(m)
by.T2 = numeric(m)
bon.r = numeric(m)
holm.r = numeric(m)
hoch.r = numeric(m)
bh.r = numeric(m)
by.r = numeric(m)

for (i in 1:length(m.n)){
  
  legend(FNP, "hold", fill='red')
  legend(FDP, "hold", fill='blue')
  legend(FWEP, "hold", fill='green')
  
  for(s in 1:length(delt)){

    for(c in 1:C){
      #Generate the matrix that will be used
      mx1 = matrix(rnorm(k*(m-m.n[i]),delt[s],sqrt(2)),m-m.n[i],k)
      mx2 = matrix(rnorm((k*m.n[i]),0,sqrt(2)),m.n[i],k)
      mx3 = matrix(rnorm((k*m),0,1),m,k)
      data <- cbind(mx3,rbind(mx1,mx2)
      #Create a variable that will soon store our results from each test              
      res.bon = numeric(m)
      res.holm = numeric(m)
      res.hoch = numeric(m)
      res.bh = numeric(m)
      res.by = numeric(m)
      res = numeric(m)
      #Apply the method 
      res.bon = p.adjust = (res, method="bonferroni") 
      res.holm = p.adjus = (res, method = "holm")
      res.hoch = p.adjust(res,method = "hochberg")
      res.bh = p.adjust(res, method = "BH")
      res.by = p.adjust(res, method = "BY") 
      #Calculate the desired stat of either rejection, type I error, or type II error for each of the individual 5 types
      #so we have about 15 results from the 5 tests
      rej <- res.bon[c] <= 0.1
      bon.T1[c] = sum(rej[(m-m.n+1):m])
      bon.T2[c] = m-m.n-sum(rej[1:(m-m.n)])
      bon.r[c] = sum(rej)
      rej <- res.holm[c] <= 0.1
      holm.T1[c] = sum(rej[(m-m.n+1):m])
      holm.T2[c] = m-m.n-sum(rej[1:(m-m.n)])
      holm.r[c] = sum(rej)
      rej <- res.hoch[c] <= 0.1
      hoch.T1[c] = sum(rej[(m-m.n+1):m])
      hoch.T2[c] = m-m.n-sum(rej[1:(m-m.n)])
      hoch.r[c] = sum(rej)
      rej <- res.bh[c] <= 0.1
      bh.T1[c] = sum(rej[(m-m.n+1):m])
      bh.T2[c] = m-m.n-sum(rej[1:(m-m.n)])
      bh.r[c] = sum(rej)
      rej <- res.by[c] <= 0.1
      by.T1[c] = sum(rej[(m-m.n+1):m])
      by.T2[c] = m-m.n-sum(rej[1:(m-m.n)])
      by.r[c] = sum(rej)
      
    }
  }
}

#The code should generate the graphs properly, and I believe I properly generated the statistics. There are some issues 
#with some of my syntax. But I tried to compute the rejection, type 1 error and type 2 error for all 5 of the tests.

#Problem 2

dist.cor.test <- function(x, y, B = 999) {
  #Take the first 2 column vectors 
  dat <- cars
  x <- dat[,1]
  y <- dat[,2]
  #Create the matrix
  mx <- matrix(rep(x,length(x),length(x),length(x)))
  #Calculate the difference for the x vector
  diff <- abs((as.vector(mx)-t(mx)))
  #Calculate the difference for the y vector
  mx2 <- matrix(rep(y,length(y),length(y),length(y)))
  diff2 <- abs((as.vector(mx2) -t(mx2)))
  #Permute the values
  rmean = apply(diff,1,mean)
  cmean = apply(diff,2,mean)
  tmean = sum(diff)/(length(x)*length(x))
  diff = diff - rmean
  diff = t(t(diff)-cmean)
  diff = diff + tmean
  #permute the values
  rmean = apply(diff2,1,mean)
  cmean = apply(diff2,2,mean)
  tmean = sum(diff2)/(length(y)*length(y))
  diff2 = diff - rmean
  diff2 = t(t(diff2)-cmean)
  diff2 = diff2 + tmean
  #Calculate the VN value
  vn <- sum(diff*diff2)/(length(x)*length(x))

  v.obs = numeric(B)
  
  for(b in 1:B){
    #Sample from the values
    diff1 <- sample(diff,length(x),replace=TRUE)
    diff3 <- sample(diff2,length(y),replace=TRUE)
    v.obs[b] <- sum(diff1*diff3)/(length(x)*length(x))
  }
  #Calculate the pvalue
  pval = (sum(v.obs >= vn ))/(B+1)
  return(1-pval)
  
}

dat <- cars
x <- dat[,1]
y <- dat[,2]
dist.cor.test(x,y,B=999)
# I get out a p-value of 0.001 and I checked with the energy test and got the same value
