# PROBLEM 1
# Part A
#There could be some rounding error, thus it results in the total not being 100. As we would have exp
#Part B
c <- c(47,52,2)
barplot(c,main = "Survey On Who Should Fund Pre-Kindergarten Programs",col = c("yellow","maroon","pink"),xlab="Federal or State", 
names = c("Federal Government","Each State","Unsure/No Answer"),ylab="Approval Rate(%)",legend = c("Federal Government","Each State",
"Unsure/No Answer"),args.legend = list(title = "SES", x = "topright", cex = .85), ylim = c(0, 100))
a <- c(47,52,2)
color <- c("cyan","black","red")
portions <- round(a/sum(a)*100,1)
portions <- paste(portions, "%", spe=" ")
pie(a,main = "Survey On Who Should Fund Pre-Kindergarten Programs", col=color, labels = portions, cex=0.85)
legend("topright", c("Federal Government","Each State","Unsure/No Answer"), cex=0.95, fill=color)
#Part C
#Formalize the hypothesis test:
#H0: p1 = p2 = 1/2
#H1: p1 != p2
chisq.test(c)
#Due to having a p-value that is significantly less than 0.05, we are not able to accept the null hypothesis
#therefore, we are able to conclude that there is a relationship between people's preference.
# Part D
#The assumptions that my conclusion relies on are that after seeing the p-value from Chi-squared test that the
#p-value is no significant enough to make a decision on rejecting the null hypothesis and from their test of 
#independence.

# PROBLEM 2 
S <- 5 
ratiovals <- function(S,N)
  {
  m = c(10,20, 50, 100, 200, 500,1000,2000,5000) # The number of samples that we have to use from 
  for(i in m) # Loop through m vetor of values
{
    for(i in 10^4) # Go through 10000 samples
      {
         LR = 0 # Intialize the vector
         c = numeric(S) # Formulate a vector column 
         size <- 1:S # Use this sample size
         x <- sample(size, N, replace =TRUE) #Take a sample from this with replacement
          for(i in 1:S) # Loop through all of the S
         {
            S <- S # intialize S
            c[i] <- sum(x==i) #Place values into this
            if(c[i] != 0) # If the values are not even equal to 0, then run the loop
            {
            LR <- LR + 2*c[i]*log(c[i]/(N/S)) #Calculate the likelihood ratio value 
          }
        }
  
  return(LR) # Return the value 
   }
  }
}

par(mfrow=c(3,3)) #Put into a plot that is 3x3
M = c(10,20, 50, 100, 200, 500,1000,2000,5000)
for (m in M)
  {
  curve(dchisq(x, m), 0, 3*m, n = 1e4, lwd=2, col=4)
  title(sprintf('d.o.f = %d', round(m)))
}

# This takes account of all of the values that must be used to make a bootstrap via a monte carlo simulation that makes
#many values. And then we also have values from that plots. 

# PROBLEM 3
# Part A
library("gmodels", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
admissions <- fread('http://media.pearsoncmg.com/aw/aw_deveaux_stats_2/activstats/text/Ch03_Magnet_schools.txt') #Grab data
tab = table(admissions) #Get all of the counts
tab/sum(tab) #Look at their percentages
admissions = as.factor(admissions)
summary(admissions)
CrossTable(admissions$Ethnicity, admissions$`Admission Decision`,prop.t=TRUE, prop.r=TRUE, prop.c=TRUE) #Make a contigency table
# From the contingency table, I discovered that Black/Hispanic had the highest admission rate, while more White 
# people got turned away from the school.
# Part B
# Formalize the hypothesis test:
# H0 = there is a relationship between admission and ethnicity 
# H1 = this is not the case (admissions and ethnicity have no correlation for determing admissions)
ctab = cbind(tab[,"Accepted"], tab[,"Turned Away"]) #Look at only Accepted and Turn Away values, ignore Wait-listed
ctab
# Use the chisquare test
chisq.test(ctab, sim=TRUE, B = 1e4) #Perform the test on it
#Here, the p-value is quite small, thus it shows that we are able to reject the null hypothesis without a doubt.

