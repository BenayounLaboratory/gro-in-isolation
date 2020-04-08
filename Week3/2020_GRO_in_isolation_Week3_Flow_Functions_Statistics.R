##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-04-08

#######   Week 3 code snippets    #######

# comment lines that are there for humans reading the script start with "#"
# All other lines will be interpreted/executed by R


####################################################### 
###########        Flow control in R        ########### 
# Conditional flow/Branching in R 

x = 15

# simple branching
if (x > 0 ) {
  print ("x is positive")
} else {
  print ("x is negative")
}

# less simple branching
if (x > 0 ) {
  print ("x is positive")
} else if  (x < 0 ){
  print ("x is negative")
} else {
  print ("x is null")
}


# Looping in R: using ‘for’ 
for(i in 1:10) {
  print(i)
}

for(i in 1:10) {
  print(i * i)
}

# Looping in R: using ‘while’ 
i <- 0

while (i<= 10) {
  
  print(i)
  i <- i +2
  print(i)
  
}


####################################################### 
###########         Functions in R         ########### 
# simple function
addNumbers <- function(x, y){
  z <- x + y
  return (z)
}

addNumbers(5,3)           # positional argument matching
addNumbers(x = 5, y = 3)  # named argument matching
addNumbers(y = 3, x = 5)  # would also work (but isn't pretty)

# simple function with default values
addNumbers <- function(x = 2, y = 3){
  z <- x + y
  return (z)
}
addNumbers(5,3)           # default values are "overidden"
addNumbers()              # use default values


# a more complex function
firstletterUp <- function(mystring) {
  mystring <- tolower(mystring)
  substr(mystring, 1, 1) <- toupper(substr(mystring, 1, 1))
  return(mystring)
}

firstletterUp("SCREAMING")
firstletterUp("sCREAMING")
firstletterUp("sCReAMinG")


####################################################### 
###########     Statistical tests in R      ###########

##### Statistical testing: the t-test
set.seed(100)

# one sample t-test
x <- rnorm(50, mean = 10, sd = 0.5)         # sample 50 numbers from a normal distribution centered on 10 wth a sd of 0.5
y <- rnorm(50, mean = 20, sd = 0.5)         # sample 50 numbers from a normal distribution centered on 20 wth a sd of 0.5

t.test(x, mu=10)                            # result displayed on console, not significant
t.test(y, mu=10)                            # result displayed on console, significant

my.test <- t.test(y, mu=10)                 # result in a variable
my.test$p.value                             # exact p-value value


# two sample t-test
t.test(x, y)                                # result displayed on console, not significant

my.test2 <- t.test(x, y)                    # result in a variable
my.test2$p.value                            # exact p-value value


##### Statistical testing: the Wilcoxon/Mann-Whitney test
x <- c(11, 5, 6, 7, 9, 16, 18, 6)
y <- c(25, 17, 34, 56, 27,56) 

wilcox.test(x, mu = 100)                      # result displayed on console, significant

wilcox.test(x, y)                             # result displayed on console, significant

my.test3 <- wilcox.test(x, y)                 # result in a variable
my.test3$p.value                              # exact p-value value

##### Statistical testing: tails
# when you have a prior hypothesis on the direction of the difference
x <- c(11, 5, 6, 7, 9, 16, 18, 6)
y <- c(25, 17, 34, 56, 27,56) 

wilcox.test(x, y, alternative = "greater")   # testing whether y is "greater" than x
wilcox.test(x, y, alternative = "less")      # testing whether y is "less" than x
wilcox.test(x, y, alternative = "two.sided") # testing whether y is "different" from x

# g for greater #=> Wilcoxon rank sum test
#=> #=> data: x and y #=> W = 35, 
p-value = 0.1272 #=> alternative hypothesis: true location shift is greater than 0


##### Statistical testing: correlation
data(iris)
cor(iris$Sepal.Length,iris$Petal.Length)    # defaults to Pearson's correlation
cor(iris$Sepal.Length,iris$Petal.Length, method = "pearson")    # same as above, but explicit
cor(iris$Sepal.Length,iris$Petal.Length, method = "spearman")   # use Spearman's rank correlation

cor.test(iris$Sepal.Length,iris$Petal.Length, method = "pearson")    # is Pearson's correlation significant
cor.test(iris$Sepal.Length,iris$Petal.Length, method = "spearman")   # is Spearman's rank correlation significant

my.test4 <- cor.test(iris$Sepal.Length,iris$Petal.Length, method = "spearman")  # save test result in variable
my.test4$p.value                                                                # get value of p-value 


##### Statistical testing: significance of overlap, or differences in proportions
fisher.test(matrix(c(200,300,500, 15500),2,2), alternative = "greater")



##### Statistical testing: ANOVA
data(iris)

# Box plot
boxplot(Sepal.Length ~ Species, 
        data = iris,
        ylab = "Sepal Length", 
        xlab = "", 
        col = c("deepskyblue", "orchid", "coral"),
        las = 2)


## classical parametric ANOVA
my.aov.fit <- aov(Sepal.Length ~ Species, data = iris) # ANOVA model

plot(my.aov.fit)                                      # diagnotic plots for the fit
summary(my.aov.fit)                                   # display Type I ANOVA table
TukeyHSD(my.aov.fit)                                  # post-hoc test

## non-parametric ANOVA
# Kruskal Wallis Test One Way Anova by Ranks
kruskal.test(Sepal.Length ~ Species, data = iris)    # response is numeric and is modeled on a factor variable


##### Statistical testing: survival
install.packages(c('survminer','survminer',"BiocManager"))  ### package installation commands (necessary for survival analysis)
BiocManager::install("RTCGA.clinical") # data for examples  ### package installation commands (necessary for survival analysis)

library(RTCGA.clinical)  # Load necessary packages (more on that next week)
library(survival)        # Load necessary packages (more on that next week)
library(survminer)       # Load necessary packages (more on that next week)

BRCAOV.survInfo <- survivalTCGA(BRCA.clinical, OV.clinical, extract.cols = "admin.disease_code") # get cancer survival data
head(BRCAOV.survInfo)

# perform Kaplan Meier curve calculation
cancer.fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code, data = BRCAOV.survInfo)

# plot Kaplan Meier curve
plot(cancer.fit, col = c("red","blue"), xlab = "Time (days)", ylab = "Survival")
legend("topright",c("brca","ov"), col = c("red","blue"), pch = "_", bty = 'n')

# log-rank test for differences in survival
surv_diff <- survdiff(Surv(times, patient.vital_status) ~ admin.disease_code, data = BRCAOV.survInfo)
surv_diff

# A 1-step visualization and test survminer
ggsurvplot(cancer.fit, data = BRCAOV.survInfo, pval = TRUE)   ## show p-value of log-rank test.


####################################################### 
###########       Power Analysis in R       ###########
install.packages('pwr')  ### package installation commands (necessary for survival analysis)
library(pwr)             # Load necessary packages (more on that next week)

# Power analysis - t-test
pwr.t.test(n           = 25,                    # calculate  power
           d           = 0.75, 
           sig.level   = 0.01, 
           alternative = "greater")  

pwr.t.test(power       = 0.6,                    # calculate  n
           d           = 0.75, 
           sig.level   = 0.01, 
           alternative = "greater")  


# Power analysis - ANOVA
pwr.anova.test(k         = 5,                    # calculate  n
               f         = 0.25,
               sig.level = 0.05,
               power     = 0.8)


