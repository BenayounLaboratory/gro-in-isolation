##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-03-25

#######   Week 1 code snippets    #######

# comment lines that are there for humans reading the script start with "#"
# All other lines will be interpreted/executed by R


####### getting help
?boxplot       # call help page for the function 'boxplot'
??boxplot      # search help pages for the word 'boxplot'


####### Assigning values to variables
# assign value to a numerical variable
x <- 12.6      # Best practice in R
x = 12.6       # accepted, but not optimal
12.6 -> x      # equivalent, but clunky

# assign a vector of values
y <- c(3, 7, 9, 11)

# ‘:’ means the series of integers between
a <- 1:6

# create a vector of equally spaced values
b <- seq(0.5, 0, -0.1) 


####### R data types
# How to verify the type of data in R:
class(x)

# How to change the type of data in R:
as.numeric(x)
as.character(x)

# Factors in R
data(iris)                   # load the 'iris' dataset
class(iris$Species)          # "factor"
iris$Species                 # the levels of the factor variable are also displayed on screen


# Missing values in R
1 + NA                      # NA
max(c(NA, 4, 7))            # NA
max(c(NA, 4, 7), na.rm=T)   # 7



####### R data structures
# Vectors and matrices in R
a = c(1,2,3)
a*2

# data frame
data(iris)       # load bundled data into memory, part of the "datasets" base R package
head(iris)       # by default, head shows the first 10 lines of a dataframe or matrix
View (iris)

# lists and accessing lists
doe <- list(name="john", age=28, married=F)   # mixed type of data

doe$name                                    
doe$age
doe$married                                    


####### Subsetting data structures in R
# subsetting vectors
my.vector <- c(6,2,25,3.4,5)
my.vector[3]

# Subsetting a matrix or data frame in R
head(iris)                  # by default, head shows the first 10 lines of a dataframe or matrix
iris[3, 2]                  # shows the element on line 3, column 2
iris["3", "Sepal.Width"]    # shows the element on line named "3", column named "Sepal.Width"
iris[,"Sepal.Width"]        # shows the column named "Sepal.Width"

my.iris <- head(iris)       # by default, head shows the first 10 lines of a dataframe or matrix
my.iris[c(1,3),]            # extract lines 1 and 3 of the iris data
my.iris[c(T,F,T,F,F,T),]    # use a TRUE/FALSE to subset: TRUE is included, FALSE is not included


my.iris$Petal.Length                    # numeric vector: the column with data on Petal length
my.iris$Petal.Length > 1.4              # logical vector: a vector of TRUE or FALSE depending on whether Petal length
my.iris[my.iris$Petal.Length > 1.4, ]   # subsetted dataframe, where only the lines with Petal length > 1.4 are included


####### Apply an operation to several lines or columns in R
# operations on a vector
my.vector <- c(6,2,25,3.4,5)
my.vector * 3
my.vector + 1
mean(my.vector)
max(my.vector)

# operations on a matrix, data frame or list
my.iris <- head(iris)[,1:4] # do not take species column (in the 5th column)

# sum per row
apply(my.iris, 1, sum)

# sum per column
apply(my.iris, 2, sum)


####### Reading and writing data in R
# read data from a file (here, read the Puromycin data that I previously exported myself)
read.table("R_Puromycin_dataset_export.txt",  # name of the input file (and relative position of the file)
           header = TRUE,                     # is there header information
           sep = "\t")                        # print column names (variable names)


# good option at the top of a script to always include to avoid bad surprises when reading in text files in R
options(stringsAsFactors = F)

# write data to a text file
write.table(iris, 
            file = "iris_dataset_export.txt",      # name of the output file
            quote = FALSE,                         # name of the output file
            sep = "\t",                            # separate fields with tabulation
            row.names = FALSE,                     # don't print row names (here, they are just numbers)
            col.names = TRUE)                      # print column names (variable names)

# create a computer generated time stamp for the file name
Sys.Date()                                               # date from your computer as year-month-day
paste(Sys.Date(), "iris_dataset_export.txt", sep = "_")  # save the date of your analysis
paste0(Sys.Date(), "_iris_dataset_export.txt")           # alternative syntax


# save data to a R-formatted file
save(iris, file="2020-03-25_iris.Rdata")   # RData or Rda are the accepted extensions

# load back R data in memory from file
load("2020-03-25_iris.Rdata")


######## At home, try testing functions with other R datasets instead of "iris"
data(chickwts)
data(DNase)
data(esoph)
data(InsectSprays)
data(Puromycin)
data(ToothGrowth)


