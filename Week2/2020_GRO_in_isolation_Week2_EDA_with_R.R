##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-03-31

#######   Week 2 code snippets    #######

# comment lines that are there for humans reading the script start with "#"
# All other lines will be interpreted/executed by R


################# Data cleanup in R #################
####### Dealing with NA values
x <- c(3, NA, 4, NA, NA)
is.na(x[2])
is.na(x)

data(airquality)
is.na(airquality$Ozone)     # logical vector
!is.na(airquality$Ozone)    # ! means "not", so it flips the logical vector

####### Cleaning data in R
data(iris)

duplicated(iris)              # Identify duplicate rows
iris[!duplicated(iris), ]     # Remove duplicate rows

as.numeric()
as.character()


####### Data manipulation in R
# create/engineer new variables
data(iris)
iris$sepal.Surface <- iris$Sepal.Length * iris$Sepal.Width

# Summarize information over a common denominator (e.g. species, gene, sex, etc.)
iris.per.species <- aggregate(iris[,1:4],  by = list(iris$Species), FUN = mean)
iris.per.species

# Aggregating data
table(iris$Species)
with(iris, table(Species))

# Transposing a matrix
head(iris)
head( t(iris) )


################# Exploratory data analysis and base graphics in R #################
####### Overview of your data
data(iris)
summary(iris)             # information by columns

mean(iris$Sepal.Length)
var(iris$Sepal.Length)

#### Quantiles and percentiles
data(iris)
quantile(iris$Sepal.Length)


####### Descriptive statistics: barplot
data(state) 
state.x77 <- data.frame(state.x77)                                        # convert to data frame for ease of manipulation
barplot(state.x77$Population, names.arg = rownames(state.x77), las = 2)   # barplot of values


####### Descriptive statistics: boxplot
set.seed(100)
x <- rnorm(100, mean=0, sd=1)   # sample from the normal distribution, mean = 0 and sd is 1
boxplot(x)

data(iris)
boxplot(iris[,-5])             # plot only numerical columns (not species)

####### Descriptive statistics: histogram
data(iris)
hist(iris$Sepal.Length)                        # Basic histogram (on frequency)
hist(iris$Sepal.Length, breaks = 15)           # set number of bins

hist(iris$Sepal.Length, breaks = 15, freq = F) # Basic histogram (on density ratehr than frequency)
rug(jitter(iris$Sepal.Length))                 # add rug and density overlayed
lines(density(iris$Sepal.Length))


####### Descriptive statistics: plotting a data frame
# pairwise scatterplots
data(iris)
plot(iris)

data(mtcars)
plot(mtcars)

####### Descriptive statistics: scatterplots
data(iris)
plot(iris$Sepal.Length, iris$Petal.Length,
     xlab = "Sepal Length",
     ylab = "Petal Length",
     pch = 16,
     cex = 2,
     col = "tomato"
     )

# is there a correlation between values?
cor.test(iris$Sepal.Length,
          iris$Petal.Length,
          method = 'pearson')

cor.test(iris$Sepal.Length,
	  iris$Petal.Length,
	  method = 'spearman')

####### Descriptive statistics: linear models
data(iris)
my.fit <- lm(Sepal.Length ~ Petal.Length, data = iris) # regress sepal length on petal length
summary(my.fit)      # coefficients and significance of coefficients
plot(my.fit)         # diagnotics plots, identify points with leverage. Most plots should be linear. Outliers are noted


data(mtcars)
my.fit <- lm(mpg ~ cyl, data = mtcars) # regress sepal length on petal length
summary(my.fit)      # coefficients and significance of coefficients
plot(my.fit)         # diagnotics plots, identify points with leverage. Most plots should be linear. Outliers are noted


####### Advanced descriptive statistics: PCA
data(iris)
my.pca <- prcomp(t(iris[,-5]),scale = TRUE)    # calculate principal components
x <- my.pca$x[,1]                              # 1st coordinates
y <- my.pca$x[,2]                              # 2nd coordinates

my.summary <- summary(my.pca)
my.summary                                     # Importance of components:


plot(x,y, pch = 16, cex=2,                     # plot PCA resuts
	col= "tomato",
	xlab = 'PC1‘,
	ylab = 'PC2') 

plot(x,y, pch = 16, cex=2,                     # plot PCA resuts, automatically color by species
     col= iris$Species,
     xlab = 'PC1',
     ylab = 'PC2') 


#################  R base graphics: Make it pretty(-er) #################
####### histogram beautification
data(iris)
hist(iris$Sepal.Length)                                 # basic mode
            
hist(iris$Sepal.Length,                                 # add axes names, title
     main = 'Iris Sepal Length',            
     xlab = 'Iris Sepal Length',            
     ylab = 'Frenquency')            
            
hist(iris$Sepal.Length,                                 # increase number of bins, add a color
     main = 'Iris Sepal Length',
     xlab = 'Iris Sepal Length',
     ylab = 'Frenquency',
     col = "pink",
     breaks = 40)
     
####### boxplot beautification
data(state)
state.x77 <- data.frame(state.x77)

boxplot(state.x77$Income)                                # basic mode

boxplot(state.x77$Income,                                # add axes names and ranges, title, color
        main = 'Income per capita in US states (1977)',
        ylab = 'Dollars',
        col = "pink",
        ylim = c(0,8000))

boxplot(state.x77$Income,                                # do not display outliers
        main = 'Income per capita in US states (1977)',
        ylab = 'Dollars',
        col = "pink",
        ylim = c(0,8000),
        outline = F)

# Group by factor variable for display
data(iris)
boxplot(Sepal.Length ~ Species, data = iris,
        main = 'Sepal length by species',
        ylab = 'cm',
        col = c('deeppink', 'deepskyblue','gold'),
        ylim = c(0,10))

####### pie chart beautification
my.iris.data <- table(iris$Species)

# Simple Pie Chart
pie(my.iris.data, labels = names(my.iris.data), 
    main="Pie Chart of species in iris data")

pie(my.iris.data, labels = names(my.iris.data),        # customize colors
	main="Pie Chart of species in iris data",
	col=c('deeppink', 'deepskyblue','gold'))


####### scatterplot beautification
data(iris)
plot(iris$Sepal.Length, iris$Petal.Length)

plot(iris$Sepal.Length, iris$Petal.Length,              # add axes names and ranges, title, color
     col = "deeppink",
     xlim = c(0,10), ylim = c(0,10),
     xlab = "Sepal Length (cm)",
     ylab = "Petal Length (cm)",
     pch = 16
)

smoothScatter(iris$Sepal.Length, iris$Petal.Length,    # give information about data density on the scatter
     xlim = c(0,10), ylim = c(0,10),
     xlab = "Sepal Length (cm)",
     ylab = "Petal Length (cm)")


####### Base Graphics: line plots
data(airquality)
plot(1:dim(airquality)[1], airquality$Temp,
     col = "deeppink",
     xlim = c(0,dim(airquality)[1]), ylim = c(50,100),
     xlab = "Data point", ylab = "Temperature",
     type = 'l')


legend("topleft", c("temperature"),                    # add a legend
       col = "deeppink", pch = "_", bty = 'n')


####### Base Graphics: creating panels
par(mfrow=c(2,2))
hist(iris$Sepal.Length, 
     main = 'Iris Sepal Length',
     xlab = 'Iris Sepal Length',
     ylab = 'Frenquency',
     col = "pink",
     breaks = 40)
boxplot(Sepal.Length ~ Species, data = iris,
        main = 'Sepal length by species',
        ylab = 'cm',
        col = c('deeppink', 'deepskyblue','gold'),
        ylim = c(0,10))
plot(iris$Sepal.Length, iris$Petal.Length,
     col = "deeppink",
     xlim = c(0,10), ylim = c(0,10),
     xlab = "Sepal Length (cm)", ylab = "Petal Length (cm)",
     main = "Scatterplot analysis",
     pch = 16)
par(mfrow=c(1,1))



####### Save graphical outputs to files
data(iris)
my.fit <- lm(Sepal.Length ~ Petal.Length, data = iris)

pdf("My_output_file_name.pdf", width = 5, height = 5)
plot(iris$Sepal.Length, iris$Petal.Length,
     col = "deeppink",
     xlim = c(0,10), ylim = c(0,10),
     xlab = "Sepal Length (cm)",
     ylab = "Petal Length (cm)",
     pch = 16)
abline(my.fit, col = "red", lty = "dashed")     #  regression line
dev.off()

# Best practices: make meaningful and informative file names
my.out.pdf <- paste(Sys.Date(),'correlation_plot.pdf',sep = "_")

pdf(my.out.pdf, width = 5, height = 5)

plot(iris$Sepal.Length, iris$Petal.Length,
     col = "deeppink",
     xlim = c(0,10), ylim = c(0,10),
     xlab = "Sepal Length (cm)",
     ylab = "Petal Length (cm)",
     pch = 16)
abline(my.fit, col = "red", lty = "dashed")     #  regression line

dev.off()






    