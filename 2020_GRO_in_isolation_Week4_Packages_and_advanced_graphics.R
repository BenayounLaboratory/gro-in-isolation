##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-04-15

#######   Week 4 code snippets    #######

# comment lines that are there for humans reading the script start with "#"
# All other lines will be interpreted/executed by R


####################################################### 
###########        R software packages        ########### 

# get all the relevant info on loaded packages
sessionInfo()

# one way to save that info for when you write up your paper
# with correct package versions
sink(file = paste(Sys.Date(),"_session_Info.txt", sep =""))
sessionInfo()
sink()


# Installing packages
# CRAN packages:
install.packages()

# Bioconductor packages:
BiocManager::install()

# github packages:
remotes::install_github()
devtools::install_github()

# If you NEED to install from source
# you will need to download the tar.gz code
# and then run "R CMD INSTALL packagename.tar.gz" in the console


# find out which packages are ALREADY installed:
installed.packages()


####################################################### 
###########        Advanced R graphics      ########### 

########## Lattice graphs
library(lattice) 

# kernel density plots of sepal length by species
data(iris)
densityplot(~ Sepal.Length | Species, 
            main="Density Plot by Species",  
            xlab="Sepal Length (cm)", layout=c(1,3),data = iris)

# boxplots of muder rates in high vs low educated states, as a function of freezing days
data(state)
state.x77 <- data.frame(state.x77)

state.x77$education <- "Low Education"
state.x77$education[state.x77$HS.Grad > 55] <- "High Education"

state.x77$Freeze <- "Less_than_60d"
state.x77$Freeze[state.x77$Frost > 60] <- "60_to_100d"
state.x77$Freeze[state.x77$Frost > 100] <- "More_than_100d"

bwplot(Freeze ~ Murder | education,
       ylab="Freezing days", xlab="Murders", 
       main="Murder by education and days with frost in US states", 
       layout=c(1,2),
       data = state.x77)


########## ggplot2 graphs
library(ggplot2) 

# observations (points) are overlayed and jittered
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")

theme_set(theme_bw())
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")

theme_set(theme_light())
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")


theme_set(theme_dark())
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")


theme_set(theme_classic())
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")

theme_set(theme_minimal())
qplot(Species, Sepal.Length, data=iris,
      geom=c("boxplot", "jitter"), 
      fill=Species, main="Sepal Length by Species",
      xlab="", ylab="Sepal Length (cm)")

# more info of ggplot2 themes
# https://ggplot2.tidyverse.org/reference/ggtheme.html



# Separate regressions of Sepal Length on petal Length for each species
qplot(Petal.Length, Sepal.Length,
      method="lm", formula= y ~ x,
      data=iris,
      geom=c("point", "smooth"), 
      color=Species, 
      main="Regression of petal vs sepal length by species", 
      xlab="Petal Length", ylab="Sepal Length")


####################################################### 
######   Useful R ‘advanced’ plotting packages   ###### 

##### pretty heatmaps using pheatmap
library(pheatmap)
data(iris)

pheatmap(iris[,1:4],show_rownames = F)

pheatmap(t(iris[,1:4]),show_colnames = F, scale = "row")

# change color palette
pheatmap(t(iris[,1:4]),show_colnames = F, scale = "row", 
         color = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))



##### bootstrapped clustering using pvclust
library(pvclust)
data(state); state.x77 <- data.frame(state.x77)

# do hierarchical clustering
my.pv <- pvclust(t(state.x77),nboot=10); # 10 bottstraps for speed when webcasting
## my.pv <- pvclust(t(state.x77),nboot=1000); # use at least 100-1000 if for publication

plot(my.pv)


##### Venn Diagrams using Vennerable
library(Vennerable)
data(state); state.x77 <- data.frame(state.x77)

high.income <- rownames(state.x77)[state.x77$Income > 4500]
high.life <- rownames(state.x77)[state.x77$Life.Exp > 70]
high.grad <- rownames(state.x77)[state.x77$HS.Grad > 50]

my.criteria <- list("Income_over_4500"    = high.income,
                    "Life_expectancy_over_70"    = high.life,
                    "High_School-Grad_over_50%"  = high.grad)
my.Venn <- Venn(my.criteria)

# with scaling
plot(my.Venn, doWeights=T)

# without scaling
plot(my.Venn, doWeights=F)
# equivalent to Vennerable::plot(my.Venn, doWeights=F)
doEuler


# shapes rather than circles
plot(my.Venn, doWeights=F, type="squares")
plot(my.Venn, doWeights=F, type="triangles")
plot(my.Venn, doWeights=F, type="AWFE")
plot(my.Venn, doWeights=F, type="ChowRuskey")

# More advanced plotting
# with or without class labels
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=FALSE,Faces=FALSE,DarkMatter=FALSE))
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=TRUE,Faces=FALSE,DarkMatter=FALSE))

# with or without color filling
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=FALSE,Faces=FALSE,DarkMatter=FALSE))
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=FALSE,Faces=TRUE,DarkMatter=FALSE))

# with or without universe size
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=FALSE,Faces=FALSE,DarkMatter=FALSE))
plot(my.Venn, doWeights=F, show=list(FaceText="signature",SetLabels=FALSE,Faces=FALSE,DarkMatter=TRUE))


# check old vignette:
#https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/Vennerable/inst/doc/Venn.pdf?revision=58&root=vennerable
