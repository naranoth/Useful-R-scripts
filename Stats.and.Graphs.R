
#To summarize, random forests are much simpler to train for a practitioner; it's easier to find a good, robust model.

# Libraries ---------------------------------------------------------------

# * * STATS packages ------------------------------------------------------
library(validate) #check_that function!
library(vcd) #Visualizing Categorical Data - cool categorical data graphs (but not ggplot-like)
library(PredictABEL) #Assessment of risk prediction models: ROC, NRI!!!!
library(epiR) #has epi.tests package for 2x2 tables
#
#
library(OneR) #Machine learning categorization package
library(pastecs) #Descriptive statistics: stat.desc(df) = table of descriptive stats
library(Hmisc) #describe(mydf) #Cs(so, it, goes)    #data analysis, has describe function:# missing, distinct, Mean, and cutoff values (0.05, 0.25, etc)
library(doBy)
library(rms) #validate function, etc
library(pROC) #HAS ROC CONFIDENCE INTERVALS! bunch of awesome ROC analyses
library(mlr) #machine learning
library(ROCR) #ROC Analysis
library(psych) # summary statistics by group
library(tigerstats) #Elementary statistics package...mainly for teaching
library(MASS) #supports "Modern Applied Statistics with S" Lots of regression modelling, etc:  functions for estimating linear models through Generalized Least Squares, fitting negative binomial linear models, robust fitting of linear models, and Kruskal’s non-metric multidimensional scaling. 
library(agricolae) #Tukey HSD, STATS for Agricultural research
library(Design)
library(crrstep) #Cox regression stepwise fine and gray method
library(aod)
library(outliers) #Detect outliers, Grubbs, etc
library(CPE) #For calculating concordance probability estimates for Cox models
library(hdnom) #The hd nomogram package
library(arules) #Mining association rules; has discretize function
library(mice) #For imputation, always check the histogram afterwards to make sure it didn't change!
library(FactoMineR)
library(factoextra)
#ggplot plotting packages:
library(GGally) #Another ggplot survival package, more than just survival though...matrices, etc
library(ggtree) #using ggplot on phylogenetic trees
library(DescTools) #Tools for descriptive statistics
library(survminer) #ggplot survival - the main one that I like
library(plotly) #making interactive plots
library(shiny) #Making webpages with plots
library(irr) #has cohen's kappa for 2 raters


# PLOTTING packages -------------------------------------------------------
library(gplots) #for heatmap.2()
library(corrplot) #plot correlation matrix
library(mvtnorm)
library(survminer)
library(ggbio)
library(ggthemes)
library(ggsci) #scientific paper color schemes
library(ggpubr) #Make scientific paper graphics from SURVMINER author
library(ggvis) #like ggplots2 but using pipes
library(cowplot) #Simple add on to ggplot to make scientific ready graphics
library(extrafont) #extrafonts like arial, etc
library(fpc) #kmeans cluster plotting of fpc

# ML ----------------------------------------------------------------------

###Factor analysis
library(factoextra)
library(FactoMineR)
library(ggfortify) ##Unified plotting tools for statistics commonly used, such as GLM, time series, PCA families, clustering and survival analysis. The package offers a single plotting interface for these analysis results and plots in a unified style using 'ggplot2'.
    ##AUTOPLOT ROCKS!
library(cluster)
inst
install.packages("ggvis")
install.packages("ggfortify")
install.packages("easyGgplot2")
install.packages("cowplot")
install.packages("extrafont")

# STATISTICS --------------------------------------------------------------

#Correspondence Analysis is Principle Components Analysis (PCA) for categorical variables as opposed to continuous ones. Data must be non-negative


# ONE LINERS --------------------------------------------------------------
correlate.columns <- function(vars, dat) sapply(vars, function(y) sapply(vars, function(x) assocstats(table(dat[,x], dat[,y]))$cramer))



# * Discretize variables --------------------------------------------------

#Trying discretize function (also discrete in Hmisc package)
#"discretize only implements unsupervised discretization. See packages discretization or RWeka for supervised discretization."
hcc$vim_ctc_lvl <- discretize(hcc$VIM_CTC, method = "cluster", categories = 3, labels = c("low", "med", "high"))
table(hcc$vim_ctc_lvl)
low  med high 
62   25   13 


hcc$ck_ctc_lvl <- discretize(hcc$CK_CTC, method = "cluster", categories = 3, labels = c("low", "med", "high"))
table(hcc$ck_ctc_lvl)
low  med high 
55   29   16


#Stuff from interwebs of interest:
Two Categorical Variables

Checking if two categorical variables are independent can be done with Chi-Squared test of independence.

This is a typical Chi-Square test: if we assume that two variables are independent, then the values of the contingency table for these variables should be distributed uniformly. And then we check how far away from uniform the actual values are.

There also exists a Crammers V that is a measure of correlation that follows from this test



# * * Descriptive Statistics ----------------------------------------------

#Getting them for all columns
#pastecs package
x <- stat.desc(hcc.only) #This actually works well

#Psych package describeBy for grouping variables
describeBy(hcc.study$CK_CTC, group = hcc.study$AJCC_5)


#Also this for correlation
correlate.columns <- function(vars, dat) sapply(vars, function(y) sapply(vars, function(x) assocstats(table(dat[,x], dat[,y]))$cramer))

#DescTools package

Desc(hcc.study$CK_CTC, plotit = TRUE)

#hcc.study$CK_CTC (integer)

#length         n       NAs    unique        0s      mean    meanCI
8e+01     8e+01         0     2e+01     1e+01  5.38e+00  4.21e+00
100.0%      0.0%               13.8%            6.54e+00

.05       .10       .25    median       .75       .90       .95
0.00      0.00  1.00e+00  3.50e+00  9.00e+00  1.31e+01  1.40e+01

#range        sd     vcoef       mad       IQR      skew      kurt
2.30e+01  5.24e+00  9.74e-01  4.45e+00  8.00e+00  1.14e+00  9.38e-01

#lowest : 0 (1e+01), 1e+00 (1e+01), 2e+00 (1e+01), 3e+00 (6e+00), 4e+00 (7e+00)
#highest: 1e+01, 1e+01 (5e+00), 2e+01, 2e+01, 2e+01

#Then get 3 plots of frequency, values, etc etc (histogram, density plot and empiric dist function plots)....looks nice

#Pairwise boxplots etc...
Desc(Stage ~ Vim_CTC, hcc.study, digits = 2, plotit = TRUE)

# Stage ~ Vim_CTC

#Summary: 
n pairs: 8e+01, valid: 8e+01 (98.8%), missings: 1e+00 (1.2%), groups: 4


#Early Stage  Locally Advanced        Metastatic           non-HCC
mean            8.97e-01          3.36e+00          9.44e+00          1.05e-01
median              0.00          3.00e+00          8.00e+00              0.00
sd              2.48e+00          2.65e+00          6.23e+00          4.59e-01
IQR                 0.00          4.75e+00          6.00e+00              0.00
n                  3e+01             2e+01             9e+00             2e+01
np                36.71%            27.85%            11.39%            24.05%
  NAs                    0                 0                 0                 0
0s                 2e+01             5e+00                 0             2e+01

Kruskal-Wallis rank sum test:
  Kruskal-Wallis chi-squared = 43.213, df = 3, p-value = 2.217e-09
Warning:
  Grouping variable contains 1 NAs (1.25%).



Proportions of Stage in the quantiles of Vim_CTC:
  
  0   (0,4]   (4,20]
Early Stage        51.1%   23.5%     6.7%
Locally Advanced   10.6%   58.8%    46.7%
Metastatic          0.0%   11.8%    46.7%
non-HCC            38.3%    5.9%     0.0%


----------
  
  #Longer example to look at all columns
  x <- statsBy(hcc.study, hcc.study$TX_stage)
x <- Desc(hcc.study)
x$Samples_to_use
x$Vim_CTC
> x$Vim_CTC
$xname
#[1] "Vim_CTC"

#$label
#NULL

#$class
#[1] "integer"

#$classlabel
[1] "integer"

$length
[1] 80

$n
[1] 80

$NAs
[1] 0

$main
[1] "14 - Vim_CTC (integer)"

$unique
[1] 11

$`0s`
[1] 48

$mean
[1] 2.3375

$meanSE
[1] 0.4486068

$quant
min    .05    .10    .25 median    .75    .90    .95    max 
0.0    0.0    0.0    0.0    0.0    4.0    8.0    8.2   20.0 

$range
[1] 20

$sd
[1] 4.012461

$vcoef
[1] 1.716561

$mad
[1] 0

$IQR
[1] 4

$skew
[1] 2.185399

$kurt
[1] 5.098594

$small
val freq
1   0   48
2   1    2
3   2    6
4   3    3
5   4    6

$large
val freq
1  20    1
2  16    1
3  14    1
4  12    1
5   8    6

$freq
level   freq   perc  cumfreq  cumperc
1       0  5e+01  60.0%    5e+01    60.0%
  2       1  2e+00   2.5%    5e+01    62.5%
  3       2  6e+00   7.5%    6e+01    70.0%
  4       3  3e+00   3.8%    6e+01    73.8%
  5       4  6e+00   7.5%    6e+01    81.2%
  6       6  5e+00   6.2%    7e+01    87.5%
  7       8  6e+00   7.5%    8e+01    95.0%
  8      12  1e+00   1.2%    8e+01    96.2%
  9      14  1e+00   1.2%    8e+01    97.5%
  10     16  1e+00   1.2%    8e+01    98.8%
  11     20  1e+00   1.2%    8e+01   100.0%
  
  $maxrows
[1] 12

$x
[1]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
[44]  0  0  0  0  0  1  1  2  2  2  2  2  2  3  3  3  4  4  4  4  4  4  6  6  6  6  6  8  8  8  8  8  8 12 14 16 20


# * * * Detecting outliers ------------------------------------------------
X <- hcc.study$CK_CTC[hcc.study$Cancer == 0]
#For HCC study healthy/NMLD patients, what are the outliers!

library(outliers)
#grubbs.flag <- function(x) {
outliers <- NULL
test <- x
grubbs.result <- grubbs.test(test)
pv <- grubbs.result$p.value
while(pv < 0.05) {
  outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
  test <- x[!x %in% outliers]
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
}
return(data.frame(X=x,Outlier=(x %in% outliers)))
}

ggplot(grubbs.flag(X),aes(x=X,color=Outlier,fill=Outlier))+
  geom_histogram(binwidth=diff(range(X))/30)+
  theme_bw()


###With error code inserted:
grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  # throw an error if there are too few values for the Grubb's test
  if (length(test) < 3 ) stop("Grubb's test requires > 2 input values")
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    # stop if all but two values are flagged as outliers
    if (length(test) < 3 ) {
      warning("All but two values flagged as outliers")
      break
    }
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
  }
  return(data.frame(X=x,Outlier=(x %in% outliers)))
}

# * * * * Creating Quartiles ----------------------------------------------

This should do it:
  
  tableOne <- within(tableOne, quartile <- as.integer(cut(salesPrice, quantile(salesPrice, probs=0:4/4), include.lowest=TRUE)))
#...Some details:

#  The within function is great for calculating new columns. You don't have to refer to columns as tableOne$salesPrice etc.

tableOne <- within(tableOne, quartile <- <<<some expression>>>)
The quantile function calculates the quantiles (or in your case, quartiles). 0:4/4 evaluates to c(0, 0.25, 0.50, 0.75, 1).

Finally the cut function splits your data into those quartiles. But you get a factor with weird names, so as.integer turns it into groups 1,2,3,4.

Try ?within etc to learn more about the functions mentioned here...

# * * *  T test ------------------------------------------------------------------

###T test
t.test(muts.clin$mutation_load[muts.clin$tail_only=="Head"], muts.clin$mutation_load[muts.clin$tail_only=="Tail"])


###Chi goodness of fit for poisson
install.packages("vcd")
library(vcd)
set.seed(2014);y=rpois(200,5)
set.seed(2014);y=rnorm(100, 5, 0.3) # goodfit asks for non-negative values
# output the results
gf = goodfit(y,type= "poisson",method= "ML")
plot(gf,main="Count data vs Poisson distribution")
summary(gf)


# * * ANOVA ---------------------------------------------------------------
x <- aov(formula = log_muts ~ tumor_loc, data=crc.comp)
summary(x)


# * * * * Turkey HSD for multiple comparisons -----------------------------

TukeyHSD(x) #x from aov() above


# * * * Regression analysis --------------------------------------------------------------

#OS:
xyplot(OS~mutation_load, data=muts.clin,type=c("p","r"))
lmGC(OS~mutation_load, data=muts.clin)
#Linear Regression

#Correlation coefficient r =  -0.03874 

#Equation of Regression Line:

#  OS = 564.9443 + -0.0025 * mutation_load 

#Residual Standard Error:	s   = 426.16 
#R^2 (unadjusted):		R^2 = 0.0015 

#PFS
xyplot(PFS~mutation_load, data=muts.clin,type=c("p","r"))
lmGC(PFS~mutation_load, data=muts.clin)
#Linear Regression

#Correlation coefficient r =  -0.1449 

#Equation of Regression Line:

#  PFS = 548.255 + -0.0123 * mutation_load 

#Residual Standard Error:	s   = 418.9478 
#R^2 (unadjusted):		R^2 = 0.021 

PFS.muts.lm <- lm(formula = PFS ~ mutation_load, data = muts.clin)
summary(PFS.muts.lm)
#Call:
lm(formula = PFS ~ mutation_load, data = muts.clin)

#Residuals:
Min     1Q Median     3Q    Max 
-533.4 -287.6 -117.2  182.4 2157.3 

#Coefficients:
Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   548.255013  29.787759  18.405   <2e-16 ***
mutation_load  -0.012251   0.004781  -2.562   0.0109 *  
  ---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 418.9 on 306 degrees of freedom
(206 observations deleted due to missingness)
#Multiple R-squared:  0.021,	Adjusted R-squared:  0.0178 
#F-statistic: 6.565 on 1 and 306 DF,  p-value: 0.01088


#* * * Chi Square --------------------------------------------------------------
set.seed(1234)
X <- rep(c("A","B"),20)
Y <- sample(c("C","D"),40,replace=T)

table(X,Y)
chisq.test(table(X,Y),correct=F)
# I don't use Yates continuity correction

#Let's make a matrix with tons of columns

Data <- as.data.frame(
  matrix(
    sample(letters[1:3],2000,replace=T),
    ncol=25
  )
)

# You want to select which columns to use
columns <- c(3,7,11,24)
vars <- names(Data)[columns]

# say you need to know which ones are associated with each other.
out <-  apply( combn(columns,2),2,function(x){
  chisq.test(table(Data[,x[1]],Data[,x[2]]),correct=F)$p.value
})
out
out <- cbind(as.data.frame(t(combn(vars,2))),out)


# * * * *multiple chi2 on DF ----------------------------------------------

apply(hcc.pred.only[,1:39], 2, function(x) summary(xtabs(~hcc.pred.only$OS+x)))


# * * * *  multiple regression by COXPH -----------------------------------

apply(hcc.pred.only[,1:39], 2, function(x) summary(coxph(Surv(OS, Mortality) ~ x, data = hcc.pred.only))$logtest["pvalue"])
#Age                            Sex                           Race                         CK_CTC 
#2.987643e-01                   8.110581e-01                   7.689362e-01                   8.551028e-02 
#Vim_CTC                       PDL1_CTC                      HCC_cause                          Stage 
#2.798673e-02                   6.190101e-02                   6.241415e-01                   1.425349e-03 
#Milan                     UCSF_Stage                       Region_5                      Cirrhosis 
#1.065285e-02                   1.120606e-04                   9.352803e-04                   7.869889e-01 
#number_of_lesions                  Lesion_1_size Multifocal_or_outside_criteria                Cumulitave_size 
#4.372673e-03                   9.954941e-01                   1.267951e-04                   9.574462e-01 
#Bilobar           Vascular_involvement                           Mets                      Bilirubin 
#5.156470e-04                   8.615879e-03                   9.519186e-02                   2.657777e-04 
#Cr                            INR                    MELD_Physio                        Albumin 
#3.481763e-05                   3.943513e-06                   1.075331e-04                   7.199202e-02 
#Encephalopathy                        Ascites                    Child_Class                           ECOG 
#1.479163e-03                   2.939822e-01                   3.122008e-02                   2.840703e-07 
#BCLC                AFP_most_recent                        AFP_MAX                        Surgery 
#1.025920e-03                   1.735016e-08                   1.552584e-06                   3.522264e-03 
#Chemo_Study                    Progression                      Mortality                            RFS 
#6.864266e-02                   9.842990e-04                   3.108624e-15                   7.638403e-01 
#PFS                             OS                        OS_DODx 
#2.341395e-03                   1.053746e-08                   3.202095e-08 


# PCA ---------------------------------------------------------------------

#From CNV Final Figures, *** PCA PLOTS
#Input above from "c" as matrix
c <- genes[,c(4:14)] #DOUBLE CHECK THE NUMBER OF COLUMNS
c.pca <- PCA(c, quali.sup = 1)
plot(c.pca, habillage = 1, 
     col.hab = c("green", "blue", "red"), 
     title = "Dataset projected onto PC1-2 Subspace")


#Manual PCA
c <- genes[,c(4:10,12:14)] #REMOVING CTC 4 AS OUTLIER
row.names(c) <- c$cytoband
c <- c[,c(2:10)]
c <- as.matrix(c)
head(cov(c))

#By using the function eigen the eigenvalues and eigenvectors of the covariance matrix are computed
Eigenvalues <- eigen(cov(c))$values
Eigenvectors <- eigen(cov(c))$vectors

#Now, the Principal Components can be estimated via a matrix multiplication
PC <- as.matrix(c) %*% Eigenvectors
cov(PC)

#In a next step we calculate the proportions of the variation explained by the various components:
#
print(round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2))
round(cumsum(Eigenvalues)/sum(Eigenvalues) * 100, digits = 2)

#[1] 95.26 99.05 99.88 99.95 99.97 99.99 99.99 100.00 100.00 100.00
#[11] 100.00
#
#
#Thus, the first component round(Eigenvalues[1]/sum(Eigenvalues)*100, digits=2) explains 95.26,
#and the first three eigenvectors of the covariance matrix explain 99.88 of the total variation in the data.
#This suggest that the effective dimension of the space of yield curves could be three and any of the yield
#curves from our data set can be described by a linear combination of the first three loadings, while the
#relative error being very small.

###PCA using prcomp
#The best way to do PCA with R is to use the function prcomp from the package stats. prcomp uses
#as arguments simply a data matrix. Furthermore, with the argument scale = TRUE (default: scale = FALSE) the variables can be scaled to a unit variance before the analysis takes place.
#
#To run PCA on this data we use
R> c.pca <- prcomp(c)
#Note: To reproduce our previous calculation we use the default case (scale = FALSE). The PrintOutput of c.pca gives us the estimated standard deviations as well as the rotations (loadings).

#PCA and factor analysis: 
##http://www.statmethods.net/advstats/factor.html
##http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining
mydata = c
# entering raw data and extracting PCs 
# from the correlation matrix 
fit <- princomp(mydata, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

###Exploratory Factor Analysis
#The factanal( ) function produces maximum likelihood factor analysis.

# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation 
fit <- factanal(mydata, 3, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(mydata),cex=.7) # add variable names

## Principal Axis Factor Analysis
library(psych)
fit <- factor.pa(mydata, nfactors=3, rotation="varimax")
fit # print results

## Determine Number of Factors to Extract
## A crucial decision in exploratory factor analysis is how many factors to extract. The nFactors package offer a suite of functions to aid in this decision. 
library(nFactors)
ev <- eigen(cor(mydata)) # get eigenvalues
ap <- parallel(subject=nrow(mydata),var=ncol(mydata),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

# PCA Variable Factor Map 
# The FactoMineR package offers a large number of additional functions for exploratory factor analysis. This includes the use of both quantitative and qualitative variables, as well as the inclusion of supplimentary variables and observations. 
library(FactoMineR)
result <- PCA(mydata) # graphs generated automatically

#Thye GPARotation package offers a wealth of rotation options beyond varimax and promax.
# GRAPHICS and PLOTTING ---------------------------------------------------


#examples of good plots using lattice:
#http://motioninsocial.com/tufte/



# Plotting descriptive statistics -----------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#Density plot and boxplot for head vs tail and mutation burden
xtabs(~tail_only, data=muts.clin)
#tail_only
#Head Tail 
#416   50 
rowPerc(xtabs(~tail_only, data=muts.clin))
#tail_only  Head  Tail Total
#           89.27 10.73   100
favstats(mutation_load ~ tail_only, data=muts.clin)
#tail_only min     Q1 median      Q3    max     mean       sd   n missing
#1      Head   1  142.0 2454.0 4772.00 109362 3521.892 6797.270 416       0
#2      Tail  12 1982.5 3028.5 9374.75  23389 6315.580 6956.621  50       0

#Ploting this
densityplot(~mutation_load|tail_only, data = muts.clin)
boxplot(mutation_load~tail_only, data = muts.clin)


# GGSCI Palattes ----------------------------------------------------------

mypal = pal_npg()(10)
mypal
#GGSCI figure out what colors I want

mypal = pal_uchicago(palette = "light")(5) #first 5 colors of this palette
mypal
library(scales)
show_col(mypal) #Display as a graph
mycols <- rev(mypal)
mycols
mycols2 <- c("#800000FF", "#D6D6CEFF")

# GGPUBR ------------------------------------------------------------------

library(ggpubr)
#Interesting package with lots of features
desc_statby(data = pdac.only, measure.var = "CTC_Count_4mL", grps = "Stage", ci = 0.95)
ggstripchart(data = pdac.only, x = "Strat.variable" , y = "CTC_Count_4mL",palette = "npg")

# SURVMINER ---------------------------------------------------------------
http://www.sthda.com/english/rpkgs/survminer/articles/Playing_with_fonts_and_texts.html
http://www.sthda.com/english/wiki/survminer-0-2-4
#Also categorize variables .. see example in webpage


# Adjusted survival curves for the variable "sex"
ggcoxadjustedcurves(res.cox, data = lung,
                    variable  = lung[, "sex"],   # Variable of interest
                    legend.title = "Sex",        # Change legend title
                    palette = "npg",             # nature publishing group color palettes
                    curv.size = 2                # Change line size
)
# * * DETERMINE CUTPOINTS IN SURVMINER! -----------------------------------

# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(myeloma, time = "time", event = "event",
                         variables = c("DEPDC1", "WHSC1", "CRIM1"))
summary(res.cut)
##        cutpoint statistic
## DEPDC1    279.8  4.275452
## WHSC1    3205.6  3.361330
## CRIM1      82.3  1.968317
# 2. Plot cutpoint for DEPDC1
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "DEPDC1", palette = "npg")



# * * * ggplot survival examples for checking data ------------------------------------------

#Checking event times, etc
ggplot(clin.segs, 
       aes(x = OS_MONTHS,
           group = OS_STATUS,
           colour = OS_STATUS,
           fill = OS_STATUS
       )) +
  geom_density(alpha = 0.5) +
  scale_x_discrete(breaks = seq(0, 60, by = 6))

mle.surv <- 
  survfit(
    Surv(OS_MONTHS,os_deceased) ~ 1,
    data = clin.segs %>%
      dplyr::mutate(os_deceased = OS_STATUS == 'DECEASED')
  )
ggsurvplot(mle.surv)
