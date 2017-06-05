#List of awesome libraries: http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html#tk.drr_mlt

# LIBRARIES ---------------------------------------------------------------
update.packages(checkBuilt = T)
# Useful libraries
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
inject.dots <- function(df) {names(df) <- sub(" ", ".", names(df));df}
library(plyr) #Need to load plyr before dplyr OR ELSE!
library(tidyverse) #Loads ggplot2, tibble, tidyr, readr, purrr, dplyr
library(dtplyr) #data.table and dplyr code together!
library(stringr) #manipulate strings... str_sub!

library(survival) #Main survival analysis package...really good
library(survMisc) #Extension of survival package, has some good plotting functions
library(Epi) #Epidemiology stuff...Lexis diagrams, etc. Powerful but different format than survival package
library(Hmisc) # BEST FUNCTIONS: describe() and Cs() #data analysis, has describe function:# missing, distinct, Mean, and cutoff values (0.05, 0.25, etc)
library(DescTools) #Tools for descriptive statistics
library(survminer) #ggplot survival - the main one that I like

library(Biobase) #The base bioconductor package
library(BiocInstaller) #gets biocLite
library(Biostrings) #defines containers
library(BSgenome) #full genome sequences for many species
library(GenomicRanges) #genomic interval sets
library(GenomicFeatures) #retrieve and manage genomic features from public databases
library(maftools) #For manipulating maf files in LOTS of ways
library(biomaRt) #Downloading MART objects
library(Hmisc) #describe(mydf) #Cs(so, it, goes)    #data analysis, has describe function:# missing, distinct, Mean, and cutoff values (0.05, 0.25, etc). check out: http://aliquote.org/articles/tech/hmisc/hmisc.html
library(xlsx) #make excel sheets
library(corrplot)
#For data preparation stuff
#library(plyr)
# Tidy the dataset.
#library(ggplot2) # Visualise data.
#
library(magrittr) #PIPES!!! Get it? magritt!
library(reshape)
library(reshape2)
library(janitor)# Basic data cleaning...find duplicates, remove empty columns, etc #tabyl(mydf, sort = TRUE) %>% add_totals_row()
library(car) #recode()!!!! turn numeric / continuous into categorical: recode(x, "1:3='Low'; 4:7='Mid'; 8:hi='High'")
library(scales) #comma(mynumvec) #format data for graphing
library(diffobj) # diffObj(x,y) #Tells you how two objects are different!
library(vcd)
library(car) #ScatterplotMatrix
library(ggdendro)
library(hexbin)
library(tidytext) #text-mining
# data.table and dplyr loaded in tidyverse
#library(data.table)
#library(dplyr) # Data preparation and pipes %>%.
library(ggsci) #ggplot scientific palettes
library(ggpubr) #Similar, getting ggplot ready for publication
#library(tidyverse) #The tidyverse is a set of packages that work in harmony because they share common data representations and API design. The tidyverse package is designed to make it easy to install and load core packages from the tidyverse in a single command.
#library(purrr) #spliting dataframes and feeding to functions...but super complicated
library(gdata) #read in EXCEL FILES! also, 
  #rename variables in columns with: data <- rename.vars(data, c("x","y","z"), c("first"        ,"second","third"))
  #reorder factor levels: trt4 <- reorder(trt, new.order=c("PLACEBO", "300 MG", "600 MG",        "1200 MG"))
library(bit64) #Allows long and large numbers

# * * STATS packages ------------------------------------------------------
library(validate) #check_that function!
library(vcd) #Visualizing Categorical Data - cool categorical data graphs (but not ggplot-like)
library(PredictABEL) #Assessment of risk prediction models: ROC, NRI!!!!
library(epiR) #has epi.tests package for 2x2 tables
#
#
library(OneR) #Machine learning categorization package
library(arules)
library(pastecs) #Descriptive statistics: stat.desc(df) = table of descriptive stats
library(Hmisc) #data analysis, has describe function:# missing, distinct, Mean, and cutoff values (0.05, 0.25, etc)
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
library(gmodels) #Basic data cleaning made easy, such as finding duplicates by multiple columns, making R-friendly column names and removing empty columns. It also has some nice tabulating tools, like adding a total row and doing tables with percentages and easy crosstabs.

# * * PLOTTING ------------------------------------------------------------
library(ggplot2) #duh
library(GGally) #Another ggplot survival package, more than just survival though...matrices, etc
library(ggtree) #using ggplot on phylogenetic trees
library(DescTools) #Tools for descriptive statistics
library(survminer) #ggplot survival - the main one that I like
library(plotly) #making interactive plots
library(shiny) #Making webpages with plots
library(ggiraph) #create ggplot interactive graphs
library(DT) #create searchable tables
library(ggsci) #Scientific journal themes for ggplot2
library(ggthemes) #Extra themes for ggplot2



# * * Bioconductor packages ---------------------------------------------
library(maftools)
##
##
library(qtlcharts)
library(Rpressa)
library(ComplexHeatmap) #For making heatmaps (esp of expression association stuff)
library(cgdsr) #CBioPortal download package
library(R.matlab) #reading and writing MATLAB mat files
library(phyloseq) #does ggplots for phylogenetic trees
library(ConsensusClusterPlus) #Build RNASeq heatmaps
library(DESeq2) #differential expression analysis
library(vcfR) #reading VCF files
library(myvariant) #VCF file annotation
library(htSeqTools) #toolset for NGS, has giniCoverage and ssdCoverage to look at coverage uniformity
library(TEQC) #NGS QC tools

#Copy number packages
library(QDNAseq) #Awesome CGH calling package...interacts w lots of other packages too, best one yet for HCC CNV data
library(QDNAseq.hg19) #get bin data for hg19
library(CGHcall) #GOOD PLOTS, takes data from QDNAseq
library(CGHregions) #good plots from CGHcall stuff
library(cn.mops) #good for low level CNV?
library(DNAcopy)
library(copynumber) #Good for plotting copy number stuff I think
library(ape) #phylogenetics and evolution package
library(snapCGH) #aCGH for agilent etc
#biocLite("CGHbase")
library(CGHbase)

#Somatic Signatures Stuff
library(VariantAnnotation)
library(SomaticCancerAlterations)
library(SomaticSignatures)
library(sca)
library(sva)


###Annotation packages:
library(BSgenome.Hsapiens.UCSC.hg19)
##
##
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#reload for maftools
library(NMF)

####
##
.
#

#RNASeq Packages:
library(GOplot) #For plotting GO enrichment, circle graphs, etc
library(cummeRbund)
library(ballgown)
library(CancerSubtypes)

# * * * * * Other/random packages -----------------------------------------

library(reinstallr) #reinstall packages lost during R upgrade
library(googlesheets) #read in googlesheets 
library(splitstackshape) #like excel text to columns

# * * creating documents --------------------------------------------------

#creating powerpoint excel and word documents
library(xlsx) #make excel sheets
library(ReporteRs) #make word and powerpoint documents
#Both are from: http://www.sthda.com/english/wiki/exporting-data-from-r
library(officer)
library(rio) #easily import almost any file type, from spread sheets to SAS, XML, google sheets, and CLIPBOARD. Has Import(), Export() and Convert()
#Really cool stuff
library(DiagrammeR) #Make awesome interconnected diagrams


library(devtools) #Installing github packages, etc
#Install stuff:
install.packages("splitstackshape")
install_github("seandavi/Rpressa")
biocLite("Rpressa")
biocLite("BiocUpgrade")
biocLite("CancerSubtypes")
install.packages("irr")
install.packages("officer")
library(rJava)
install.packages("rJava")
install.packages("xlsx")
install.packages("ReporteRs")
biocLite("DESeq2")
biocLite("htSeqTools")
biocLite("SomaticSignatures")
install.packages("ggiraph")
install.packages("validate")
install.packages("shiny")
install.packages("plotly")
install.packages("factoextra")
install.packages("FactoMineR")
install.packages("xlsx")
install.packages("ReporteRs")
install.packages("MASS")
biocLite("sva")
biocLite("sca")
biocLite("SomaticSignatures")
biocLite("BSgenBSgenome.Hsapiens.1000genomes.hs37d5")
biocLite("mygene")
install.packages(c("rattle", "randomForest", "lubridate", "FSelector"))
biocLite("FSelector")
install.packages("rJava",type='source')
install.packages("RWeka")
install.packages('cairoDevice') 
install.packages("tibble")
biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
biocLite("signeR")
install.packages("sos")
install.packages("hexbin")
biocLite("survminer")
biocLite("lazyeval")
biocLite("cummeRbund")
biocLite("dtplyr")
biocLite("ballgown")
biocLite("GOplot")
install.packages("purrr")
biocLite("TEQC")
install.packages("pastecs")
install.packages("epiR")
install.packages("outliers")
install.packages("gdata")
biocLite("copynumber")
require(vcd)

#Other packages
library(sos) #Search the help pages of R packages...I think google is better
library(mygene)
library(cgdsr)
library(FSelector) #Feature Selection
library(cairoDevice)
library(lubridate) # Handle dates.
library(tigerstats)
library(psych) #Also has describe but gives skew, kurtosis, SE, median, SD, etc instead
library(lubridate)
library(jsonlite)
library(signeR) #Somatic Sigs associated package
library(lazyeval) #Nonstandard formula evaluation (e.g. NSE) - relates to programming mainly
library(rattle) # GUI for R! - The weather dataset and normVarNames().
library(randomForest) # Impute missing values using na.roughfix().



#Detaching unused packages:
detach("TxDb.Hsapiens.UCSC.hg19.knownGene", unload = TRUE)
detach("MASS", unload = TRUE)
detach("rattle", unload = TRUE) # The weather dataset and normVarNames().
detach("cgdsr", unload = TRUE) #CBioPortal download package
detach("cairoDevice", unload = TRUE)
detach("jsonlite", unload = TRUE)

detach("stringi", unload = T)
library(stringi)




# Install packages --------------------------------------------------------
install.packages("rio")
install.packages("reinstallr")
install.packages("googlesheets")
install.packages("arules")
install.packages("mice")
install.packages("doBy")
install.packages("tigerstats")
install.packages("crrstep")
install.packages("Design")
install.packages("R.matlab")
biocLite("copynumber")
install.packages("PredictABEL")
install.packages("agricolae")
install.packages("OneR")
install.packages("ape")
biocLite("phyloseq")
biocLite("ggtree")
biocLite("snapCGH")
# Functions to add --------------------------------------------------------

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

#Repair excel dates:
repairExcelDates <- function(x, yearcol=3, fmt="%m/%d/%Y") {
x <-  do.call(rbind, lapply(strsplit(x, "/"), as.numeric))
year <- x[,yearcol]
if(any(year>99)) stop("dont'know what to do")
x[,yearcol] <- ifelse(year <= as.numeric(format(Sys.Date(), "%Y")), year+2000, year + 1900) 
# if year <= current year then add 2000, otherwise add 1900
x <- apply(x, 1, paste, collapse="/")
as.Date(x, format=fmt)
}
# Making names work -------------------------------------------------------

#Getting rid of spaces...need to run multiple times for multiple spaces
inject.dots <- function(df) {names(df) <- sub(" ", ".", names(df));df}

#Or use this:
names(ctm2) <- gsub(" ", "_", names(ctm2))

#or, as mentioned in the first answer (though not in a way that would fix all spaces):

spaceless <- function(x) {colnames(x) <- gsub(" ", "_", colnames(x));x}
newDF <- spaceless(ctm2)

#using base R
impute.med <- function(x) {
  z <- median(x, na.rm = TRUE)
  x[is.na(x)] <- z
  return(x)
}

#Using plyr
impute.med <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
dat2 <- sapply(dat, function(x){
  if(is.numeric(x)){
    impute.med(x)
  } else {
    x
  }
}
)

max.values <- function(x){
  if(is.numeric(x)){max(x,na.rm = 1)}
  else{max(as.character(x),na.rm=1)}
}  

values.incolumn.greaterthan.0.5 <- function(data){
  new_df <- data[sapply(data,is.numeric)]
  sapply(new_df, function(x) sum(abs(x) > 0.5))
} #Computes the number of values that are greater than 0.5...e.g copy number changes by chromosome arm
col.great.5(arm)
#arm.1p  arm.1q  arm.2p  arm.2q  arm.3p  arm.3q  arm.4p  arm.4q  arm.5p  arm.5q  arm.6p  arm.6q  arm.7p  arm.7q  arm.8p  arm.8q  arm.9p  arm.9q 
#52     171      10      14      20      15      53      73      65      55      68      77      54      62     180     159      58      46 


# REFERENCES AND BLOGS ----------------------------------------------------

#http://simplystatistics.org/ #Great blog posts about stats and data science
#www.r-bloggers.org

#http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html#tk.drr_mlt
# Excellent list of good packages and descriptions.

#https://bcbio-nextgen.readthedocs.io/en/latest/
#Best practices and a AWS server that can be used directly!

http://core-genomics.blogspot.com/2016/01/flowcnv-seq-almost-novel-metthod-for_29.html
#Good blog, this post is about scCNV Seq
# Other stuff -------------------------------------------------------------

https://github.com/chrisalbon/code_r #some easy to use R codes

#Correlation of nominal and continuous variables:
http://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618

#BEST PRACTICES FOR PREPARING DATA:
#http://www.sthda.com/english/wiki/best-practices-in-preparing-data-files-for-importing-into-r

#Exporting data from R
#http://www.sthda.com/english/wiki/exporting-data-from-r
### INCLUDES EDITABLE R GRAPHS TO POWERPOINT!!!!!

#Hmisc! has: Two of my favorites: describe, a more robust summary function, and Cs, which creates a vector of quoted character strings from unquoted comma-separated text. Cs(so, it, goes) creates c("so", "it", "goes"). 

# Notes -------------------------------------------------------------------

My number of signatures results in a huge anount of the variance not being explained...need to figure out whats going on there.



# READING IN DATA ------------------------------------------------------------
#fread creates a data.table which is NOT A DATA FRAME
#readr's read_table etc creates a data frame wrapped as a tbl_df

###Scan: allows you to scan in data
more.data = scan()
#then just copy and paste it into the window that pops up. enter space to end 


### fread - READING IN LARGE DATA TABLES
#Best way to read in large data tables:
library(data.table)
mydata <- fread("mylargefile.txt")

#Alt:
#Using read_table in readr

#Trying to store as much data as you can in databases rather than flat files. (As well as being a better permanent storage medium, data is passed to and from R in a binary format, which is faster.) read.csv.sql in the sqldf package, as described in JD Long's answer, imports data into a temporary SQLite database and then reads it into R.
#http://stackoverflow.com/questions/1727772/quickly-reading-very-large-tables-as-dataframes-in-r?rq=1

#Readr package:
#http://www.sthda.com/english/wiki/fast-reading-of-data-from-txt-csv-files-into-r-readr-package



# * * Paste CLIPBOARD -----------------------------------------------------

data <- read.table(pipe("pbpaste"), sep = "\t", header = F)
rnas.current <- rownames(groups)
rnas.current <- as.factor(rnas.current)
rnas.current <- as.data.frame(rnas.current)
data$samples <- data$V1
colnames(data) <- c("V1", "rnas.current")
x <- merge(x = rnas.current, y = data, by = "rnas.current", all.x = T)
no.rnas <- data %>% filter(V1 %in% rnas.current)


# DATABASE PREPARATION ----------------------------------------------------
#Convert blanks to NAs
convert_blank_to_na <- function(x) {
  if (!purrr::is_character(x)) {
    warning('input vector is not character - returning original input')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}

#Make Tidy Names
tidy.name.vector <- make.names(name.vector, unique=TRUE)

#make.names() makes syntactically valid names out of character vectors. A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number.
#Additionally, flag unique=TRUE allows you to avoid possible dublicates in new column names.


#From
crc.comp <- crc.comp %>% dplyr::mutate_each(funs = funs(convert_blank_to_na), everything()) #runs convert_blank_to_na on all columns
crc.comp <- crc.comp %>% mutate_if(is.character, as.factor) #Makes everything a factor
crc.comp <- as.data.frame(crc.comp)



# VALIDATING DATA ---------------------------------------------------------



# * * * Figure out who is missing values ----------------------------------------
clin.segs %>%
  dplyr::filter(is.na(OS_STATUS) | OS_STATUS == '') %>%
  dplyr::select(OS_STATUS, OS_MONTHS) %>%
  str()
#3 ppl are missing OS_STATUS

## confirm 4 fewer observations than original
assert_that(nrow(clin_data) == nrow(clinical_data) - 4)


# * * * VALIDATE PACKAGE --------------------------------------------------
data(women)
cf <- check_that(women, height > 0, weight > 0, height/weight > 0.5)
summary(cf)
#Alternative with pipes:
women %>% check_that(height > 0, weight > 0, height/weight > 0.5) %>% summary()
#check that returns the number of items checked and how many passed the rule as well as how many NAs
barplot(cf,main="Checks on the women data set") #plot this!

#Validator objects
v <- validator(height > 0, weight > 0, height/weight > 0)
v
#The validator object has stored the rule and assigned names to them for future reference. To check this, we confront the data set with the validation rules we’ve just defined:
cf <- confront(women,v)
cf

#more complex with BMI
v <- validator(
  BMI := (weight*0.45359)/(height*0.0254)^2
  , height > 0
  , weight > 0
  , BMI < 23
  , mean(BMI) > 22 & mean(BMI) < 22.5
)
v
cf <- confront(women,v)
summary(cf)

#Validator objects:


# CLEANING DATA -----------------------------------------------------------

#https://www.r-bloggers.com/a-data-cleaning-example/



# STRINGS -----------------------------------------------------------------



#* * * * Print string one per line -----------------------------------------------
wd <- getwd()
cat("Current working dir: \n", wd)
#   Current working dir:
#   /Users/colincourt/NGS/aSomaticSignature/NEW.combine.dataset

sprintf("Current working dir: ", wd)
#[1] "Current working dir: "

sprintf("Current working dir: %s", wd) #%s = 
#[1] "Current working dir: /Users/colincourt/NGS/aSomaticSignature/NEW.combine.dataset"
cat(sprintf("Current working dir: %s\n", wd))
#Current working dir: /Users/colincourt/NGS/aSomaticSignature/NEW.combine.dataset 
message(sprintf("Current working dir: %s\n", wd))
* Current working dir: /Users/colincourt/NGS/aSomaticSignature/NEW.combine.dataset*

###Not working
rs <- rowSums(H.tcga)
rs
print(cat(" \n", rs))


# * * * * Replacing values in a string ------------------------------------

#Replace a comma with a space:
z <- gsub(",", " ", sig.genes)


# * * * * Split elements in a string based on a character -----------------
x <- c(as = "asfef", qu = "qwerty", "yuiop[", "b", "stuff.blah.yech")
x
strsplit(x, "e")
strsplit("a.b.c", ".")

#Had a string where some of the values had multiple gene names in each item... this worked to seperate them into individual genes:
sig.genes <- x.sig$gene #Gene names from cuffdiff, some w/ multiple genes per XLOC value
head(sig.genes)
z <- unlist(strsplit(sig.genes, ","))
table(z)


###WORKED!!!

#Other option
#Sample data
myDat <- read.table(text=
                      "pages|count
                    [page 1, page 2, page 3]|23
                    [page 2, page 4]|4
                    [page 1, page 3, page 4]|12", header=TRUE, sep="|") 

#Then
# if factors, convert to characters
pages <- as.character(myDat$page)

# remove brackets.  Note the double-escape's in R
pages <- gsub("(\\[|\\])", "", pages)

# split on comma
pages <- strsplit(pages, ",")

# find the largest element
maxLen <- max(sapply(pages, length))

# fill in any blanks. The t() is to transpose the return from sapply
pages <- 
  t(sapply(pages, function(x)
    # append to x, NA's.  Note that if (0 == (maxLen - length(x))), then no NA's are appended 
    c(x, rep(NA, maxLen - length(x)))
  ))

# add column names as necessary
colnames(pages) <- paste(c("First", "Second", "Third"), "Page")

# Put it all back together
data.frame(pages, Count=myDat$count)

#
mydf <- data.frame(
  ID = c("A01", "A02", "A03", "A04", "B01", "B02"),
  Name = c("Boyd", "Rufus", "Dana",
           "Carole", "Ramona", "Kelley"),
  Likes = c("1,2,4,5,6", "1,2,5,6", "1,3,4",
            "2,3,6", "1,2,3,5", "1,4,6"),
  Siblings = c("Reynolds , Albert , Ortega",
               "Cohen , Bert , Montgomery",
               "Pierce", "Colon , Michelle , Ballard",
               "Snyder , Joann ,", "James , Roxanne ,"))
mydf
#There are two typical R approaches to split data like this up. First is strsplit, which will create a list of your split data:

X <- lapply(mydf, function(y) {
  strsplit(as.character(y), ",")
})
X

Y <- lapply(mydf, function(y) {
  read.csv(text = as.character(y), header = FALSE,
           strip.white = TRUE, blank.lines.skip = FALSE,
           fill = TRUE, stringsAsFactors = FALSE)
})
Y

# DATAFRAMES -------------------------------------------------------------
  #Columns = vectors, rows = values for an instance
  # [x, y] = x as rows, y as columns ALWAYS
dspath <- "http://rattle.togaware.com/weather.csv"
weather <- read.csv(dspath)

#* *  query DF structure, factors etc -----------------------------------------

dim(weather)
#[1] 366  24
names(weather)
#Lists column names
str(weather)
#'data.frame':	366 obs. of  24 variables:
#  $ Date         : Factor w/ 366 levels "2007-11-01","2007-11-02",..: 1 2 3 4 5 6 7 8 9 10 ...
#$ Location     : Factor w/ 1 level "Canberra": 1 1 1 1 1 1 1 1 1 1 ...
#$ MinTemp      : num  8 14 13.7 13.3 7.6 6.2 6.1 8.3 8.8 8.4 ...
#$ MaxTemp      : num  24.3 26.9 23.4 15.5 16.1 16.9 18.2 17 19.5 22.8 ...
summary(weather)
#Date         Location      MinTemp          MaxTemp         Rainfall       Evaporation    
#2007-11-01:  1   Canberra:366   Min.   :-5.300   Min.   : 7.60   Min.   : 0.000   Min.   : 0.200  
#2007-11-02:  1                  1st Qu.: 2.300   1st Qu.:15.03   1st Qu.: 0.000   1st Qu.: 2.200  
#2007-11-03:  1                  Median : 7.450   Median :19.65   Median : 0.000   Median : 4.200  
#2007-11-04:  1                  Mean   : 7.266   Mean   :20.55   Mean   : 1.428   Mean   : 4.522  
#2007-11-05:  1                  3rd Qu.:12.500   3rd Qu.:25.50   3rd Qu.: 0.200   3rd Qu.: 6.400  
#2007-11-06:  1                  Max.   :20.900   Max.   :35.80   Max.   :39.800   Max.   :13.800  
(Other)   :360  
glimpse(data.frame or table) #shows columns as rows
head(clinical) #shows first 6 rows of each column
summary(clinical) #shows summary of responses in each column (no #, yes #, etc)
str(clinical) #shows columns w/ how many factors they have. 
dim(clinical) #dimensions or nrow() or ncol() or length() = nrows
attach() #Makes the data frame columns directly usable (no need for $)
df[x,y] = NULL #removes a column or row
names(clinical) #lists column names
unique(loc)
describe(tcga.maf$n_depth)
#tcga.maf$n_depth 
#n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90 
#1173160        0     2940    0.999    73.58    102.4        3        4        7       25       76      167 
#.95 
#271 

#lowest :    0    1    2    3    4, highest: 6080 6117 6124 6146 6802 
table(comb.maf$Variant_Type)

                     DEL      deletion of <=200bp                      INS     insertion of <=200bp 
                  223856                   214936                    28177                    91750 
single base substitution                      SNP 
                 2813606                   611115 
              
comb.maf %>% count(Variant_Type) #Or count(comb.maf$Variant_Type)
# A tibble: 6 × 2
              Variant_Type       n
                    <fctr>   <int>
1                      DEL  223856
2      deletion of <=200bp  214936
3                      INS   28177
4     insertion of <=200bp   91750
5 single base substitution 2813606
6                      SNP  611115

#Using xtabs
xtabs(~tail_only, data=muts.clin)
#tail_only
#Head Tail 
#416   50 
#and with tigerstats to get the percentages
rowPerc(xtabs(~tail_only, data=muts.clin))
#tail_only  Head  Tail Total
#           89.27 10.73   100


# * * * Other DF queries and stats ----------------------------------------
#also 
colPerc()!!!!

quantile(clin.t.mut$maximum_tumor_dimension, na.rm = TRUE)
  #0%  25%  50%  75% 100% 
  #1.5  3.0  3.5  4.5 12.0 
#Plotting empiric cumulative frequency distribution of data points
ecdf(clin.t.mut$maximum_tumor_dimension)
  #Empirical CDF 
  #Call: ecdf(clin.t.mut$maximum_tumor_dimension)
  #x[1:38] =    1.5,    1.8,      2,  ...,      9,     12
plot(ecdf(clin.t.mut$maximum_tumor_dimension))
    ###SUPER USEFUL!

# * DATAFRAME manipulation --------------------------------------------------
######Dataframe maniplation stuff


#* * *  Column as row names  ----------------------------------------------------

alex.sigs3 <- as.data.frame(alex.sigs3)
rownames(alex.sigs3) <- alex.sigs3[,1] #column 1 as row names
alex.sigs3 <- alex.sigs3[,-1] #delete column 1


#subsettting rows
  rows_to_keep <- c(TRUE, FALSE, TRUE, FALSE)
  limited_writers_df <- writers_df[rows_to_keep,]
  limited_writers_df
#Adding rows
  new_row <- c(50, 22, "Roberto", "Bolano", "MALE", "2003-07-15")
  writers_df_large <- rbind(writers_df, new_row)

#Deleting rows
gene.panel.genes.gsea <- gene.panel.genes[-c(2,5,7,11,12),]

### Adding column between other columns
dat$C <- NA
dat <- dat[, c("A", "C", "B")]

#A  C          B
#1  0.596068 NA -0.7783724
#2 -1.464656 NA -0.8425972

###You can also use append
dat <- data.frame(A = rnorm(2), B = rnorm(2))
as.data.frame(append(dat, list(C = NA), after = 1))

#A   C          B
#1 -0.7046408  NA  0.2117638
#2  0.8402680  NA -2.0109721


#Rows to list format
xy.list <- split(xy.df, seq(nrow(xy.df)))
#And if you want the rownames of xy.df to be the names of the output list, you can do:
xy.list <- setNames(split(xy.df, seq(nrow(xy.df))), rownames(xy.df))
  
#Going long to wide with stack and unstack
Subject <- c(1,2)
Gender <- c("M", "F")
Read <- c(10, 7)
Write <-c(8, 4)
Listen <- c(7, 6)
observations_wide <- data.frame(Subject, Gender, Read, Write, Listen)
observations_wide
long_format <- stack(observations_wide,select=c(Read,Write,Listen))
long_format

wide_format <- unstack(long_format,values ~ ind)
wide_format

#Merge data frames -- do it with merge() based on column names
data2 <- data.frame(Age.At.Death=c(22,40,72,41), Location=5:8)
new_writers_df <- merge(writers_df, data2)
data2 <- data.frame(x.Age.At.Death=c(21,39,71,40), Location=5:8)
  #This assumes the columns are the same even if the names of them are not


#Merge data.frames with the superior rbindlist – data.table we meet again. As the community knows: “rbindlist is an optimized version of do.call(rbind, list(...)), which is known for being slow when using rbind.data.frame.”

I(factors) # will "insulate" the items, making them class "AsIs"


#* * * Join/Merge Dataframes ---------------------------------------------------

#Trying rbind
all.comb.motifs <- rbind(tcga.comb.motif, icgc.comb.motif)
tcga = 611115
icgc = 1443448
combined = 2054563
nrow(tcga.comb.motif) + nrow(icgc.comb.motif) #2054563
#Works

comb <- merge(icgc.comb2, tcga.comb, by = "specimen_id", all = TRUE)
1443448 + 36840
> 1443448 + 36840
[1] 1480288
#Worked!

#Inner join: merge(df1, df2) will work for these examples because R automatically joins the frames by common variable names, but you would most likely want to specify merge(df1, df2, by = "CustomerId") to make sure that you were matching on only the fields you desired. You can also use the by.x and by.y parameters if the matching variables have different names in the different data frames.

#Outer join: merge(x = df1, y = df2, by = "CustomerId", all = TRUE)

#Left outer: merge(x = df1, y = df2, by = "CustomerId", all.x = TRUE)

#Right outer: merge(x = df1, y = df2, by = "CustomerId", all.y = TRUE)

#Cross join: merge(x = df1, y = df2, by = NULL)




# * * * * MERGE MULTIPLE DATABASES ----------------------------------------

#Using reduce:http://stackoverflow.com/questions/32526889/merge-multiple-data-tables-with-the-same-column-names
Reduce((function() {counter = 0
function(x, y) {
  counter <<- counter + 1
  d = merge(x, y, all = T, by = 'x')
  setnames(d, c(head(names(d), -1), paste0('y.', counter)))
}})(), list(DT1, DT2, DT3, DT4, DT5))

#OR:

list(DT1 = DT1, DT2 = DT2, DT3 = DT3, DT4 = DT4, DT5 = DT5) %>%
  bind_rows(.id = "source") %>%
  mutate(source = paste(source, "y", sep = ".")) %>%
  spread(source, y)


# * * * Filtering all NA cols or mostly NA rows ---------------------------

###Using ColSums:
new.tcga2 <- new.tcga
new.tcga2 <- new.tcga2[, colSums(is.na(new.tcga2)) != nrow(new.tcga2)] #Removed 22 columns
new.tcga2 <- new.tcga2[, colSums(is.na(new.tcga2)) < nrow(new.tcga2) * 0.75] ###REMOVE ALL COLUMNS WITH >75% NAS!!
colnames(new.tcga2)
new.tcga <- new.tcga2[, c(343, 1:342)]

#Find out which cols have only NAs
mvc <- sapply(tcga.maf, function(x) sum(is.na(x)))
mvn <- names(which(mvc == nrow(tcga.maf)))
mvn
mvn2 <- names(which(mvc >= 0.7*nrow(tcga.maf)))
mvn2
tcga.maf <- tcga.maf[, !(colnames(tcga.maf) %in% mvn2)] #Get the boulean output but doesn't drop them
tcga.maf[ , -which(names(tcga.maf) %in% mvn2)] 
#[1] -16 -21 -22 -47 -68 -71 -72 -73 -74 -75 -76 -77 -78 -79




#* * FILTER/REMOVE ROWS -------------------------------------------------------------
###Using %not in%

#  Filterning out useless Variant_Classifications
useless.vars <- c("3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron")
T3.mut.maf <- filter(T2.mut.maf, Variant_Classification %not in% useless.vars)
#866635 ==> 324028 mutations

"Filtering rows"

# base R
crime.ny.2005 <- crime.by.state[crime.by.state$Year==2005 &
                                  crime.by.state$State=="New York", ]

# dplyr
crime.ny.2005 <- filter(crime.by.state, State=="New York", Year==2005)


#Remove all letters from a column
pretxbx$pT_stage <- gsub('[a-z]', '', pretxbx$pT_stage, ignore.case = TRUE)
comb.maf2$Chromosome <- gsub('[[:lower:]]', '', comb.maf2$Chromosome)

# Reading in TSV files
array <- read.table("exp_array.tsv", sep = "\t", header = TRUE)
exp <- read.table("exp_seq.tsv", sep = "\t", header = TRUE)

# Turning Stacked TSV file into expression type file/matrix etc
exp2 <- exp[ , c(1, 8, 9)]
explist <- split(exp2, exp$icgc_donor_id)
expmat <- as.data.frame(explist)

# better option:
xlist <- llply(explist, subset, select = c(normalized_read_count)) #create a new list with just the column I want from
#each of the dataframes in the table
patients <- names(xlist)            #create vector of the names of the patients
expmap <- as.data.frame(xlist, colnames(patients))   #rework the list as a data frame with colnames
colnames(rawmatrix) <- patients #change the column names

#* * * Adding or deleting values (all alphanumeric, etc) within a columns values-----------------------------------------------
###Adding
#Using regex pattern
icgc.comb$Chromosome <- sub("^", "chr", icgc.comb$Chromosome)

#sprintf is a lot more powerful than plain concatenation.
dat$V1 <- sprintf('chr%i', dat$V1)

#We can also use interaction:
  
df$V1 <- interaction( "chr", df$V1, sep = "")
df

#Or using sqldf:
library(sqldf)    
df$V1 <- as.character(df$V1)
df$V1 <- sqldf("select 'chr'|| V1 as V1 from df") 

### Deleting
#You can use gsub for this:
gsub('[[:digit:]]+', '', x)

#or
gsub('[0-9]+', '', x)

# unlist, so that you have a vector
x <- unlist(x)

#How to remove non-alphanumeric characters
preg_replace("/[^A-Za-z0-9 ]/", '', $string);


#* * * Deleting duplicate rows -------------------------------------------------

# A sample data frame:
df <- read.table(header=TRUE, text='
                 label value
                 A     4
                 B     3
                 C     6
                 B     3
                 B     1
                 A     2
                 A     4
                 A     4
                 ')


# Is each row a repeat?
duplicated(df)
#> [1] FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE  TRUE

# Show the repeat entries
df[duplicated(df),]
#>   label value
#> 4     B     3
#> 7     A     4
#> 8     A     4

# Show unique repeat entries (row names may differ, but values are the same)
unique(df[duplicated(df),])
#>   label value
#> 4     B     3
#> 7     A     4

# Original data with repeats removed. These do the same:
unique(df)
#>   label value
#> 1     A     4
#> 2     B     3
#> 3     C     6
#> 5     B     1
#> 6     A     2
df[!duplicated(df),]
#>   label value
#> 1     A     4
#> 2     B     3
#> 3     C     6
#> 5     B     1
#> 6     A     2

#* * FILTER/REMOVE COLUMNS ---------------------------------------------------

#Removing columns
badcols <- c("TUMOR_TISSUE_SITE", "SAMPLE_TYPE_ID", "SAMPLE_TYPE", "INFORMED_CONSENT_VERIFIED", "OTHER_PATIENT_ID", "OTHER_SAMPLE_ID", "PATHOLOGY_REPORT_FILE_NAME", "PATHOLOGY_REPORT_UUID") #These are the columns I don't want, just FYI
clin.h2 <- subset(clin.h, select=-c(TUMOR_TISSUE_SITE, SAMPLE_TYPE_ID, SAMPLE_TYPE, INFORMED_CONSENT_VERIFIED, OTHER_PATIENT_ID, OTHER_SAMPLE_ID, PATHOLOGY_REPORT_FILE_NAME, PATHOLOGY_REPORT_UUID))
#--

# Remove columns that you don’t want from a data frame:
class(expmat)
[1] "data.frame"
colsdontwant <- c(1,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,41,43,44,46,47,49,50,52,53,55,56,58,59,61,62,64,65,67,68,70,71,73,74,76,77,79,80,82,83,85,86,88,89,91,92,84,95) 
expmat2 <- expmat[, ! names(expmat) %in% colsdontwant, drop = FALSE]

#   Rename columns
# df = dataframe
# old.var.name = The name you don't like anymore
# new.var.name = The name you want to get

names(df)[names(df) == 'old.var.name'] <- 'new.var.name'
colnames(df)[colnum] <- "new.var.name"
#This code pretty much does the following:
  
#   names(df) looks into all the names in the df
#   [names(df) == old.var.name] extracts the variable name you want to check
#   <- 'new.var.name' assigns the new variable name.
# I think it's worth specify that [names(df) == old.var.name] actually returns a vector with true/false values. So it has the potential to change multiple column names if, for example, regular expressions are used.
#   For regular expression results, use something like names(df) = sub('pattern', 'replacement', names(df))


#* * REORDER/RENAME COLUMNS --------------------------------------------------


### This is an old question, but it is worth noting that you can now use setnames from the data.table package.
library(data.table)
setnames(DF, "oldName", "newName")

# or since the data.frame in question is just one column: 
setnames(DF, "newName")

# And for reference's sake, in general (more than once column)
nms <- c("col1.name", "col2.name", etc...)
setnames(DF, nms)


###   Reorder columns - this is tough in R
names(trading)
[1] "OpenDate"   "CloseDate"  "Symbol"     "Action"     "Lots"       "SL"         "TP"         "OpenPrice"
[9] "ClosePrice" "Commission" "Swap"       "Pips"       "Profit"     "Gain"       "Duration"   "Trader"   
[17] "System"

#old method:
> trading = trading[, c(16:17, 1:15)]
> names(trading)
[1] "Trader"     "System"     "OpenDate"   "CloseDate"  "Symbol"     "Action"     "Lots"       "SL"       
[9] "TP"         "OpenPrice"  "ClosePrice" "Commission" "Swap"       "Pips"       "Profit"     "Gain"     
[17] "Duration"

#using setdiff
> refcols <- c("Trader", "System")
> #
  > trading <- trading[, c(refcols, setdiff(names(trading), refcols))]
> names(trading)
[1] "Trader"     "System"     "OpenDate"   "CloseDate"  "Symbol"     "Action"     "Lots"       "SL"       
[9] "TP"         "OpenPrice"  "ClosePrice" "Commission" "Swap"       "Pips"       "Profit"     "Gain"     
[17] "Duration"

#Meta Data Cleansing
#Make all names lowercase
names(weather)
#[1] "Date"          "Location"      "MinTemp"       "MaxTemp"       "Rainfall"      "Evaporation"   "Sunshine"   
names(weather) <- normVarNames(weather)



#* * * * * plyr and dplyr df manipulation ------------------------------------------

"Filtering rows"

# base R
crime.ny.2005 <- crime.by.state[crime.by.state$Year==2005 &
                                  crime.by.state$State=="New York", ]

# dplyr
crime.ny.2005 <- filter(crime.by.state, State=="New York", Year==2005)

"Arranging and ordering"

# base R
crime.ny.2005 <- crime.ny.2005[order(crime.ny.2005$Count, 
                                     decreasing=TRUE), ]

# dplyr
crime.ny.2005 <- arrange(crime.ny.2005, desc(Count))

#   The base R solution ranks each row by value of "Count" in decreasing order, and uses the rank vector to        subset the "crime.ny.2005" data frame. The dplyr solution appears to be about 20% faster.

"Selecting columns"

# base R
crime.ny.2005 <- crime.ny.2005[, c("Type.of.Crime", "Count")]

# dplyr
crime.ny.2005 <- select(crime.ny.2005, Type.of.Crime, Count)

#   This example is relatively self-explanatory. Here the base R solution appears to be faster, by about 30%.

"Creating new columns"
# base R
crime.ny.2005$Proportion <- crime.ny.2005$Count /
  sum(crime.ny.2005$Count)

# dplyr
crime.ny.2005 <- mutate(crime.ny.2005, 
                        Proportion=Count/sum(Count))

"Aggregation and summarization"
# base R
summary1 <- aggregate(Count ~ Type.of.Crime,
                      data=crime.ny.2005,
                      FUN=sum)
summary2 <- aggregate(Count ~ Type.of.Crime,
                      data=crime.ny.2005,
                      FUN=length)
summary.crime.ny.2005 <- merge(summary1, summary2,
                               by="Type.of.Crime")

# dplyr
by.type <- group_by(crime.ny.2005, Type.of.Crime)
summary.crime.ny.2005 <- summarise(by.type,
                                   num.types = n(),
                                   counts = sum(Count))


"All together now"
# We haven't showcased the best part of dplyr yet... it presents itself when combining all of these statements:

# base R
crime.by.state <- read.csv("CrimeStatebyState.csv")
crime.ny.2005 <- crime.by.state[crime.by.state$Year==2005 &
                                  crime.by.state$State=="New York", 
                                c("Type.of.Crime", "Count")]
crime.ny.2005 <- crime.ny.2005[order(crime.ny.2005$Count, 
                                     decreasing=TRUE), ]
crime.ny.2005$Proportion <- crime.ny.2005$Count /
  sum(crime.ny.2005$Count)
summary1 <- aggregate(Count ~ Type.of.Crime,
                      data=crime.ny.2005,
                      FUN=sum)
summary2 <- aggregate(Count ~ Type.of.Crime,
                      data=crime.ny.2005,
                      FUN=length)
final <- merge(summary1, summary2,
               by="Type.of.Crime")


# dplyr
crime.by.state <- read.csv("CrimeStatebyState.csv")
final <- crime.by.state %>%
  filter(State=="New York", Year==2005) %>%
  arrange(desc(Count)) %>%
  select(Type.of.Crime, Count) %>%
  mutate(Proportion=Count/sum(Count)) %>%
  group_by(Type.of.Crime) %>%
  summarise(num.types = n(), counts = sum(Count))


# * * * MERGE COLUMNS to new column ---------------------------------------

#Using tidyr "unite"
###NB: remember that if using strings need to use unite_ instead..unite() is only for numeric!

#* * *  Editing Chr strings in DF --------------------------------------------------------

#stringr - str_sub
hw <- "Hadley Wickham"
str_sub(hw, 1, 1)
str_sub(hw, end = 6)
str_sub(hw, 8, 14)
str_sub(hw, -1)
str_sub(hw, -7)
str_sub(hw, end = -7)



#* * * Finding Values in Data Frames -------------------------------------------

### Find values that match between columns
use %in% as follows
A$C %in% B$C
# [1]  TRUE FALSE  TRUE  TRUE

# Find which column names match between data frames:
colnames(Q.mut.maf) %in% colnames(mut.maf)
[1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
[19]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[37] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE



# * * CREATE COLUMNS FROM OTHER COLUMN VALUES -----------------------------

#https://www.analyticsvidhya.com/blog/2015/11/8-ways-deal-continuous-variables-predictive-modeling/
#8 Ways to deal with Continuous Variables in Predictive Modeling
# * * * * Using cut -------------------------------------------------------
#How can we do this in R? There’s a great function in R called cut() that does everything at once.  It takes in a continuous variable and returns a factor (which is an ordered or unordered categorical variable).  Factor variables are extremely useful for regression because they can be treated as dummy variables.  
mydata$Agecat2<-cut(mydata$Age, seq(0,30,5))

mydata$Agecat3<-cut(mydata$Age, seq(0,30,5), right=FALSE)

mydata$Agecat4<-cut(mydata$Age, seq(0,30,5), right=FALSE, labels=c(1:6))

# * * * * MUTATE() ----------------------------------------------------------
other=data.frame(name=c("a","b","a","c","d"),result=c("Y","N","Y","Y","N"))
data.table::setDT(other)
other[ , table(result), by = name]
df.diamonds_ideal <- mutate(df.diamonds_ideal, price_per_carat = price/carat)

head(df.diamonds_ideal)
# carat   cut     color price clarity   price_per_carat
# 0.23    Ideal     E   326     SI2        1417.391
# 0.23    Ideal     J   340     VS1        1478.261
# 0.31    Ideal     J   344     SI2        1109.677
# 0.30    Ideal     I   348     SI2        1160.000
# 0.33    Ideal     I   403     SI2        1221.212
# 0.33    Ideal     I   403     SI2        1221.212


# * * *get counts as columns ----------------------------------------------

http://stackoverflow.com/questions/40454138/spread-columns-by-count-in-r-dplyr

car <- c("a","b","b","b","c","c","a","b","b","b","c","c")
type <- c("good", "regular", "bad","good", "regular", "bad","good", "regular", "bad","good", "regular", "bad")
car_type <- data.frame(car,type)


# * * NESTED LIST for stats by col value ----------------------------------

nested <- iris %>%
  group_by(Species) %>%
  nest()

## # A tibble: 3 × 2
##      Species              data
##       <fctr>            <list>
## 1     setosa <tibble [50 × 4]>
## 2 versicolor <tibble [50 × 4]>
## 3  virginica <tibble [50 × 4]>

means <- map(nested$data, colMeans)

## [[1]]
## Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
##        5.006        3.428        1.462        0.246 
## 
## [[2]]
## Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
##        5.936        2.770        4.260        1.326 
## 
## [[3]]
## Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
##        6.588        2.974        5.552        2.026
# TIDYVERSE ---------------------------------------------------------------


# * CHAINING INTRO --------------------------------------------------------

# Subset the diamonds dataset to 'Ideal' cut diamonds
# THEN, keep (select) the variables: carat, cut, color, price, clarity
# THEN, add new variable called price_per_carat (mutate)

df.diamonds_ideal_chained <- diamonds %>%
  filter(cut=="Ideal") %>%
  select(carat, cut, color, price, clarity) %>%
  mutate(price_per_carat = price/carat)

head(df.diamonds_ideal_chained)
# carat   cut     color price clarity   price_per_carat
# 0.23    Ideal     E   326     SI2        1417.391
# 0.23    Ideal     J   340     VS1        1478.261
# 0.31    Ideal     J   344     SI2        1109.677
# 0.30    Ideal     I   348     SI2        1160.000
# 0.33    Ideal     I   403     SI2        1221.212
# 0.33    Ideal     I   403     SI2        1221.212

### AND FOR GRAPHS:
# dplyr + ggplot
# HISTOGRAM of price, ideal cut diamonds
diamonds %>%                                        # Start with the 'diamonds' dataset
  filter(cut == "Ideal") %>%                        # Then, filter down to rows where cut == Ideal
  ggplot(aes(price)) +                            # Then, plot using ggplot
  geom_histogram() +                              # and plot histograms
  facet_wrap(~ color)                             # in a 'small multiple' plot, broken out by 'color' 

# LISTS -------------------------------------------------------------------


#* *  Lists and referencing list items ----------------------------------------

n = c(2, 3, 5) 
s = c("aa", "bb", "cc", "dd", "ee") 
b = c(TRUE, FALSE, TRUE, FALSE, FALSE) 
x = list(n, s, b, 3)   # x contains copies of n, s, b
x[2]
#[[1]] 
#[1] "aa" "bb" "cc" "dd" "ee"

x[[2]]
#[1] "aa" "bb" "cc" "dd" "ee"


# R CODING ----------------------------------------------------------------

#* *  basic loops -------------------------------------------------------------

#Check if file is available
if(file.exists(myfile)){
  data <- read.delim(myfile)
  ...rest of your code here...
  ...
} else{
  stop("Could not find input file")
}

#for loop of histograms
for(i in 1:3){
  hist(data[,i])
}

#example multiline code
for(i in 1:3){
  ...process the data....
  hist(data[,i])
  ...customise the plot...
  ...export...
  ...etc...
  
}


#Biobase is the package to display microarray structure
#Each column is a sample and each row is a gene

#The MA plot
M <- log2(evals[,1]) - log2(evals[,2])
A <- 0.5*(log2(evals[,1]) + log2(evals[,2]))
plot(A,M)

#The Tilde (~) is the R version of a formula
boxplot(mygene~myfactor) #plots mygene against myfactor


#Call a specific row based on a value in a column --FAILED
specimens[which(specimens$icgc_specimen_id == SP110836)] #nope
v <- c(SP110836)
out <- specimens[specimens$icgc_specimen_id %in% v] #nope


#* Important Functions, Programming, etc -----------------------------------

###   ifelse  ###
#Basically returns TRUE or FALSE based on the results of the function
ifelse(function, "the test parameter is true", "the test parameter is false")

ifelse(TRUE, "the test parameter is true", "the test parameter is false")
[1] "the test parameter is true"
ifelse(FALSE, "the test parameter is true", "the test parameter is false")
[1] "the test parameter is false"

x <- c(T, F, T)
ifelse( x, "it's true", "it's false")
[1] "it's true"  "it's false" "it's true" 

mtcars
mtcars[ , "efficiency" ] <-
  ifelse(
    (mtcars[ , "mpg" ] >= 20) ,
    "at least 20 mpg"  ,
    "less than 20 mpg"
  )
mtcars

changetoNA <- function(colnum,df) {
  col <- df[,colnum]
  if (is.numeric(col)) {  #edit: verifying column is numeric
    col[col == -1 & is.numeric(col)] <- NA
  }
  return(col)
}
df <- data.frame(sapply(1:5, changetoNA, df))

changetobinary <- function(colnum,df) {
  for(i in colnum){
    df$colnum[df$colnum>=1] <- 1
  }
  return(col)
}
#for loop of histograms
for(i in 1:3){
  hist(data[,i])
}

x<-1:9
length(x)
if (length(x)<=10) {
x<-c(x,10:20);print(x)}
printx <- function(x) {
  if (length(x)<=10) {
    x<-c(x,10:20);print(x)}
}
printx(5)

mt

# R learning modules ------------------------------------------------------


#R learning modules, classes, etc
#http://watson.nci.nih.gov/~sdavis/assets/tutorials/rintro/RIntro.pdf
#Factors are a special class that look like characters but ARE NOT
    #Can be recognized by lack of " " around them.
  #Stored in R as numbers with a key name

#Most data structures are vectors with "attributes" and "classes"
citizen<-c("uk","us","no","au","uk","us","us","no","au")
citizenf<-factor(citizen)
citizenf #no " "
as.numeric(citizenf)
#Most data structures are vectors with "attributes" and "classes"
attributes(citizenf) #levels and factors are attributes
unclass(citizenf)
table(citizenf)
#Matrix are numeric, data frames are whatever the vector of that column is
###ALWAYS [ROW,COLUMN]
x <- 1:10
y <- rnorm(10)
mat <- cbind(x,y)
mat
mat2 <- rbind(x,y)
mat2
z = paste0('a',1:10)
tab <-cbind(x,y,z)
tab<-data.frame(x,y,z)
rownames(tab)<-paste0("row",1:10)
colnames(tab) = paste0('col',1:3)
tab[,1]>7 #Logical, True if >7 for column 1
tab[1:3,] #shows rows 1-3
tab[tab[,1]>7,] #displays rows where column 1 >7
tab[tab[,1]>7,3] #Creates a factor with levels = a1-a10, returns the three levels that match (a8-10)
tab[tab[,1]>7,2:3]
tab[tab$x>7,3] #returns nothing
tab$z[tab$x>3] #NULL

#List is a collection of objects of any type
    #Indexed either by name (my.list$name3) or component number (my.list[[3]])
      #A data frame is a list of matched column vectors


#Programming / Control structures in R
#for (x in set) {operations}
#while (x in condition){operations}
# if (condition) {
#     some operations
# } else { other operations }

x <- 1:9
if (length(x) <=10) {
  x <-c(x,10:20);print(x)
}

if (length(x)<5) {
  print(x)
} else {
  print(x[5:20])
}

for (i in x) print(i)
for (i in x) i
i

#loop over character vector
y<-c('a','b','hi there')
for (i in y) print(i)
#While loop
j<-1
while(j<10) { # do this while j<10
  print(j)
  j<-j+2
  } # at each iteration, increase j by 2

###Apply - better than a loop in R
apply(mat,1,fun) # over rows--second argument is 1
apply(mat,2,fun) # over columns--second argument is 2
  # In either case, the output is a vector.

lapply(list, function) #applies the function to every element of list
sapply(list or vector, function) 
  #applies the function to every element of list or vector, and returns a vector, when possible (easier to process)
tapply(x, factor, fun) #uses the factor to split vector x into groups, and then applies fun to each group

# create a list
my.list <- list(a=1:3,b=5:10,c=11:20)
my.list
# Get the mean for each member of the list return a vector
sapply( my.list, mean)
# Get the full summary for each member of the list, returned as a list
lapply( my.list, summary)
# Find the mean for each group defined by a factor
my.vector <- 1:10
my.factor <- factor(c(1,1,1,2,2,2,3,3,3,3))
tapply(my.vector, my.factor, mean)

#Functions are objects and are assigned to names, just like data.
myFunction = function(argument1,argument2) {
  expression1
  expression2
}

add1 = function(x) { # this function adds one to the first argument and returns it
  x + 1
  }
add1(17)
add1(c(17,18,19,20)) #Apply function to each member of a vector
#Use edit to edit a function in a new window
add2 = edit(add1)

            #### ##
#http://www.bioinformatics.babraham.ac.uk/training/Advanced_R_Course/Advanced%20R%20course%20booklet.pdf

a.vector <- c(4,8,2,3,8,4)
names(a.vector) <- c("bob", "sam", "sue", "eve", "don", "jon")
a.vector
order(a.vector)
a.vector[order(a.vector)]

#Lists
to.test %in% days #Sees which indices are in both (returns boolean)


### Looping ###
apply (data.frame, rows(1)/cols(2), function name, function arguments)
  #The vectors produced by apply must always be the same type (boolean, etc)
#Test what class is returned:
apply(example.data,1,class)

tapply(vector of values, vector of categories, function to apply)
  #e.g: give mean for different chromosome numbers in a data.frame containing all chromosome numbers

#lapply() = for lists --> returns a list
#sapply() = for lists --> returns a vector


# Segmentation plotting ---------------------------------------------------

plot_segments = function(PAAD.snp6wgl_hg19.txt){

  seg = read.delim("PAAD.snp6wgl_hg19.txt", sep="\t")
  
  seg.spl = split(seg,as.factor(as.character(seg$Chromosome)))
  
  
  
  pdf(file=paste(segments_file,"pdf",sep="."),paper="special",width=12,onefile=T,pointsize=8)
  
  
  
  par(mfrow=c(4,1))
  
  
  
  for(i in 1:length(seg.spl)){
    
    
    
    x = seg.spl[[i]]    
    
    
    
    plot(x$Start,x$Segment_Mean,xlim = c(x[1,3],x[nrow(x),4]),pch = "",ylim = c(-3,3),xlab = paste("chr",names(seg.spl[i]),sep="_"),ylab = "log2 ratio")
    
    
    
    points(x$End,x$Segment_Mean,pch = "")
    
    
    
    segments(x0 = x$Start , y0 = x$Segment_Mean , x1 = x$End, y1 = x$Segment_Mean, lwd = 2, col = "maroon")
    
    
    
    abline(h = 0, lty = 1,lwd = 0.5)    
    
  }
  
  dev.off()
  
}

# GENERIC Variables -------------------------------------------------------
#We will store the dataset as the generic variable ds (short for dataset). This will make the following steps somewhat generic and often we can just load a different dataset into ds and these steps can simply be re-run without change.
dsname <- "weather"
ds <- get(dsname)
dim(ds)
#[1] 366  24

# Make Granges object -----------------------------------------------------



cnv <- makeGRangesFromDataFrame(paadCNV.nogl) ##########
summary(cnv)



# GRanges subsetting gene CNVs ------------------------------------------------------
#This one worked for getting Gene CNVs: https://support.bioconductor.org/p/67118/
geneRanges <- function(db, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }

splitColumnByOverlap <-function(query, subject, column="ENTREZID", ...){
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

setwd("~/NGS/HCC.NGS.CNV/clustering/CTC.primary.only")
segcopy <- read_delim("SegCopy.txt", delim = "\t")
df <- makeGRangesFromDataFrame(segfix, keep.extra.columns = F)
df
genes <- cosmic[,c(20:22,1)]
good <- c(1:22, "X", "Y", "M")
genes <- genes %>% filter(chromosome_name %in% good)
colnames(genes) <- Cs(chromosome, start, end, gene.symbol)
genes$loc <- round((genes$start + genes$end) / 2, 0)
genes <- genes[, c(1,5,4)]
genes$start <- genes$loc
colnames(genes)[2] <- "end"
genes <- genes[,c(1,4,2,3)]
genes$chromosome <- paste('chr', genes$chromosome, sep = '')
#
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns = T)
genes.gr

#Getting the genes into the format of the Segcopy/fixed/norm, etc so that there are 6194 rows and only saving the gene symbol names
symInCnv <- splitColumnByOverlap(genes.gr,df,column = "gene.symbol")
y <- unstrsplit(symInCnv, sep = ", ")
y <- as.data.frame(y) #6194 row and 1 character string of gene names in the appropriate places!
colnames(y) <- "gene.symbol"

seg.genes <- cbind(segcopy, y)
seg.genes <- seg.genes[,c(1:3,29,4:28)]


#Alternative GRanges merging:
##Example from: https://support.bioconductor.org/p/54470/
gr2 <- GRanges("chr1", IRanges(c(5, 100, 1500), c(5, 100, 1500),
                               names=paste0("rsid:", letters[1:3])), score=1:3)
gr1 <- GRanges("chr1", IRanges(c(4,1000,99), c(98,120000,999),
                               names=paste0("dnase:", letters[1:3])), score=4:6)
gr1
gr2
ranges <- subsetByOverlaps(gr2,gr1)
ranges
hits <- findOverlaps(gr2, gr1)
hits
rsid <- CharacterList(split(names(gr1)[subjectHits(hits)],
                            queryHits(hits)))
rsid
snpscore <- CharacterList(split(gr1$score[subjectHits(hits)],
                                queryHits(hits)))
snpscore
mcols(ranges) <- DataFrame(mcols(ranges), rsid, snpscore)
mcols
ranges
# Somatic Signatures ------------------------------------------------------

library(SomaticSignatures)
library(SomaticCancerAlterations)
biocLite("BSgenome.Hsapiens.1000genomes.hs37d5")
library(BSgenome.Hsapiens.1000genomes.hs37d5)
sca_metadata = scaMetadata()
sca_metadata
sca_data = unlist(scaLoadDatasets())
sca_data$study = factor(gsub("(.*)_(.*)", "\\1", toupper(names(sca_data))))
sca_data = unname(subset(sca_data, Variant_Type %in% "SNP"))
sca_data = keepSeqlevels(sca_data, hsAutosomes())
sca_vr = VRanges(
  seqnames = seqnames(sca_data),
  ranges = ranges(sca_data),
  ref = sca_data$Reference_Allele,
  alt = sca_data$Tumor_Seq_Allele2,
  sampleNames = sca_data$Patient_ID,
  seqinfo = seqinfo(sca_data),
  study = sca_data$study)
sca_vr
sort(table(sca_vr$study), decreasing = TRUE)

sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
head(sca_motifs)
sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

head(round(sca_mm, 4))
plotMutationSpectrum(sca_motifs, "study")


# GTOM Horvath ------------------------------------------------------------

# This function computes the GTOMm DISSIMILARITY
#
# Input:
# - adjmat1, a symmetric adjacency matrix with binary entries
# - m, the order of GTOM
#
# Output:
# - The GTOMm dissimilarity matrix
#
# Andy M. Yip and Steve Horvath
# January 2005

if(exists("GTOMmdist1")) rm(GTOMmdist1);
GTOMmdist1 = function(adjmat1,m=1){
  if (m!=round(abs(m))){
    stop("m must be a positive integer!!!", call.=TRUE);}
  if (any(adjmat1!=0 & adjmat1!=1)){
    stop("The adjacency matrix must be binary!!!", call.=TRUE);}
  
  B <- adjmat1;
  if (m>=2) {
    for (i in 2:m) {
      diag(B) <- diag(B) + 1;
      B = B %*% adjmat1;}}   # number of paths with length at most m connecting each pair
  B <- (B>0);                    # m-step reachability matrix
  diag(B) <- 0;                  # exclude each node being its own neighbor
  B <- B %*% B;                  # number of common k-step neighbors that each pair of nodes share
  
  Nk <- diag(B);                 # number of common k-step neighbors that each node possesses
  B <- B +adjmat1;
  diag(B) <- 1;
  denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
  diag(denomTOM) <- 1;
  1 - B/denomTOM                 # turn the GTOM matrix into a dissimilarity
}

setwd("~/NGS/HCC.NGS.CNV/example")
#https://labs.genetics.ucla.edu/horvath/GTOM/old/
#steve horvath!
#
# Read in the expression levels
datExpr = read.csv("YEASTCellCycle4000.csv", header=T, row.names=1)
datExpr = t(datExpr)

# Read in the gene information
datInfo = read.csv("YEASTCellCycle4000_geneinfo.csv", header=T, row.names=1)

# Correlation matrix
cormat = cor(datExpr,use="p")

# Adjacency matrix
tau = 0.7
AdjMat1 <- as.matrix(abs(cormat)>=tau)
diag(AdjMat1)=0

# Restrict our attention to 1000 highly connected genes
Degree = apply(AdjMat1,1,sum)
DegreeRank = rank(-Degree)
keep.node = DegreeRank <= 1000
AdjMat1.1000 = AdjMat1[keep.node, keep.node]
datExpr.1000 = datExpr[, keep.node]
datInfo.1000 = datInfo[keep.node, ]

# Compute the GTOMm dissimilarity (m=2 in this example)
dissGTOM2 = GTOMmdist1(AdjMat1.1000, m=2)

# Hierarchical Clustering
hierGTOM2 = hclust(as.dist(dissGTOM2),method="average")
plot(hierGTOM2, main="GTOM2 Dissimilarity Measure", labels=F, xlab="", sub="")

# BioMART Example ---------------------------------------------------------

# Code snippet to get chromosomal position and cytoband from gene symbols using
# biomaRt
# Kenneth Daily, 2014

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Only use standard human chromosomes
normal.chroms <- c(1:22, "X", "Y", "M")

# Filter on HGNC symbol and chromosome, retrieve genomic location and band
my.symbols <- c("RB1", "TP53", "AKT3")

my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                    mart = ensembl)


# RANDOM STUFF ------------------------------------------------------------

##  Print items in a list 1 per line in the Rstudio output
cat(print(QT.LN.genes.in.both.filtered.and.unfiltered, row.names = FALSE, quote = FALSE), sep="\n")

cat(1:5,sep="\n")
1
2
3
4
5

#This is particularly useful when you need a line-by-line list of the variables in a data frame, which you can get with:
  
  cat(names(dataframe),sep="\n")
  
  
  -----------
    

getGenes(head(entrezID), fields="ensembl.gene")
#Doesn't result in a lot of objects :(

#Trying biomart
listMarts()
#hmmmm...

require(org.Hs.eg.db)
x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

my_genes <- c("LRRK2","PINK1")

#symbol and OMIM
a <- select(org.Hs.eg.db,
            keys = my_genes,
            columns=c("ENTREZID", "SYMBOL","OMIM"),
            keytype="SYMBOL")
a
b <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
            keys = a$ENTREZID,
            columns=c('GENEID', 'TXCHROM', 'TXSTART', 'TXEND', 'TXID', 'ENSEMBL'),
            keytype="GENEID")
b


#Use a random sample for your exploratory analysis or to test code. Here is a code snippet that will give you a convenient function: row.sample(yourdata, 1000) will reduce your massive file to a random sample of 1,000 observations.

biocLite("tigerstats")


# MACHINE LEARNING --------------------------------------------------------


# TUTORIAL FROM BIOSTARS --------------------------------------------------

#https://www.biostars.org/p/87580/


# OneR package ------------------------------------------------------------
#https://cran.r-project.org/web/packages/OneR/vignettes/OneR.html
iris <- iris
data <- optbin(iris)
model <- OneR(data, verbose = T)
#Attribute    Accuracy
#1 * Petal.Width  96%     
#2   Petal.Length 95.33%  
#3   Sepal.Length 74.67%  
#4   Sepal.Width  55.33%  
#---
#  Chosen attribute due to accuracy
#and ties method (if applicable): '*'
summary(model)
prediction <- predict(model, data)
eval_model(prediction, data)

#Breast cancer classification
data("breastcancer")
data <- breastcancer
set.seed(12) # for reproducibility
random <- sample(1:nrow(data), 0.8 * nrow(data))
data_train <- optbin(data[random, ], method = "infogain")
data_test <- data[-random, ]
model_train <- OneR(data_train, verbose = TRUE)
summary(model_train)

# IMPUTATION --------------------------------------------------------------

#Median from mlr package
imputeMedian()

#From: https://www.kaggle.com/mrisdal/titanic/exploring-survival-on-the-titanic
#Replace missing values with median
# Replace missing fare value with median fare for class/embarkment
full$Fare[1044] <- median(full[full$Pclass == '3' & full$Embarked == 'S', ]$Fare, na.rm = TRUE)


# Multiplot function ------------------------------------------------------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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



# CBioPortal R code -------------------------------------------------------

biocLite("cgdsr")
library(cgdsr)

#From Website
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)
    #getCancerStudies...  OK
    #getCaseLists (1/2) ...  OK
    #getCaseLists (2/2) ...  OK
    # etc, etc, etc... IT WORKS!

# Get list of cancer studies at server
studs <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[2,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list
getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

# documentation
help('cgdsr')
help('CGDS')


#Trying it out:
# Get list of cancer studies at server
studs <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = "lihc_tcga" #Older method
mycancerstudy = getCancerStudies(mycgds)[70,1] #Better method
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list
getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[8,1]
####Does not have synonymous mutation data!

#Figuring out what other info is out there:
mycaselist = getCaseLists(mycgds, mycancerstudy)
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)

lihc.linear.CNA <- getGeneticProfiles(mycgds, mycancerstudy)[2,1]
cna <- getProfileData(mycgds,


# FirebrowseR -------------------------------------------------------------

biocLite("FirebrowseR")

# COSMIC Data Files -------------------------------------------------------

#Reading in my custom cosmic data file
cosmic.muts <- fread("~/NGS/references/COSMIC/cosmic.mutations.by.site.txt")
  #26384 genes, with # of mutations total and # mutations per site.


# MAFTOOLS ----------------------------------------------------------------


# * *Somatic Sigs in maftools ----------------------------------------

###Important note from NMF package:
  #2.4 Multiple runs
  #When the seeding method is stochastic, multiple runs are usually required to achieve stability or a
  #resonable result. This can be done by setting argument nrun to the desired value. For performance
  #reason we use nrun=5 here, but a typical choice would lies between 100 and 200:

###Different approaches to deciding on factorization rank r (= # of sigs):
#Several approaches have then been proposed to choose the optimal value of r. For example,
#(Brunet et al. 2004) proposed to take the first value of r for which the cophenetic coefficient starts
#decreasing, (Hutchins et al. 2008) suggested to choose the first value where the RSS curve presents
#an inflection point, and (Frigyesi et al. 2008) considered the smallest value at which the decrease
#in the RSS is lower than the decrease of the RSS obtained from random data.

###From https://confluence.broadinstitute.org/display/GDAC/Analysis+Run+Release+Notes
#How NMF clustering is done at BROAD
#The cophenetic correlation coefficients and average silhouette values are used to determine the k with the most robust #clusterings. From the plot of cophenetic correlation versus k, we select modes and the point preceding the greatest decrease #in cophenetic correlation coefficient, and from these choose the k with the highest average silhouette value.


#Example Code:
write.table(hcc.maf, file = 'mut.maf', quote = FALSE, sep = '\t', row.names = FALSE)
maf = read.maf(maf = "mut.maf", removeSilent = FALSE, useAll = FALSE)
#Creating trinucleotide Matrix
hcc.tnm = trinucleotideMatrix(maf = maf, ref_genome = '/Users/colincourt/NGS/references/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa', prefix = 'chr', add = TRUE, ignoreChr = 'chr23', useSyn = TRUE)
#Creating Signatures
hcc.sig = extractSignatures(mat = hcc.tnm, nTry = 10, plotBestFitRes = TRUE)
#Plotting Signatures
maftools::plotSignatures(hcc.sig) 

#Loading in the Alexandrov Signatures:
alex.sigs <- fread("~/NGS/references/Alexandrov.sigs/Alexandrov_signatures.txt")
#Cleaning it up:
alex.sig.all <- alex.sigs
alex.sig.all <- inject.dots(alex.sig.all) #Run twice
colnames(alex.sig.all)
alex.sig.all$sub.type <- gsub(">", "", alex.sig.all$Substitution.Type)
alex.sig.all$first <- str_sub(alex.sig.all$Trinucleotide, 1, 1)
alex.sig.all$middle <- "."
alex.sig.all$last <- str_sub(alex.sig.all$Trinucleotide, -1, -1)
alex.sig.all$type <- with(alex.sig.all, paste0(alex.sig.all$first, alex.sig.all$middle, alex.sig.all$last))
alex.sig.all$sigs <- with(alex.sig.all, paste0(alex.sig.all$sub.type, " ", alex.sig.all$type))
#Do this to get rid of the messy columns
alex.sig.all <- alex.sig.all[, c(36, 4:30)]
rownames(alex.sig.all) <- alex.sig.all$sigs
alex.sig.all$sigs <- NULL


#NMF PACKAGE: How MAFTOOLS reads in alex.sigs -------------------------------

sigs = data.table::fread(input = system.file('extdata', 'signatures.txt', package = 'maftools'), stringsAsFactors = FALSE, data.table = FALSE)
colnames(sigs) = gsub(pattern = ' ', replacement = '_', x = colnames(sigs))
rownames(sigs) = sigs$Somatic_Mutation_Type
sigs = sigs[,-c(1:3)]
sigs = sigs[,1:22] #use only first 21 validated sigantures
sigs = sigs[rownames(w),]

###How program works
#1) transposes maf -> mat matrix 
mat <- t(hcc.tnm) 
#mat = motifs as rows and samples as columns

#Estimating number of signatures to use
nmfTry = nmfEstimateRank(mat, seq(2,6), method='brunet', nrun=10, seed=123456)

#Summary of NMF try
nmf.sum = summary(nmfTry)

#Getting point where cophenetic correlation coefficient starts decreasing
nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing
  #bestfit doesnt work: Error in x[j] : invalid subscript type 'closure'
n = 5

#Performing NMF using the optimum number of signatures
conv.mat.nmf = NMF::nmf(x = mat, rank = n) #Using rank = 5

#Signatures
w = NMF::basis(conv.mat.nmf)
w = apply(w, 2, function(x) x/sum(x)) #Scale the signatures (basis)
colnames(w) = paste('Signature', 1:ncol(w),sep='_')

#Contribution
h = NMF::coef(conv.mat.nmf)
#Code to use, actual code below
h = apply(h, 2, function(x) x/sum(x)) #Scale contributions (coefs)
rownames(h) = paste('Signature', 1:nrow(h),sep='_')

#For single signature, contribution will be 100% per sample
if(n == 1){
  h = h/h
  rownames(h) = paste('Signature', '1', sep = '_')
}else{
  h = apply(h, 2, function(x) x/sum(x)) #Scale contributions (coefs)
  rownames(h) = paste('Signature', 1:nrow(h),sep='_')
}

  #Doing Cosine similarity one signature at a time
sig = w[,1] #96 motifs, percentage contribution per motif
a.sig = sigs[,1]
x = sigs[,1]
coSineMat = rbind(coSineMat, apply(sigs, 2, function(x){
  crossprod(sig, x)/sqrt(crossprod(x) * crossprod(sig))
}))   #This is rows = w matrix signatures, columns = Alexandrov signatures
?crossprod
#This is formally equivalent to (but usually slightly faster than) the call t(x) %*% y (crossprod) or x %*% t(y) (tcrossprod).
#Multiplies the equivalent number in each signature (e.g. A[C>A]A)
a <- t(sig) %*% a.sig # == 0.0113
b <- sqrt(t(sig) %*% a.sig) # == 0.106
c <- t(sig) %*% sig # == 0.0168
a / (b * c) # == 6.35???
a1 = crossprod(sig, x) #= 0.0113
b1 = sqrt(crossprod(x)) # = 0.248
c1 = crossprod(sig) # =0.0168
a1/(b1*c1) # = 2.72???


###Actual MAFTOOLS code
#corMat = c()
coSineMat = c()
for(i in 1:ncol(w)){
  sig = w[,i]
  coSineMat = rbind(coSineMat, apply(sigs, 2, function(x){
    crossprod(sig, x)/sqrt(crossprod(x) * crossprod(sig)) #Estimate cosine similarity against all 21 signatures
  }))
  #corMat = rbind(corMat, apply(sigs, 2, function(x) cor.test(x, sig)$estimate[[1]])) #Calulate correlation coeff.
}
#rownames(corMat) = colnames(w)
rownames(coSineMat) = colnames(w)

# for(i in 1:nrow(corMat)){
#   message('Found ',rownames(corMat)[i], ' most similar to validated ',names(which(corMat[i,] == max(corMat[i,]))), '. Correlation coeff: ', max(corMat[i,]), sep=' ')
# }

for(i in 1:nrow(coSineMat)){
  message('Found ',rownames(coSineMat)[i], ' most similar to validated ',names(which(coSineMat[i,] == max(coSineMat[i,]))), '. CoSine-Similarity: ', max(coSineMat[i,]), sep=' ')
}

####THIS IS WHAT BECOMES THE EXTRACT SIGNATURES OBJECT!!!
return(list(signatures = w, contributions = h, coSineSimMat = coSineMat, nmfObj = conv.mat.nmf))
}


    ###   Better way of figuring out how many signatures to use:
# perform 10 runs for each value of r in range 2:6
estim.r <- nmf(esGolub, 2:6, nrun = 10, seed = 123456)
plot(estim.r) #Plots the different points to see what the optimum should be


# Bad code that doesnt work -----------------------------------------------

#Trying with factors to get the order right -- NOT WORKING
w_df <- melt(alex.sig.all)
w_df$X2 <- str_sub(w_df$X2, 11)
w_df$X3 <- factor(w_df$X2, levels = x, is.ordered(w_df$X2))
x <- unique(w_df$X2)
class(w_df$X2)
p <- ggplot(w_df)
p = p + geom_bar(aes_string(x = "X1", y = "value", fill = "X2"), stat = "identity", position = "identity")
p = p + facet_grid(X2 ~ .)
p = p + xlab("Motif") + ylab("Contribution")
p



###Trying to undo strings as factors .... MUCH better to just import it!
#Creating Alex.sigs.all file for all 21 signatures
alex.sig.all <- alex.sig[, c(1:25)]
alex.sig.all$sub.type <- gsub(">", "", alex.sig.all$`Substitution Type`)
alex.sig.all$first <- str_sub(alex.sig.all$Trinucleotide, 1, 1)
alex.sig.all$middle <- "."
alex.sig.all$last <- str_sub(alex.sig.all$Trinucleotide, -1, -1)
alex.sig.all$type <- with(alex.sig.all, paste0(alex.sig.all$first, alex.sig.all$middle, alex.sig.all$last))
alex.sig.all$sigs <- with(alex.sig.all, paste0(alex.sig.all$sub.type, " ", alex.sig.all$type))

alex.sigs.all <- alex.sig.all[, c(31, 4:25)]
rownames(alex.sigs.all) <- alex.sigs.all$sigs
alex.sigs.all$sigs <- NULL
#alex.sigs.all = motifs as row names and signatures as columns = Works!
alex.sigs.all2 <- as.data.frame(alex.sigs.all)
colnames(alex.sigs.all2) <- str_sub(colnames(alex.sigs.all2), 11)
colnames(alex.sigs.all2) <- as.character(colnames(alex.sigs.all2))
colnames(alex.sigs.all2) <- gsub("^", "S", colnames(alex.sigs.all2))
alex.sigs.all2 <- transform(alex.sigs.all2, class=as.numeric(as.character(alex.
#Try re-reading in the signatures with stringsasfactors = FALSE!
alex.sigs.all2 <- as.data.frame(alex.sigs.all2)
alex.sigs.all2$1A
alex.sigs.all2 <- as.data.frame(sapply(alex.sigs.all, as.character))
alex.sigs.all2 <- data.matrix(alex.sigs.all2)
colSums(alex.sigs.all2)
class(alex.sigs.all2$`1A`
colnames(alex.sigs.all2) <- as.character(colnames(alex.sigs.all2))



# printing out list of rowsums --------------------------------------------

rowSums(H.tcga) %>% order(decreasing = TRUE) %>% print()
#Same as:
print(order(rowSums(H.tcga)))
rowSums(H.tcga) %>% cat() %>% order() #Nope
