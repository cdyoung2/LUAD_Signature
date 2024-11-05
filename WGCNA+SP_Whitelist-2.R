##########################################################################################################
## Eric Dammer - adapted code for WGCNA from Neelroop Parikshak, Vivek Swarup, and Divya Nandakumar
## SeyfriedLab&ProteomicsCorePipeline.R
##
## Shared with Morehouse students Corey, Kaylin, Tiara, and Courtney, Shakayla in James Lillard Lab, Morehouse School of Medicine
##  -applied to Tiara's curation of RNA-seq on CLL leukemia Feb 20, 2019
##
## Goal: This code will prepare FPKM data for coexpression network analysis and systems biology.
##########################################################################################################

## CODE BLOCKS:
# Set up environment, load data,
# Clean up abundance matrix of rows/gene products with too many 0 FPKM/TPM or NA protein abundance values
# Perform TAMPOR robust median polish removal of unwanted batch covariance (preserving wanted biological variance)


#########################################################################################################
# Prepare the workspace.

rm(list=ls(all=TRUE))
dev.off()
f <- "\f"
cat(f) # cat("\014") #alt= > cat("\f")
options(stringsAsFactors = FALSE)

# (Always run the below code after loading a saved.image .Rdata)
#########################################################################################################
# Set up environment

## Set folders and libraries
options(stringsAsFactors=FALSE)
rootdir <- "/Users/coyoung/Projects/WGCNA/" # This is the folder containing all of the analysis scripts, input, and output for this project
functiondir <- "Code"
datadir <- "Input/FPKM"
outputfigs <- "Figures"
outputtabs <- "Tables"
rdata <- "RData"
pca <- "PCA"
iGraphs <- "Input/iGraphs"

# Set the working directory and subfolders for input/output.
setwd(rootdir)

functiondir = paste0(rootdir,functiondir,"/")
datadir = paste0(rootdir,datadir,"/")
outputfigs = paste0(rootdir,outputfigs,"/")
outputtabs = paste0(rootdir,outputtabs,"/")
rdata = paste0(rootdir,rdata,"/")
pca = paste0(rootdir,pca,"/")
iGraphs = paste0(rootdir,iGraphs,"/")
dir.create(file.path(outputfigs)) #harmless warning given if already exists
dir.create(file.path(outputtabs))
dir.create(file.path(datadir))
dir.create(file.path(rdata))
dir.create(file.path(iGraphs))

## Load any libraries used in the following code blocks
library(WGCNA) # Network analysis package
library(NMF) # this package has a great annotated heatmap function - aheatmap
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Cairo) # nicer graphics, anti-aliased, etc. --text from windows output PDFs using CairoPDF() function may not load in Illustrator, though -- so also use the pdf() standard output function when generating PDF figures
##Only for macs:
#CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic="Arial:style=Bold Italic,BoldItalic",symbol="Symbol")

#note: other libraries are loaded as needed in certain blocks but also listed below for completeness
library(reshape2)     #TAMPOR block
library(gridExtra)    #TAMPOR block
library(ggpubr)       #TAMPOR block
library("doParallel") #Bootstrap Regression block and GlobalNetworkPlots
library("biomaRt")    #Fisher Exact Test/list overlap block (enables cross-species lookup)
library("NMF")        #WGCNA GlobalNetworkPlots - eigengene heatmap and clustering
library("plotly")     #Volcano Figure Generation
library("stringr")    #GO-Elite block
library(boot)         #bootstrap regression block
library(tictoc)

## Load the sample info.

#Check and change your input filenames below
inputAbundanceFile="InitialFPKM_533.csv"        #abundance.CSV: rows are genes/proteins, columns are samples, this may change based on input (e.g. entire LUAD project, or focus on primary diagnoses)
inputTraitsFile="final-clinical.csv"     #traits.CSV: rows are samples, columns are traits (as many as possible should be numerically coded), the tcga ids are ordred by stage
#IMPORTANT: Make sure your inputTraitsFile has a "Group" column calling out (expected) different subsets of case-samples. Groups can be numerically coded if a single severity scale applies, but each group (number) should be used for at least 3 samples.

##Part of most output filenames specific to this project
projectFilesOutputTag="TCGA_LuCa"


## Load and clean the data
#########################################################################################################

dat <- read.csv(file=paste0(datadir,inputAbundanceFile),header=TRUE,row.names=1)
dat.original = dat
range(dat) #(0.0 - 43812.4) #get a sense of whether these are all positive >0, doesn't work?
#dat<-log2(dat) #only if input protein abundance data not already log2 transformed
#colnames(dat)<-gsub("\\.","-",colnames(dat)) # one can replace . with -
dim(dat) #60483, 502

## Match up metadata and plot it for visualization - we want to know what variables are correlated with other variables
numericMeta <- read.csv(paste0(datadir,inputTraitsFile),header=TRUE,row.names=1)
dim(numericMeta) #522, 157
rownames(numericMeta)==colnames(dat) #sanity check -- are sample names in same order?

#dat<-dat[,na.omit(match(rownames(numericMeta),colnames(dat)))] #cull to keep only samples we have traits for, and match the column (sample) order of dat to the row order of numericMeta
numericMeta <- numericMeta[match(colnames(dat),rownames(numericMeta)),] #use this line instead of above if you have more samples in your traits file than you do in abundance data; trait sample (row) order will be matched to column names of dat; or use both lines if matching needs to be enforced due to different missing samples in both traits and abundance data
rownames(numericMeta)==colnames(dat) #sanity check -- are sample names in same order?
dim(dat) #60483, 502

traits<-numericMeta   #used after normalization loop completes, only.

## Replacing all zeros (0's) w/ NA in cleadDat
datZeroPositions <- dat
dat[dat==0] <- NA

## Removal of samples with NA's in traits that are using in regression code later (i.e.regvars) 
# No NAs in regvars, they will be prograted throughout the sample, others should not be effected. Code below removes samples w/NA in Group column (Group column can & may be changed)
numericMeta2<-numericMeta[!is.na(numericMeta$Group),] # removing rows with NA in for condition
numericMeta = numericMeta2
dat<-dat[,na.omit(match(rownames(numericMeta),colnames(dat)))] # rows in numericaMeta == columns in cleanDat
rownames(numericMeta)==colnames(dat) #sanity check -- are sample names in same order?
dim(dat) #60483, 494

numericMeta2<-numericMeta[!is.na(numericMeta$age_at_initial_pathologic_diagnosis),] 
numericMeta = numericMeta2
dat<-dat[,na.omit(match(rownames(numericMeta),colnames(dat)))] 
rownames(numericMeta)==colnames(dat)
dim(dat) #60483, 475

numericMeta2<-numericMeta[!is.na(numericMeta$pathologic_t.1),] 
numericMeta = numericMeta2
dat<-dat[,na.omit(match(rownames(numericMeta),colnames(dat)))] 
rownames(numericMeta)==colnames(dat)
dim(dat) #60483, 472

#########################################################################################################
#=======================================================================================================#
                          #  TAMPOR --> followed by 50% missingness code  #
#=======================================================================================================#
#########################################################################################################
#DID NOT RUN TAMPOR: Run TAMPOR on all biological samples as denominator (noGIS=TRUE option; use case #4)
## Obtain current TAMPOR.R source from https://www.github.com/edammer/TAMPOR
# ====================================================================================================== 

# source('/Users/coyoung/Projects/WGCNA/Code/TAMPOR.R')
# TAMPORlist <- TAMPOR(dat, traits, noGIS=TRUE, batchPrefixInSampleNames=FALSE, parallelThreads=8, minimumBatchSize=4)
# 
# ## Save TAMPOR output and full session image
# saveRDS(TAMPORlist, file="TAMPORlist.rds")
# save.image(paste0(rdata,"TAMPORcompleted_",projectFilesOutputTag,".Rdata"))

#TAMPORlist<-readRDS("TAMPORlist.rds")  #you can load an existing TAMPORlist output if saved previously, instead of running the above 4 lines of code.

# ## Output PDF of 3 captured graphical output pages
# pdf(paste0(outputfigs,"TAMPOR_graphicsOut_",projectFilesOutputTag,".pdf"),width=18,height=12)
# TAMPORlist$meanSDplots
# TAMPORlist$MDSplots
# TAMPORlist$convergencePlots
# dev.off()

#Prepare to plot histograms comparing log2(FPKMs) before and after TAMPOR
## Dat original FPKM data of 0 values and log2 transform
dat.Original.0asNA<-dat.original
dat.Original.0asNA[dat.Original.0asNA==0]<-NA
dat.Original.0asNA<-log2(dat.Original.0asNA)

# ## Set cleanDat to post-TAMPOR normalized FPKM
# cleanRelAbun <- TAMPORlist$cleanRelAbun
# numericMeta <- TAMPORlist$traits
# traits <- TAMPORlist$traits

cleanRelAbun <- dat

## Enforce <50% missingness (1 less than half of cleanRelAbun columns (or round down half if odd number of columns))
dim(cleanRelAbun) #(60483, 472)
LThalfSamples<-length(colnames(cleanRelAbun))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanRelAbun)) %% 2)==1) { 0.5 } else { 1.0 }

# If operating on log2(FPKM) data, remove rows with >=50% originally 0 FPKM values (only if there are rows to be removed)
# note- TAMPOR has already done this so dim(cleanRelAbun) before running this should match dim(cleanRelAbun) after running.
# note 2 - only applies here if samples were dropped after TAMPOR, in which case some rows could have >50% missing values again
dim(cleanRelAbun) #60483, 472
temp2<-data.frame(ThrowOut=apply(cleanRelAbun,1,function(x) length(x[x<=0])>LThalfSamples))
if (nrow(temp2) > 1) cleanRelAbun<-cleanRelAbun[!temp2$ThrowOut,]  #run only if there are >1 rows to be removed; error if temp2 is a vector (1 row to remove)
dim(cleanRelAbun)  #31548, 472

save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#########################################################################################################
#=======================================================================================================#
                            #  Check Histograms & Distribution of Data  #
#=======================================================================================================#

## Output PDF of histograms
# Nov-9-2019 Courtney's code via Dammer to print out all of the plots associated with TAMPOR, check Courtney's presentation for explanation

# run page 1 code before opening PDF, look at comment for page 2 code, and adjust pdf page 2 variables according to the data being plotted

# page 1

# Histogram original FPKM data, without 0 values, log2 transformed
hist(as.double(unlist(dat.Original.0asNA)),breaks=200,main=paste0("Pre-Missingness log2(non-zero FPKM data)\n",nrow(dat.Original.0asNA)," transcripts"),xlab="log2(FPKM)")
# Histogram of original FPKM data, without 0 values, and selecting only which rows that survived missingness control
hist(as.double(unlist(dat.Original.0asNA[match(rownames(cleanRelAbun),rownames(dat.Original.0asNA)),match(colnames(cleanRelAbun),colnames(dat.Original.0asNA))])),breaks=200,main=paste0("Post-Missingness Matched log2(non-zero FPKM data)\n",nrow(cleanRelAbun)," transcripts only present post-50% Zero control"),xlab="log2(FPKM)")
# Histogram of log2 transformed values, exactly the same as previous plot (both sets of values are log2 transformed, only cleanRelAbun data is shown)
hist(log2(as.double(unlist(cleanRelAbun))),breaks=200,main=paste0("Post-Missingness log2(non-zero FPKM data)\n",nrow(cleanRelAbun)," transcripts"),xlab="log2(relative abundance)")

# Histogram of after all modficications, last histogram. Use this histogram to inform how to get your cleanDat data normal density peak centered at 0. Follow directions below:

# set the last hist() -> variable, then max(variable$density) will show you the x-value bin midpoint 
# subtract that decimal bin center from all cleanDat matrix values
# max(blah$density) = 0.2238862

hist(log2(as.double(unlist(cleanRelAbun.FINAL))),
     breaks=200,main=paste0("Post-Missingness log2(non-zero Normalized FPKM data)\nONLY ROWS with >50% values above 1\n",nrow(cleanRelAbun.FINAL)," transcripts"),
     xlab="log2(relative abundance)")

pdf(file=paste0(outputfigs,"NoTAMPOR-Histograms-Pre&PostMissingness-60kVs30kTranscripts.pdf"), width = 24, height = 8)

#PDF row/column count and margin settings
par(mfrow=c(1,3))
par(mar=c(5,3,5,1))

# page 1 (copy of above code)

# Histogram original FPKM data, without 0 values, log2 transformed
hist(as.double(unlist(dat.Original.0asNA)),breaks=200,main=paste0("Pre-Missingness log2(non-zero FPKM data)\n",nrow(dat.Original.0asNA)," transcripts"),xlab="log2(FPKM)")
# Histogram of original FPKM data, without 0 values, and selecting only which rows that survived missingness control
hist(as.double(unlist(dat.Original.0asNA[match(rownames(cleanRelAbun),rownames(dat.Original.0asNA)),match(colnames(cleanRelAbun),colnames(dat.Original.0asNA))])),breaks=200,main=paste0("Post-Missingness Matched log2(non-zero FPKM data)\n",nrow(cleanRelAbun)," transcripts only present post-50% Zero control"),xlab="log2(FPKM)")
# Histogram of log2 transformed values, exactly the same as previous plot (both sets of values are log2 transformed, only cleanRelAbun data is shown)
hist(log2(as.double(unlist(cleanRelAbun))),breaks=200,main=paste0("Post-Missingness log2(non-zero FPKM data)\n",nrow(cleanRelAbun)," transcripts"),xlab="log2(relative abundance)")

# page 2

# We will look at partial data histograms to answer whether your data has unusual properties that would prevent us from removing rows with 50%+ values < x FPKM; x=1 typically
# Are values of transcripts (rows) with a median near the left peak split between left and right peaks? 
# (Check output of page1 code to screen without PDF open), and adjust below 2 values...)
valueLeftOfLeftPeak=-2.0
valueNearRightPeak=2.6

# Histogram of transcripts with FPKM primarily contributing to the left tail of left (noise) peak
hist(log2(as.double(unlist(cleanRelAbun[which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))< valueLeftOfLeftPeak),]))),
     breaks=200,main=paste0("Are values of transcripts with a median value near the left peak\nsplit between left and right peaks?\npost-Missingness log2(non-zero FPKM data)\n",length(which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))< valueLeftOfLeftPeak))," transcripts"),
     xlab="log2(relative abundance)")
# Histogram of transcripts with FPKM primarily contributing to the right tail of right (true signal) peak
hist(log2(as.double(unlist(cleanRelAbun[which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))> valueNearRightPeak),]))),
     breaks=200,main=paste0("Are values of transcripts with a median near the right peak\nsplit between left and right peaks?\npost-Missingness log2(non-zero FPKM data)\n",length(which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))> valueNearRightPeak))," transcripts"),
     xlab="log2(relative abundance)")
# Histogram of transcripts with median FPKM values between the left and right peaks
hist(log2(as.double(unlist(cleanRelAbun[which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))> valueLeftOfLeftPeak & apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))< valueNearRightPeak),]))),
     breaks=200,main=paste0("Are values with a median in between the two peaks\nin the full data split between left and right peaks?\npost-Missingness log2(non-zero FPKM data)\n",length(which(apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))> valueLeftOfLeftPeak & apply(log2(cleanRelAbun),1,function(x) median(x,na.rm=TRUE))< valueNearRightPeak))," transcripts"),
     xlab="log2(relative abundance)")

# page 3
# Histogram of all non-zero normalized FPKM data. Same as plot 3 on page 1.
hist(log2(as.double(unlist(cleanRelAbun))),
     breaks=200,main=paste0("Post-Missingness log2(non-zero Normalized FPKM data)\n",nrow(cleanRelAbun)," transcripts"),
     xlab="log2(relative abundance)")
# Reprint of last histogram, selecting only values for Normalized FPKM > (greater than) 1, even in rows with fewer than 50% values <=1 (psot figure generation: had to had GT1 to nrow(cleanRelAbun) b/c numbers were incorrect)
cleanRelAbun.GT1 <- cleanRelAbun
cleanRelAbun.GT1[cleanRelAbun<=1]<- NA
hist(log2(as.double(unlist(cleanRelAbun.GT1))),
     breaks=200,main=paste0("Post-Missingness log2(non-zero Normalized FPKM data)\nONLY VALUES > 1\n",nrow(cleanRelAbun.GT1)," transcripts"),
     xlab="log2(relative abundance)")
# Final histogram, selecting only rows of cleanRelAbun post-TAMPOR (if ran) with <50% values <=1
# ...first prepare data to enforce <50% of row values meeting below criteria
# (<50% is 1 less than half of cleanRelAbun columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanRelAbun))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanRelAbun)) %% 2)==1) { 0.5 } else { 1.0 }
# Operating on Normalized FPKM data which has NAs instead of 0, and we remove rows with >=50% FPKM values <=1 (only if there are rows to be removed)
cleanRelAbun.FINAL <- cleanRelAbun
cleanRelAbun.FINAL[is.na(cleanRelAbun)]<-0
temp3<-data.frame(ThrowOut=apply(cleanRelAbun.FINAL,1,function(x) length(x[x<=1])>LThalfSamples))
### Run in NoTAMPOR code: 
if (nrow(temp3) > 1) cleanRelAbun.FINAL<-cleanRelAbun.FINAL[!temp3$ThrowOut,]   #run only if there are >1 rows to be removed; error if temp3 is a vector (1 row to remove)
#put back NAs for zeroes
cleanRelAbun.FINAL[cleanRelAbun.FINAL==0] <- NA
#plot histogram
hist(log2(as.double(unlist(cleanRelAbun.FINAL))),
     breaks=200,main=paste0("Post-Missingness log2(non-zero Normalized FPKM data)\nONLY ROWS with >50% values above 1\n",nrow(cleanRelAbun.FINAL)," transcripts"),
     xlab="log2(relative abundance)")


dev.off()

dim(cleanRelAbun.FINAL) #13486, 472
cleanDat <- log2(cleanRelAbun.FINAL) # does not need to be run b/c log2(cleanRelAbun.FINAL) will happen in sourced code below, unless commented out

# Only run if you want you have FPKM centered around 0 for histograms 
#source('/Users/coyoung/Projects/WGCNA/Code/Sourced/Centered-FPKMs.R') 
#cleanDat.centered <- cleanRelAbun.centered
save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#########################################################################################################
#=======================================================================================================#
                                      #  Check and Remove Outliers  #
#=======================================================================================================#

if(!exists("numericMeta")) numericMeta<-traits
library(WGCNA)

sdout=3 #Z.k SD fold for outlier threshold
outliers.noOLremoval<-outliers.All<-vector()
cleanDat.noOLremoval<-cleanDat

for (repeated in 1:5) {
  normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)
  
  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))
  targets.All=numericMeta
  
  cleanDat <- cleanDat[,!outliers] 
  numericMeta <- targets <- targets.All[!outliers,]
  outliers.All<-c(outliers.All,outliers)
} #repeat 5 times

#All outliers removed
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
#"There are 7 total outlier samples removed in 5 iterations:"
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
dim(cleanDat) #13486, 465
save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))

#########################################################################################################
#=======================================================================================================#
                                  #  Parallel Bootstrap Regression  #
#=======================================================================================================#

# Changing gender binary coding from 0/1 to 1/2
gender.2<-rep(NA,nrow(numericMeta))
gender.2[numericMeta$gender.1==0]<-1
gender.2[numericMeta$gender.1==1]<-2
#View(gender.2)
numericMeta$gender.2 <- gender.2
#View(numericMeta$gender.2)

cleanDat.unreg<-cleanDat
dim(numericMeta) #492, 158

library("doParallel")
#when Eric is running at Emory (requires RSA public key-ssh, & manual run of shell script from command prompt to start server backend):  
#parallelThreads=30
#clusterLocal <- makeCluster(c(rep("haplotein.biochem.emory.edu",parallelThreads)), type = "SOCK", port=10191, user="edammer", rscript="/usr/bin/Rscript",rscript_args="OUT=/dev/null SNOWLIB=/usr/lib64/R/library",manual=FALSE)
##OR to run parallel processing threads for regression: :
parallelThreads=15 #max is number of processes that can run on your computer at one time
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")

registerDoParallel(clusterLocal) 

##OR for no parallel processing skip all doParallel functions (much slower): #parallelThreads=1

library(boot)
boot <- TRUE
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  

# # No NAs in regvars, they will be prograted throughout the sample, others should not be effected. Code below removes samples w/NA in Group column (Group column can & may be changed)
# numericMeta2<-numericMeta[!is.na(numericMeta$Group),] # removing rows with NA in for condition
# numericMeta = numericMeta2
# cleanDat<-cleanDat[,na.omit(match(rownames(numericMeta),colnames(cleanDat)))] # rows in numericaMeta == columns in cleanDat
# rownames(numericMeta)==colnames(cleanDat) #sanity check -- are sample names in same order?
# 
# numericMeta2<-numericMeta[!is.na(numericMeta$age_at_initial_pathologic_diagnosis),] 
# numericMeta = numericMeta2
# cleanDat<-cleanDat[,na.omit(match(rownames(numericMeta),colnames(cleanDat)))] 
# rownames(numericMeta)==colnames(cleanDat)
# 
# numericMeta2<-numericMeta[!is.na(numericMeta$pathologic_t.1),] 
# numericMeta = numericMeta2
# cleanDat<-cleanDat[,na.omit(match(rownames(numericMeta),colnames(cleanDat)))] 
# rownames(numericMeta)==colnames(cleanDat)

dim(cleanDat) #13470, 463
dim(numericMeta) #463, 158

condition <- as.factor(numericMeta$Group)
condition2 <- as.factor(numericMeta$pathologic_t.1)

age=as.numeric(numericMeta$age_at_initial_pathologic_diagnosis)
sex=as.factor(numericMeta$gender.2) # 0 becomes 1 & 1 becomes 2
#PMI=as.numeric(numericMeta$PMI)

regvars <- as.data.frame(cbind(condition, condition2, sex, age)) #,PMI))
#EBD: note the order of columns for sex and age have been swapped from what I do, this is fine, and gives preference to first modelling sex, then age

## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
rownames(normExpr.reg) <- rownames(cleanDat)
colnames(normExpr.reg) <- colnames(cleanDat)
coefmat <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(regvars)+2) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

#another RNG seed set for reproducibility
set.seed(8675309); #R's random number generator, seed is a number

if (parallelThreads > 1) {
  
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    set.seed(8675309)
    cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(cleanDat), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(cleanDat)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
    coefmat <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% {
      set.seed(8675309)
      options(stringsAsFactors=FALSE)
      library(boot)
      thisexp <- as.numeric(cleanDat[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~ condition +condition2+sex+age)#+PMI)  ## run 1000 resamplings
      ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      bs.stats
      #cat('[bootstrap] Done for Protein ',i,'\n') #will not be visible
    }
    #    normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
    #EBD: you had: coefmat[i,5]*regvars[,"age"]) -- but coefmat only had 4 columns, one for the main term, and one for each of your modelled covariates, plus a 4th column for residuals
    #EBD: I changed the column position for sex back (since I put condition back), and corrected the position for age
    normExpr.reg <-foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { (cleanDat[i,]- coefmat[i,4]*regvars[,"sex"] - coefmat[i,5]*regvars[,"age"]) }
    #EBD: condition is column 2 in coefmat
    #EBD: sex is column 4 in coefmat
    #EBD: age is column 5 in coefmat
    ### All regression code below deos not apply to conditionS set when parallelThreads=15 & boot=TRUE ###
  } else { #linear model regression; faster but incomplete regression of Age, Sex, PMI effects, SO NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~condition +condition2+sex+age,data=regvars) #+PMI
      #EBD: I swapped age and sex here as well in case you ever use lm regression (BOOT=FALSE)
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      #cat('Done for Protein ',i,'\n')
    }
  } #end parallel option -- Average run time estimate printed in console in minutes is calculated based on benchmark of a 2+ GHz intel Xeon with 30 threads 
} else {
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(cleanDat[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~ condition +condition2+sex+age)  #+PMI) ## run 1000 resamplings
      #                       R=numboot, formula=thisexp~ condition +age+sex+PMI) ## run 1000 resamplings
      ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      coefmat[i,] <- bs.stats
      normExpr.reg[i,] <- thisexp - bs.stats[4]*regvars[,"sex"]- bs.stats[5]*regvars[,"age"]  #- bs.stats[5]*regvars[,"PMI"]
      #EBD: again, I swapped age and sex so you would get the same result if using BOOT=TRUE and parallelThreads=1
      cat('[bootstrap] Done for Protein ',i,'\n')
    }
  } else { #linear model regression; faster but NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~condition +condition2+sex+age,data=regvars)  #+PMI
      #EBD: swapped age and sex for lm with single thread
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      cat('Done for Protein ',i,'\n')
    }
  } #end nonparallel option
}

write.csv(cleanDat,file=paste0(outputtabs,"cleanDat.NoTAMPOR.output[log2(FPKMoverCentralTendency)]_OutliersRemoved.csv"))
write.table(cleanDat,file=paste0(outputtabs,"/cleanDat.NoTAMPOR.Regressed_",projectFilesOutputTag,"_",nrow(cleanDat),"x",ncol(cleanDat),".txt"),sep="\t") #check and apply changes to static tail of filename if necessary

##Overwrite cleanDat with regressed data
##(DO NOT RERUN OUT OF CONTEXT)
############################################
cleanDat<-normExpr.reg
rownames(cleanDat)<-rownames(cleanDat.unreg)
############################################

save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))

#########################################################################################################
#=======================================================================================================#
# StarProtocols: Removal of samples (TCGA IDs of miRNA & mRNA platform) that are lableled "do not used" in whitelist downloaded from StarProtocols webpage: https://star-protocols.cell.com/protocols/600#key-resources-table 
#=======================================================================================================#
cleanDat2 <- cleanDat[,which(!colnames(cleanDat)=="TCGA.05.5420")]
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.38.4629")]
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.38.4630")] 
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.44.2661")]
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.44.6146")] 
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.44.6147")] 
cleanDat2 <- cleanDat2[,which(!colnames(cleanDat2)=="TCGA.44.6775")]
cleanDat <- cleanDat2
rm(cleanDat2)
numericMeta <- numericMeta[match(colnames(cleanDat),rownames(numericMeta)),]
# Only 4 of the TCGA ID where in the cleanDat document after outlier removal, regeression steps etc. 
dim(cleanDat) #13486, 461
dim(numericMeta) #461, 201

#########################################################################################################
#=======================================================================================================#
                                    #  WGCNA blockwiseModules()  #
#=======================================================================================================#

enableWGCNAThreads()
library("doParallel")
stopCluster(clusterLocal)
parallelThreads=15 #set to # of threads on your computer
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
tic("timing")

#Check power and connectivity
powers <- seq(8,12,by=0.25)
sft <- pickSoftThreshold(t(cleanDat), powerVector=powers, blockSize=nrow(cleanDat)+1000, verbose=3, corFnc="bicor", networkType="signed")

#10 10.20    0.804 -2.58          0.982    48.9      40.3    194
#11 10.50    0.818 -2.57          0.982    44.0      35.8    182
#12 10.80    0.835 -2.53          0.984    39.6      31.8    172
#13 11.00    0.844 -2.52          0.983    35.8      28.3    163

#blockSize > total number of rows!
toc() # timing 127.646 sec elapsed, time taken to run the sft thresholding code, depending on RAM size, parallelThreads, sequence of power values (i.e. 5-12) & increments in which the analysis moves in. This could take upwards of 2-3 hrs

# plotting SFT R^2 by power, chose value at the elbow with value R.sq above 0.8
plot(sft[[2]][,1],sft[[2]][,2],xlab="Power (Beta)",ylab="SFT R^2")
tic("sleeping")
print("Running power blockwiseModules...")
power=10.2 #choose power at elbow of SFT R≤ curve approaching asymptote near or ideally above 0.80
enableWGCNAThreads()
getDoParWorkers() #6
getDoParRegistered() #TRUE
getDoParName() # doParallelSNOW"
getDoParVersion() # "1.0.15"

## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=100,
                        mergeCutHeight=0.07,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets
print("...done running power blockwiseModules")
toc() #1801 sec
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))
table(net$colors)["grey"]

#turquoise      M1    turquoise 1810
#blue           M2         blue 1278
#brown          M3        brown 1142
#yellow         M4       yellow  822
#green          M5        green  793
#red            M6          red  651
#black          M7        black  629
#pink           M8         pink  518
#magenta        M9      magenta  410
#purple        M10       purple  335
#greenyellow   M11  greenyellow  281
#tan           M12          tan  220
#salmon        M13       salmon  220
#cyan          M14         cyan  208
#midnightblue  M15 midnightblue  198
#lightcyan     M16    lightcyan  197
#grey60        M17       grey60  192
#lightgreen    M18   lightgreen  136
# grey 3446

save.image(paste0(rdata,"NoTAMPOR+Liu+SP-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))

#########################################################################################################
#=======================================================================================================#
            # ANOVA grouping/GNP and kME table, see code of "T1-T4" for ANOVA of T stage  #
            # Adding/Merging Liu Endpoints into NumericaMeta to Re-run GlobalNetworkPlots #
            # see T-Staging ANOVA to adapt GNPs code for other grouping columns #
#======================================================================================================#
pre.Liu.numericMeta = numericMeta # dim(pre.Liu.numericMeta) = 465, 158
Liu.numericMeta = numericMeta
Liu.data <- read.csv(paste0(datadir,"Liu_Data.csv"),header=TRUE,row.names=1)
rownames(Liu.data)==rownames(Liu.numericMeta)
dim(Liu.data) #522, 43
dim(Liu.numericMeta) #465, 158
Liu.data <- Liu.data[match(rownames(Liu.numericMeta),rownames(Liu.data)),]
rownames(Liu.data)==rownames(Liu.numericMeta)
Liu.numericMeta2 = merge(Liu.numericMeta,Liu.data,by = 0,all = TRUE)
rownames(Liu.numericMeta2) <- Liu.numericMeta2[,1]
Liu.numericMeta2$Row.names <- NULL
Liu.numericMeta2 <- Liu.numericMeta2[match(colnames(cleanDat),rownames(Liu.numericMeta2)),] # must order numericMeta to cleanDat again
numericMeta = Liu.numericMeta2

#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)
Grouping<-rep(NA,nrow(numericMeta)) #numericMeta$Group  #typically there is a column "Group" loaded as a column in the traits.csv file
Grouping = numericMeta$Group #numericMeta$Group <- Grouping
Grouping <- paste0("Stg",Grouping)
Grouping[Grouping=="StgNA"] <- NA

## Output GlobalNetworkPlots and kMEtable
FileBaseName=paste0(projectFilesOutputTag,"_pwr",power,"_PAMstageTRUE")

pdf(file=paste0(outputfigs,"2Courtney-Staging+Survival-Darker-SP-Limited-NoTAMPOR-GlobalNetworkPlots-",FileBaseName,".pdf"),width=16,height=12)
#pdf(file=paste0("1.GlobalNetworkPlots-",FileBaseName,".pdf"),width=16,height=12)

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)


#=======================================================================================================#
    # numericIndices can be edited depending on how many traits you will have & want to display in GNPs              (numericIndices2 has been split to fit on pages, then GNP will be combined)
#=======================================================================================================#
# numericIndices2 <- list(1, 2, 21, 35, 37:40, 42:49, 51:57, 66:68, 73:77, 104, 156, 163, 165, 174, 177, 180, 181, 184, 187, 188, 192, 196) #limited GNPs + lastfudate & days2death(like Courtney's paper)
# numericIndices2 <- list(21, 35, 37:40, 42:49, 51:57, 66:68, 73:77, 104, 156, 163, 165, 174, 177, 180, 181, 184, 187, 188, 192, 196) #limited GNPs
# numericIndices2 <- list(21, 35, 37:40, 42:49, 51:57, 66:68, 73:77, 79, 104, 156, 160:201) #limited (full_Liu) GNPs
# numericIndices2 <- list(1:3, 5:11, 13:15, 17, 19, 21:22, 24:26, 28:31, 33, 35, 37:40, 42:49, 51:58, 60, 62, 64, 66:68, 70, 72:79, 81:86, 88, 95, 97, 99, 104:105, 107:111, 156,158, 160:167)
# numericIndices2 <- list(112:155)
# numericIndices2 <- list(1:3, 5:11, 13:15, 17, 19, 21:22, 24:26, 28:31, 33, 35, 37:40, 42:49, 51:58, 60, 62, 64, 66:68, 70, 72:79, 88, 99, 104, 105, 108:111) # supplemental pt1
# numericIndices2 <- list(155:158, 160:201) # supplemental pt2
# numericIndices2 <- list(1, 2, 156, 163, 165, 174, 177, 178, 180, 181, 183:185, 187, 188, 192, 196) # SurvivalBlock
#  numericIndices2 <- list(1, 2, 156, 163, 165, 192, 196, 174, 178, 181, 177, 181, 180, 185, 183, 184, 188)  # Ordered SurvivalBlock
# numericIndices2 <- list(1, 163, 164, 165, 173, 174, 176, 177, 183, 184) # Ordered SurvivalBlock, removed DFI, DSS, overall survival year & days to last followup & CRs
# numericIndices2 <- list(35, 37, 42, 47, 48, 51, 56, 38:40, 43:46, 49, 52:55, 57) # StageBlock ordered, used sorted unlist() line below
# numericIndices2 <- list(1, 2, 163, 164, 173, 174, 176, 177, 183, 184, 35, 37, 42, 47, 48, 51, 56, 38:40, 43:46, 49, 52:55, 57) # survival & stage GNP together
numericIndices2 <- list(1, 163, 160, 51, 56, 52:54, 57, 55) # fig 1
# numericIndices <- sort(unlist(numericIndices2))
 numericIndices <- unlist(numericIndices2) # for ordered list
 # numericIndices <- unlist(numericIndices2) # remove sort to order anyway you like 
# numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
#Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))

## Plot eigengene dendrogram/heatmap - using bicor
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))


######################
## Find differences between Groups (as defined in Traits input file); Finalize Grouping of Samples for ANOVA

# This gets ANOVA (Kruskal-Wallis) nonparametric p-values for groupwise comparison of interest.
# look at numericMeta (traits data) and choose traits to use for linear model-determination of p value
head(numericMeta)
# Change below line to point to a factored trait, which will define groups for ANOVA
regvars <- data.frame(as.factor( numericMeta$Group ), as.numeric(numericMeta$age_at_initial_pathologic_diagnosis), as.numeric(numericMeta$gender.2))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) #sex and age effects are removed by the linear model

pvec <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec) <- colnames(MEs)


######################
## Get sigend kME values
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
cexy <- if(nModules>75) { 0.8 } else { 1 }
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
               xLabels = colnames(numericMeta)[numericIndices],
               yLabels = paste0("ME",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y= cexy,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)

######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
xlabAngle <- if(nModules>75) { 90 } else { 45 }

par(mar = c(5, 15, 3, 3)) #par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = cexy,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-0.5,0.5),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)

## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
# This is where we will first use the Grouping vector of string group descriptions we set above.
toplot <- MEs

MEs.copy<-MEs
MEs.copy[MEs.copy>= 0.15]<-0.15
MEs.copy[MEs.copy<= -0.15]<-0.15
toplot <- MEs.copy

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)

pvec <- pvec[match(names(pvec),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  K-W p=",signif(pvec,2),sep="") # you can rename your heatmap rows to include not only the module rank with color, but also the nonparametric P-value.

# add any traits of interest you want to be in the legend
Gender=as.numeric(numericMeta$gender.2)
Gender[Gender==1]<-"Female" #confirm this is accurate for your data!
Gender[Gender==2]<-"Male"
metdat=data.frame(Group=Grouping,Age=as.numeric(numericMeta$age_at_initial_pathologic_diagnosis), Gender=Gender)

# set colors for the traits in the legend
heatmapLegendColors=list('Group'=c("purple","dodgerblue","goldenrod","red"), #alphabetical order Stg0-Stg4
                         'Age'=c("white","darkgreen"), #young to old
                         'Gender'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))
#Possible Supervised clustering: subset toplot for rows with modules of interest, look into modules & modules2 for possible edits, make edits in clustering alogirthms, only cluster cases using eigenegenes of interest
library(NMF)
par(mfrow=c(1,1))
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=TRUE) ## 

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=NA, Colv=TRUE) ## Cluster columns

# Only use the next 3 aheatmap blocks if you make edits to clustering or set toplot2 (order cases in order of t-stg or pathological stage)
# aheatmap(x=na.omit(toplot2), ## Numeric Matrix
#          main="Plot of Eigengene-Trait Relationships-SAMPLES in T-Stage Order, no clustering",
#          annCol=metdat2,
#          #annRow=data.frame(Modules=colnames(MEs)),
#          annRow =data.frame(Modules2),
#          annColors=heatmapLegendColors2,
#          border=list(matrix = TRUE),
#          scale="row",
#          distfun="correlation",hclustfun="average", ## Clustering options
#          cexRow=0.8, ## Character sizes
#          cexCol=0.8,
#          col=blueWhiteRed(100), ## Color map scheme
#          treeheight=80,
#          Rowv=NA, Colv=NA) ## NA: Do not cluster columns - keep given order
# 
# aheatmap(x=na.omit(toplot2), ## Numeric Matrix
#          main="Plot of Eigengene-Trait Relationships-SAMPLES in T-Stage Order, Column Clustered",
#          annCol=metdat2,
#          #annRow=data.frame(Modules=colnames(MEs)),
#          annRow =data.frame(Modules2),
#          annColors=heatmapLegendColors2,
#          border=list(matrix = TRUE),
#          scale="row",
#          distfun="correlation",hclustfun="average", ## Clustering options
#          cexRow=0.8, ## Character sizes
#          cexCol=0.8,
#          col=blueWhiteRed(100), ## Color map scheme
#          treeheight=80,
#          Rowv=NA, Colv=TRUE) ## NA: Do not cluster columns - keep given order
# 
# aheatmap(x=na.omit(toplot2), ## Numeric Matrix
#          main="Plot of Eigengene-Trait Relationships-SAMPLES in T-Stage Order - Edits",
#          annCol=metdat2,
#          #annRow=data.frame(Modules=colnames(MEs)),
#          annRow =data.frame(Modules2),
#          annColors=heatmapLegendColors2,
#          border=list(matrix = TRUE),
#          scale="row",
#          #distfun="correlation",hclustfun="average", ## Clustering options
#          cexRow=0.8, ## Character sizes
#          cexCol=0.8,
#          col=blueWhiteRed(100), ## Color map scheme
#          treeheight=80,
#          Rowv=NA, Colv=TRUE) ## NA: Do not cluster columns - keep given order

########## see T_Staging ANOVA & "" for edits to heatmap, no ordering, order by T-stage etc.  ###########

## Change the below code in the for loop using the following session output

#These are your numerically coded traits:
colnames(numericMeta)[numericIndices] #choose traits for correlation scatterplots (verboseScatterplot functions below)

#These are your ANOVA sample groups and the number of samples in each
table(Grouping) #alphabetically ordered, you choose the order of groups in the boxplot function by typing them in

## Make changes after checking output on console for the above 2 lines
par(mfrow=c(4,5)) #rows,columns on each page (try to choose to keep all plots for each module on one page or row(s), without squeezing too much in)
par(mar=c(5,6,4,2))

for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(Grouping,c("Stg1","Stg2","Stg3","Stg4")),col=colnames(MEs)[i],ylab="Eigengene Value",main=rownames(toplot)[i],xlab=NULL,las=2)  #you choose the order of groups for boxplots
  verboseScatterplot(x=numericMeta[,"age_at_initial_pathologic_diagnosis"],y=toplot[i,],xlab="# Age",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_T1"],y=toplot[i,],xlab="# TNM (tumor staging T1)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #verboseScatterplot(x=numericMeta[,"years_smoked"],y=toplot[i,],xlab="Years Smoked",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_T2"],y=toplot[i,],xlab="TNM (tumor staging T2)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_i"],y=toplot[i,],xlab="Stage_i",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_ii"],y=toplot[i,],xlab="Stage_ii",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_iii"],y=toplot[i,],xlab="Stage_iii",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_iv"],y=toplot[i,],xlab="Stage_iv",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Stage_iii.5"],y=toplot[i,],xlab="Stage_iii.5",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"pathologic_stage.1"],y=toplot[i,],xlab="Pathologic_Stage",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"pathologic_stage_3.5"],y=toplot[i,],xlab="Pathologic_Stage_3.5",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)  
  verboseScatterplot(x=numericMeta[,"vital_status"],y=toplot[i,],xlab="Alive/Dead",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"pathologic_t_breakdown"],y=toplot[i,],xlab="pathologic_t_breakdown",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"pathologic_t_breakdown_3.5"],y=toplot[i,],xlab="pathologic_t_breakdown_3.5",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"pathologic_t.1"],y=toplot[i,],xlab="pathologic_t.1",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"cigarettes_per_day"],y=toplot[i,],xlab="cigarettes_per_day",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Black"],y=toplot[i,],xlab="Black",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #  verboseScatterplot(x=numericMeta[,"R1"],y=toplot[i,],xlab="R1",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"White"],y=toplot[i,],xlab="White",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #  verboseScatterplot(x=numericMeta[,"R0"],y=toplot[i,],xlab="R0",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #  verboseScatterplot(x=numericMeta[,"R2"],y=toplot[i,],xlab="R2",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"Overall_Survival_Year"],y=toplot[i,],xlab="Overall_Survival",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  verboseScatterplot(x=numericMeta[,"SH.GOE.15.yrs"],y=toplot[i,],xlab="R0",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #verboseScatterplot(x=numericMeta[,"days_to_new_tumor_event_after_initial_treatment"],y=toplot[i,],xlab="Days to New Tumor Event",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21)
  #  frame()
}

dev.off()
save.image(paste0(rdata,"NoTAMPOR+Liu+SP-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#=======================================================================================================#
                            # Providing ENSG-IDs with gene names #
#=======================================================================================================#
#Does no harm to make sure 0's are out, columns matched to traits (actually we can use numericMeta, or set numericMeta to traits on line 204)
#dat.original[dat.original==0]<-NA # makes all 0's to NA's
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]
rownames(numericMeta)==colnames(cleanDat)
source('/Users/coyoung/Projects/WGCNA/Code/Sourced/Biomart.r')
write.csv(cleanDat,file=paste0(outputtabs,"cleanDat.Regressed_Abundance.symbolsBiomaRt.csv"))


#======================================================================================================#
                          # Writing Module Membership/kME table #
#======================================================================================================#
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
#write.table(kMEtable,file=paste0("2.ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
#write.table(kMEtable,file=paste0(outputtabs,"ModuleAssignments-NoTAMPOR-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
write.csv(kMEtable,file=paste0(outputtabs,"ModuleAssignments-NoTAMPOR-",FileBaseName,".csv"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel.

#======================================================================================================#
                    # NetworkScreening for Different clincial traits #
#======================================================================================================#
# FindReplace function to change '-- in cigs per day column into NA so I can run NetworkScreening function 
# Link: https://www.rdocumentation.org/packages/DataCombine/versions/0.2.21/topics/FindReplace 
# BiocManager::install("DataCombine")
dim(numericMeta) #461, 201
library(DataCombine)
# Create replacements data frame
#Replaces <- data.frame(from = c("'--"), to = c("NA"))
#newNum2 <- FindReplace(data = numericMeta, Var = "cigarettes_per_day", replaceData = Replaces, from = "from", to = "to", vector = FALSE)
#newNum2$cigarettes_per_day <- as.numeric(newNum2$cigarettes_per_day)
#numericMeta=newNum2

# networkScreening: gives you 2 sets of columns (weighted correlation [normally disregard] & unweighted correlation, bicor, q-score, p-value & z-score). Not to different from bicorAndPvale fxn. Bascially, sorts through cleanDat & finds what is most correlated to a clincal trait (this can be a gene), then you can sort by correlation to see what is most correlated.
# Dammer: Networkscreening works best w/continous traits
#add pack-years #########

myDataframe<-numericMeta # setting numericMeta to new df
myColnameVector<-c("Overall_Survival_Year", "pathologic_n.1", "Stage_N0", "Stage_N2", "pathologic_t.1", "Stage_T1", "Stage_T2", "pathologic_t_breakdown_3.5", "pathologic_stage.1", "Stage_i", "Stage_ii", "Stage_iii",  "pathologic_stage_3.5",  "Stage_iii.5", "Stage_iv",  "race.1", "Black", "White",  "histological_type.1", "vital_status", "pathologic_t_breakdown", "cigarettes_per_day", "number_pack_years_smoked", "karnofsky_performance_score", "days_to_death", "OS", "OS.time", "DFI", "DFI.time") # column names of interest, all columns must be typeof() = "integer"
# for (col in myColnameVector) {print(typeof(numericMeta[,col]))  } # for loop to give typeof for each column name

netScreen.lists <- list() # saving each dataframe output by networkScreening in a list, saves output into list so your output will not only be from the last loop.
for (iterator in myColnameVector) { #for every iterator in the myColnameVector group
  #cat(myDataframe[,iterator]) 
  #print(myDataframe[,iterator]) # lines just to make sure all values are caputred 
  netScreen = paste0(outputtabs,"netScreens")
  dir.create(file.path(netScreen))
  temp.netScreen.traits <- networkScreening(myDataframe[,iterator], MEs, t(cleanDat), #can change y=clinical trait given as numeric vector (1 valure per sample)
                                            corFnc = "cor", corOptions = "use = 'p'",
                                            oddPower = 3,
                                            blockSize=nrow(cleanDat)+1000,
                                            minimumSampleSize = 4,
                                            addMEy = TRUE, removeDiag = FALSE,
                                            weightESy = 0.5, getQValues = TRUE)
  write.csv(temp.netScreen.traits,file=paste0(netScreen,"NetworkScreening-NoTAMPOR-",iterator,".csv"))
  netScreen.lists[[iterator]]<- temp.netScreen.traits # #list stores all iterations' output with separate elements named the same as your loop iterator variable.
}

#======================================================================================================#
                              # GO-ELite for networkScreening() output #
#======================================================================================================#
# Make susre all netScreen outputs are in: "/Users/coyoung/Projects/WGCNA/Tables/netScreens"
# write.csv as mentioned in code (output from netScreenforGO.Final.R) and open in excel as detailed in the code to remove all formatting changes 
source('/Users/coyoung/Projects/WGCNA/Code/Sourced/netScreenforGO.Final.R') #output result
  # add in background row, take top genes with z.standard of less than -3 (bottom) and greater than 3 for each trait (top), take the postive z.standard values (top) & use a column of sequential numbers in another sheet to reverse the order. Now you have the top/bottom (3/-3) signifcant valuesby z.standard. Bottom/top= negative/postive z.standard
# Open a new Excel sheet, select the Data tab, then click 'From Text' in the Get External Data group. Browse to the CSV file and select 'Import'....
# use file to build GO-elite input. Using all signficant values that are postively correlated as top & all values that are negatively correlated as bottom 
# save final file as "inputForGoElite.csv" in '/Users/coyoung/'
# run edited GO-Elite code: 
source('/Users/coyoung/Projects/WGCNA/Code/Sourced/Scale-Edit-GoElite.R') #output result
# source('/Users/coyoung/Projects/WGCNA/Code/Sourced/FixedAxis-Scale-Edit-GoElite.R') #output result, has fixed GO-ELite axes
# just open the file and use if you are using db other than the orginial and C3 dbs

#======================================================================================================#
  # ANOVA / DiffEx -- list of ANOVAout dataframes for different subgroup comparisons, if necessary #
#======================================================================================================#
groupsToTest=c("ALL") #unique(numericMeta$Group)

ANOVAoutList<-list()
#for (caseSubset in groupsToTest) {
#  cleanDat.sub<-cleanDat[,match(rownames(numericMeta)[numericMeta$Group==caseSubset],colnames(cleanDat))]
#  Grouping.sub <- numericMeta$Region.Treatment[numericMeta$Group==caseSubset] #character strings
caseSubset="ALL"
cleanDat.sub<-cleanDat
Grouping.sub<-Grouping

outfile = paste0("ANOVA_diffEx-",caseSubset,"-",FileBaseName)
data = as.data.frame(cbind(colnames(cleanDat.sub), Grouping.sub, t(cleanDat.sub)))
colnames(data)[1:2]<-c("CODE","SampleType") # just changes the 1st 2 rows to CODE & SampleType
#test run gets column headers for output
i=3
aov<-aov(data[,i]~SampleType, data=data) #SampleType is a factor 
#data2 <- as.numeric(as.character(colnames(data[3:31393])))
#dffinalSep$PC1 <- as.numeric(as.character(dffinalSep$PC1))
#aov<-aov(data2[,i]~SampleType, data2=data2)

anovaresult<-anova(aov)
tuk <- TukeyHSD(aov)
tukresult1<-data.frame(tuk$SampleType) #ASSUMES NO PAIRWISE COMPARISONS ARE MISSING FOR FIRST cleanDat.sub protein (ROWS OF THIS data frame)--this is the template for comparison columns in ANOVAoutList[[caseSubset]]
j=length(rownames(tukresult1))
comparisonList<-rownames(tukresult1)

line = c(paste("Protein", "F-Value", "Pr(>F)", sep=","))
for (a in 1:length(comparisonList)) {
  line=c(paste(line,comparisonList[a],sep=","))
}
for (a in 1:length(comparisonList)) {
  line=c(paste(line,paste0("diff ",comparisonList[a]),sep=","))
}

## Fast code with apply
dataclipped <- data[, 3:ncol(data)] # as.matrix(as.numeric(data[,3:ncol(data)]))
SampleType <- data$SampleType
ANOVAoutList[[caseSubset]] <- apply(dataclipped, 2, function(x) {
  if(length(unique(SampleType[which(!is.na(x))]))<2) { #<2, handles "only one contrast level..." 
    x=data.frame(x=rep(1,length(x)), SampleType = SampleType)  #handles too many missing values
  } else {
    x <- data.frame(x = as.double(x), SampleType = SampleType) #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
  }
  aov <- aov(x ~ SampleType, data = x)
  anovaresult <- anova(aov)
  tuk <- TukeyHSD(aov)
  tukresult <- data.frame(tuk$SampleType)
  if (length(rownames(tukresult)) == length(rownames(tukresult1))) {
    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]))
  } else {
    tukresult <- tukresult[match(rownames(tukresult1), rownames(tukresult)), ]
    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]))
  }
})
ANOVAoutList[[caseSubset]] <- t(ANOVAoutList[[caseSubset]])
ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
ANOVAcols <- ANOVAcols[2:length(ANOVAcols)]
if (length(unique(SampleType))==2) { #for single pairwise comparison, essentially Pr(>F) from ANOVA is equivalent to T-test result (except 1 vs 2 non-NA measurements are allowed)
  ANOVAoutList[[caseSubset]][,3] <- p.adjust(ANOVAoutList[[caseSubset]][,2],method="BH") #get BH FDR for ANOVA/T-test p values
  ANOVAoutList[[caseSubset]][,2:3]<-ANOVAoutList[[caseSubset]][,c(3,2)]
  ANOVAcols[2] <- "FDR (BH)"
}
colnames(ANOVAoutList[[caseSubset]]) <- ANOVAcols
ANOVAoutList[[caseSubset]] <- as.data.frame(ANOVAoutList[[caseSubset]])

## Add network module colors
ANOVAoutList[[caseSubset]]$NETcolors <- net$colors

#*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
write.csv(ANOVAoutList[[caseSubset]],file=paste0(outputtabs,"/",outfile,".csv"))
#*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#  cat(paste0("Finished ANOVA for subset ",caseSubset,". [",length(groupsToTest)," total groups to test]\n"))
#}
ANOVAout<-ANOVAoutList[["ALL"]]
# Saved image on 12/4/21 after loading NoTAMPOR+Liu+SP... image & running ANOVA code
save.image(paste0(rdata,"NoTAMPOR+Liu+SP-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))

#======================================================================================================#
                          # Generation of orignal Dammer Volcano plots #
#======================================================================================================#
flip <- NULL
library(ggplot2)

cutoff <- log2(1) # NO change minimum to be counted in the volcano bookends; log2(1.25) for 25% FC min.
cutoff

sigCutoff <- 0.05 # p value for Volcano significant hit counting; dashed line at -log10(sigCutoff)

#** for (caseSubset in groupsToTest) { #*+*+*+*+*
#** baseNameVolcanoes <- paste0(caseSubset,".",FileBaseName)
#** ANOVAout<-ANOVAoutList[[caseSubset]]
#** instead of above loop for different ANOVA groupings, use the non-list ANOVAout already set for the single comparison of 2 groups (collapsed to T-test)
#** using below code instead #**
baseNameVolcanoes <- paste0(FileBaseName) #**

# We want to note the positions of the important genes in each volcano
BIGspots <- c()
BIGspots <- unique(c('CXCL13', 'CXCR5', 'ARRB1', 'ARRB2', 'GRK1', 'GRK2', 'GRK3', 'GRK4', 'GRK6', 'GRK7', 'GNAQ', 'GNA11', 'GNA12', 'GNA13', 'GNA14', 'GNA15', 'GNAI1', 'GNAI2', 'GNAI3', 'GNB3', 'GNB3s',  'GNG8", GNG9', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIK3R4', 'PIK3R5', 'PIK3R6', 'PIK3CA', 'PIK3CB', 'PIK3CG', 'PIK3CD',  'DOCK2', 'RAC1', 'RAC2', 'RAC3', 'RHOA', 'RHOB', 'RHOC', 'RHOG', 'RHOH', 'RHOBTB1', 'RHOBTB2', 'RHOBTB3', 'JUN', 'VAV1', 'CTNNB1', 'BRCA1', 'CREB1', 'ELK1',  'HDAC8', 'RB1', 'Actin', 'Calpain', 'CAV1' ,'CFL1', 'Cofilin', 'CTTN', 'Dynamin', 'Erm', 'EZR', 'FAK', 'Src', 'G3BP1', 'KRT18', 'MAP2K1','MAP2K1', 'NF2', 'NTRK2', 'Pak', 'phosphatase', 'PTEN', 'PTK2', 'PXN', 'Rac', 'Ras', 'Rock', 'SRC', 'Talin', 'VASP')) #unique( c('TARDBP','C9ORF72','SCA36','SCA3','SQSTM1','SFPQ','MSN','PLEC','HEPACAM') )
colnames(ANOVAout) # choose a column for -log10 p (change column number to point to testIndex below)
n <- nrow(ANOVAout)
#***********************************
testIndexMasterList <- c(3,4,5,6,7,8)
#flip <- as.vector(c(3,4,5,6,7,8)) # edit this vector, to flip sign of this/these comparison(s); e.g. put WT as denominator for log2(ratio) if you see it first in colnames(ANOVAout)
useNETcolors<-TRUE
splitColors<-TRUE
#***************************** (stop here, edit above)
dexComps <- list()
iter <- length(testIndexMasterList) + 1
comparisonIDs <- data.frame(dfVariable = rep(NA, length(testIndexMasterList)), Comparison = rep(NA, length(testIndexMasterList)))
numberOfNonComparisonColumns=length(colnames(ANOVAout)) - length(which(grepl("diff ",colnames(ANOVAout))))*2
numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons
for (i in testIndexMasterList) {
  iter <- iter - 1
  # dexRows<-which(ANOVAout[,i]<sigCutoff) #choose rows where the DEX p<sigCutoff
  comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("-", ".", colnames(ANOVAout)[i])), paste0(as.character(gsub("-", " vs ", colnames(ANOVAout)[i])))))
  dexComps[[comparisonIDs[iter, 1]]] <- ANOVAout
  if (!is.na(match(i, flip))) {
    dexComps[[comparisonIDs[iter, 1]]][, i + numComp] <- -1 * as.numeric(dexComps[[comparisonIDs[iter, 1]]][, i + numComp])
    comparisonIDs[iter, 2] <- gsub("(*.*) vs (*.*)", "\\2 vs \\1", comparisonIDs[iter, 2]) # flip label "vs" in comParisonIDs$Comparison[iter]
  }
}
comparisonIDs # list element names and Logical comparisons for those retrievable Dex measurements in the list elements
ls(dexComps) # list elements are dataframes with the DEX entries for that comparison

## volcano plots with (module or other) colors, SPLITTABLE on COLOR

pointColorsVectorListForPlots <- list()
volcListModColorsWeb <- list()
volcListModColors <- list()
dfListModColors <- list()

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  df <- eval(parse(text = "dexComps[[comparisonIDs$dfVariable[iter]]]"))
  cat(paste0("Processing ANOVA column ", testIndex, " (", comparisonIDs$Comparison[iter], ") for volcano...\n"))
  # correct 0 Tukey pValues to ANOVA p (in column 2); it's better than taking -log10 of 0 in the next step
  df[which(df[, testIndex] == 0), testIndex] <- as.numeric(df[which(df[, testIndex] == 0), 2])
  
  ## Check if ANOVA pVal is Significant and above FC cutoff defined above. Thresholds are used to set volcano point colors
  df[, testIndex][is.na(df[, testIndex])] <- 1 # p=0.9999 instead of NA
  df[, testIndex + numComp][is.na(df[, testIndex + numComp])] <- 0 # log2(difference)=0 instead of NA
  
  df$negLogP <- -log10(as.numeric(df[, testIndex]))
  
  df$threshold1 <- as.numeric(rep(0, n))
  ## Any COMPARISON SIGNIFICANT (uses ANOVA p in column 2 of df instead of Tukey p): # for (i in 1:n) { if (abs(as.numeric(df[i,testIndex+numComp]))<cutoff | df[i,2]>sigCutoff ) {df$threshold1[i]=3} else { if (df[i,testIndex+numComp]<cutoff) {df$threshold1[i]=2} else {df$threshold1[i]=1}} }
  for (i in 1:n) {
    if (abs(as.numeric(df[i, testIndex + numComp])) < cutoff | as.numeric(df[i, testIndex]) > sigCutoff) {
      df$threshold1[i] <- 3
    } else {
      if (as.numeric(df[i, testIndex + numComp]) < cutoff) {
        df$threshold1[i] <- 2
      } else {
        df$threshold1[i] <- 1
      }
    }
  }
  df$threshold1 <- as.factor(df$threshold1)
  
  df$Symbol <- rownames(df) #for symbol only:  do.call("rbind", strsplit(as.character(rownames(df)), "[|]"))[, 1]
  
  ## Color Interesting Gene Product Spots DIFFERENTLY as 4th color if doing blue/red/green (no module colors) -- (4=gold1 below)
  if(useNETcolors) { 
    df$color1<-df$NETcolors
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  } else {
    df$color1 <- as.numeric(df$threshold1)
    df$color1[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  }
  df$color1 <- as.factor(df$color1)
  
  if(useNETcolors) {
    df$size1 <- as.numeric(df$threshold2)
  } else {
    df$size1 <- as.numeric(df$color1)
  }
  df$size1[df$size1 < 4] <- 2.5
  df$size1[df$size1 == 4] <- 6
  
  #df$color2: actual spot color as.character() ; df$color1 is a factorable dummy
  if(useNETcolors) {
    df$color2<-as.character(df$color1)
    df$color2[df$threshold2 == 4] <- "gold1" #for BIGspots
  } else {
    df$color2 <- as.numeric(df$color1)
    df$color2[df$color2 == 1] <- "darkred"
    df$color2[df$color2 == 2] <- "darkgreen"
    df$color2[df$color2 == 3] <- "dodgerblue"
    df$color2[df$color2 == 4] <- "gold1" #for BIGspots
  }                              #gold1 is also the 435th, last WGCNA unique color
  
  #df$color3 is outline color, where outlined pch symbols (21) used
  df$color3 <- df$color2
  df$color3[df$color3 == "gold1"] <- "black" #for BIGspots outline
  df$pch <- as.numeric(df$threshold2)
  df$pch[df$pch < 4] <- 16 # unfilled circles (use color2)
  df$pch[df$pch == 4] <- 21 # filled, outlined circles (border uses color3)
  
  #put gold1 back to module color for fill of BIGspots if (useNETcolors)
  if (useNETcolors) { df$color2[df$color2=="gold1"] <- as.character(df$color1[df$color2=="gold1"]) }
  
  df <- df[order(df$size1, decreasing = FALSE), ] # puts larger dots on top (at bottom of df)
  
  
  #splitColors TRUE/FALSE: make one volcano with all colors (FALSE), or make volcanoes for each color (TRUE)
  # SPLIT DATA FRAME FOR VOLCANO PLOT BY COLORS (if multiple eachColorSplit items)
  df.AllColors <- df
  eachColorSplit <- if (splitColors) {
    unique(df.AllColors$NETcolors)
  } else {
    c("allcolors")
  }
  for (eachColor in eachColorSplit) {
    if (splitColors) {
      df.oneColor <- df.AllColors[which(df.AllColors$NETcolors == eachColor), ]
      df.oneColor$color1 <- factor(df.oneColor$NETcolors)
    }
    else {
      df.oneColor <- df.AllColors
    } # end if (splitColors)
    
    names(df.oneColor)[testIndex + numComp] <- "xdata" # x=as.numeric(df.oneColor[,testIndex+numComp])
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".") # colnames(df)[testIndex]
    pointColorsVectorListForPlots[[list_element]] <- data.frame(color1 = factor(as.integer(df.oneColor$color1)), color2 = as.character(df.oneColor$color2), color3 = as.character(df.oneColor$color3), size = as.numeric(df.oneColor$size), pch = as.numeric(df.oneColor$pch)) #*** df.oneColor$netColors
    
    
    volcano1 <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      # scale_colour_manual(values = unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[order(unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[,1]),2] )+ #THIS COLOR(S) IS LOOKED UP ACTIVELY DURING GRAPHICS OUTPUT IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      geom_point(aes(fill = pointColorsVectorListForPlots[[list_element]][, "color3"]), alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = pointColorsVectorListForPlots[[list_element]][, "pch"], color = pointColorsVectorListForPlots[[list_element]][, "color3"]) +
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(as.expression(bquote("Difference, log"[2] ~ .(comparisonIDs$Comparison[iter])))) + # colnames(df.oneColor)[testIndex]
      ylab(as.expression(bquote("-log"[10] ~ "p value"))) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      # geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    # web version doesn't use as.expression! plotly fails with those, so we rebuild the volcano for the web.
    volcanoweb <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      scale_colour_manual(values = unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[order(unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[, 1]), 2]) + # THIS COLOR(S) IS LOOKED UP ACTIVELY BY PLOTLY IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      # scale_y_continuous(breaks = seq(0, 8, by = 1))+
      # geom_point(aes(fill=pointColorsVectorListForPlots[[list_element]][,"color2"]), alpha=0.66, size=pointColorsVectorListForPlots[[list_element]]$size, pch=16, color=pointColorsVectorListForPlots[[list_element]][,"color3"]) + #pch=pointColorsVectorListForPlots[[list_element]][,"pch"] not using variable symbol types (21 for outlined circle only)
      geom_point(alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = 16) + # pch=pointColorsVectorListForPlots[[list_element]][,"pch"] just uses the higher pch code in the web render.
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(paste0("Difference, log2 ", comparisonIDs$Comparison[iter])) +
      ylab(paste0("-log10 p value")) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      # geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    volcListModColors[[list_element]] <- volcano1
    volcListModColorsWeb[[list_element]] <- volcanoweb
    
    print(volcano1) # prints to active output (separate page)
    rm(volcano1)
    rm(volcanoweb)
    dfListModColors[[list_element]] <- df.oneColor
  } # closes for(eachColor...
} # closes for(testIndex...

if(splitColors) { dir.create(file.path(outputfigs, "/SplitVolcano/")) }

# Print to PDFs, one per color (per comparison, if multiple)
iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    if(splitColors) {
      file <- paste0(outputfigs,"/SplitVolcano/","Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".pdf")
    } else {
      file <- paste0(outputfigs,"/","Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".pdf")
    }
    pdf(file = file, colormodel="srgb", height = 8, width = 8)
    par(mfrow = c(1, 1))
    par(mar = c(6, 8.5, 3, 3))
    print(volcListModColorsWeb[[list_element]])
    #     volcListModColors[[list_element]])    #bigspots colors not handled correctly for PDF output of modColored spots:
    dev.off()
  }
}

## Print to html volcano plots (one per module colors and per comparison, if applicable)
library(plotly)

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    webPlot <- ggplotly(volcListModColorsWeb[[list_element]])
    if(splitColors) {
      tempfilename=paste0(outputfigs,"/SplitVolcano/","HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".html")
    } else {
      tempfilename=paste0(outputfigs,"/","HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".html")
    }
    htmlwidgets::saveWidget(webPlot, tempfilename, selfcontained = TRUE, libdir = "delete.me")
  }
}
#** } #ends for (caseSubset in groupsToTest) { #*+*+*+*+*

# For Volcano plots that spand all modules per comparion per PDF sheet run:
# tHis can also be done by having sourced code of Volcano plots without edits ran again and just setting splitcolors = FALSE and running the code again.
source('/Users/coyoung/Projects/WGCNA/Code/Sourced/Volcano_AllModules.R')
#save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#########################################################################################################
#=======================================================================================================#
      # iGRAPHs (Multiple Toggle Options, e.g. BioGRID interactome overlap) // CONNECTIVITY PLOT #
#=======================================================================================================#
#(re)Calculate MEs (no grey)
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL


#specify filename
baseName=FileBaseName #inputExprMat #PDF filename base
#check these settings
#########################################################################################################
PPIedges=FALSE  #TRUE will be a lot slower...
vertexsize=16   #8 for regular, 16 for large balls
species="human" #current option "mouse" will convert bioGRID to mouse symbols before drawing PPI edges.
#CAIRO=TRUE      #Open/Write PDF using CairoPDF or pdf functions; now both file versions are created
#########################################################################################################

#(re)Establish M# for module colors table
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))

#load interactome [full bioGRID 3.4 (01/25/2016)]
bioGrid <- read.table(file=paste0(iGraphs,"BIOGRID-Homo_sapiens-3.4.133.SIMPLEsymbols_editpad.txt"),sep="\t",header=TRUE)
#------

if (species=="mouse") {
  library(biomaRt)
  human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
  mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
  
  ## If converting bioGrid from one species to another
  genelist.convert<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=unique(c(as.vector(bioGrid[,"From"]),as.vector(bioGrid[,"To"]))), mart=human, attributesL=c("mgi_symbol","external_gene_name"),martL = mouse)
  
  bioGrid.mouse.convert<-cbind(genelist.convert[match(bioGrid[,1],genelist.convert[,1]),3],genelist.convert[match(bioGrid[,2],genelist.convert[,1]),3])
  bioGrid.human<-bioGrid
  bioGrid<-bioGrid.mouse.convert
  dim(bioGrid)
  bioGrid<-as.data.frame(bioGrid[-which(is.na(bioGrid)),])
  colnames(bioGrid)<-colnames(bioGrid.human)
  dim(bioGrid)
}

if (PPIedges==FALSE) { bioGrid<-data.frame(From=c(1),To=c(1)) } #if you want to skip bolding edges for PPIs from BioGrid

#Symbol to add to iGraphs that don't include it, if checking for protein-protein interaction edges (add 4, one for each corner; can edit output later to remove unwanted ones.)
symbols2Add=c() #c("PiB","APP","SGIP1","APOE")

#ADmagmaList<-read.csv(file="c:/Users/Eric Dammer/Documents/wgcnatest/magma/input/AD_Mean_Zstat_genePValues_1234genesFINAL.csv",header=TRUE)
RNAbindingGOlist<- c()

#Get KMEs
# KME.vis=signedKME(cleanDat, kMEdat,corFnc="bicor",outputColumnName = "");
KME.vis=kMEdat[,] # -ncol(kMEdat) #minus grey column if calculating signedKME on the fly
annot=as.data.frame(do.call("rbind",strsplit(as.character(rownames(cleanDat)),"[|]")))
KME.vis$Symbol=annot$V1
geneInfo.vis=as.data.frame(cbind(KME.vis$Symbol,net$colors, KME.vis))  

geneInfo.vis=geneInfo.vis[,-ncol(geneInfo.vis)] # check if last column is Ensembl gene id
colnames(geneInfo.vis)[1]= "Symbol"
colnames(geneInfo.vis)[2]= "Module.Color"
colnames(geneInfo.vis)<-gsub("kMEME","kME",colnames(geneInfo.vis))

library(igraph);
library(RColorBrewer);
library(WGCNA);
library(Cairo);

softPower = power
adjacency = adjacency(cleanDat, power=softPower, type="signed",corFnc="bicor")
TOM = TOMsimilarity(adjacency)

TOM.matrix = as.matrix(TOM);
#Get the top connected genes in the module
uniquemodcolors = gsub("kME","",gsub("kMEME","",colnames(kMEdat[,]))); #-ncol(kMEdat) last column if KMEdat.vis was calculated on the fly... moduleColors
#OR SELECT MODULES INSTEAD OF ALL
#uniquemodcolors = c("orange")

for (CAIRO in c(TRUE,FALSE)) {
  
  if (CAIRO) { CairoPDF(file=paste0(outputfigs,"/iGraph_Modules-",FileBaseName,"-CAIRO-LargeNodes.pdf"),width=16,height=12) } else {
    pdf(paste0(outputfigs,"/iGraph_Modules-",FileBaseName,"-nonCAIRO-LargeNodes.pdf"),height=9,width=10)
  }
  # for (i in 1:length(sigmodcolors))  {
  # mod=sigmodcolors[i];	
  # numgenesingraph = 50;
  # numconnections2keep = 1500;
  for (mod in uniquemodcolors)  {
    #mod="darkslateblue"
    numgenesingraph = 100;
    numconnections2keep = 700;
    cat('module:',mod,'\n');
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="NA",]
    geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="",]
    
    colind = which(colnames(geneInfo.vis)== paste("kME",mod,sep=""));
    rowind = which(geneInfo.vis[,2]==mod);
    cat(' ',length(rowind),'probes in module\n');
    submatrix = geneInfo.vis[rowind,];
    orderind = order(submatrix[,colind],decreasing=TRUE);
    if (length(rowind) < numgenesingraph) {
      numgenesingraph = length(rowind);
      numconnections2keep = numgenesingraph/2 * (numgenesingraph/6 - 1); #added /2 and /6 9/14/2015
    }
    if (length(rowind)<10) { innercircleNum=length(rowind) } else { innercircleNum=10 }
    cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
    submatrix = submatrix[orderind[1:numgenesingraph],];
    #Identify the columns in the TOM that correspond to these hub probes
    matchind = match(submatrix$Symbol,annot$Symbol);
    reducedTOM = TOM.matrix[matchind,matchind];
    
    orderind = order(reducedTOM,decreasing=TRUE);
    connections2keep = orderind[1:numconnections2keep];
    reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] = 1;
    
    g0 <- graph.adjacency(as.matrix(reducedTOM[1:innercircleNum,1:innercircleNum]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMata <- layout.circle(g0)
    
    if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) {   
      g0 <- graph.adjacency(as.matrix(reducedTOM[11:ncol(reducedTOM),11:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatb <- layout.circle(g0)
      
      g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75)
    } else { if (ncol(reducedTOM) > 10) { #****
      
      g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatb <- layout.circle(g0)
      
      g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMatc <- layout.circle(g0)
      g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
      layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc)
    }
      else { layoutMat <- layoutMata } #****
    }
    
    #PLOT DONE BELOW WITH ADDED PHYSICAL INTERACTION EDGES (OR USE THIS SINGLE LINE FOR NO PHYSICAL INT HIGHLIGHTING)
    #  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(rbind(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*16,main=paste(mod,"module"))
    
    #Add symbols2Add
    iter=as.numeric(0)
    
    if (length(symbols2Add)==4) {
      for (symbol2Add in symbols2Add) {
        iter=iter+1
        if (iter==1) { position=cbind(-0.7071068,-0.7071068) } 
        if (iter==2) { position=rbind(position,cbind(0.7071068,-0.7071068)) }
        if (iter==3) { position=rbind(position,cbind(-0.7071068,0.7071068)) } 
        if (iter==4) { position=rbind(position,cbind(0.7071068,0.7071068)) }
        #ADD APP to graph g3 (to be used if APP not already in graph)
        g2 <- vertex(1) #,size=40,color="green", label=symbol2Add
        #layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,cbind(-0.7071068,-0.7071068)))
        g3 <- g1+g2
        V(g3)$name <- c(1:(numgenesingraph+1))
        #WORKS# plot(g3,edge.color="grey",vertex.color=c(rep(mod,nrow(layoutMat)-1),"green"),vertex.label=as.character(c(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=c(submatrix[,colind]^2*16,15),main=paste(mod,"module"))
        
        g1<-g3
        numgenesingraph=numgenesingraph+1
      }
    }
    
    #x***moved out of loop
    if (length(symbols2Add)==4) {  if (ncol(reducedTOM) < 51) { layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,position)) } else { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, rbind(layoutMatc,position*1.33)) }
    } else { if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75) } else { if(ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc) } } #do not add 4 node positions if symbols2Add is empty (length 0)
    }
    symbolList<-c(as.character(submatrix$Symbol),symbols2Add)
    if (length(symbols2Add)==4) { vertexColorVec<-c(rep(mod,numgenesingraph)-4,"green","darkgreen","steelblue","darkslateblue"); vertexSizeMat<-c(submatrix[,colind]^2*vertexsize,rep(15,4)); } else { vertexColorVec<-rep(mod,numgenesingraph); vertexSizeMat<-submatrix[,colind]^2*vertexsize; }
    
    #FIND EDGES OVERLAPPING WITH BIOGRID PAIRS
    listboldEdges<-matrix(ncol=2,nrow=0)
    if(nrow(bioGrid)>0) {
      for(i in 1:numgenesingraph) {
        for(j in i:numgenesingraph) {
          if(!(length(which(bioGrid$From==symbolList[i]&bioGrid$To==symbolList[j]))+length(which(bioGrid$From==symbolList[j]&bioGrid$To==symbolList[i])))==0) { listboldEdges <- rbind(listboldEdges,c(i,j)) }
        } }
      
      ##Remove self-loops
      if (PPIedges) { if (length(nrow(listboldEdges))>0 & nrow(listboldEdges)>0) { 
        for(i in nrow(listboldEdges):1) {
          if (i<=nrow(listboldEdges)) {
            if(listboldEdges[i,1]==listboldEdges[i,2]) {
              listboldEdges <-listboldEdges[-i,]
              if (!length(nrow(listboldEdges))>0) { 
                cat('NO PHYSICAL INTERACTIONS FOUND FOR THIS MODULE.\n')
                break #HANDLE NO BOLD EDGES CASE
              }
            }
          }
        }
      }
      }
      
    }
    
    if (is.vector(listboldEdges)) { listboldEdges<-t(data.frame(listboldEdges)) }
    
    newgraph <- g1 %>%
      #  delete_edges(listboldEdges) %>%
      set_edge_attr("color",value="lightgrey") %>%
      set_edge_attr("width",value=1) %>%
      set_edge_attr("curved",value=0) 
    
    
    if (length(symbols2Add)==4) {
      if (is.vector(listboldEdges)) { 
        cat('ONLY 1 SINGLE INTERACTION FOUND.\n')
        newgraph <- newgraph %>%
          add_edges(c(listboldEdges[1],listboldEdges[2]),color="steelblue",width=2,curved=0)
      } else {
        if (dim(listboldEdges)[1]>0) { 
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[1]))==0)&(listboldEdges[k,1]==(numgenesingraph-3)|listboldEdges[k,2]==(numgenesingraph-3))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#33BB33",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[2]))==0)&(listboldEdges[k,1]==(numgenesingraph-2)|listboldEdges[k,2]==(numgenesingraph-2))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#338833",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[3]))==0)&(listboldEdges[k,1]==(numgenesingraph-1)|listboldEdges[k,2]==(numgenesingraph-1))) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="steelblue",width=3,curved=0)
            }
          }
          for(k in 1:nrow(listboldEdges)) {
            if((length(which(submatrix$Symbol==symbols2Add[4]))==0)&(listboldEdges[k,1]==numgenesingraph|listboldEdges[k,2]==numgenesingraph)) { 
              newgraph <- newgraph %>%
                add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="darkslateblue",width=3,curved=0)
            }
          }
          
        }
      }
    } else {
      if (dim(listboldEdges)[1]>0) { 
        for(k in 1:nrow(listboldEdges)) {
          newgraph <- newgraph %>%
            add_edges(c(listboldEdges[k,1],listboldEdges[k,2]),color="#33BB33",width=3,curved=0)
        }
      }
    }
    
    highlightNodes<-match(RNAbindingGOlist,symbolList) #ADmagmaList$Gene.Symbol,symbolList
    highlightNodes<-sort(na.omit(highlightNodes))
    if(length(highlightNodes)>0) { if(mod=="yellow") {vertexColorVec[highlightNodes]<-"cyan" } else {vertexColorVec[highlightNodes]<-"yellow" } }
    
    plot(newgraph,vertex.color=vertexColorVec,vertex.label=as.character(symbolList),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout=layoutMat,vertex.size=vertexSizeMat,main=paste0(orderedModules[which(orderedModules[,"Color"]==mod),"Mnum"]," ",mod," module"))
  }
  dev.off();
  
} # for (CAIRO in c(TRUE,FALSE)) loop... <repeat>
#save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#########################################################################################################
#=======================================================================================================#
            # One-Step GO-ELITE WITH USER PARAMETERS - by Eric Dammer, Divya Nandakumar #
              #- performs GO-Elite v1.2.5 with Fisher Exact Test for enrichment p<0.05 #
                              # and 5 minimum genes per ontology #
#=======================================================================================================#
# GO-Elite is a python package to perform ontology analysis on gene sets. 
# The python script for GO-Elite can be downloaded from http://www.genmapp.org/go_elite/
# Alternatively, there is also a GUI which can be downlaoded from the same website.
# Custom databases can be downloaded from http://software.broadinstitute.org/gsea/msigdb/
# including TRANSFAC and curated pathways collected from various sources, e.g. the C2 DB.
# GO-Elite requires python 2.7 be installed, and FET with command-line requires
# that a bugfix be applied; copy GO_Elite.py included to GO-Elite program subfolder:
# GOeliteFolder/GO-Elite_v.1.2.5-Py/
# just open the file (paly-GOelit) or copy and paste to a different R script sheet and use if you are using db other than the orginial and C3 dbs
#########################################################################################################

dev.off()
options(stringsAsFactors=FALSE)

#EDIT THESE VARIABLES (USER PARAMETERS)
#########################################################################################################
#fileName <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#fileName <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv" #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#INPUT CSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background for GO-Elite.
#  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.
#filePath <- gsub("//","/",outputfigs)
filePath <- '/Users/coyoung/'
#Folder that will contains the input file specified above, and which contain the outFilename project Folder.
#for system call from R to python to work, avoid folders with spaces; 
#but in Windows, spaces are handled by the script, making them passable.
outFilename <- "Modules-GO_C3_custom_SymbolsFile"
#SUBFOLDER WITH THIS NAME WILL BE CREATED
GOeliteFolder <- "/Users/coyoung/GOelite/"#"c:/Users/Ericda~1/Documents/wgcnatest/GOelite/"
#has subfolder and python script for GO-Elite v1.2.5 (script should be edited per authors' instructions to Divya Nandakumar to correct bug for Fisher Exact Test
# i.e.,          GO-Elite_v.1.2.5-Py/GO_Elite.py
#for system call from R to python to work, avoid folders with spaces; 
#but in Windows, spaces are handled by the script, making them passable.
maxBarsPerOntology=5
#Ontologies per ontology type, used for generating the PDF report; does not limit GO-Elite output
speciesCode="Hs"
#Hs for homo sapiens, Dm for fly; Mm, mouse; Rn, rat... (must have database downloaded via command line)
#if you use the GUI for GO-Elite, create a copy of the folders for command line, and delete databases downloaded via GUI.
downloadDB=FALSE
#If TRUE, the database files for speciesCode species will be downloaded "(be patient)" from Ensembl (v62 preferred)
pythonPath <- "/usr/bin/"
#python.exe for python v2.7 is here.
panelDimensions=c(1,1)
#rows, columns for each page of PDF output
color=c("darkseagreen3","lightsteelblue1","lightpink4")
#colors respectively for ontologyTypes:
#"Biological Process","Molecular Function","Cellular Component"
#must be valid R colors
modulesInMemory=TRUE
#uses cleanDat, net, and kMEdat from pipeline already in memory
ANOVAgroups=FALSE
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!

#ADVANCED OPTION
#########################################################################################################
#customDBcmd=paste0("--customSet ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/Databases/EnsMart62Plus/C2/ --dataToAnalyze all ")
#customDBcmd=paste0("--customSet ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/Custom_Databases/EnsMart62Plus/C2/ --dataToAnalyze all ")
customDBcmd=paste0("--customSet ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/Custom_Databases/EnsMart62Plus/C3/ --dataToAnalyze all ")
#set to "" if you have no custom database. Change C2/C3 to other custom db's (C4/C6/C7/H).
customPanelDimensions=c(1,1)
customReportMaxBars=20

#REST OF THE CODE IS AUTOMATIC
#########################################################################################################
## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
# this run in terminal window
filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
pythonPath=paste0(paste( sapply(do.call(c,strsplit(pythonPath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
GOeliteFolder=paste0(paste( sapply(do.call(c,strsplit(GOeliteFolder,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")


## The input files for GO-Elite are text files with the gene list as the 1st column, a symbol identified (gene symbol, uniprot etc) as the 2nd column
## Different accepted inputs are given in the tutorial
## Commonly used symbols - Gene Symbol - Sy (example of input file below)
### GeneSymbol		SystemCode (Symbol format)
###	  GFAP		Sy
###	  APOE		Sy
## All input files are placed in one folder

## The background file is prepared similarly and is placed in a separate folder
## The initial part of the code prepares files for GO-Elite. This can be skipped if the files are being made manually as described above.
## The second part of the code runs GO-ELite either from R (using the system command) or can be run using the terminal (in mac)
## The second part requires GO-Elite to be installed and path to the GO-Elite installation site indicated following python
## The 3rd part of the code plots the results from the GO-Elite results folder. When using the GUI the 1st 2 parts can be skipped and only the 3rd part can be used for plotting

##-------------------------------##
## Preparing files for GO-Elite ##
## Takes in the module assignment file as input with 1st column having gene names, 2nd column having color assignments followed by kME values

dir.create(file.path(filePath, outFilename))

##1a. GO Elite of the significant up and down (p<0.05) proteins in the current cleanDat
### GO Elite analysis for ANOVA-defined categories ###
if (ANOVAgroups) {
  sigThresh=0.05;
  
  #colnames(ANOVAout)
  #********************
  ##This code relies on a pre-existing data frame ANOVAout already in the environment, processed through volcano output for selected pairwise comparisons of interest.
  #numComp=6 #number of pairwise comparisons for ANOVA+Tukey p value columns, which are followed by the same-order log2(mean difference) columns
  #********************
  
  ANOVAout$Symbol <- do.call("rbind", strsplit(as.character(rownames(ANOVAout)), "[|]"))[, 1]
  
  DEXlistsForGO<-list()
  iter=0
  for (i in testIndexMasterList) {
    iter=iter+1;
    j=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".down")
    k=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".up")
    if (length(intersect(i,flip))==1) {
      #flipped sign (all diffs >0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
    } else {
      #do not flip sign (all diffs <0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
    }
  }
  
  #write lists to GOElite input files, and also the background file
  for (i in names(DEXlistsForGO)) { 
    dfGO<-data.frame(GeneSymbol=DEXlistsForGO[[i]],SystemCode=rep("Sy",length(DEXlistsForGO[[i]])))
    write.table(unique(dfGO),file=paste(filePath,outFilename,"/",i,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
  }
  #write background
  background <- unique(ANOVAout$Symbol)
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  dir.create(file.path(paste0(filePath,outFilename),"background"))
  write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  nModules=length(names(DEXlistsForGO))
  WGCNAinput=FALSE
  
} else { #NOT creating lists from ANOVA/volcano up & down groups
  
  
  ##1b. GO Elite of the WGCNA modules from specified input file or in the current cleanDat, net, and kME table
  ### GO Elite analysis for WGCNA-defined modules ###
  
  ##use data structures in memory if modulesInMemory=TRUE; otherwise read csv following the template which is written in an earlier R session and edited in excel to produce module membership table saved as .csv:
  if (modulesInMemory) {
    modulesData <- cbind(rownames(cleanDat),net$colors,kMEdat)
    WGCNAinput=TRUE
  } else {
    modulesData <- read.csv(paste(filePath,fileName,sep=""),header=TRUE, sep=",");
    # check if this is a WGCNA modules/kME table or simple list input
    if(length(na.omit(match("net.colors",colnames(modulesData))))>0) {
      WGCNAinput=TRUE
      # Remove first and last column if it contains sort information ("Original Order" and "color|kMEin")
      modulesData <- modulesData[,-c(1,ncol(modulesData))];
    } else {
      WGCNAinput=FALSE
    }
  }
  
  if (WGCNAinput) {
    library(WGCNA) #for labels2colors
    
    # Include column with Symbol (if it is gene symbol, if not use appropriate code as given in GO-Elite manual)
    modulesData$SystemCode <- rep("Sy",nrow(modulesData)) 
    
    # Assign Names of First columns, in case they are non standard
    colnames(modulesData)[1]<-"Unique.ID" #This should have Symbol|UniprotID
    colnames(modulesData)[2]<-"net.colors" #This should have colors
    
    #Split out symbols from UniprotIDs, keep symbols in column 1
    rownames(modulesData)<-modulesData$Unique.ID
    modulesData$Unique.ID<-do.call("rbind",strsplit(as.character(modulesData$Unique.ID), "[|]"))[,1]
    
    ## Creating background file for GO Elite analysis
    background <- unique(modulesData[,"Unique.ID"])
    background <- cbind(background,rep("Sy",length=length(background)))
    colnames(background) <- c("GeneSymbol","SystemCode")
    dir.create(file.path(paste0(filePath,outFilename),"background"))
    write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
    
    # Separate into independent module txt files for analysis by GO-Elite (CREATE INPUT FILES)
    greySubtractor=if(length(which(modulesData$net.colors=="grey"))>0) { 1 } else { 0 } #remove grey from count of modules
    nModules <- length(unique(modulesData$net.colors))-greySubtractor
    moduleColors <- uniquemodcolors <- labels2colors(c(1:nModules)) 
    for (i in 1:length(moduleColors)) {
      moduleName <- moduleColors[i]
      ind <- which(colnames(modulesData) == gsub("kMEME","kME",paste("kME",moduleName,sep="")))
      moduleInfo <- modulesData[modulesData$net.colors == gsub("ME","",moduleName), c(1,ncol(modulesData),ind)]
      colnames(moduleInfo) <- c("GeneSymbol","SystemCode","kME")
      if (moduleName == "blue" | moduleName == "brown" | moduleName == "green" | moduleName == "cyan") { write.table(moduleInfo,file=paste(filePath,outFilename,"/",moduleName,"_2_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
      } else {
        write.table(unique(moduleInfo),file=paste(filePath,outFilename,"/",moduleName,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
      }
    }
  } else { #input is not WGCNA kME table format
    
    ##1c. GO Elite of the WGCNA modules from specified input file, which must be in a column-wise list format, and including longest such list as background.
    # We process the input file as simple lists by column in the CSV (largest list used as background)
    
    #reread the file to a list of gene symbol (or UniqueID) lists
    modulesData <- as.list(read.csv(paste(filePath,fileName, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 
    
    nModules <- length(names(modulesData))
    for (a in 1:nModules) {
      modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
      modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
      modulesData[[a]] <- do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1]
    }
    ## Creating background file for GO Elite analysis
    background <- modulesData[order(sapply(modulesData,length),decreasing=TRUE)][[1]]
    background <- unique(background)
    background <- cbind(background,rep("Sy",length=length(background)))
    colnames(background) <- c("GeneSymbol","SystemCode")
    dir.create(file.path(paste0(filePath,outFilename),"background"))
    write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
    
    # Separate Symbol Lists into independent module txt files for analysis by GO-Elite (CREATE INPUT FILES)
    modulesData[[ names(modulesData[order(sapply(modulesData,length),decreasing=TRUE)])[1] ]] <- NULL
    nModules = nModules -1 #no background
    listNames <- uniquemodcolors <- names(modulesData)
    for (i in listNames) {
      listName <- i
      listInfo <- cbind(modulesData[[listName]],rep("Sy",length=length(modulesData[[listName]])))
      colnames(listInfo) <- c("GeneSymbol","SystemCode")
      write.table(unique(listInfo),file=paste(filePath,outFilename,"/",listName,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    }
  } #end if (WGCNAinput)
} #end else for if (ANOVAgroups)


##2. GO Elite Python Call
####----------------------- GO-ELite Analysis in Command Prompt/ Terminal or from R ------------------------------------#####

# Input and Denominator Files were prepared in the specified format (3 columns - genelist, format (Sy- for gene symbol) and kME)
# Each module is a separate txt file. Denominator file contains background list whcih is all identified proteins in this case
# All input files are placed in input or geneInfo folder and background file in denominator folder

if (downloadDB) {
  cat(paste0("Downloading Ensembl v62 database for species code ",speciesCode,". (Be patient.)\n\n"))
  system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --update Official --species ",speciesCode," --mod Ensembl --version EnsMart62Plus")) )
}

commandLine=gsub("'",'"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/"))
cat(paste0("NOW RUNNING THE FOLLOWING COMMAND:\n\n", commandLine,"\n\n(Estimated time for ", nModules, " lists to complete: ",round((30*nModules)/60,1)," minutes)\n Start time: ",Sys.time(),"\n"))
system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/")) )


##3. Output Report of Z-Score Barplots, processing all GO-Elite output files
#########################################################################################################
#-----------------------------------------Plotting for modules------------------------------------------#
# this script plots the top 3 ontologies for biological process, mol function and cell component for each module
#color scheme 4 ontology type key/legend (can be changed in user parameters, editing the "color" vector)
ontologyTypes=c("Biological Process","Molecular Function","Cellular Component")

if(ANOVAgroups) {
  xlabels <- names(DEXlistsForGO)
  xlabels.frame <- data.frame(Colors=rep(NA,length(xlabels)),Labels=xlabels)
  uniquemodcolors <- names(DEXlistsForGO) #not set above
} else {
  xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
  xlabels1 <- paste("M",seq(1:nModules),sep="")
  xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
}

setwd(paste0(filePath,outFilename,"/"))
pdf(paste0("GO-Elite_",outFilename,".pdf"),height=10,width=8)
op <- par(mfrow=panelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
#frame()
GOEliteOUTfileTrailer<-if(ANOVAgroups | !WGCNAinput) { c("-GO_z-score_elite.txt"); } else { c("_Module-GO_z-score_elite.txt"); }
summary <- list()
for(i in c(1:(length(uniquemodcolors)))){
  thismod=uniquemodcolors[i]	
  if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
    if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
      tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
    } else {
      next
    } 
  } else {
    tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t") #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
  }
  if (length(tmp[,2]) == 0) next
  tmp = tmp[,c(2,3,9,10,11)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
  tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
  tmp2 = tmp1[order(tmp1$Ontology.Type,decreasing=T),] #was tmp2
  tmp3 = tmp2[tmp2$Ontology.Type == "biological_process",][c(1:maxBarsPerOntology),]
  tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "molecular_function",][c(1:maxBarsPerOntology),] )
  tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "cellular_component",][c(1:maxBarsPerOntology),] )
  tmp3 <- na.omit(tmp3)
  #	tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),] #added this row, if you want to mix ontology types and sort by Z.Score only
  tmp3 <- tmp3[rev(rownames(tmp3)),]
  
  summary[[i]] <- tmp3
  
  ### To color bars by mol function, cell component or biological process
  for (j in 1:nrow(tmp3)){
    if (tmp3$Ontology.Type[j] == "molecular_function"){
      tmp3$color[j] <- color[2]
    } else if (tmp3$Ontology.Type[j] == "cellular_component"){
      tmp3$color[j] <- color[3]
    } else if (tmp3$Ontology.Type[j] == "biological_process"){
      tmp3$color[j] <- color[1]
    }
    # tmp3$color[j] <- uniquemodcolors[i] #module color for all bars, instead of different colors by ontology type
  } 
  
  if (tmp3$Z.Score == F) next
  par(mar=c(4,15,4,3))
  xlim <- c(0,1.1*max(tmp3$Z.Score))	
  moduleTitle <- xlabels.frame[i,"Labels"]
  xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=1.5,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
  abline(v=1.96,col="red", cex.axis = 0.5)
  axis(2, at=xh, labels = tmp3$Ontology.Name, tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
}

par(op) # Leaves the last plot
dev.off()


if(!customDBcmd=="") {
  library(stringr)
  pdf(paste0("GO-Elite_",outFilename,"-CUSTOM_db.pdf"),height=10,width=8)
  op <- par(mfrow=customPanelDimensions,oma=c(0,0,3,0))
  frame()
  #legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
  #legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
  
  GOEliteOUTfileTrailer<-if(ANOVAgroups | !WGCNAinput) { c("-UserSuppliedAssociations_z-score_elite.txt"); } else { c("_Module-UserSuppliedAssociations_z-score_elite.txt"); }
  summary <- list()
  for(i in c(1:(length(uniquemodcolors)))){
    thismod=uniquemodcolors[i]	
    if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
      if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
        tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
      } else {
        next
      } 
    } else {
      tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t")
    }
    if (length(tmp[,2]) == 0) next
    tmp = tmp[,c(1,1,7,8,12)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
    tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
    tmp2 = tmp1
    tmp3 = tmp2
    tmp3 <- na.omit(tmp3)
    tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),][c(1:customReportMaxBars),] #added this row, if you want to mix ontology types and sort by Z.Score only
    tmp3 <- tmp3[rev(rownames(tmp3)),]
    
    summary[[i]] <- tmp3
    
    ### To color bars by mol function, cell component or biological process
    for (j in 1:nrow(tmp3)){
      tmp3$color[j] <- color[1]
      if(WGCNAinput) { tmp3$color[j] <- uniquemodcolors[i] } #module color for all bars, instead of different colors by ontology type
    } 
    
    #	if (tmp3$Z.Score == F) next
    if (is.na(max(tmp3$Z.Score))) tmp3<-na.omit(tmp3)
    
    par(mar=c(4,15,4,3))
    xlim <- c(0,1.1*max(tmp3$Z.Score))	
    moduleTitle <- xlabels.frame[i,"Labels"]
    xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=1.5,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
    abline(v=1.96,col="red", cex.axis = 0.5)
    axis(2, at=xh, labels = str_to_title(gsub("_"," ",tmp3$Gene.Set.Name)), tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
  }
  
  par(op) # Leaves the last plot
  dev.off()
} #end if(!customDBcmd=="")

setwd(rootdir)
#########################################################################################################
          # For GO-Elite run on networkScreening output see "networkScreenforG0" codes #
#########################################################################################################

save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))


#########################################################################################################
#=======================================================================================================#
              # Secondary ANOVA on T-Staging OR highly correlated values in GNP #
#=======================================================================================================#
#code derived from: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/ 

############################################
#column of interest 
Grouping.2nd <- numericMeta$pathologic_t.1
Grouping.2nd <- paste0("TStg", Grouping.2nd)
cleanDat.sub<-cleanDat
############################################

ANOVA2data = as.data.frame(cbind(colnames(cleanDat.sub), Grouping.2nd, t(cleanDat.sub)))
colnames(ANOVA2data)[1:2]<-c("TCGA_ID","Grouping.2nd") # just changes the 1st 2 rows to CODE & SampleType

# improved ANOVA code from datasheet/data table code immediately below

datalist2 = list()
for(gene in colnames(ANOVA2data[3:13488]))
{
  print(gene)
  b = aov(ANOVA2data[[gene]] ~ Grouping.2nd, data = ANOVA2data)
  p = TukeyHSD(b) # bascially the adjusted p-value, upper & lower quartiles 
  tukresult1 <- data.frame(p$Grouping.2nd) 
  datalist = list() # creates object that contains lists, empty 
  for (i in rownames(tukresult1)) # i is the comparision name 
  {
    beta=tukresult1[i,] # isolate row for each comparison of tukey results
    rownames(beta)=NULL
    #colnames(beta)<-paste0(i,colnames(beta),sep=".")
    beta=beta[-c(2:3)] # removes columns that have upper (upr) & lower (lwr), upper & lower quatiles?
    datalist[[i]] <- beta #gathers all differences in p-values across all genes in all comparisons 
    big_data = do.call(cbind, datalist) # appending cols onto big_data for each gene
    datalist2[[gene]] <- big_data
  }
  big_data2 = do.call(rbind, datalist2)
}
write.table(big_data2,file=paste0(outputtabs,"2nd-ANOVA.txt"), append=TRUE,sep=",") # follow by opening readng data into excel (i.e. opening excel sheet and importanting data w/comma separators)


# # To develop a datasheet of ANOVA results
# for(gene in colnames(ANOVA2data[2:13488]))
# {
#   print(gene)
#   b = aov(ANOVA2data[[gene]] ~ Grouping.2nd, data = ANOVA2data)
#   p = TukeyHSD(b)
#   tukresult1 <- data.frame(p$Grouping.2nd)
#   g = cbind(gene, tukresult1)
#   write.table(g,file=paste0(outputtabs,"Dectukeyresults.csv"), append=TRUE,sep=",") 
#   # write.table(g,file=paste0(outputtabs,"ANOVA.2nd.Tukeyresults.csv") append=TRUE, sep=",")
# }

#=======================================================================================================#
            # SKIP: Boxplot for select gens based on ANOVA, skip if no boxplots are warranted #
#=======================================================================================================#

pdf(paste0(outputfigs,"SecondANOVA",".pdf"),height=10,width=8)
#pdf(paste0(outputfigs,"/iGraph_Modules-",FileBaseName,"-nonCAIRO-LargeNodes.pdf"),height=9,width=10)

for(gene in colnames(ANOVA2data[2:13488]))
{
  p <- ggboxplot(ANOVA2data, x = "Grouping.2nd", y = gene, color = "Grouping.2nd", ylab = "Expression", xlab = "Pathologic_T_Staging", legend = "none", palette = c("red", "gold", "mediumblue", "orange"), main = gene, add = "jitter") +rotate_x_text(angle = 45) + stat_compare_means(method = "kruskal.test", label.y = 6) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = NULL)
  print(p)
}
dev.off()

# the first stat_compare_means does and anova on all groups and displays the p-value at y axis value "6". I did this to ensure that all of my boxplots would have the same axis length
## the second stat_compare means does individual comparisons between all groups and specifies that the control group is NULL, this could be changed to a specific control group in which each of the levels in the group would be compared to the reference (i.e. control group) group.
### the third stat_compare means in this code adds pvalues and brackets for individual comparisons


# # Tiara: Alternative code just in case you want to display pvalues between select groups in your figure
# my_comparisons <- list( c("A", "B"), c("A", "C"), c("A", "D"),c("A", "E"),c("A", "Normal CD34+ Control")) # subsititue "Normal CD34+ Control"for your control, this should alos match your ref.group in the below code
# 
# # pdf(paste0(outputfigs,"SecondANOVA",".pdf"),height=10,width=8)
# for(gene in colnames(ANOVA2data[2:13488]))
# {
#   p <- ggboxplot(ANOVA2data, x = "Grouping.2nd", y = gene, color = "Grouping.2nd", ylab = "Expression", xlab = "Pathologic_T_Staging", legend = "none", palette = c("red", "gold", "mediumblue", "orange"), main = gene, add = "jitter") +rotate_x_text(angle = 45) + stat_compare_means(method = "kruskal.test", label.y = 6) + stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = NULL) + stat_compare_means(comparisons = my_comparisons)
#   print(p)
# }
# dev.off()

#save.image(paste0(rdata,"NoTAMPOR-updated-complete.saved.image.",projectFilesOutputTag,".Rdata"))

#########################################################################################################
#=======================================================================================================#
                            # Identify distinct clusters in patietns #
#=======================================================================================================#
# Tiara: I’m providing you with the code I used to identify distinct clusters in patients, KM and chi-square statistics for the clusters, and an alternative heatmap visualization method(pheatmap).
# code derived from: https://compgenomr.github.io/book/clustering-grouping-samples-based-on-their-similarity.html 
# distance & hierarchal clustering fxn from NMF heatmap were used to cluster patients.See lin for information on lcuster part of code (mulitple ways to calculate distance & cluster patients)




