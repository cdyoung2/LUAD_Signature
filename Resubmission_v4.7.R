###############################################################################################################
# Project: WGCNA Manuscript Resubmission
# Parts 1-4: Pulling count data, running DESeq2, regression on sex & age
# Part: Addressing generalizabilty comments
# 
###############################################################################################################

rm(list=ls(all=TRUE))
dev.off()
#options(stringsAsFactors = FALSE)
options(stringsAsFactors=FALSE)
rootdir="/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Resubmission/"
setwd(rootdir)

# Download TCGA Data Transfer Tool
# https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

# Downloading Data Using a Manifest File 
# https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/

# Run in terminal: /Users/youngcd/Downloads/gdc-client download -m /Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Signature_Comparison/2023wrapup/STAR\ Counts\ -\ Gene\ Expression\ Quant/gdc_manifest.2024-09-19.txt
###############################################################################################################

# Part 1: Downloading & setting up Counts Data
setwd("./count_files/")
caseDirs<-list.dirs()
caseDirs<-caseDirs[!grepl("logs",caseDirs)]  #remove ./logs/ subfolders from vector of folder names.
caseDirs<-caseDirs[!grepl("\\.$",caseDirs)]
length(caseDirs) #509 directories for 509 cases (original Dammer run was 557)

# BiocManager::install("TCGAutils")
library(TCGAutils)
UUID.Filename.df<-read.delim(file="gdc_manifest.2024-09-19.txt",sep="\t",header=TRUE)
caseBarcodes.df<-TCGAutils::UUIDtoBarcode(UUID.Filename.df$id, from_type = "file_id")

library(doParallel)
num_cores <- detectCores(logical = TRUE) 
num_cores <- num_cores - 1
clusterLocal <- makeCluster(c(rep("localhost",num_cores)),type="SOCK")
registerDoParallel(clusterLocal)
parallel::clusterExport(cl=clusterLocal, varlist=c("UUID.Filename.df"), env=globalenv())

loadSTARtables <- function(caseDir) {
  thisFile <- list.files(path=caseDir, pattern=".tsv")
  fileUUID = UUID.Filename.df$id[match(thisFile,UUID.Filename.df$filename)]
  caseSTARtable<-read.csv(file=paste0(caseDir,"/",thisFile),sep="\t",header=FALSE)
  caseSTARtable<-caseSTARtable[-1,]
  colnames(caseSTARtable)<-caseSTARtable[1,]
  caseSTARtable<-caseSTARtable[-1,]
  return(caseSTARtable)
}

#caseSTARtables<-parLapply(cl=clusterLocal, caseDirs, loadSTARtables)
caseSTARtables<-list()
caseSTARtables<-foreach::foreach (caseDir=caseDirs) %dopar% loadSTARtables(caseDir)
names(caseSTARtables) <- caseBarcodes.df$associated_entities.entity_submitter_id[match(gsub("^\\./","",caseDirs),caseBarcodes.df$file_id)]
geneIDs.allTables<-unique(unlist( parLapply(cl=clusterLocal, caseSTARtables,function(x) x$gene_id) ))
length(geneIDs.allTables)
#[1] 60664
#output not in df format
#geneNamesLookupTable<-unique(unlist( parLapply(cl=clusterLocal, caseSTARtables,function(x) return(c(x$gene_name,x$gene_id))) ))
#slow
#options(future.globals.maxSize= 10000*1024^2)  #10,000 megabyte memory limit
#geneNamesLookupTable<-foreach::foreach (caseBarcode=names(caseSTARtables), .combine=rbind) %dopar% caseSTARtables[[caseBarcode]][,c("gene_name","gene_id")]
#geneNamesLookupTable.unique.df<-data.frame(EnsemblID=geneIDs.allTables, geneName=geneNamesLookupTable[match(geneIDs.allTables,geneNamesLookupTable[,1]),2])
# more efficient:
geneNamesLookupTable2 <- unique( lapply( lapply(caseSTARtables,function(x) return(x[,1:2])), rbind) )[[1]]
geneNamesLookupTable2$UniqueID <- paste0(geneNamesLookupTable2$gene_name, "|", geneNamesLookupTable2$gene_id)
geneNamesLookupTable2$UniqueID.noVer <- paste0(geneNamesLookupTable2$gene_name, "|", as.data.frame(do.call(rbind,strsplit(geneNamesLookupTable2$gene_id,"[.]")))[,1])  #without .version numbers (not unique, yet)

#Finalize COUNTS table for all genes, all LUAD samples (includes paired biospecimens (primary tumor, blood, normal tissue), and some technical replicate aliquots)
parallel::clusterExport(cl=clusterLocal, varlist=c("geneIDs.allTables"), env=globalenv())
#for FPKM: FPKMlist.allCases.allGeneIDs<- parLapply(cl=clusterLocal, caseSTARtables,function(x) x$fpkm_unstranded[match(geneIDs.allTables,x$gene_id)])
COUNTSlist.allCases.allGeneIDs<- parLapply(cl=clusterLocal, caseSTARtables,function(x) x$unstranded[match(geneIDs.allTables,x$gene_id)])
COUNTStable.allCases.allGeneIDs<-as.data.frame(COUNTSlist.allCases.allGeneIDs,check.names=FALSE)
COUNTStable.allCases.allGeneIDs<-apply(COUNTStable.allCases.allGeneIDs,2,as.numeric)
rownames(COUNTStable.allCases.allGeneIDs)<-geneNamesLookupTable2$UniqueID[match(geneIDs.allTables,geneNamesLookupTable2$gene_id)]
###############################################################################################################

# Part 2: Building, editing & matching count matrix to FPKM matrix
# Loading in numericMeta & matching to cleanDat samples (TCGA-IDs); see Resubmission_v2.R for Dammer full scale code with addtional steps
numericMeta <- readRDS('/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/RData/numericMeta.rds')
cleanDat.org <- readRDS('/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/RData/cleanDat.rds')
cleanDat <- cleanDat.org

# Filtering counts data to 461 samples
#cleanDat <- COUNTStable.allCases.allGeneIDs[,colnames(tumor.exp)] (ended just using tumor.exp)
dim(COUNTStable.allCases.allGeneIDs) # 60664, 509
newCounts <- COUNTStable.allCases.allGeneIDs
id_cleaned <- substr(colnames(newCounts), 1, 12)
id_cleaned <- gsub("-", ".", id_cleaned)
colnames(newCounts) <- id_cleaned
newCounts <- newCounts[, match(rownames(numericMeta), colnames(newCounts))]

# Filtering counts data to 13486 genes
newCounts <- newCounts[-c(1:4), ] # Remove the first 4 rows from newCounts
rownames(newCounts) <- gsub("\\.\\d+$", "", rownames(newCounts))
# Extract the part after the "|" in both newCounts and cleanDat rownames without modifying them
newCounts_suffix <- sub(".*\\|", "", rownames(newCounts))
cleanDat_suffix <- sub(".*\\|", "", rownames(cleanDat))
# Match based on the suffix and reorder newCounts based on cleanDat
newCounts_matched <- newCounts[match(cleanDat_suffix, newCounts_suffix),]
# Keep the full rownames from cleanDat in the matched newCounts
rownames(newCounts_matched) <- rownames(cleanDat)
newCounts <- newCounts_matched

# Adding batch # & condition to numericMeta (grouping/staging)
numericMeta$Grouping<-rep(NA,nrow(numericMeta))
numericMeta$Grouping = numericMeta$Group
numericMeta$Grouping <- paste0("Stg",numericMeta$Grouping)
#numericMeta$Grouping[Grouping=="StgNA"] <- NA
head(numericMeta$Grouping)
# [1] "Stg2" "Stg2" "Stg3" "Stg1" "Stg1" "Stg3"
numericMeta$condition <- numericMeta[,202]
numericMeta$Batch <- factor(numericMeta$Batch)
numericMeta$condition <- factor(numericMeta$condition)

# Check for rows with any NAs (DESeq2 does not run w/NAs in countData)
rows_with_na <- rownames(newCounts)[apply(newCounts, 1, function(x) any(is.na(x)))]
length(rows_with_na) #114
dim(newCounts) #13486, 461
newCounts <- newCounts[!(rownames(newCounts) %in% rows_with_na), ]
dim(newCounts) #13372, 461
###############################################################################################################

# Part 3: Running simple DESeq2 workflow (normalization & diffEx)
library(DESeq2)
library(BiocParallel)
library(parallel)
detectCores() #10
multicoreWorkers() #8
register(MulticoreParam(workers = 8))  

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = newCounts,
                              colData = numericMeta,
                              design = ~ Batch + condition)
dds  <- DESeq(dds, parallel = TRUE)
res <- results(dds)

# Comparing other conditions without changing the reference
res_Stg3_vs_Stg2 <- results(dds, contrast = c("condition", "Stg3", "Stg2")) # Compare Stg3 vs Stg2
res_Stg4_vs_Stg2 <- results(dds, contrast = c("condition", "Stg4", "Stg2")) # Compare Stg4 vs Stg2
res_Stg4_vs_Stg3 <- results(dds, contrast = c("condition", "Stg4", "Stg3")) # Compare Stg4 vs Stg3

res_all <- list(
  Stg2_vs_Stg1 = results(dds, name="condition_Stg2_vs_Stg1"),
  Stg3_vs_Stg1 = results(dds, name="condition_Stg3_vs_Stg1"),
  Stg4_vs_Stg1 = results(dds, name="condition_Stg4_vs_Stg1"),
  Stg3_vs_Stg2 = results(dds, contrast=c("condition", "Stg3", "Stg2")),
  Stg4_vs_Stg2 = results(dds, contrast=c("condition", "Stg4", "Stg2")),
  Stg4_vs_Stg3 = results(dds, contrast=c("condition", "Stg4", "Stg3")))

saveRDS(res_all, file = '/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/res_all.rds')
###############################################################################################################

# Part 4: Regression code
normalized_counts <- counts(dds, normalized=TRUE)
dim(normalized_counts) #13372, 461

# Replace all 0 values with 0.5 in the normalized_counts matrix as log2() 0 cause Inf in df & cause an error in regression 
normalized_counts[normalized_counts == 0] <- 0.5

normalized_counts <- log2(normalized_counts)
dim(normalized_counts) #13372, 461

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

normalized_counts.unreg<-normalized_counts
dim(numericMeta) #461, 203

library("doParallel")
#when Eric is running at Emory (requires RSA public key-ssh, & manual run of shell script from command prompt to start server backend):  
#parallelThreads=30
#clusterLocal <- makeCluster(c(rep("haplotein.biochem.emory.edu",parallelThreads)), type = "SOCK", port=10191, user="edammer", rscript="/usr/bin/Rscript",rscript_args="OUT=/dev/null SNOWLIB=/usr/lib64/R/library",manual=FALSE)
##OR to run parallel processing threads for regression: :
parallelThreads=8 #max is number of processes that can run on your computer at one time
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
# normalized_counts<-normalized_counts[,na.omit(match(rownames(numericMeta),colnames(normalized_counts)))] # rows in numericaMeta == columns in normalized_counts
# rownames(numericMeta)==colnames(normalized_counts) #sanity check -- are sample names in same order?
# 
# numericMeta2<-numericMeta[!is.na(numericMeta$age_at_initial_pathologic_diagnosis),] 
# numericMeta = numericMeta2
# normalized_counts<-normalized_counts[,na.omit(match(rownames(numericMeta),colnames(normalized_counts)))] 
# rownames(numericMeta)==colnames(normalized_counts)
# 
# numericMeta2<-numericMeta[!is.na(numericMeta$pathologic_t.1),] 
# numericMeta = numericMeta2
# normalized_counts<-normalized_counts[,na.omit(match(rownames(numericMeta),colnames(normalized_counts)))] 
# rownames(numericMeta)==colnames(normalized_counts)

dim(normalized_counts) #13470, 463
dim(numericMeta) #463, 158

condition <- as.factor(numericMeta$Group)
condition2 <- as.factor(numericMeta$pathologic_t.1)

age=as.numeric(numericMeta$age_at_initial_pathologic_diagnosis)
sex=as.factor(numericMeta$gender.2) # 0 becomes 1 & 1 becomes 2
#PMI=as.numeric(numericMeta$PMI)

regvars <- as.data.frame(cbind(condition, condition2, sex, age)) #,PMI))
#EBD: note the order of columns for sex and age have been swapped from what I do, this is fine, and gives preference to first modelling sex, then age

## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(normalized_counts),ncol=ncol(normalized_counts))
rownames(normExpr.reg) <- rownames(normalized_counts)
colnames(normExpr.reg) <- colnames(normalized_counts)
coefmat <- matrix(NA,nrow=nrow(normalized_counts),ncol=ncol(regvars)+2) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

#another RNG seed set for reproducibility
set.seed(8675309); #R's random number generator, seed is a number

if (parallelThreads > 1) {
  
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    set.seed(8675309)
    cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(normalized_counts), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(normalized_counts)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
    coefmat <- foreach (i=1:nrow(normalized_counts), .combine=rbind) %dopar% {
      set.seed(8675309)
      options(stringsAsFactors=FALSE)
      library(boot)
      thisexp <- as.numeric(normalized_counts[i,])
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
    #    normExpr.reg <- matrix(NA,nrow=nrow(normalized_counts),ncol=ncol(normalized_counts))
    #EBD: you had: coefmat[i,5]*regvars[,"age"]) -- but coefmat only had 4 columns, one for the main term, and one for each of your modelled covariates, plus a 4th column for residuals
    #EBD: I changed the column position for sex back (since I put condition back), and corrected the position for age
    normExpr.reg <-foreach (i=1:nrow(normalized_counts), .combine=rbind) %dopar% { (normalized_counts[i,]- coefmat[i,4]*regvars[,"sex"] - coefmat[i,5]*regvars[,"age"]) }
    #EBD: condition is column 2 in coefmat
    #EBD: sex is column 4 in coefmat
    #EBD: age is column 5 in coefmat
    ### All regression code below deos not apply to conditionS set when parallelThreads=15 & boot=TRUE ###
  } else { #linear model regression; faster but incomplete regression of Age, Sex, PMI effects, SO NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(normalized_counts)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(normalized_counts[i,])~condition +condition2+sex+age,data=regvars) #+PMI
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
    for (i in 1:nrow(normalized_counts)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(normalized_counts[i,])
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
    for (i in 1:nrow(normalized_counts)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(normalized_counts[i,])~condition +condition2+sex+age,data=regvars)  #+PMI
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

#write.csv(normalized_counts,file=,"normalized_counts.csv")
write.table(normalized_counts,file="/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/normalized_counts.Regressed.txt",sep="\t") #check and apply changes to static tail of filename if necessary

##Overwrite normalized_counts with regressed data
##(DO NOT RERUN OUT OF CONTEXT)
############################################
normalized_counts_regressed<-normExpr.reg
rownames(normalized_counts_regressed)<-rownames(normalized_counts.unreg)
#########################################################################################################

# Over-writing cleanDat with normalized_counts_regressed
cleanDat <- normalized_counts_regressed
dim(cleanDat) #13372, 461

# WGCNA
enableWGCNAThreads()
allowWGCNAThreads()
#library("doParallel")
#stopCluster(clusterLocal)
#parallelThreads=8#parallelThreads=15 #set to # of threads on your computer
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
#tic("timing")

#Check power and connectivity
powers <- seq(8,12,by=0.25)
sft <- pickSoftThreshold(t(cleanDat), powerVector=powers, blockSize=nrow(cleanDat)+1000, verbose=3, corFnc="bicor", networkType="signed")
#toc() #300 sec

#5   9.00    0.789 -3.16          0.939    70.2      62.5    219
#6   9.25    0.802 -3.13          0.940    62.6      55.2    206
#7   9.50    0.819 -3.07          0.944    55.9      48.7    194

# plotting SFT R^2 by power, chose value at the elbow with value R.sq above 0.8
plot(sft[[2]][,1],sft[[2]][,2],xlab="Power (Beta)",ylab="SFT R^2")
#tic("timing")
print("Running power blockwiseModules...")
power=9.25 #choose power at elbow of SFT R≤ curve approaching asymptote near or ideally above 0.80
enableWGCNAThreads()
allowWGCNAThreads()
getDoParWorkers() #6
getDoParRegistered() #TRUE
getDoParName() # doParallelSNOW"
getDoParVersion() # "1.0.15"

# Issue w/this cleanDat in modulePreservation, will need to run goodSamplesGenes() to remove genes with too many missing value later anyway why not run now and get a better net
cleanDat_badRemoved <- goodSamplesGenes(t(cleanDat))
cleanDat_badRemoved$allOK #must be TRUE, if FALSE removed
target.cleanDat_badremoved <- cleanDat[which(cleanDat_badRemoved$goodGenes!=FALSE),]
dim(target.cleanDat_badremoved) #13308, 461
cleanDat <- target.cleanDat_badremoved

## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=100,
                        mergeCutHeight=0.07,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets
print("...done running power blockwiseModules")
#toc() #1801 sec
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))
table(net$colors)["grey"]

# turquoise      M1    turquoise 1456
# blue           M2         blue  940
# brown          M3        brown  929
# yellow         M4       yellow  815
# green          M5        green  777
# red            M6          red  776
# black          M7        black  725
# pink           M8         pink  458
# magenta        M9      magenta  430
# purple        M10       purple  392
# greenyellow   M11  greenyellow  332
# tan           M12          tan  301
# salmon        M13       salmon  298
# cyan          M14         cyan  297
# midnightblue  M15 midnightblue  285
# lightcyan     M16    lightcyan  269
# grey60        M17       grey60  211
# lightgreen    M18   lightgreen  170
# lightyellow   M19  lightyellow  170
# royalblue     M20    royalblue  144
# grey 3133

#save.image("/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/counts_WGCNA.Rdata")
save.image("/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/counts_WGCNAv2.Rdata")
#########################################################################################################

# MODULE PRESERVATION

options(stringsAsFactors=FALSE)

# Configuration Section
rootdir <- "/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/CountsvFPKM_Preservations/" ##path to working dir, ending with /
template.Rdata.path <- "/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/NoTAMPOR+Liu+SP-updated-complete.saved.image.TCGA_LuCa.Rdata" #template network; full or 'forORA' RData relative or absolute file path
target.Rdata.path <- "/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/counts_WGCNAv2.Rdata" #target network (usually ...for ORA.RData)

template.net.label<- "Org"
target.net.label<-   "Count"

nPermute=500	#number of permutations (2+ hours for 500, recommended for publication.) 200 minimum for reliable calculation.

calculateOrLoadRDS <- "calculate"  #or [anything else] to load   #RDS file must exist with above labels and permutations in the paste0 concatenate format used when calculating.
#--------------------------------


template.cleanDat.variable<- "cleanDat"
template.net.variable<-      "net"

target.cleanDat.variable<-   "cleanDat"
target.net.variable<-      "net"
#########################################################################################################

setwd(rootdir)
workingDir<- rootdir

# target data
template.cleanDat.variable<- "cleanDat"
template.net.variable<-      "net"

load(template.Rdata.path) #template
setwd(workingDir)
rootdir<-workingDir
template.net<-eval(parse(text=template.net.variable))
template.cleanDat<-eval(parse(text=template.cleanDat.variable))
dim(template.cleanDat) #13486, 461

# template data
target.cleanDat.variable<-   "cleanDat"
target.net.variable<-      "net"

load(target.Rdata.path)  #target
setwd(workingDir)
target.net<-eval(parse(text=target.net.variable))
target.cleanDat<-eval(parse(text=target.cleanDat.variable))
dim(target.cleanDat) #13308, 461
#########################################################################################################

## Run an automated network analysis
table(template.net$colors) 
#get module color information
colorsign=template.net$colors

# setup matched information for rep1 and rep2
setLabels = c(template.net.label,target.net.label)
multiExpr = list(rep1 = list(data=t(template.cleanDat)), rep2 = list(data=t(target.cleanDat)))
multiColor = list(rep1=colorsign)

library(WGCNA)

# The number of permutations drives the computation time
# of the module preservation function. For a publication use 200 permutations.
# But for brevity, let's use a small number
nPermutations1=nPermute
# Set it to a low number (e.g. 3) if only the medianRank statistic
# and other observed statistics are needed.
# Permutations are only needed for calculating Zsummary
# and other permutation test statistics.
# set the random seed of the permutation test analysis
set.seed(1)


if(calculateOrLoadRDS=="calculate") {
  system.time({
    mp = modulePreservation(multiExpr, multiColor,
                            networkType="signed",corFnc="bicor",
                            referenceNetworks=1, nPermutations = nPermutations1,
                            randomSeed=1, quickCor=0, verbose=3)
  })
  
  ## Save the results of the module preservation analysis 
  saveRDS(mp, file=paste0("Module.Preservation-",template.net.label,"_presIn_",target.net.label,"-",nPermute,"permutations.RDS"))
} else {
  ## If needed, reload the data:
  if(file.exists(paste0(rootdir,paste0("Module.Preservation-",template.net.label,"_presIn_",target.net.label,"-",nPermute,"permutations.RDS")))) {
    mp<-readRDS(file=paste0(rootdir,paste0("Module.Preservation-",template.net.label,"_presIn_",target.net.label,"-",nPermute,"permutations.RDS")))
  } else {
    print(paste0("Configured .RDS file does not exist:\n","Module.Preservation-",template.net.label,"_presIn_",target.net.label,"-",nPermute,"permutations.RDS\nin path/folder:  ",rootdir,"\n\n"))
    break;
  }
}


# specify the reference and the test networks
ref=1; test=2

Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats

# View Z statistics from the permutation test analysis
Z.PreservationStats


# Let us now visualize the data.
#library(Cairo)
# Switch between Cairo and non-Cairo versions

for (version in c("nonCairo")) {
  # Use the standard pdf() function to save the output
  pdf(file=paste0("Module.Preservation-", template.net.label, "_presIn_", target.net.label, "-", nPermute, "permutations-nonCairo.pdf"), width=16, height=12)
  
  modColors = rownames(Obs.PreservationStats)
  moduleSize = Obs.PreservationStats$moduleSize
  
  # We will omit the grey module (background genes) and the gold module (random sample of genes)
  selectModules = !(modColors %in% c("grey", "gold"))
  
  # Text labels for points
  point.label = modColors[selectModules]
  
  # Composite preservation statistics
  medianRank = Obs.PreservationStats$medianRank.pres
  Zsummary = Z.PreservationStats$Zsummary.pres
  
  # Set up 1x2 plot layout
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  
  # Plot medianRank versus module size
  plot(moduleSize[selectModules], medianRank[selectModules], col = 1, bg = modColors[selectModules], 
       pch = 21, main = "medianRank Preservation", cex = 2, ylab = "medianRank", xlab = "Module size", log = "x")
  labelPoints(moduleSize[selectModules], medianRank[selectModules], point.label, cex = 1, offs = 0.03)
  
  # Plot Zsummary versus module size
  plot(moduleSize[selectModules], Zsummary[selectModules], col = 1, bg = modColors[selectModules], 
       pch = 21, main = "Zsummary Preservation", cex = 2, ylab = "Zsummary", 
       ylim = c(-1, max(Zsummary[selectModules], na.rm = TRUE) + 3), xlab = "Module size", log = "x")
  labelPoints(moduleSize[selectModules], Zsummary[selectModules], point.label, cex = 1, offs = 0.03)
  
  # Add threshold lines for Zsummary
  abline(h = 0)
  abline(h = 2, col = "blue", lty = 2)
  abline(h = 10, col = "red", lty = 2)
  
  # Close the PDF device
  dev.off()
}
save.image("/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/counts+ModulePresWGCNA.Rdata")
#########################################################################################################

#options(stringsAsFactors=FALSE)
#<...> run pipeline code or load RData with background network (cleanDat and net) objects

#DATA LOADING SECTION (SWAP IN YOUR DATA)
#########################################
#load("/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Rdata/counts+ModulePresWGCNA.Rdata")
# This RData contains a second uniquely-named net list object and cleanDat.transposed data frame or matrix

#Original WGCNA
# template.cleanDat
# template.net
#Counts WGCNA
# target.cleanDat
# target.net

rootdir="/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/CountsvFPKM+Preservations+ORA//"
setwd(rootdir)
outputfigs<-outputtabs<-rootdir

library(WGCNA)

##ORA (overrepresentation analysis) of network modules' gene symbols across 2 networks

#CONFIGURATION SECTION
##################################
inputDat1=template.cleanDat           #cleanDat for first network (rows are gene products, columns are samples)
inputDat2=target.cleanDat             #cleanDat for second network (''   '');     used for background
net1=template.net                     #network list object output by WGCNA::blockwiseModules() function for network 1
net2=target.net                       #''                                   ''                          for network 2

species1="human"                      # this would be your fly species if fly-human
species2="human"                      # no symbol conversion happens if species1==species2

# test="p" #"p" or "FDR" #both are handled below with a for loop

outfiletitle="Count_vs_ogWGCNA_modules"                             #e.g. "Net2_vs_Net1_modules"
heatmapTitle="Count Network Modules vs. ogWGCNA Network Modules"    #e.g. "Net2 Network Protein Modules vs. Net1 Network Protein Modules"

HeatmapNet1moduleType="ogWGCNA"         #Short phrase-nickname for network 1 module labels
HeatmapNet2moduleType="Counts"         #Short phrase-nickname for network 2 module labels

oneOrTwoTailed=1                     #Odds Ratio greater only (1-tailed) or both greater and less than tails (2-tailed)
orderByRelatedness=TRUE              #keep net$ME order (TRUE)? ...or size-rank order (FALSE)?

#If you want to manually tweak column order in the first-pass reordered data, call out column order of the AUTO-REORDERD (first-pass) ORA heatmap.
columnReorderOfAutoReorderedOutput=c()  #c(2,3,1,23,4,9,10,5,6,7,27,8,28,11,13,14,18,12,15,16,17,19:22,24:26,29:44)
#if not, leave this variable set to an empty vector, i.e.: c()
##################################


## NO EDITING NECESSARY BELOW -- RUN AS ONE CODE CHUNK

library(WGCNA)

if (!species1==species2) {
  fly = useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
  #	fly = useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
  human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
}

ORoneTailed <- function(q,k,m,t) {
  q #<-  ## Intersection of test list and reference list, aka number of white balls drawn
  m #<-  ## All genes in reference list, aka number of draws
  k #<-  ## All genes in test list, aka number white balls
  t #<-  ## Total number of genes assessed, aka black plus white balls
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE,alternative="greater")
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

ORtwoTailed <- function(q,k,m,t) {
  q #<-  ## Intersection of test list and reference list, aka number of white balls drawn
  m #<-  ## All genes in reference list, aka number of draws
  k #<-  ## All genes in test list, aka number white balls
  t #<-  ## Total number of genes assessed, aka black plus white balls
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

## count overlaps and run the analysis
ORA <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  if(oneOrTwoTailed==1) { empvals <- ORoneTailed(q,k,m,t) } else { empvals <- ORtwoTailed(q,k,m,t) }
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}


#***
geneInfo.BLSA<-inputDat1
geneInfo.BLSA$Symbol=do.call("rbind",strsplit(as.character(rownames(inputDat1)),"[|]"))[,1]

geneInfo.Emory<-inputDat2
geneInfo.Emory$Symbol=do.call("rbind",strsplit(as.character(rownames(inputDat2)),"[|]"))[,1]

##Convert to same species lists
genelist1 <- data.frame(unlist(geneInfo.BLSA$Symbol),ncol=1)
colnames(genelist1)<-"species1"
genelist1<-data.frame(genelist1[which(!genelist1[,1]==""),"species1"],ncol=1)
colnames(genelist1)<-"species1"
genelist1<-data.frame(genelist1$species1[!grepl("'", genelist1$species1)],ncol=1)
colnames(genelist1)<-c("species1","species2")


### Must rerun whole chunk of code from *** above or $Symbol will lookup previously lookedup values!
if (!species1==species2) {
  if (species1=="fly") {
    genelist.clean<-getLDS(attributes="external_gene_name", filters="external_gene_name", values=genelist1$species1, mart=fly, attributesL="hgnc_symbol",martL = human)
    #		write.table(listHumanFly[,2],file=paste(outputDir,outputFile,"_AllGenes_human.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    #		write.table(modulesData$GI,file=paste(outputDir,outputFile,"_AllGenes_fly.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  }
  if (species1=="human") {
    genelist.clean<-getLDS(attributes="hgnc_symbol", filters="hgnc_symbol", values=genelist1$species1, mart=human, attributesL="external_gene_name", martL=fly)
    #		write.table(listHumanFly[,2],file=paste(outputDir,outputFile,"_AllGenes_fly.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    #		write.table(modulesData$GI,file=paste(outputDir,outputFile,"_AllGenes_human.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  }
  geneInfo.BLSA$SymbolOriginal<-geneInfo.BLSA$Symbol
  geneInfo.BLSA$Symbol<-genelist.clean[match(geneInfo.BLSA$SymbolOriginal,genelist.clean[,1]),2]
}
background.BLSA=as.character(unique(geneInfo.BLSA$Symbol))
background.Emory=as.character(unique(geneInfo.Emory$Symbol)) #[geneInfo.Emory$Symbol!=''])) #Ensembl.Gene.ID))
###***



#\/----------------------------REPEAT FOR BOTH FDR AND P-VALUE CALCULATIONS--------------------------------------------\/
for (test in c("FDR","p")) {
  
  
  if (orderByRelatedness) {
    if (length(which(gsub("ME","",colnames(net1$MEs))=="grey")) >0) { net1.colorOrder=gsub("ME","",colnames(net1$MEs))[-which(gsub("ME","",colnames(net1$MEs))=="grey")] } else { net1.colorOrder=gsub("ME","",colnames(net1$MEs)) }
    uniquemodcolors.BLSA=net1.colorOrder
    if (length(which(gsub("ME","",colnames(net2$MEs))=="grey")) >0) { net2.colorOrder=gsub("ME","",colnames(net2$MEs))[-which(gsub("ME","",colnames(net2$MEs))=="grey")] } else { net2.colorOrder=gsub("ME","",colnames(net2$MEs)) }
    uniquemodcolors.Emory=net2.colorOrder
    
    nModules<-max(length(net1.colorOrder),length(net2.colorOrder))
    #modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
    orderedModules<-data.frame(Mnum=seq(1:nModules),Color=labels2colors(c(1:nModules)))
    #modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
    
    Emory.Order<-as.numeric(orderedModules$Mnum[match(uniquemodcolors.Emory,orderedModules$Color)])
    BLSA.Order<-as.numeric(orderedModules$Mnum[match(uniquemodcolors.BLSA,orderedModules$Color)])
  } else {
    net1.colorCount<-length(unique(net1$colors)) - length(which(unique(net1$colors)=="grey")) #-1 remove GREY#*** In rare case where net$colors has no grey, should be -0.
    net2.colorCount<-length(unique(net2$colors)) - length(which(unique(net2$colors)=="grey"))
    
    Emory.Order <- c(1:net2.colorCount)
    BLSA.Order <- c(1:net1.colorCount) 
    #Emory.Order <- c(21,3,4,19,11,9,20,10,18,12,2,6,1,15,7,5,13,8,14,16,17,22,23) #You can manually define order
    
    uniquemodcolors.BLSA=labels2colors(BLSA.Order)
    uniquemodcolors.Emory=labels2colors(Emory.Order)
  }
  
  ORmat=matrix(NA,nrow=length(uniquemodcolors.Emory),ncol=1)
  Pmat=matrix(NA,nrow=length(uniquemodcolors.Emory),ncol=1)
  overlapList=matrix(NA,ncol=length(uniquemodcolors.BLSA),nrow=length(uniquemodcolors.Emory))
  
  for (i in 1:length(uniquemodcolors.BLSA)){
    thismod= uniquemodcolors.BLSA[i]
    thisGene= geneInfo.BLSA$Symbol[net1$colors==thismod] #[net$colors==thismod,"Symbol"]
    testpath <- as.character(unique(thisGene)) ## Module Genes
    oraMat1=matrix(NA,ncol=9,nrow=length(uniquemodcolors.Emory))
    
    for(j in 1:length(uniquemodcolors.Emory)){
      thismod1= uniquemodcolors.Emory[j]
      thisGene1= geneInfo.Emory$Symbol[net2$colors==thismod1]
      refpath <- as.character(unique(thisGene1)) ## Module Genes
      
      testbackground <- background.BLSA
      refbackground <- background.Emory
      
      oraout=ORA(testpath,refpath,testbackground,refbackground)
      oraMat1[j,]=oraout
      
      ## GENERATE LIST OF ACTUAL OVERLAPPING SYMBOLS
      overlapList1 <- intersect(refpath,testpath)
      if (length(as.character(unlist(overlapList1)))>0) { overlapList[j,i] <- paste(unlist(overlapList1),collapse=",") }
    }
    
    ORmat=cbind(ORmat,as.numeric(oraMat1[,1]))
    Pmat=cbind(Pmat,as.numeric(oraMat1[,2]))
  }
  
  ORmat.Array =ORmat[,-1]
  Pmat.Array =Pmat[,-1]
  FDRmat.Array <- matrix(p.adjust(Pmat.Array,method="BH"),nrow=nrow(Pmat.Array),ncol=ncol(Pmat.Array))
  
  colnames(overlapList) <- colnames(ORmat.Array) <- colnames(Pmat.Array) <- colnames(FDRmat.Array) <- paste(HeatmapNet1moduleType,uniquemodcolors.BLSA,sep='.')
  rownames(overlapList) <- rownames(ORmat.Array) <- rownames(Pmat.Array) <- rownames(FDRmat.Array) <- paste(HeatmapNet2moduleType,uniquemodcolors.Emory,sep='.')
  
  ###########
  if (test=="p") {
    ##P-value-based
    dispMat <- -log10(Pmat.Array)*sign(log2(ORmat.Array))
    FDRmat.Array<-Pmat.Array
    testtext="-log10(P Value)"
    testfileext="pval"
  } else {
    ##FDR-based: 
    outputData=rbind("FET pValue", Pmat.Array, "FDR (BH) corrected", FDRmat.Array, "Overlap (Symbols)", overlapList)
    write.csv(outputData,file=paste0(rootdir,"ORA-",oneOrTwoTailed,"tailed-",outfiletitle,".FULL.csv"))
    
    dispMat <- -log10(FDRmat.Array)*sign(log2(ORmat.Array)) ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios
    testtext="-log10(FDR, BH)"
    testfileext="FDR"
  }
  ###########
  
  
  ## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g. * 
  txtMat <- dispMat #ORmat.Array
  txtMat[FDRmat.Array>=0.05] <- ""
  txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
  txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
  txtMat[FDRmat.Array <0.005] <- "***"
  
  txtMat1 <- signif(dispMat,2) #ORmat.Array
  txtMat1[txtMat1<1.5] <- ""
  
  
  textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '');
  textMatrix1= matrix(textMatrix1,ncol=ncol(Pmat.Array),nrow=nrow(Pmat.Array))
  
  #if (rownames(net2)[100]==rownames(speakeasyModules)[100] && rownames(net2)[200]==rownames(speakeasyModules)[200]) { #check if net2 is speakeasyModules; if so, reassign the SE module numbers to the modules!
  #	ordinalConversion<-unique(speakeasyModules[-nrow(speakeasyModules),c(3:4)])
  #	Emory.Order2<-as.vector(ordinalConversion[Emory.Order,1])
  #} else { Emory.Order2<-Emory.Order }
  Emory.Order2<-Emory.Order
  
  bw<-colorRampPalette(c("#0058CC", "white"))
  wr<-colorRampPalette(c("white", "#CC3300"))
  
  colvec<-if(oneOrTwoTailed==1) { wr(100) } else { c(bw(50),wr(50)) }
  zvec=if(oneOrTwoTailed==1) { c(0,10) } else { c(-10,10) }
  
  pdf(paste(rootdir,"/ORA-",oneOrTwoTailed,"tailed-",outfiletitle,".FULL.",testfileext,"-021221.pdf",sep=""), width=16,height=8)
  par( mar = c(8, 12, 3, 3) );
  par(mfrow=c(1,1))
  
  labeledHeatmap(Matrix=dispMat,
                 colorLabels=TRUE,
                 setStdMargins = FALSE,
                 yLabels= paste("ME",uniquemodcolors.Emory,sep=""),
                 ySymbols=as.vector(as.character(paste(HeatmapNet2moduleType,"-M",Emory.Order2,sep=""))),
                 xLabels= paste("ME",uniquemodcolors.BLSA,sep=""),
                 xSymbols=as.vector(as.character(paste(HeatmapNet1moduleType,"-M",BLSA.Order,sep=""))),
                 colors=colvec,
                 textMatrix = textMatrix1,
                 cex.text=0.55,
                 cex.lab.x=1,
                 zlim=zvec,
                 main=paste0(heatmapTitle," ORA/FET Signed ",testtext))
  
  dev.off()
  
  ###########################REORDER MATRIX####################################
  #find significantly correlated modules and find the order to subset them in a diagonal-positive overlap pattern, then repeat above code
  
  orderVecEmory=matrix(NA,ncol=1,nrow=1) #length(uniquemodcolors.Emory))
  orderVecBLSA=matrix(NA,ncol=1,nrow=1) #length(uniquemodcolors.Emory))
  
  BetterThanCutoff=1 #(1, include all modules in reordered PDFs; 0.5, only keep modules with a p-significance <0.5...)
  
  #transposed - works well
  for (j in 1:length(uniquemodcolors.BLSA)) {
    bestMatch= which(dispMat[,j]==max(dispMat[,j])[1])[1]  #*** outer [1] added to insure singular value returned.
    if (dispMat[bestMatch,j] >BetterThanCutoff) {
      orderVecBLSA=rbind(orderVecBLSA,j)
      truncatedRow=dispMat[-bestMatch,j]
      if (is.na(orderVecEmory [1,1])) { orderVecEmory=matrix(bestMatch,ncol=1)
      } else {
        if (is.na(match(bestMatch,orderVecEmory[,1]))) {
          orderVecEmory=rbind(as.matrix(orderVecEmory,ncol=1),bestMatch)
        }
      }
      bestMatch1= which(dispMat[,j]==max(truncatedRow)[1])[1]  #*** outer [1] added to insure singular value returned.
      truncatedRow1=truncatedRow[-match(max(truncatedRow),truncatedRow)]
      if (is.na(match(bestMatch1,orderVecEmory)) & abs(dispMat[bestMatch1,j]) >BetterThanCutoff) {
        orderVecEmory=rbind(as.matrix(orderVecEmory,ncol=1),bestMatch1)
        bestMatch2= which(dispMat[,j]==max(truncatedRow1)[1])
        truncatedRow2=truncatedRow1[-which(truncatedRow1==bestMatch2[1])]
        if (is.na(match(bestMatch2[1],orderVecEmory)) & abs(dispMat[bestMatch2[1],j]) >BetterThanCutoff) {  # *** added [1] for case where multiple e.g. 0 matches to max=0 exist.
          orderVecEmory=rbind(as.matrix(orderVecEmory,ncol=1),bestMatch2)
        }
      }
    }
  }
  
  if (!dim(orderVecBLSA)[1]==1) { #handle case where there is ZERO overlap
    #*** *** *** *** ***
    Emory.ReOrder<-as.vector(orderVecEmory)	#*** as.vector(orderVecEmory) may have worked when in sequential, not relatedness, order
    #Defined by convoluted code above ##
    BLSA.ReOrder<-as.vector(orderVecBLSA[2:nrow(orderVecBLSA),1])  #*** as.vector(orderVecBLSA[2:nrow(orderVecBLSA),1]) may have worked when in sequential, not relatedness, order
    #^Can be manually defined (as below) to include all, even completely unmatched, modules
    
    
    
    ## Add back 'BLSA' (x-axis) removed modules below threshold at right side of heatmap (lacking red/heat)
    BLSA.ReOrder<-c(BLSA.ReOrder,c(1:(length(table(net1$colors)) -length(which(unique(net1$colors)=="grey"))))[which(!c(1:(length(table(net1$colors)) -length(which(unique(net1$colors)=="grey")))) %in% BLSA.ReOrder)])
    #***(new)
    ## Add back y-axis removed modules below threshold at bottom of heatmap (lacking red/heat)
    Emory.ReOrder<-c(Emory.ReOrder,c(1:(length(table(net2$colors)) -length(which(unique(net2$colors)=="grey"))))[which(!c(1:(length(table(net2$colors)) -length(which(unique(net2$colors)=="grey")))) %in% Emory.ReOrder)])
    
    
    
    ##02/12/21 Further Apply MANUAL reorder of x (columns) using variable set in settings
    if(length(columnReorderOfAutoReorderedOutput)==length(BLSA.ReOrder)) BLSA.ReOrder<-as.numeric(as.character(BLSA.ReOrder)[columnReorderOfAutoReorderedOutput])
    #BLSA.ReOrder
    
    
    
    uniquemodcolors.ordered.Emory=labels2colors(Emory.Order[Emory.ReOrder])
    uniquemodcolors.ordered.BLSA=labels2colors(BLSA.Order[BLSA.ReOrder])
    
    dispMat.ordered<-dispMat[Emory.ReOrder,BLSA.ReOrder]
    FDRmat.Array.ordered<-FDRmat.Array[Emory.ReOrder,BLSA.ReOrder]
    
    txtMat <- dispMat.ordered
    txtMat[FDRmat.Array.ordered>=0.05] <- ""
    txtMat[FDRmat.Array.ordered <0.05&FDRmat.Array.ordered >0.01] <- "*"
    txtMat[FDRmat.Array.ordered <0.01&FDRmat.Array.ordered >0.005] <- "**"
    txtMat[FDRmat.Array.ordered <0.005] <- "***"
    
    txtMat1 <- signif(dispMat.ordered,2)
    txtMat1[txtMat1<1.5] <- ""
    
    textMatrix.ordered = paste( txtMat1, '\n', txtMat , sep = '');
    textMatrix.ordered= matrix(textMatrix.ordered,ncol=length(BLSA.ReOrder),nrow=length(Emory.ReOrder))
    
    #Transpose back to original rows and columns, but rows still reordered and culled
    dispMat.ordered.untransposed<-t(dispMat.ordered)
    textMatrix.ordered.untransposed<-t(textMatrix.ordered)
    
    #if (rownames(net2)[100]==rownames(speakeasyModules)[100] && rownames(net2)[200]==rownames(speakeasyModules)[200]) { #check if net2 is speakeasyModules; if so, reassign the SE module numbers to the modules!
    #	ordinalConversion<-unique(speakeasyModules[-nrow(speakeasyModules),c(3:4)])
    #	Emory.ReOrder2<-as.vector(ordinalConversion[Emory.ReOrder,1])
    #} else { Emory.ReOrder2<-Emory.ReOrder }
    Emory.ReOrder2<-Emory.ReOrder
    
    pdf(paste(rootdir,"/ORA-",oneOrTwoTailed,"tailed-",outfiletitle,".reordered.","GTcutoff",BetterThanCutoff,".",testfileext,"-021221.pdf",sep=""), width=16,height=8)
    par( mar = c(8, 12, 3, 3) );
    par(mfrow=c(1,1))
    
    labeledHeatmap(Matrix=t(dispMat.ordered.untransposed),
                   colorLabels=TRUE,
                   setStdMargins = FALSE,
                   yLabels= paste("ME",uniquemodcolors.ordered.Emory,sep=""),
                   ySymbols=as.vector(as.character(paste(HeatmapNet2moduleType,"-M",Emory.Order[Emory.ReOrder2],sep=""))),  #***
                   xLabels= paste("ME",uniquemodcolors.ordered.BLSA,sep=""),
                   xSymbols=as.vector(as.character(paste(HeatmapNet1moduleType,"-M",BLSA.Order[BLSA.ReOrder],sep=""))),     #***
                   colors=colvec,
                   textMatrix = t(textMatrix.ordered.untransposed),
                   cex.text=0.55,
                   cex.lab.x=1,
                   zlim=zvec,
                   main=paste0(heatmapTitle," ORA/FET Signed ",testtext))
    
    
    labeledHeatmap(Matrix=t(dispMat.ordered),
                   colorLabels=TRUE,
                   setStdMargins = FALSE,
                   xLabels= paste("ME",uniquemodcolors.ordered.Emory,sep=""),
                   xSymbols=as.vector(as.character(paste(HeatmapNet2moduleType,"-M",Emory.Order[Emory.ReOrder2],sep=""))),  #***
                   yLabels= paste("ME",uniquemodcolors.ordered.BLSA,sep=""),
                   ySymbols=as.vector(as.character(paste(HeatmapNet1moduleType,"-M",BLSA.Order[BLSA.ReOrder],sep=""))),     #***
                   colors=colvec,
                   textMatrix = t(textMatrix.ordered),
                   cex.text=0.55,
                   cex.lab.x=1,
                   zlim=zvec,
                   main=paste0(heatmapTitle," (transposed) ORA/FET Signed ",testtext))
    
    dev.off()
    
  } 		#*** *** *** *** *** ZERO Overlap case handled
  #/\----------------------------REPEATS FOR BOTH FDR AND P-VALUE CALCULATIONS-------(all warnings follow)----------------/\
}
#########################################################################################################

# Part
# Addressing generalizabilty comments #
setwd("/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Generalizabililty-SignatureInOtherDataset/")

cptac_fpkm <- read.table(file = '/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Generalizabililty-SignatureInOtherDataset/CPTAC, GDC (LUAD) cbioportal/data_mrna_seq_fpkm.txt', header = TRUE, row.names = 1) # genes=rows, sample=cols
og.cptac_fpkm <- cptac_fpkm
cptac_clinical <- read.table(file = '/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Generalizabililty-SignatureInOtherDataset/CPTAC, GDC (LUAD) cbioportal/data_clinical_patient-CDY_edited.txt', sep = "\t", header = TRUE, row.names = 1, quote = "\"")
og.cptac_clinical <- cptac_clinical

# Editing row & colnames of CPTAC files
#Do you need to split your ENSG/ENST from .# (version number)? If yes, keep everything before the "."
if (length(which(grepl("\\.",rownames(cptac_fpkm)[1])))>0) rownames(cptac_fpkm)<- as.data.frame(do.call("rbind",strsplit(rownames(cptac_fpkm),"[.]")))[,1]

# Find the Second Period and Keep Text Before It
# Keep only the first 9 characters of each column name
colnames(cptac_fpkm) <- substring(colnames(cptac_fpkm), 1, 9)

# Replace all "-" with "." in the rownames of cptac_clinical
rownames(cptac_clinical) <- gsub("-", ".", rownames(cptac_clinical))

library(biomaRt)
#First connect to ensembl database/mart
mart=useMart('ensembl')

#Find out what Datasets/species are available
listDatasets(mart)

#Find human mart
listDatasets(mart)[which(grepl("sapiens",listDatasets(mart)))]

#set up connection to ensembl mart for human
human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

tableSymbols<-getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters="entrezgene_id", values=rownames(cptac_fpkm),mart=human)

#sanity checks
dim(cptac_fpkm)
dim(tableSymbols)
#put symbols in cptac_fpkm order
symbols<-tableSymbols$hgnc_symbol[match(rownames(cptac_fpkm),tableSymbols$entrezgene_id)]
length(symbols)
tail(symbols,100)

#replace NA with empty string where no lookup was possible.
symbols[is.na(symbols)]<-""
tail(symbols,100)

uniqueIDs<-paste0(symbols,"|",rownames(cptac_fpkm))
tail(uniqueIDs)
rownames(cptac_fpkm)<-uniqueIDs

# Find the intersection of row names in cptac_clinical and column names in cptac_fpkm
matching_names <- intersect(rownames(cptac_clinical), colnames(cptac_fpkm))

# Filter cptac_clinical and cptac_fpkm based on matching names
cptac_clinical <- cptac_clinical[matching_names, , drop = FALSE]
cptac_fpkm <- cptac_fpkm[, matching_names, drop = FALSE]

# Check the dimensions to ensure both are 175
print(dim(cptac_clinical))  #175 rows
print(dim(cptac_fpkm))      #175 columns

rownames(cptac_clinical)==colnames(cptac_fpkm)
write.csv(cptac_fpkm,file=paste0("cptac_fpkm.symbolsBiomaRt.csv"))

# just follow "FPKM_4signature_comparison_byCorr&ROC&Clustering-EBD.R" no need to do other stuff. Just follow what is done here and modify for different dataset, only look at 1 year across different signatures (ROC/AUC figure)
table(cptac_clinical$PATH_STAGE)
# removing those w/o stage
dim(cptac_clinical) #175,  59
cptac_clinical <- cptac_clinical[cptac_clinical$PATH_STAGE != "", ]
dim(cptac_clinical) #172, 59

matching_names <- intersect(rownames(cptac_clinical), colnames(cptac_fpkm))
# Filter cptac_clinical and cptac_fpkm based on matching names
cptac_clinical <- cptac_clinical[matching_names, , drop = FALSE]
cptac_fpkm <- cptac_fpkm[, matching_names, drop = FALSE]

# Replace zero values with NA to prepare for log2 transformation
cptac_fpkm[cptac_fpkm == 0] <- NA

# Remove rows with more than 50% missing values
cleanDat.removeRows.50percent.idx <- which(apply(cptac_fpkm, 1, function(x) length(which(is.na(x)))) > ncol(cptac_fpkm) / 2)
cleanDat.upTo50percentNA <- cptac_fpkm[-cleanDat.removeRows.50percent.idx, ]
cleanDat.upTo50percentNA <- log2(cleanDat.upTo50percentNA)  # Log-transform

# Remove rows with more than 90% missing values
cleanDat.removeRows.90percent.idx <- which(apply(cptac_fpkm, 1, function(x) length(which(is.na(x)))) > ncol(cptac_fpkm) * 0.90)
cleanDat.upTo90percentNA <- cptac_fpkm[-cleanDat.removeRows.90percent.idx, ]
cleanDat.upTo90percentNA <- log2(cleanDat.upTo90percentNA)  # Log-transform

# Check dimensions
print(dim(cleanDat.upTo50percentNA))
print(dim(cleanDat.upTo90percentNA))

# Replace missing values with a conservative value for log2-transformed data
cleanDat.upTo90percentNA[is.na(cleanDat.upTo90percentNA)] <- min(cleanDat.upTo90percentNA, na.rm = TRUE) - 1
#cleanDat.upTo50percentNA[is.na(cleanDat.upTo50percentNA)] <- min(cleanDat.upTo50percentNA, na.rm = TRUE) - 1

# Cluster 3.0 is equivalent to:
# Normalization step 1: Median subtraction for each sample column
norm.upTo90percentNA <- apply(cleanDat.upTo90percentNA, 2, function(x) x - median(x, na.rm = TRUE))
# Normalization step 2: Scale each column so the sum of the squares is 1
norm.scaled.upTo90percentNA <- apply(norm.upTo90percentNA, 2, function(x) x / sqrt(sum(x^2, na.rm = TRUE)))

#tumor.exp.unlog<-2^(norm.scaled.upTo90percentNA)  # returns all values positive, centering total population near 1 (may be multimodal/multi-peak population histogram with noise peaks to left)
#range(tumor.exp.unlog - tumor.exp.scaled.unlog)

# Create PDF to visualize histograms
pdf("Histograms-Normalized-Scaled-FPKM.pdf", width = 16, height = 7)
par(mfrow = c(1, 2))  # Side-by-side plots

# Histogram for up to 90% missing data
h<-hist(norm.scaled.upTo90percentNA, breaks = 250, main = "Up to 90% missing log2(FPKM) values", xlab = "Expression Value", col = "lightblue")
max(h$counts[40:326]) #69399
which(h$counts==69399) #159
h$mids[159] #0.00245
h$mids[162] #0.00275
abline(v=0.00275,col="red")

nonNoisePeak.shift <- 0.00275
#norm.tumor.exp.scaled.recentered<-2^(norm.scaled.upTo90percentNA - nonNoisePeak.shift) 
norm.tumor.exp.scaled.recentered<-2^(norm.scaled.upTo90percentNA) 

hist(norm.tumor.exp.scaled.recentered, breaks = 250, main = "Up to 90% missing unlogged shifted log2(FPKM) values", xlab = "Expression Value", col = "lightblue") #this is a histogram of an distribution of an unlogged non-noise peak center at 1

# Histogram for up to 50% missing data
# hist(norm.scaled.upTo50percentNA, breaks = 250, main = "Up to 50% missing log2(FPKM) values", xlab = "Expression Value", col = "lightgreen")

dev.off()

# Calculate Prognostic/predictor ratios based on signature gene lists from 4 studies 
# The data from which ratios are calculated
dim(norm.tumor.exp.scaled.recentered) #35589, 172
str(norm.tumor.exp.scaled.recentered)

# Sample traits
numericMeta <- cptac_clinical
dim(numericMeta) #172, 59

#Curated LUAD Signature Gene Lists (as R list)
hiExpr.hiSurv<-loExpr.hiSurv<-list()

#Song et al, Sci Rep (2022) Per  https://www.nature.com/articles/s41598-022-14323-6.pdf  #Hazard ratios of 11 genes in  IRG signature (Figure 2, panel D)
loExpr.hiSurv[["Song"]]<-c('GPC3','IL7R','NMUR1','PTPRE','SEMA4D')
hiExpr.hiSurv[["Song"]]<-c('ADM','NMI','PSEN1','PVR','SERPINE1','SHK1')
#Paste the following signature into https://xenabrowser.net/heatmap/
# =ADM+NMI+PSEN1+PVR+SERPINE1+SHK1-ADM-NMI-PSEN1-PVR-SERPINE1-SHK1

#Shedden et al, Nat Med (2009) per https://www.nature.com/articles/nm.1790
#459 genes covering 545 affymetrix array probes curated from http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6.html
loExpr.hiSurv[["Shedden"]]<-c('RFC2','PSMB2','SMNDC1','DSP','AP2B1','CALM3','WARS1','PGK1','RAN','STMN1','SCD','H2AZ1','DNAJA1','FKBP4','PPM1G','ACTR3','UBE2V1','PAICS','EIF5B','LDHB','PFKP','CSE1L','GNAI3','SLC7A5','PCNA','SLC2A1','PSMA5','ADRM1','TOP2A','MOB1A','EIF4A3','YME1L1','GSS','BUB3','RUVBL2','PPIF','SRM','UBE2N','SAR1A','NME1','DARS1','ACP1','PNP','MYBL2','TUBG1','TWF1','MCM5','MTHFD2','NCAPD2','HDAC2','RRM2','PSRC1','MCM6','BIRC5','UBAC1','NUP62','SLC16A1','TOMM40','TK1','UBE2K','PSMD12','PCLAF','ATG5','GMFB','FOXM1','TYMS','PLOD2','NRAS','RNMT','CCNB2','CAD','LSM4','R3HDM1','UBE2S','HEXIM1','SLC16A3','CDC20','NBN','ADM','HK2','UBE2C','PYGL','SSX2IP','KIF2A','MTIF2','CDYL','DNM1L','HAT1','SPAG5','KRR1','CDK1','THOP1','DTYMK','LMNB1','DCK','MTF2','EZH2','MAD2L1','CCNA2','STC2','MRPL19','CIAO1','PTTG1','GGH','ARL4D','RTCA','CPT1A','HCCS','BUB1B','DLGAP5','FANCA','DGUOK','IGF2BP3','MRPL12','CDC6','PAWR','RFC4','PDCD2','ZWINT','TRIP13','AURKA','CDC45','RFC3','RAD51AP1','CDKN2C','NDC80','WASF1','CKS2','PPID','GCH1','SMC2','DBF4','ELOVL6','TEAD4','PMAIP1','GTSE1','TRMU','SEC23A','AK4','TTF2','KIF11','CDC7','FAM216A','RAD54L','MMP12','STC1','EXO1','ACD','NEK2','TROAP','SRD5A1','AFP','CDC25A','KIF23','WDHD1','TTK','MELK','PLK4','SAP30','CENPA','PFN2','CDK5R1','CCNE2','CENPE','PRIM1','ENPP1','ORC1','RECQL','TMEFF1','VRK2','CDC25C','TIMM8A','GPSM2','RBL1','CHEK1','H2AX','WDR76','ANGPT2','SNX7','CEP43','PRIM2','F12','DUSP9','NOLC1','POLE2','FOXG1','NMU','SLBP','SNRPA1','HMGA1','GINS1','GABPB1','MPHOSPH9','KNTC1','IL2RA','KIF14','PRMT1','CD70','CHRNA5','APOBEC3B','POLR3G','KLRC1','SEMA3A','LDHC','ERN1','IL18RAP','GLMN','HMMR','GPR19','TUBA4B','CENPI','DR1','PRKAA2','PKP2','CENPF','CSNK1D','ACOT7','PPARD','AP2S1','CCL7','MYO7A','UGT8','BRCA2','ESM1','PRPS1','PITX1','PTTG3P','H2BC9','PA2G4','CCT5','TCP1','HMGB2','SNRPB','PCMT1','SEPHS1','NDUFA9','XPNPEP1','NSD2','UBE2V2','PLIN2','HBS1L','EIF4E2','KIF2C','PPAT','PKM','AURKB','NCBP1','HDGFL3','EXOSC2','RIPK2','GPR37','BUB1','CDKN2A','CDKN3','CBX5','ITCH','CDT1','PSME3','SNCG','SPC25','SS18','GPD2','PGM3','TPX2','HAUS3','CTSV','RGS20','GCFC2','CORT','IL1RAP','PPFIA1','ARTN','CDKN2D','SMG1P5','TUBA3C','PRKDC','SLC9A2','NAB1','MSH6','LRPPRC','KPNA2','GINS4','IPO5','SZRD1','MKI67','MCM4','VEGFA','ANKLE2','TUBA4A','WAPL','NUP210','TFDP1','NEMP1','PVR','NCAPD3','CBS','MAPKAPK5','KLC1','NCAPH','FANCI','DNAJC9','CNOT9','TNNT1','ZNF248','TEX30','COL2A1','RECQL4','CCNE1','OIP5','DNA2','FCHO1','MED8','MYBL1','PAX8','FUT9','TBC1D31','GNAS','RALA','ARPP19','PTTG2','UBE2D1','MAGEA3','WDR43','PTBP3','RIF1','CCNB1','WDR62','PHTF2','ATP2B1','ARFGEF2','PSMA7','XRCC3','TUBA3E','KRT8P11','SOD2','PPIAP21','KRT8P17','EIF3JP1','OR7E12P','STEAP1B','SKA1','NAA50','JPT1','MTCH2','YKT6','SNX6','SYNCRIP','TFG','MIF','CDC27','STRN4','MRPL42','EIPR1','PRC1','NUSAP1','MRPL13','AGPAT5','CMAS','TSR1','GOLT1B','SLC38A1','CKAP2','DDA1','SLC25A10','TACC3','UBA6','ZWILCH','COMMD8','KIF4A','RAB22A','CARHSP1','SFXN1','CMC2','KCTD5','ERO1A','BPNT2','CEP55','DTL','MRGBP','NUP37','NCAPG','GINS3','HJURP','CENPM','KIF20A','TMEM38B','SMC6','ATAD2','PPP2R3C','FBXO5','CENPU','NETO2','MAGOHB','NABP2','YEATS4','CCDC59','MTPAP','PNKP','TFPT','RASAL2','TPRKB','NTAQ1','MEMO1','PHF10','RPP25','PBK','TWSG1','FBXO11','EGLN3','CENPQ','KIF15','VANGL1','LRFN4','SHCBP1','NEIL3','DSN1','ELOVL4','CENPN','NCAPG2','ERCC6L','QSER1','PRR7','DENND1A','CCDC93','ZSCAN5A','ASPM','E2F8','PARPBP','HELLS','MFSD12','MUC16','KMT5AP1','RC3H2','LRIF1','DEPDC1','TAF7L','POMT2','MCM10','TBRG4','PSAT1','KAZALD1','FIP1L1','YEATS2','KIF18A','CDCA3','CDCA8','GINS2','EIF4EBP1','NUDT3','PIMREG','AGMAT','INTS13','BRIP1','MAP6D1','C6orf120','CLPB','KIF18B','RACGAP1','GAPDH')
hiExpr.hiSurv[["Shedden"]]<-vector() #c()
#Paste the following signature into https://xenabrowser.net/heatmap/
# =RFC2+PSMB2+SMNDC1+DSP+AP2B1+CALM3+WARS1+PGK1+RAN+STMN1+SCD+H2AZ1+DNAJA1+FKBP4+PPM1G+ACTR3+UBE2V1+PAICS+EIF5B+LDHB+PFKP+CSE1L+GNAI3+SLC7A5+PCNA+SLC2A1+PSMA5+ADRM1+TOP2A+MOB1A+EIF4A3+YME1L1+GSS+BUB3+RUVBL2+PPIF+SRM+UBE2N+SAR1A+NME1+DARS1+ACP1+PNP+MYBL2+TUBG1+TWF1+MCM5+MTHFD2+NCAPD2+HDAC2+RRM2+PSRC1+MCM6+BIRC5+UBAC1+NUP62+SLC16A1+TOMM40+TK1+UBE2K+PSMD12+PCLAF+ATG5+GMFB+FOXM1+TYMS+PLOD2+NRAS+RNMT+CCNB2+CAD+LSM4+R3HDM1+UBE2S+HEXIM1+SLC16A3+CDC20+NBN+ADM+HK2+UBE2C+PYGL+SSX2IP+KIF2A+MTIF2+CDYL+DNM1L+HAT1+SPAG5+KRR1+CDK1+THOP1+DTYMK+LMNB1+DCK+MTF2+EZH2+MAD2L1+CCNA2+STC2+MRPL19+CIAO1+PTTG1+GGH+ARL4D+RTCA+CPT1A+HCCS+BUB1B+DLGAP5+FANCA+DGUOK+IGF2BP3+MRPL12+CDC6+PAWR+RFC4+PDCD2+ZWINT+TRIP13+AURKA+CDC45+RFC3+RAD51AP1+CDKN2C+NDC80+WASF1+CKS2+PPID+GCH1+SMC2+DBF4+ELOVL6+TEAD4+PMAIP1+GTSE1+TRMU+SEC23A+AK4+TTF2+KIF11+CDC7+FAM216A+RAD54L+MMP12+STC1+EXO1+ACD+NEK2+TROAP+SRD5A1+AFP+CDC25A+KIF23+WDHD1+TTK+MELK+PLK4+SAP30+CENPA+PFN2+CDK5R1+CCNE2+CENPE+PRIM1+ENPP1+ORC1+RECQL+TMEFF1+VRK2+CDC25C+TIMM8A+GPSM2+RBL1+CHEK1+H2AX+WDR76+ANGPT2+SNX7+CEP43+PRIM2+F12+DUSP9+NOLC1+POLE2+FOXG1+NMU+SLBP+SNRPA1+HMGA1+GINS1+GABPB1+MPHOSPH9+KNTC1+IL2RA+KIF14+PRMT1+CD70+CHRNA5+APOBEC3B+POLR3G+KLRC1+SEMA3A+LDHC+ERN1+IL18RAP+GLMN+HMMR+GPR19+TUBA4B+CENPI+DR1+PRKAA2+PKP2+CENPF+CSNK1D+ACOT7+PPARD+AP2S1+CCL7+MYO7A+UGT8+BRCA2+ESM1+PRPS1+PITX1+PTTG3P+H2BC9+PA2G4+CCT5+TCP1+HMGB2+SNRPB+PCMT1+SEPHS1+NDUFA9+XPNPEP1+NSD2+UBE2V2+PLIN2+HBS1L+EIF4E2+KIF2C+PPAT+PKM+AURKB+NCBP1+HDGFL3+EXOSC2+RIPK2+GPR37+BUB1+CDKN2A+CDKN3+CBX5+ITCH+CDT1+PSME3+SNCG+SPC25+SS18+GPD2+PGM3+TPX2+HAUS3+CTSV+RGS20+GCFC2+CORT+IL1RAP+PPFIA1+ARTN+CDKN2D+SMG1P5+TUBA3C+PRKDC+SLC9A2+NAB1+MSH6+LRPPRC+KPNA2+GINS4+IPO5+SZRD1+MKI67+MCM4+VEGFA+ANKLE2+TUBA4A+WAPL+NUP210+TFDP1+NEMP1+PVR+NCAPD3+CBS+MAPKAPK5+KLC1+NCAPH+FANCI+DNAJC9+CNOT9+TNNT1+ZNF248+TEX30+COL2A1+RECQL4+CCNE1+OIP5+DNA2+FCHO1+MED8+MYBL1+PAX8+FUT9+TBC1D31+GNAS+RALA+ARPP19+PTTG2+UBE2D1+MAGEA3+WDR43+PTBP3+RIF1+CCNB1+WDR62+PHTF2+ATP2B1+ARFGEF2+PSMA7+XRCC3+TUBA3E+KRT8P11+SOD2+PPIAP21+KRT8P17+EIF3JP1+OR7E12P+STEAP1B+SKA1+NAA50+JPT1+MTCH2+YKT6+SNX6+SYNCRIP+TFG+MIF+CDC27+STRN4+MRPL42+EIPR1+PRC1+NUSAP1+MRPL13+AGPAT5+CMAS+TSR1+GOLT1B+SLC38A1+CKAP2+DDA1+SLC25A10+TACC3+UBA6+ZWILCH+COMMD8+KIF4A+RAB22A+CARHSP1+SFXN1+CMC2+KCTD5+ERO1A+BPNT2+CEP55+DTL+MRGBP+NUP37+NCAPG+GINS3+HJURP+CENPM+KIF20A+TMEM38B+SMC6+ATAD2+PPP2R3C+FBXO5+CENPU+NETO2+MAGOHB+NABP2+YEATS4+CCDC59+MTPAP+PNKP+TFPT+RASAL2+TPRKB+NTAQ1+MEMO1+PHF10+RPP25+PBK+TWSG1+FBXO11+EGLN3+CENPQ+KIF15+VANGL1+LRFN4+SHCBP1+NEIL3+DSN1+ELOVL4+CENPN+NCAPG2+ERCC6L+QSER1+PRR7+DENND1A+CCDC93+ZSCAN5A+ASPM+E2F8+PARPBP+HELLS+MFSD12+MUC16+KMT5AP1+RC3H2+LRIF1+DEPDC1+TAF7L+POMT2+MCM10+TBRG4+PSAT1+KAZALD1+FIP1L1+YEATS2+KIF18A+CDCA3+CDCA8+GINS2+EIF4EBP1+NUDT3+PIMREG+AGMAT+INTS13+BRIP1+MAP6D1+C6orf120+CLPB+KIF18B+RACGAP1+GAPDH

#Soltis et al, (2022) 155 RNA-based metastasis-free survival (MFS) associated gene list [from Supplemental data Table 2C] in https://www.sciencedirect.com/science/article/pii/S2666379122003780
loExpr.hiSurv[["Soltis"]]<-vector() #c()
hiExpr.hiSurv[["Soltis"]]<-c('ABCD3','ACSS1','ADAMTS4','ADAMTS5','ADCY9','ADGRF5','AHR','ALDH6A1','ANKLE2','ANKMY2','AQP4','ARSD','ATP13A4','ATPAF1','AVEN','B3GLCT','B4GAT1','BDH2','BTD','C16orf89','CA4','CAB39L','CALU','CAMK2D','CARD14','CAT','CD302','CLIP1','CLUAP1','CNPY3','COL11A1','COL3A1','COL5A2','COMMD1','COMMD10','CPE','CTSH','CYB5A','CYP24A1','CYP4B1','CYP4X1','DAPK2','ECHDC2','EGFR','EPHX2','FAAH','FAM76A','FAM83A','FBXO18','FDXR','FGA','FGG','FLNC','FMO5','FNDC1','FOSL2','FREM2','FURIN','GAPDH','GDF10','GDF15','GFPT2','GPD1L','HIRIP3','HMGCL','HSD17B8','IFIT1','IFT57','IFT80','IRS2','KIF16B','KIF26B','KSR1','KYNU','LBP','LILRB3','LZIC','MAMDC2','MATR3','MCCC1','MISP','MMP28','MTHFD2','MTX3','MUC4','MUC5B','MYO1B','MYO1E','NAMPT','NDRG2','NFIX','NNMT','NRP2','OSCP1','PDE3A','PEX3','PIGK','PLA2G4F','PLEKHG2','PLLP','POSTN','PPIF','PPP2R5A','PRSS12','PUS1','PUS10','RAB29','RAB33B','RHOF','RIPK2','RPS6KA5','SACM1L','SBNO2','SCD5','SCP2','SDS','SELENBP1','SEMA4B','SERPINB4','SERPINE2','SETMAR','SHROOM2','SLC16A3','SLC27A3','SLC29A1','SLC2A3','SLC7A5','SND1','SRPRA','ST3GAL4','STX3','SULF1','SUSD2','TAPT1','TCEAL4','THBS1','THBS2','THSD4','THUMPD1','THY1','TMEM231','TMX4','TOE1','TOR1AIP1','TPPP3','TRAF2','TUBB3','UBAC1','UPRT','UTRN','VCAN','VMP1','VPS41','WLS','ZNF219')
#Paste the following signature into https://xenabrowser.net/heatmap/
#=ABCD3+ACSS1+ADAMTS4+ADAMTS5+ADCY9+ADGRF5+AHR+ALDH6A1+ANKLE2+ANKMY2+AQP4+ARSD+ATP13A4+ATPAF1+AVEN+B3GLCT+B4GAT1+BDH2+BTD+C16orf89+CA4+CAB39L+CALU+CAMK2D+CARD14+CAT+CD302+CLIP1+CLUAP1+CNPY3+COL11A1+COL3A1+COL5A2+COMMD1+COMMD10+CPE+CTSH+CYB5A+CYP24A1+CYP4B1+CYP4X1+DAPK2+ECHDC2+EGFR+EPHX2+FAAH+FAM76A+FAM83A+FBXO18+FDXR+FGA+FGG+FLNC+FMO5+FNDC1+FOSL2+FREM2+FURIN+GAPDH+GDF10+GDF15+GFPT2+GPD1L+HIRIP3+HMGCL+HSD17B8+IFIT1+IFT57+IFT80+IRS2+KIF16B+KIF26B+KSR1+KYNU+LBP+LILRB3+LZIC+MAMDC2+MATR3+MCCC1+MISP+MMP28+MTHFD2+MTX3+MUC4+MUC5B+MYO1B+MYO1E+NAMPT+NDRG2+NFIX+NNMT+NRP2+OSCP1+PDE3A+PEX3+PIGK+PLA2G4F+PLEKHG2+PLLP+POSTN+PPIF+PPP2R5A+PRSS12+PUS1+PUS10+RAB29+RAB33B+RHOF+RIPK2+RPS6KA5+SACM1L+SBNO2+SCD5+SCP2+SDS+SELENBP1+SEMA4B+SERPINB4+SERPINE2+SETMAR+SHROOM2+SLC16A3+SLC27A3+SLC29A1+SLC2A3+SLC7A5+SND1+SRPRA+ST3GAL4+STX3+SULF1+SUSD2+TAPT1+TCEAL4+THBS1+THBS2+THSD4+THUMPD1+THY1+TMEM231+TMX4+TOE1+TOR1AIP1+TPPP3+TRAF2+TUBB3+UBAC1+UPRT+UTRN+VCAN+VMP1+VPS41+WLS+ZNF219

#Current LUAD Study ROC Optimization 8 gene list:  ATP6V0E1+SVBP+HSDL1+UBTD1 / GNPNAT1+XRCC2+TFAP2A+PPP1R13L
loExpr.hiSurv[["thisStudy"]]<-c('GNPNAT1','XRCC2','TFAP2A','PPP1R13L')
hiExpr.hiSurv[["thisStudy"]]<-c('ATP6V0E1','SVBP','HSDL1','UBTD1')
#The above 8 genes are based on AUC outcomes of ratioed combinations of genes tested in Round 3 (17 denom; 15 numerator)
#loExpr.hiSurv <- c('hsa-miR-143-3p','hsa-miR-30d-3p','hsa-miR-99a-3p','hsa-miR-497-5p','hsa-miR-1254','hsa-miR-193a-5p','hsa-miR-141-3p','hsa-miR-141-5p','hsa-miR-7704','hsa-miR-221-5p','hsa-miR-550a-5p','TFAP2A','PPP1R13L','GNPNAT1','GRK3','CITED2','CDK1')
#hiExpr.hiSurv <- c('hsa-miR-143-3p','hsa-miR-30c-2-3p','hsa-miR-29b-2-5p','hsa-miR-4709-3p','hsa-miR-3065-3p','hsa-miR-135b-5p','hsa-miR-664a-3p','hsa-miR-22-5p','hsa-miR-200c-3p','BEND5','HSDL1','ATP6V0E1','ATAT1','SVBP','GADD45GIP1')
#note: Xena TCGA portal does not have SVBP curated RNA-Seq norm_counts.
#Paste the following signature into https://xenabrowser.net/heatmap/
# =ATP6V0E1+SVBP+HSDL1+UBTD1-GNPNAT1-XRCC2-TFAP2A-PPP1R13L



geneSymbols=as.data.frame(do.call(rbind,strsplit(rownames(norm.scaled.upTo90percentNA),"[|]")))[,1]
length(geneSymbols) #35589
length(unique(geneSymbols)) #22872 (high # of "" in gene list)

# Cleaning numericMeta & creating numericMeta$OS.time & OS 
numericMeta$OS.time <- numericMeta$OS_MONTHS
#0:living, 1:dead (same as og.numericMeta)
numericMeta$OS <- ifelse(grepl("0:LIVING", numericMeta$OS_STATUS), 0,
                         ifelse(grepl("1:DECEASED", numericMeta$OS_STATUS), 1, NA))
dim(numericMeta) #172, 61
na.idx <- which(is.na(numericMeta$OS.time))
if (length(na.idx)>0) {
  numericMeta <- numericMeta[-na.idx,]
  dim(numericMeta) #154, 61
  norm.tumor.exp.scaled.recentered <- norm.tumor.exp.scaled.recentered[,match(rownames(numericMeta),colnames(norm.tumor.exp.scaled.recentered))]
}
dim(numericMeta) #154, 61
dim(norm.tumor.exp.scaled.recentered) #35589, 154

######################################
#Calculate 4 gene signature lists' LUAD equal gene weight prognostic ratios  ...for Overall Survival ROC AUC
rootdir="/Users/youngcd/Library/CloudStorage/Box-Box/Coreys_Folder/WGCNA/Final/Generalizabililty-SignatureInOtherDataset/" 


ratioList<-list()
#ratioNames<-vector()
iter=0 #placeholder for list elements in ratioList

which(!names(hiExpr.hiSurv)==names(loExpr.hiSurv))
#same list element names & order in both lists
studyNames=names(hiExpr.hiSurv)
studyNames
#[1] "Song"      "Shedden"   "Soltis"    "thisStudy"

numerator.geneCount<-denominator.geneCount<-placeholder.check.numerator<-placeholder.check.denominator<-list()
for(study in studyNames) {
  iter=iter+1
  numerator.df<-norm.tumor.exp.scaled.recentered[geneSymbols %in% hiExpr.hiSurv[[study]],]
  denominator.df<-norm.tumor.exp.scaled.recentered[geneSymbols %in% loExpr.hiSurv[[study]],]
  if(sum(numerator.df)==0) { numerator.df=as.data.frame(matrix(rep(1,ncol(norm.tumor.exp.scaled.recentered)),nrow=1)); numerator.geneCount[[study]]=1; placeholder.check.numerator[[study]]=TRUE; } else { numerator.geneCount[[study]]=length(which(geneSymbols %in% hiExpr.hiSurv[[study]])); placeholder.check.numerator[[study]]=FALSE; }
  if(sum(denominator.df)==0) { denominator.df=as.data.frame(matrix(rep(1,ncol(norm.tumor.exp.scaled.recentered)),nrow=1)); denominator.geneCount[[study]]=1; placeholder.check.denominator[[study]]=TRUE; } else { denominator.geneCount[[study]]=length(which(geneSymbols %in% loExpr.hiSurv[[study]])); placeholder.check.denominator[[study]]=FALSE; }
  names(numerator.df)<-names(denominator.df)<-colnames(norm.tumor.exp.scaled.recentered)  #in case numerator.df or denominator.df are de novo generated without sample names; the generic names ("V##") would propagate to ratioList.
  
  ratioList[[study]]<-colSums(numerator.df) / colSums(denominator.df)
  #adjust ratio for gene counts in numerator, denominator
  ratioList[[study]]<- ratioList[[study]] / numerator.geneCount[[study]] * denominator.geneCount[[study]]
}

#list element names must match studyNames
ratioDescription<-list(Song=c('Song et al, Sci Rep (2022) - Inflammatory Response Gene (IRG) Signature'),
                       Shedden=c('Shedden et al, Nat Med (2008) - Poor Survival 459 Gene Signature'),
                       Soltis=c('Soltis et al, Cell Rep Med (2022) - Metastasis Free Survival 155 Gene Signature'),
                       thisStudy=c('This Study - ROC AUC optimization 8 Gene Signature'))

ratioNames <- sapply(studyNames,function(x) paste0(ratioDescription[[x]], "\nGenes In TCGA curated data... ",if(placeholder.check.numerator[[x]]) { "0" } else { numerator.geneCount[[x]] }," for numerator; and ",if(placeholder.check.denominator[[x]]) { "0" } else { denominator.geneCount[[x]] }," for denominator"))
ratioNames
#                                                                                                                                                     Song 
#          "Song et al, Sci Rep (2022) - Inflammatory Response Gene (IRG) Signature\nGenes In TCGA curated data... 5 for numerator; and 5 for denominator" 
#                                                                                                                                                  Shedden 
#               "Shedden et al, Nat Med (2008) - Poor Survival 459 Gene Signature\nGenes In TCGA curated data... 0 for numerator; and 459 for denominator" 
#                                                                                                                                                   Soltis 
#"Soltis et al, Cell Rep Med (2022) - Metastasis Free Survival 155 Gene Signature\nGenes In TCGA curated data... 155 for numerator; and 0 for denominator" 
#                                                                                                                                                thisStudy 
#                               "This Study - ROC AUC optimization 8 Gene Signature\nGenes In TCGA curated data... 4 for numerator; and 4 for denominator

lapply(ratioList,range)
#$Song
#[1] 0.9996238 1.0027183
#
#$Shedden
#[1] 0.9994075 1.0011988
#
#$Soltis
#[1] 0.9996493 1.0007339
#
#$thisStudy
#[1] 1.000554 1.003106

#unrelated examples of other data ranges for: v7c3 final 9 v7c with MPL & v7b
#0.992103 1.0135701 with collapse method="maxRowVariance"
#0.992134 1.0123890 with Collapse method="Average"


iter # how many pages total
#4


## Begin trying ratios computed above as predictors in ROC analysis.
blockSize=4  #Number of pages max per PDF of ROC curves
filenamePrefix="ROCcurves-SignatureSet_comparison.RatiosForLuAD-01-16.2023"


    
outcomeNumber=0
for(outcome in c("OS")) {
  cat(paste0("\nPerforming ROC analysis for outcome: ",outcome,".\n"))
  
  outcomeNumber=outcomeNumber+1
  AUClist<-ROClist<-plist<-CIlist95<-ACClist<-SENlist<-SPElist<-list()
  pdf(file=paste0(rootdir,"/3e",outcomeNumber,".Predicting_",outcome,"_",filenamePrefix,"(1to",blockSize,").pdf"),width=10,height=10)
  for (page in 1:iter) {
    prognostic.ratio=ratioList[[page]]
    
    # # subgroup data for 12 mo, 18 mo, 3 year, 5 year, and 10 year survival binary trait (some cases will be censored)
    # Grouping.12moRFS<-rep(NA,length(prognostic.ratio))
    # Grouping.12moRFS[numericMeta[,paste0(outcome,".time")] >= 12] <- "Survived"
    # Grouping.12moRFS[numericMeta[,outcome]==1 & numericMeta[,paste0(outcome,".time")] <= 12] <- "Died"
    # length(na.omit(Grouping.12moRFS))
    # #?/480 casesamples,    (was 411 before removal of 4 missing data cases)
    
    Grouping.18moRFS<-rep(NA,length(prognostic.ratio))
    Grouping.18moRFS[numericMeta[,paste0(outcome,".time")] >= 18] <- "Survived"
    Grouping.18moRFS[numericMeta[,outcome]==1 & numericMeta[,paste0(outcome,".time")] <= 18] <- "Died"
    length(na.omit(Grouping.18moRFS))
    #?/480 case samples
    
    # Grouping.3yrRFS<-rep(NA,length(prognostic.ratio))
    # Grouping.3yrRFS[numericMeta[,paste0(outcome,".time")] >= 36] <- "Survived"
    # Grouping.3yrRFS[numericMeta[,outcome]==1 & numericMeta[,paste0(outcome,".time")] <= 36] <- "Died"
    # length(na.omit(Grouping.3yrRFS))
    # #?/480 case samples,
    
    # Grouping.5yrRFS<-rep(NA,length(prognostic.ratio))
    # Grouping.5yrRFS[numericMeta[,paste0(outcome,".time")] >= 365*5] <- "Survived"
    # Grouping.5yrRFS[numericMeta[,outcome]==1 & numericMeta[,paste0(outcome,".time")] <= 365*5] <- "Died"
    # length(na.omit(Grouping.5yrRFS))
    #?/480,
    
    #  Grouping.7.5yrRFS<-rep(NA,length(prognostic.ratio))
    #  Grouping.7.5yrRFS[numericMeta$OS.time >= 365.25*7.5] <- "Survived"
    #  Grouping.7.5yrRFS[numericMeta$OS==1 & numericMeta$OS.time <= 365.25*7.5] <- "Died"
    #  length(na.omit(Grouping.7.5yrRFS))
    #  #?/480 case samples
    
    # Grouping.10yrRFS<-rep(NA,length(prognostic.ratio))
    # Grouping.10yrRFS[numericMeta[,paste0(outcome,".time")] >= 365.25*10] <- "Survived"
    # Grouping.10yrRFS[numericMeta[,outcome]==1 & numericMeta[,paste0(outcome,".time")] <= 365.25*10] <- "Died"
    # length(na.omit(Grouping.10yrRFS))
    #?/480,
    
    #install.packages("pROC")
    library(pROC)
    library(verification)
    library(reportROC) # for accuracy, sensitivity, specificity reporting.
    ## Assemble list of survival data, prognostic ratio
    survivalList<- list(Grouping.18moRFS=Grouping.18moRFS)
    #survivalList<- list(Grouping.12moRFS=Grouping.12moRFS,Grouping.18moRFS=Grouping.18moRFS,Grouping.3yrRFS=Grouping.3yrRFS,Grouping.5yrRFS=Grouping.5yrRFS,Grouping.10yrRFS=Grouping.10yrRFS)
    descriptionVec<-c("18 Months") # Changed 15 to 12 months throughout
    #descriptionVec<-c("12 Months","18 Months","3 Years","5 Years","10 Years") # Changed 15 to 12 months throughout
    colvec<-c("green","magenta","red","purple","blue") # Changed to Corey's colors + kept magenta for 18 mo
    par(mar=c(0,0,3,0))
    par(oma=c(0,0,3,0))
    par(mfrow=c(1,1)) #changed from 2 columns to 1 (from 1,2)
    roc.info<-roc.info.p<-roc.info.95ci<- roc.info.accuracy<-roc.info.specificity<-roc.info.sensitivity <- vector()
    for (i in 1:length(names(survivalList))) {
      Survived1 <- ifelse(test=survivalList[[i]][-which(is.na(survivalList[[i]]))]=="Survived", yes=1, no=0)
      # predictor.log2 <- log2(#peptides[,names(bestTpVals)[i]])
      # predictor.log2[!is.finite(predictor.log2)] <- min(predictor.log2[is.finite(predictor.log2)])-1
      # predictor <- predictor.log2 -mean(predictor.log2,na.rm=TRUE)
      predictor <- unlist(prognostic.ratio[which(!is.na(survivalList[[i]]))]) # changed from -which( to which(! ; also added unlist to force vector
      reorder<-order(predictor,decreasing=FALSE)
      Survived<-Survived1[reorder]
      predictor<-predictor[reorder]
      # plot(x=predictor, y=Survived)
      glm.fit=glm(Survived~predictor, family=binomial)
      # lines(predictor, glm.fit$fitted.values)
      # roc(Survived, glm.fit$fitted.values, plot=TRUE)
      par(pty="s")
      if (i==1) {
        suppressMessages( roc(Survived, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col=colvec[i], lwd=2, quiet=TRUE, main=paste0("Prognostic Ratio as Predictor of Outcome ",outcome,"\n12 Months to 10 Years after Measurement (all LUAD)")) )
      } else {  # Add to existing plot
        suppressMessages( roc(Survived, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col=colvec[i], lwd=2,add=TRUE,quiet=TRUE) )
      }
      suppressMessages( assign(paste0("myroc.",i), roc(Survived, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Positive Percentage", col=colvec[i],lwd=1,add=TRUE,quiet=TRUE)) ) #changed col= from "red" to colvec[i]
      roc.info <- c(roc.info, as.numeric(gsub("%","",gsub("Area under the curve: ","",eval(parse(text=paste0("myroc.",i)))$auc)))/100)
      roc.info.p <- c(roc.info.p, roc.area(obs=Survived, pred=glm.fit$fitted.values)$p.value) # significance of ROC auc from verification pkg function
      #roc.area
      roc.info.95ci <- c(roc.info.95ci, paste0(round(summary(ci(eval(parse(text=paste0("myroc.",i)))))[c(1,6)]/100,4),collapse="-")) #DeLong 95% confidence interval (CI) from pROC function ci()
      assign(paste0("myroc.more.",i), rbind( suppressMessages( reportROC(gold=Survived, predictor=glm.fit$fitted.values,plot=FALSE,positive="s") ), suppressMessages( reportROC(gold=Survived, predictor=glm.fit$fitted.values,plot=FALSE,positive="l") ) ))
      if( signif(roc.info[i],3)==signif(as.numeric(eval(parse(text=paste0("myroc.more.",i)))$AUC[1]),3) ) { #both AUCs assuming Dead/Alive swap are in 2 rows of myroc.more.{i} summing to 1; we cam match to the correct one and pull out Accuracy, Sensitivity, and Specificity.
        roc.info.accuracy <- c(roc.info.accuracy, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$ACC[1]) )
        roc.info.sensitivity <- c(roc.info.sensitivity, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$SEN[1]) )
        roc.info.specificity <- c(roc.info.specificity, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$SPE[1]) )
      } else {
        roc.info.accuracy <- c(roc.info.accuracy, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$ACC[2]) )
        roc.info.sensitivity <- c(roc.info.sensitivity, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$SEN[2]) )
        roc.info.specificity <- c(roc.info.specificity, as.numeric(eval(parse(text=paste0("myroc.more.",i)))$SPE[2]) )
      }
    }
    ROClist[[page]]<-list(myroc.1) #list of 5 ROCs becomes list element [[page]] in ROClist.TNBC
    AUClist[[page]]<-paste0(signif(roc.info,4)*100,"%")
    plist[[page]]<-roc.info.p
    CIlist95[[page]]<-roc.info.95ci
    ACClist[[page]]<-roc.info.accuracy
    
    SENlist[[page]]<-roc.info.sensitivity
    SPElist[[page]]<-roc.info.specificity
    names(ROClist[[page]])<-names(AUClist[[page]])<-names(plist[[page]])<-names(CIlist95[[page]])<-names(ACClist[[page]])<-names(SENlist[[page]])<-names(SPElist[[page]])<-descriptionVec
    
    ## Add legend to completed plot
    descriptionVecLegend<-paste0(descriptionVec," ",signif(roc.info,4)*100,"% AUC")
    legend("bottomright",legend=descriptionVecLegend,lwd=2,col=colvec)
    ## Add title to page
    #for ratio gene lists:  mtext(gsub(".","+",gsub(".."," / ",ratioNames[page],fixed=TRUE),fixed=TRUE),side=3,outer=TRUE,adj=0.5,cex=1.3)
    #for current multi-outcome loop:
    mtext(paste0(outcome,": ",ratioNames[page]),side=3,outer=TRUE,adj=0.5,cex=1.3)
    
    ## Close and Start new PDF if page==blockSize (or a multiple) 
    if (page %% blockSize==0 & !page==iter) {
      dev.off();
      blockStart=page+1
      blockEnd=if (page+blockSize>iter) { iter } else { page+blockSize }
      pdf(file = paste0(rootdir,"/3e",outcomeNumber,".Predicting_",outcome,"_",filenamePrefix,"(",blockStart,"to",blockEnd,").pdf"),width=7,height=7)
      par(mar=c(0,0,3,0))
      par(oma=c(0,0,3,0))
      par(mfrow=c(1,1)) #changed from 2 columns to 1 (from 1,2)
      #plot.new()
    }
    
    ## Progress indicator
    cat(paste0("\rDone page ",page," of ",iter," ..."))
  } #end for loop through different combinations
  #i.e., for (page in 1:iter)
  dev.off()
  
}
  # ## Arrange Tables of stats for AUC (4 signatures, 5 endpoint times, 4 outcomes)
  # AUCtable<-as.data.frame(matrix(unlist(AUClist),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(AUCtable)<-descriptionVec
  # rownames(AUCtable)<-ratioNames
  # Ptable<-as.data.frame(matrix(unlist(plist),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(Ptable)<-descriptionVec
  # CItable<-as.data.frame(matrix(unlist(CIlist95),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(CItable)<-descriptionVec
  # rownames(Ptable)<-rownames(CItable)<-ratioNames
  # pAndCItable<-sapply(1:5,function(x) paste0("p=",signif(Ptable[,x],2)," CI ",CItable[,x]))
  # rownames(pAndCItable)<-ratioNames
  # #pAndCItable # previously taken from console to excel for BrCa Paper Table (3).
  # ACCtable<-as.data.frame(matrix(unlist(ACClist),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(ACCtable)<-descriptionVec
  # SENtable<-as.data.frame(matrix(unlist(SENlist),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(SENtable)<-descriptionVec
  # SPEtable<-as.data.frame(matrix(unlist(SPElist),ncol=length(descriptionVec),byrow=TRUE))
  # colnames(SPEtable)<-descriptionVec
  # #rownames(Ptable)<-rownames(CItable)<-rownames(pAndCItable)<-
  # rownames(ACCtable)<-rownames(SENtable)<-rownames(SPEtable)<-ratioNames
  # tableOut<-cbind(c(1:nrow(AUCtable)),AUCtable,pAndCItable, Ptable,CItable, ACCtable,SENtable,SPEtable)
  # colnames(tableOut)<-c("PDF Page", paste0("AUC",paste0("(",descriptionVec,")")), paste0("P val with 95% CI",paste0("(",descriptionVec,")")), paste0("P Value",paste0("(",descriptionVec,")")), paste0("AUC 95% CI",paste0("(",descriptionVec,")")), paste0("Accuracy",paste0("(",descriptionVec,")")), paste0("Sensitivity",paste0("(",descriptionVec,")")), paste0("Specificity",paste0("(",descriptionVec,")")) )
  # write.csv(tableOut,file=paste0(rootdir,"/3e",outcomeNumber,".Predicting_",outcome,"_",filenamePrefix,"-with_ACC.SEN.SPE.csv"))
  # 
#} #end for(outcome ...)


OS.time.deadonly <- numericMeta$OS.time
OS.time.deadonly[numericMeta$OS==0] <- NA
bicor(ratioList[["thisStudy"]],OS.time.deadonly,use = "p")
# [,1]
# [1,] 0.3262091
bicorAndPvalue(ratioList[["thisStudy"]],OS.time.deadonly,use = "p")
# $bicor
# [,1]
# [1,] 0.3262091
# 
# $p
# [,1]
# [1,] 0.05581522
# 
# $Z
# [,1]
# [1,] 1.944994
# 
# $t
# [,1]
# [1,] 1.982369
# 
# $nObs
# [,1]
# [1,]   35