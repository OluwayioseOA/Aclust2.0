# Import libraries -------------------------------------------------------
library(dplyr)
library(readr)
require(geepack)
require(data.table)
library("biomaRt") 
library("ChIPpeakAnno")
library("GenomicRanges")
library(tibble)

# Step 1:Call functions & import betas -------------------------------------
##Download the code from Github into a named directory and set as working directory
##setwd()
source("aclust2.0_utils.R")

##Note: Remove NAs from betas data before proceeding!!
load("data/betasEPIC.RData") ##Sample data. Please refer to our functions to convert to mvalue if you need to covert to m-value
#load("data/betas450K.RData")
#load("data/betasMM285.RData") ##this sample data contains NAs
#betasMM285 <- na.omit(betasMM285) 

# Step 2: Get Illumina manifests for 450K and EPIC ------------------------
##make sure to download manifest from GitHub and place in the same folder as utils script
manifest2 <- get_manifest("EPIC") #to pull EPIC manifest
# manifest <- get_manifest("450K") #to pull 450K manifest
#manifest <- get_manifest("MM285") #to pull manifest for mouse array (mm10)

# Step 3: Clustering of CpGs ----------------------------------------------
##This step creates a list that accommodates three objects 
list.out <- find_cluster_list(probe.vec = rownames(betas), betas = betas, manifest = manifest, minimum.cluster.size = 2)

clusters.list <- list.out$clusters.list #required for optional GEE model in step 4
annot.betas <- list.out$annot.betas #required chromosomal annotation in step 5
cpg.clusters <- list.out$cpg.clusters #Includes only the CpG clusters which can be passed directly as clus object for gene annotations (step 5)

# Step 4: Run gee model ---------------------------------------------------
pheno <- read.csv("data/pheno_EPIC.csv") ##A dataframe with columns: sampleID, exposure and covariates. Number of observations must match the number of columns in betas
#pheno <- read.csv("data/pheno_MM285.csv")

## Define exposure, covariates (optional) and sample.id
sample.id  <- pheno$sampleID
exposure <- pheno$expo #(can also take categorical exposure variables)
covariates<- pheno %>% dplyr::select(bmi, cotinine)
identical(colnames(betas), sample.id) ##must be in the same order
#betas <- betas[, sample.id] ##to reorder betas by sample.id

###Adjusting for covariates
cluster.gee <- GEE.clusters(betas = betas, clusters.list = list.out$clusters.list, exposure = exposure,
                                 covariates= covariates, id = colnames(betas), working.cor = "ex", sample.id = sample.id) %>% 
  mutate(exposure_padjusted = p.adjust(exposure_pvalue, "BH")) ##adjusting for false positives

###No covariate adjustment 
# cluster.gee <- GEE.clusters(betas = betas, clusters.list = list.out$clusters.list, exposure = exposure,
#                                    covariates = NULL, id = colnames(betas), working.cor = "ex", sample.id = sample.id) %>%
#   mutate(exposure_padjusted = p.adjust(exposure_pvalue, "BH")) ##adjusting for false positives

# Step 5: Perform chromosomal annotations ---------------------------------
#For cluster outputs from step 3 (no phentoypes/predictors)
annot.betas <- list.out$annot.betas
clus <- list.out$cpg.clusters ## (from step 3)
annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "hsa")
#annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "mm")#mouse annotation

#For cluster outputs from gee model in step 4
annot.betas <- list.out$annot.betas
clus <- cluster.gee ## (from step 4)
annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "hsa") 
##annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "mm")#mouse annotation


