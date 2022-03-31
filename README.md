# Alcust2.0 user instructions

Step 1: Users are required to download a zipped folder from Github, which contains a sample data directory, assay-specific manifests and two R scripts- “aclust2.0.R”, and “aclust2.0_utils.R”. The data folder accommodates test phenotype and betas data for EPIC and MM285 arrays. For the betas data, the rows and columns are CpGs and samples IDs, respectively, with each row of the phenotype data representing an individual sample. Please note that complete betas and phenotype data is required to run the entire pipeline. “aclust2.0.R” script provides all the necessary guidelines for the five steps required to run the pipeline (Figure 1). All the pipeline functions are housed in “aclust2.0_utils.R”, such that calling this with a “source” function automatically loads all the functions into the R environment. These functions are executed in the following steps: 

Step 2 (“function get_manifest”). By just providing either “EPIC”, “450K” or “MM285”, this function pulls corresponding manifests downloaded from our GitHub repository and stored in the same local directory as the “aclust2.0_utils.R”. 
Step 3 (“function find_cluster_list”). This function accepts appropriate Infinium manifest, betas data and CpG IDs as inputs to produce a list that accommodates three objects- “clusters.list”,  “cpg.clusters” “annot.betas”. Both clusters.list and cpg.clusters include all identified CpG clusters completely agnostic of the exposure variables while “annot.betas” is the hg19 chromosomal annotation of each of the CpGs in the cluster. While “cpg.clusters” can be used directly for downstream analysis or passed on to step 5 for chromosomal and nearest gene annotation, “clusters.list” is required for optional GEE analysis in step 4.  Please note that this step requires complete methylation (betas) data without “NAs” and matching of its column names with sample IDs from phenotype data. 

Step 4 (“function GEE.clusters”). This function runs GEE models with the identified clusters and takes as input “clusters.list”, betas data, exposure, covariates (optional), “id” which is the column name of betas, and the correlation structure specification. The output (clus) is a dataframe of six columns including as clusters_sites, their respective exposure effect sizes, and unadjusted p-values. The significant clusters (DMRs) can be adjusted for false positives based on users’ preference. 

Step 5 (“function annot.clus.gene”). This function performs chromosomal and nearest gene annotations of clusters with “clus” and “annot.betas” as input. Specifying either “hsa” or “mm” as a model organism provides annotation to Ensembl human (hg19) or mouse (mm10) genome, respectively. The output is a dataframe that includes columns such as cluster-specific coordinates (seqnames, start, end), gene names, cluster names and distance to feature. To convert resulting annotation files from hg19 assembly to hg38, users can apply a “liftover” tool available at the UCSC genome browser.

For easy data manipulations, we also provided functions (beta2Mval & Mval2beta) to convert between beta and m-values, the latter being more statistically valid for modeling array data.

# Clustering only


setwd() ##set working directory

source("aclust2.0_utils.R")

load("data/betasEPIC.RData")

manifest <- get_manifest("EPIC")

list.out <- find_cluster_list(probe.vec = rownames(betas), betas = betas, manifest = 

manifest, minimum.cluster.size = 2)

clusters.list <- list.out$clusters.list #required for optional GEE model in step 4
annot.betas <- list.out$annot.betas #required chromosomal annotation in step 5
cpg.clusters <- list.out$cpg.clusters #Includes only the CpG clusters which can be passed directly as clus object for gene annotations (step 5)


# Optional GEE model


pheno <- read.csv("data/pheno_EPIC.csv") 
sample.id  <- pheno$sampleID
exposure <- pheno$expo #(can also take categorical exposure variables)
covariates<- pheno %>% dplyr::select(bmi, cotinine)
identical(colnames(betas), sample.id) ##must be in the same order

cluster.gee <- GEE.clusters(betas = betas, clusters.list = list.out$clusters.list, exposure = exposure, covariates= covariates, id = colnames(betas), working.cor = "ex", sample.id = sample.id) %>% mutate(exposure_padjusted = p.adjust(exposure_pvalue, "BH")) ##adjusting for false positives

annot.betas <- list.out$annot.betas
clus <- list.out$cpg.clusters ## (from step 3)
annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "hsa")

annot.betas <- list.out$annot.betas
clus <- cluster.gee ## (from step 4)
annotated.genes <- annot.clus.gene(annot.betas = annot.betas, clus = clus, model = "hsa") 



