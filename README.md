# Aclust2.0 instructions

Step 1. Users are required to download a zipped folder from Github, which contains a sample data directory, assay-specific manifests and two R scripts- “aclust2.0.R”, and “aclust2.0_utils.R”. The data folder accommodates test phenotype and betas data for EPIC and MM285 arrays. For the betas data, the rows and columns are CpGs and samples IDs, respectively, with each row of the phenotype data representing an individual sample. “aclust2.0.R” script provides all the necessary guidelines for the five steps required to run the pipeline. All the pipeline functions are housed in “aclust2.0_utils.R”, such that calling this with a “source” function automatically loads all the functions into the R environment. These functions are executed in steps 2-5 below: 

Step 2 (“function get_manifest”). By just providing either “EPICv2”, “EPICv1”, “450K” or “MM285”, this function pulls corresponding manifests downloaded from our Github repository and stored in the same local directory as the “aclust2.0_utils.R”. 

Step 3 (“function find_cluster_list”). This function accepts appropriate Infinium manifest, betas data and CpG IDs as inputs to produce a list that accommodates three objects- “clusters.list”,  “cpg.clusters” “annot.betas”. Both clusters.list and cpg.clusters include all identified CpG clusters completely agnostic of the exposure variables while “annot.betas” is the hg38 chromosomal annotation of each of the CpGs in the cluster. While “cpg.clusters” can be used directly for downstream analysis or passed on to step 5 for chromosomal and nearest gene annotation, “clusters.list” is required for optional GEE analysis in step 4.  Please note that this step filters out methylation sites missing in >20% of the samples. 

Step 4 (“function GEE.clusters”). This function runs GEE models with the identified clusters and takes as input “clusters.list”, betas data, exposure, covariates (optional), “id” which is the column name of betas, and the correlation structure specification. The output (clus) is a dataframe of six columns including as clusters_sites, their respective exposure effect sizes, and unadjusted p-values. The significant clusters (DMRs) can be adjusted for false positives based on users’ preference. 

Step 5 (“function annot.clus.gene”). This function performs chromosomal and nearest gene annotations of clusters with “clus” and “annot.betas” as input. Specifying either “hsa” or “mm” as a model organism provides annotation to Ensembl human (hg38) or mouse (mm10) genome, respectively. The output is a dataframe that includes columns such as cluster-specific coordinates (seqnames, start, end), gene names, cluster names and distance to feature.

For easy data manipulations, we also provided functions (beta2Mval & Mval2beta) to convert between beta and m-values, the latter being more statistically valid for modeling array data.

Future updates will extend functionalities to sequencing data. 


