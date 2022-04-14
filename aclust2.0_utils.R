beta2Mval<- function(x){return(log2(x/(1-x)))} #function to convert beta values to mvalues
Mval2beta <- function(x){return(2^(x) / (2^x + 1))} # function to convert beta values to mvalues

# Manifests function ----------------------------------------------
get_manifest <- function(platform = c("450K","EPIC","MM285"),...){
  if(platform == "EPIC"){
    annot <- read_csv("EPIC.hg38.manifest.csv") %>% mutate(id = Probe_ID) %>% column_to_rownames(., var = "id")
    return(annot)
  }
  else if(platform == "450K"){
    annot <- read_csv("HM450.hg38.manifest.csv") %>% mutate(id = Probe_ID) %>% column_to_rownames(., var = "id")
    return(annot)
  }
  else if (platform == "MM285") {
    annot <- read_csv("MM285.mm10.manifest.csv") %>% mutate(id = Probe_ID) %>% column_to_rownames(., var = "id")
    return(annot)
  }
  else{
    stop("You need to specify either EPIC or 450K to run")
  }
}
# Dbp.merge function -----------------------------------------------------
######################
update.clust.indicator <-
  function(which.clust, ind.1, ind.2){
    ## merge all points between ind.1 and ind.2, and previous/later point if belong to the same cluster. 
    ## which.clust is a vector of cluster numbers, so that for each point it says which cluster in belongs to. 
    ## ind.1 and ind.2 are two points that should be merged into the same vector (as are all points between them)
    
    stopifnot(ind.2 > ind.1)
    
    ind.2 <- last(which(which.clust == which.clust[ind.2]))
    which.clust[(ind.1 + 1):ind.2] <- which.clust[ind.1]
    return(which.clust)
  }

######################
Dbp.merge <- function (ordr.vec, thresh.dist, bp.thresh.dist, location.vec, 
                       dist.type = "spearman") 
{
  stopifnot(is.element(tolower(dist.type), c("spearman", "pearson", 
                                             "euclid")))
  le <- dim(ordr.vec)[2]
  dist.clust <- rep(0, le - 1)
  which.clust <- 1:le
  max.d <- 1
  for (d in 1:min(bp.thresh.dist, length(location.vec))) {
    if (min(diff(location.vec, lag = d)) <= bp.thresh.dist) 
      max.d <- d
    else break
  }
  dist.mat.by.nr <- matrix(0, ncol = le, nrow = max.d)
  for (nr in 1:max.d) {
    dist.mat.by.nr[nr, ] <- calc.dist.d.neighbor(ordr.vec, 
                                                 nr, dist.type)
  }
  loc.diff.mat <- matrix(Inf, nrow = max.d, ncol = le)
  for (i in 1:max.d) {
    diff.temp <- diff(location.vec, lag = i)
    loc.diff.mat[i, 1:length(diff.temp)] <- diff.temp
  }
  nr.number <- apply(loc.diff.mat, 2, function(x) max(0, which(x < 
                                                                 bp.thresh.dist)))
  for (k in 1:(le - 1)) {
    if (le - k > nr.number[k]) 
      d.cur <- nr.number[k]
    else d.cur <- le - k
    if (d.cur > 0) {
      for (nr in d.cur:1) {
        if (which.clust[k] == which.clust[k + nr]) 
          break
        cur.dist <- dist.mat.by.nr[nr, k]
        if (cur.dist < thresh.dist) {
          which.clust <- update.clust.indicator(which.clust, 
                                                k, k + nr)
          break
        }
      }
    }
  }
  return(which.clust)
}

# Acluster function -------------------------------------------------------
###################
calc.dist.clust.point <-
  function(clust, pnt, type = "single", dist.type = "spearman"){
    ## caculates the distance between a set of vectors and a single vector. 
    ## the distance is either the min, max or mean of distances between the vector and each of the vectors in the set
    ## the distance function is dist.type
    
    stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "pearson", "euclid")), dim(clust)[1] == length(pnt))
    
    clust <- as.matrix(clust)
    size <- ncol(clust)
    if (size == 1){
      if (is.element(dist.type , c("pearson", "spearman"))) return(1 - abs(cor(clust, pnt, use = "complete.obs", method = dist.type))) else
        if (dist.type == "euclid") {
          inds.rm <- union(which(is.na(clust)), which(is.na(pnt)))
          if (length(inds.rm) > 0) return(sqrt(sum((clust[-inds.rm] - pnt[-inds.rm])^2))) else
            return(sqrt(sum((clust - pnt)^2)))	
        }	
      
    } else{
      
      distances <- apply(clust, 2, function(x){
        if (is.element(dist.type , c("pearson", "spearman"))) return(1 - abs(cor(x, pnt, use = "complete.obs", method = dist.type))) else
          if (dist.type == "euclid") {
            inds.rm <- union(which(is.na(x)), which(is.na(pnt)))
            if (length(inds.rm) > 0) return(sqrt(sum((x[-inds.rm] - pnt[-inds.rm])^2))) else
              return(sqrt(sum((x - pnt)^2)))	
          }	
      } )
      
      if (type == "single") return(min(distances)) else
        if (type == "complete") return(max(distances)) else
          return(mean(distances)) }
  }

###################
calc.dist.clusters <-
  function(clust.1, clust.2, type = "single", dist.type = "spearman"){
    ## caculates the distance between two sets of vectors
    ## the distance is either the min, max or mean of distances between the each of the vectors in one set
    ## and each of the vectors in the other set
    ## the distance function is dist.type
    stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "pearson", "euclid")))
    
    clust.1 <- as.matrix(clust.1)
    size.1 <- ncol(clust.1)
    if (size.1 == 1) return(calc.dist.clust.point(clust.2, clust.1, type = type, dist.type = dist.type))
    
    distances <- apply(clust.1, 2, function(x){
      calc.dist.clust.point(clust.2, x, type = type, dist = dist.type)
    })
    
    if (type == "single") return(min(distances)) else
      if (type == "complete") return(max(distances)) else
        return(mean(distances))
  }


calc.dist.d.neighbor <-
  function(ordr.vec, d, dist.type = "spearman"){
    ## calucate the distances between vectors, and their d neighbor. 
    
    ### initialize variables
    le <- dim(ordr.vec)[2]
    dist.vec.d <- rep(0, le - d)
    
    ### calculate the distances vector
    if (is.element(dist.type , c("pearson", "spearman"))){
      for (i in 1:(le-d)){
        dist.vec.d[i] <- 1- abs(cor(ordr.vec[,i], ordr.vec[,i+d], use = "complete.obs", method = dist.type))
      } } else if (dist.type == "euclid"){
        inds.rm <- union(which(is.na(ordr.vec[,i])), which(is.na(ordr.vec[,i + d])))
        if (length(inds.rm) > 0) dist.vec.d[i] <- sqrt(sum((ordr.vec[-inds.rm, i] - ordr.vec[-inds.rm, i + d])^2)) else
          dist.vec.d[i] <- sqrt(sum((ordr.vec[,i] - ordr.vec[,i + d])^2))
      }
    
    dist.vec.d <- c(dist.vec.d, rep(Inf, d))		
    
    return(dist.vec.d)
  }
###################

acluster <- function (ordr.vec, thresh.dist, which.clust = NULL, location.vec = NULL, 
                      max.dist = Inf, type = "single", dist.type = "spearman") 
{
  stopifnot(is.element(type, c("single", "complete", "average")), 
            is.element(tolower(dist.type), c("spearman", "pearson", 
                                             "euclid")))
  if (!is.null(max.dist) & !all(!is.na(location.vec))) 
    stop("missing location values")
  
  # defining two functions that are used only in this function
  first <- function(x){x[1]}
  last <- function(x){x[length(x)]}
  
  le <- dim(ordr.vec)[2]
  dist.clust <- rep(0, le - 1)
  if (is.null(which.clust)) {
    which.clust <- 1:le
    dist.clust <- calc.dist.d.neighbor(ordr.vec, 1, dist.type)
  }
  else {
    dist.clust <- rep(0, le - 1)
    for (k in 1:(le - 1)) {
      if (which.clust[k] == which.clust[k + 1]) 
        dist.clust[k] <- 0
      else {
        clust.1.min <- first(which.clust[which.clust == 
                                           which.clust[k]])
        clust.1.max <- k
        clust.2.min <- k + 1
        clust.2.max <- last(which.clust[which.clust == 
                                          which.clust[k + 1]])
        dist.clust[k] <- calc.dist.clusters(ordr.vec[, 
                                                     clust.1.min:clust.1.max], ordr.vec[, clust.2.min:clust.2.max], 
                                            type = type, dist.type = dist.type)
      }
    }
  }
  no.more <- FALSE
  condition1 <- function(x) {
    return((x < thresh.dist) & is.null(max.dist))
  }
  condition2 <- function(x, loc1, loc2) {
    if (is.null(max.dist)) 
      return(FALSE)
    if (is.null(loc1) | is.null(loc2)) 
      return(FALSE)
    else return((x < thresh.dist) & !is.null(max.dist) & 
                  (max.dist >= abs(loc2 - loc1)))
  }
  while (!no.more) {
    no.more = T
    clust.1.min <- 1
    for (k in 1:(le - 1)) {
      if (which.clust[k] != which.clust[k + 1]) {
        if (condition1(dist.clust[k]) | condition2(dist.clust[k], 
                                                   location.vec[k], location.vec[k + 1])) {
          no.more = F
          clust.to.merge <- which.clust[k + 1]
          l <- k + 1
          while (which.clust[l] == clust.to.merge) {
            which.clust[l] <- which.clust[k]
            l <- l + 1
            if (l > le) 
              break
          }
          clust.1.max <- l - 1
          dist.clust[clust.1.min:(clust.1.max - 1)] <- Inf
          if (clust.1.max < le) {
            clust.2.min <- l
            clust.2.max <- l
            if (l < le) {
              l <- l + 1
              while (which.clust[l] == which.clust[clust.2.min]) {
                l <- l + 1
                if (l > le) 
                  break
              }
              clust.2.max <- l - 1
            }
            dist.clust[clust.1.max] <- calc.dist.clusters(as.matrix(ordr.vec[, 
                                                                             clust.1.min:clust.1.max]), as.matrix(ordr.vec[, 
                                                                                                                           clust.2.min:clust.2.max]), type = type, 
                                                          dist.type = dist.type)
          }
          if (clust.1.min > 1) {
            clust.0.max <- clust.1.min - 1
            clust.0.min <- clust.0.max
            if (clust.0.min > 1) {
              l <- clust.0.min - 1
              while (which.clust[l] == which.clust[clust.0.max]) {
                which.clust[l] <- which.clust[clust.0.max]
                l <- l - 1
                if (l < 1) 
                  break
              }
              clust.2.min <- l + 1
            }
            dist.clust[clust.0.max] <- calc.dist.clusters(as.matrix(ordr.vec[, 
                                                                             clust.0.min:clust.0.max]), as.matrix(ordr.vec[, 
                                                                                                                           clust.1.min:clust.1.max]), type = type, 
                                                          dist.type = dist.type)
          }
        }
        else clust.1.min <- k + 1
      }
    }
  }
  return(which.clust)
}

# Find_cluster_list function --------------------------------------------------------
find_cluster_list <- function (probe.vec, betas, manifest, minimum.cluster.size = 2, 
                               thresh.dist = 0.25, bp.thresh.dist = 999,
                               max.dist = 1000, type = "average", dist.type = "spearman", 
                               missingness_max_prop = 0.2) {
  message(paste("Checking missingness, allowing maximum missingness proportion", missingness_max_prop))
  missingness_prop <- apply(betas, 1, function(x) mean(is.na(x)))
  inds_rm_probes <- which(missingness_prop > missingness_max_prop)

  if (length(setdiff(rownames(betas), manifest$Probe_ID))!=0) {
    message(paste("Removed", length(setdiff(rownames(betas), manifest$Probe_ID)),
                  "betas CpGs missing in the manifest due to a change in Illumina manufacturing process"))
    betas <- betas %>% dplyr::filter(rownames(.)%in%manifest$Probe_ID)
    
    if (length(inds_rm_probes) > 0){
      message(paste("Removing", length(inds_rm_probes), 
                    "methylation sites due to high missingness"))
      betas <- betas[-inds_rm_probes,]
      probe.vec <- probe.vec[-inds_rm_probes]
    }
    
  }
  # if (length(inds_rm_probes) > 0){
  #   message(paste("Removing", length(inds_rm_probes), 
  #                 "methylation sites due to high missingness"))
  #   betas <- betas[-inds_rm_probes,]
  #   probe.vec <- probe.vec[-inds_rm_probes]
  # }
  
  # if(all(!is.na(betas))==F) {
  #  stop("Error: There are NAs in betas data. Please remove before proceeding!")
  #}else{
  annot.betas <- as.data.frame(manifest) %>% filter(rownames(.)%in%probe.vec) %>% 
    dplyr::select(Probe_ID, seqnames, probeTarget) %>%
    setNames(c("IlmnID", "CHR", "Coordinate")) %>% 
    mutate(IlmnID = as.character(IlmnID), CHR = as.character(CHR), Coordinate = as.numeric(Coordinate)) %>% 
    filter(!duplicated(IlmnID)) %>% data.table(., key = c("CHR","Coordinate"))
  
  CHR = annot.betas$CHR
  Coordinate = annot.betas$Coordinate
  IlmnID = annot.betas$IlmnID
  chroms <- unique(annot.betas$CHR)
  if(length(unique(CHR))==24){
    return.chroms = paste("chr",c(1:22, "X", "Y"),sep="")
  } else{
    return.chroms = paste("chr",c(1:19, "X", "Y", "M"),sep="")
  }
  chroms <- intersect(chroms, return.chroms)
  chroms <- chroms[order(as.numeric(as.character(chroms)))] 
  betas.by.chrom <- vector(mode = "list", length = length(chroms))
  sites.by.chrom <- vector(mode = "list", length = length(chroms))
  names(betas.by.chrom) <- names(sites.by.chrom) <- chroms
  for (i in 1:length(chroms)) {
    cpg.chrom <- as.character(annot.betas[CHR == chroms[i]]$IlmnID)
    betas.by.chrom[[i]] <- as.matrix(betas[cpg.chrom, ])
    if (ncol(betas.by.chrom[[i]]) == 1) {
      betas.by.chrom[[i]] <- t(betas.by.chrom[[i]])
      rownames(betas.by.chrom[[i]]) <- cpg.chrom
    }
    sites.by.chrom[[i]] <- annot.betas[CHR == chroms[i],
                                       c("IlmnID", "Coordinate"), with = F]
  }
  chrom.list = list(betas.by.chrom = betas.by.chrom, sites.locations.by.chrom = sites.by.chrom)
  clusters.by.chrom <- vector(mode = "list", length = length(chrom.list[[1]]))
  for (i in 1:length(chrom.list[[1]])) {
    betas.temp <- chrom.list[[1]][[i]]
    locations.temp <- chrom.list[[2]][[i]]
    betas.temp <- betas.temp[which(!is.na(locations.temp$Coordinate)), ]
    locations.temp <- locations.temp[!is.na(Coordinate)]
    which.clust <- Dbp.merge(t(betas.temp), 
                             thresh.dist = thresh.dist, 
                             bp.thresh.dist = bp.thresh.dist, 
                             as.numeric(locations.temp$Coordinate), 
                             dist.type = dist.type)
    clust.vec <- acluster(t(betas.temp), 
                          thresh.dist = thresh.dist, 
                          which.clust = which.clust, 
                          location.vec = chrom.list$sites.locations.by.chrom[[i]]$Coordinate, 
                          max.dist = max.dist, 
                          type = type, 
                          dist.type = dist.type)
    clusters.by.chrom[[i]] <- lapply(clust.vec, function(x) return(locations.temp$IlmnID[which(clust.vec == x)]))
    clusters.by.chrom[[i]] <- clusters.by.chrom[[i]][which(!duplicated(clusters.by.chrom[[i]]))]
    if (i == 1)
      clusters.all <- clusters.by.chrom[[1]]
    else clusters.all <- c(clusters.all, clusters.by.chrom[[i]])
  }
  
  n.sites <- unlist(lapply(clusters.all, function(x) return(length(x))))
  inds.rm <- which(n.sites < minimum.cluster.size)
  
  if (length(inds.rm) > 0) clusters.all <- clusters.all[-inds.rm]
  
  n.mod <- length(clusters.all)
  sites <- rep("", n.mod)
  n.sites <- rep(0, n.mod)
  
  for (i in 1:n.mod){
    clus.probes <- clusters.all[[i]]
    clus.betas <- betas[clus.probes,]
    
    if (length(clus.probes) == 1){
      
      ind.comp <- which(complete.cases(clus.betas))
      n.sites[i] <- 1
      sites[i] <- clus.probes
      
    } else{
      
      ind.comp <- which(complete.cases(t(clus.betas)))
      temp.betas <- clus.betas[,ind.comp]
      n.sites[i] <- length(clus.probes)
      for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
      substr(sites[i], 1, 1) <- ""
    }
    
  }
  cpg.clusters <- data.frame(n_sites_in_cluster = n.sites, cluster_sites =sites) %>%
    dplyr::mutate(cluster_name = paste0("cluster_", seq(1:nrow(.))))
  
  return(list(annot.betas = annot.betas, clusters.list = clusters.all, cpg.clusters = cpg.clusters))
  # }
}

# GEE.clusters function ---------------------------------------------------
organize.island.repeated <-
  function(X, covar.mat){
    ## function recives X matrix of betas values (columns are people) and matrix covar.mat of covariates 
    ## (columns are covariates). returns reshaped matrix. 
    data <- cbind(t(X), covar.mat)
    data.long <- reshape(data, direction = "long", varying = list(names(data)[1:dim(X)[1]]), v.names = "Beta", idvar = colnames(covar.mat), timevar="probeID", times = rownames(X))
    rownames(data.long) <- unlist(lapply(rownames(rownames(X)), function(x) rep(x,ncol(X))))    
    
    return(data.long)
    
  }

GEE.clusters <- function(betas, clusters.list, exposure, id, covariates = NULL,  working.cor = "ex", minimum.cluster.size = 2, result.file.name = NULL, sample.id){
  ##Please note the following:
  ##The ids are rownames of the betas, and must be ordered by the same ordering of the covariates matrix.
  ##This function gets a matrix of beta values, a list of clusters, exposure variable, and covariates
  ##and performs GEE. Simply supply NULL if no covariates are available. 
  ##Methylation is treated here as outcome
  #It also accepts m-value. The simple functions to convert between m-value and beta are also provided 
  ##Returns a matrix with the results of the GEE model (only the exposure effect and p-value) for each cluster.
  ##The results are ordered according to p-value
  require(geepack)
  if(!identical(colnames(betas), sample.id)) {
    stop("Error: colnames of betas must match sample IDs in pheno data and in the same order")
  }
  else{
    
    n.sites <- unlist(lapply(clusters.list, function(x) return(length(x))))
    inds.rm <- which(n.sites < minimum.cluster.size)
    
    if (length(inds.rm) > 0) clusters.list <- clusters.list[-inds.rm]
    
    
    n.mod <- length(clusters.list)
    
    effect <- rep(NA, n.mod)
    se <- rep(NA, n.mod)
    pvals <- rep(NA, n.mod)
    sites <- rep("", n.mod)
    n.sites <- rep(NA, n.mod)
    n.samp <- rep(NA, n.mod)
    
    if (!is.null(covariates)){
      
      if (is.null(dim(covariates))) covariates <- as.matrix(covariates)
      
      if (is.null(colnames(covariates)))  colnames(covariates) <- rep("", ncol(covariates))
      for (i in 1:ncol(covariates)){
        if (colnames(covariates)[i] == "") colnames(covariates)[i] <- paste("covariate", i, sep = '_')
        
      }
      
      model.expression <- "model <- geeglm(Beta ~ exposure +"
      for (j in 1:ncol(covariates)){
        model.expression <- paste(model.expression, colnames(covariates)[j], "+")  }
      
      model.expression <- paste(model.expression, "as.factor(probeID), id = as.factor(id), data = temp.long, corstr = working.cor)")
      
      model.expr.1.site <- paste("model <- geeglm(clus.betas[ind.comp] ~ exposure[ind.comp] + ")  
      for (j in 1:ncol(covariates)){
        if (j < ncol(covariates)) model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j], "+")
        else    model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j])}
      
      model.expr.1.site <- paste(model.expr.1.site, ", data = covariates[ind.comp,], id = as.factor(id)[ind.comp])")
      covariates <- as.data.frame(covariates)
      
      
      for (i in 1:n.mod){
        clus.probes <- clusters.list[[i]]
        clus.betas <- betas[clus.probes,]
      if (length(clus.probes) == 1){
        
        ind.comp <- which(complete.cases(cbind(clus.betas, exposure, id)))
        eval(parse(text = model.expr.1.site))
        
        effect[i] <- summary(model)$coef["exposure",1]
        se[i] <- summary(model)$coef["exposure",2]
        pvals[i] <- summary(model)$coef["exposure",4]
        n.sites[i] <- 1
        sites[i] <- clus.probes
        n.samp[i] <- length(ind.comp)
        
      } else{
        
        ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure,  id)))
        
        temp.betas <- clus.betas[,ind.comp]
        
        temp.covars <- data.frame(exposure[ind.comp], id[ind.comp])
        colnames(temp.covars) <- c("exposure",  "id")
        
        temp.long <- organize.island.repeated(temp.betas, temp.covars)
        temp.long <- temp.long[complete.cases(temp.long),]
        temp.long <- temp.long[order(temp.long$id),]
        
        if (length(clus.probes) == 1){
          
          ind.comp <- which(complete.cases(cbind(clus.betas, exposure, covariates, id)))
          eval(parse(text = model.expr.1.site))
          
          effect[i] <- summary(model)$coef["exposure",1]
          se[i] <- summary(model)$coef["exposure",2]
          pvals[i] <- summary(model)$coef["exposure",4]
          n.sites[i] <- 1
          sites[i] <- clus.probes
          n.samp[i] <- length(ind.comp)
          
          
        } else{
          
          ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure, covariates, id)))
          
          temp.betas <- clus.betas[,ind.comp]
          
          temp.covars <- as.data.frame(cbind(exposure[ind.comp], covariates[ind.comp,], id[ind.comp]))
          colnames(temp.covars) <- c("exposure", colnames(covariates), "id")
          
          temp.long <- organize.island.repeated(temp.betas, temp.covars)
          temp.long <- temp.long[complete.cases(temp.long),]
          temp.long <- temp.long[order(temp.long$id),]
          
          eval(parse(text = model.expression))
          
          
          effect[i] <- summary(model)$coef["exposure",1]
          se[i] <- summary(model)$coef["exposure",2]
          pvals[i] <- summary(model)$coef["exposure",4]
          n.sites[i] <- length(clus.probes)
          n.samp[i] <- length(ind.comp)
          for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
          substr(sites[i], 1, 1) <- ""
          
        }
        
      }
      
    } else{
      
      model.expression <- "model <- geeglm(Beta ~ exposure +"
      model.expression <- paste(model.expression, "as.factor(probeID), id = as.factor(id), data = temp.long, corstr = working.cor)")
      
      model.expr.1.site <- paste("model <- geeglm(clus.betas[ind.comp] ~ exposure[ind.comp] ")      
      model.expr.1.site <- paste(model.expr.1.site, ", data = covariates[ind.comp,], id = as.factor(id)[ind.comp])")
      
      for (i in 1:n.mod){
        clus.probes <- clusters.list[[i]]
        clus.betas <- betas[clus.probes,]
        
        
        if (length(clus.probes) == 1){
          
          ind.comp <- which(complete.cases(cbind(clus.betas, exposure, id)))
          eval(parse(text = model.expr.1.site))
          
          effect[i] <- summary(model)$coef["exposure",1]
          se[i] <- summary(model)$coef["exposure",2]
          pvals[i] <- summary(model)$coef["exposure",4]
          n.sites[i] <- 1
          sites[i] <- clus.probes
          n.samp[i] <- length(ind.comp)
          
        } else{
          
          ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure,  id)))
          
          temp.betas <- clus.betas[,ind.comp]
          
          temp.covars <- data.frame(exposure[ind.comp], id[ind.comp])
          colnames(temp.covars) <- c("exposure",  "id")
          
          temp.long <- organize.island.repeated(temp.betas, temp.covars)
          temp.long <- temp.long[complete.cases(temp.long),]
          temp.long <- temp.long[order(temp.long$id),]
          
          eval(parse(text = model.expression))
          
          
          effect[i] <- summary(model)$coef["exposure",1]
          se[i] <- summary(model)$coef["exposure",2]
          pvals[i] <- summary(model)$coef["exposure",4]
          n.sites[i] <- length(clus.probes)
          n.samp[i] <- length(ind.comp)
          for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
          substr(sites[i], 1, 1) <- ""
          
        }
      }
      
      
    }
  }
  result <- data.frame(exposure_effect_size = effect, exposure_effect_se = se, exposure_pvalue= pvals, n_sites_in_cluster = n.sites, n_samples = n.samp, cluster_sites =sites)  
  result <- result[order(result$exposure_pvalue),]
  
  
  if(!is.null(result.file.name)) write.table(result, file = result.file.name, col.names = T , row.names = F,  append = T)
  
  
  return(result)
  
}

# Chromosomal annotation function ------------------------------
chrom_annot <- function(x){
  prb <- unlist(strsplit(clus[x,"cluster_sites"],";"))
  prb <- prb[!(prb=="")]
  if(length(prb)==1){
    myan= annot.betas[IlmnID==prb,]
    chr=as.character(myan$CHR)
    start=as.numeric(myan$Coordinate)
    end=as.numeric(myan$Coordinate)+1
    probe=length(prb)
    inter= data.frame(chr,start,end,probe)
  }
  if(length(prb)>1){
    myan= annot.betas[IlmnID%in%prb,]
    chr=unique(as.character(myan$CHR))
    start=min(as.numeric(myan$Coordinate))
    end=max(as.numeric(myan$Coordinate))
    probe=length(prb)
    inter= data.frame(chr,start,end,probe,stringsAsFactors = FALSE)
  }
  return(inter)
}

# Nearest gene annotation --------------------------------------------

annot.clus.gene <- function(annot.betas, clus, model = c("mm", "hsa")){
  if(model=="hsa"){
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", 
                      path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  }else{
    
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                      host="https://www.ensembl.org", path="/biomart/martservice", 
                      dataset="mmusculus_gene_ensembl")
    
  }
  anno <- t(sapply(1:nrow(clus),chrom_annot)) %>% data.frame(.) %>% 
    mutate(chr = unlist(chr), start = as.numeric(start), end = as.numeric(end), 
           width = end - start, probe = unlist(probe), cluster_name = paste("cluster_",rownames(.),sep="")) %>% 
    cbind(clus) %>% GRanges(.) %>% ChIPpeakAnno:: annotatePeakInBatch(., ensembl, featureType = c("Exon"),
                                                                      output=c("nearestLocation"),multiple=c(TRUE),maxgap=0L, 
                                                                      PeakLocForDistance=c("middle"),FeatureLocForDistance=c("TSS"),
                                                                      select=c("all"),ignore.strand=TRUE) %>% data.frame(.) 
  bm <- getBM(attributes=c("ensembl_exon_id","external_gene_name"),
              filters='ensembl_exon_id', values=anno$feature, mart=ensembl) %>% 
    rename(ensembl_exon_id = "feature") %>% right_join(anno) %>% filter(!duplicated(c(cluster_name)))
  
  return(bm)
}
