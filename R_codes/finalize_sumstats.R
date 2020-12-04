#!/usr/bin/env Rscript

library(e1071)

final.sumstats <- c(Fst.s1.adm = NA, Fst.s2.adm = NA, Fst.s1.s2 = NA, mean.ASD.s1.adm = NA, var.ASD.s1.adm = NA, mean.ASD.s2.adm = NA, var.ASD.s2.adm = NA, mean.ASD.s1.s2 = NA, var.ASD.s1.s2 = NA, f3 = NA, mean.ASD.s1 = NA, var.ASD.s1 = NA, mean.ASD.s2 = NA, var.ASD.s2 = NA, mean.ASD.adm = NA, var.ASD.adm = NA, mean.het.s1 = NA, var.het.s1 = NA, mean.het.s2 = NA, var.het.s2 = NA, mean.het.adm = NA, var.het.adm = NA, mean.F.s1 = NA, var.F.s1 = NA, mean.F.s2 = NA, var.F.s2 = NA, mean.F.adm = NA, var.F.adm = NA, mean.adm.props=NA, var.adm.props=NA, skew.adm.props=NA, kurt.adm.props=NA, mode.adm.props=NA, mean.adm.angles=NA, var.adm.angles=NA, skew.adm.angles=NA, kurt.adm.angles=NA, mode.adm.angles=NA)
tmp <- sprintf('perc%d.adm.props', seq(0,100,10))
final.sumstats[tmp] <- NA
tmp <- sprintf('perc%d.adm.angles', seq(0,100,10))
final.sumstats[tmp] <- NA

## pairwise Fst computed by vcftools
get.Fsts <- function(){
    d <- read.table('result.fst')
    final.sumstats['Fst.s1.adm'] <- d[1,]
    final.sumstats['Fst.s2.adm'] <- d[2,]
    final.sumstats['Fst.s1.s2'] <- d[3,]
    final.sumstats
}


## we use 'tryCatch' tricks to prevent problems if files do not exist
tryCatch(final.sumstats <- get.Fsts(), error=function(e){invisible()}, warning=function(w){invisible()})



## F3 statistic, following Patterson 2012
compute.F3 <- function(){
    ## get allelic frequencies computed by vcftools
    p.adm <- read.table('adm.frq', skip=1)[,6]
    p.s1 <- read.table('s1.frq', skip=1)[,6]
    p.s2 <- read.table('s2.frq', skip=1)[,6]
    d <- data.frame(p.adm, p.s1, p.s2)
    ## compute numerators and denomitors of statistic
    num <- sum((d$p.adm - d$p.s1) * (d$p.adm - d$p.s2))
    denom <- sum(2 * d$p.adm * (1-d$p.adm))
    ## compute statistic
    final.sumstats['f3'] <- num/denom
    final.sumstats
}

tryCatch(final.sumstats <- compute.F3(), error=function(e){invisible()}, warning=function(w){invisible()})

compute.inbreeding.het <- function(){
    for (pop in c('s1', 's2', 'adm')){
        ## inbreeding coefficients from vcftools
        d <- read.table(sprintf('%s.het', pop), header=TRUE)
        final.sumstats[sprintf('mean.F.%s', pop)] <- mean(d$F)
        final.sumstats[sprintf('var.F.%s', pop)] <- var(d$F)
        ## we compute heterozigosity using allelic frequencies computed by vcftools
        d <- read.table(sprintf('%s.frq', pop), skip=1)
        freqs <- d[,6]
        nb.sample <- d[1,4]
        tmp <- freqs*freqs + (1-freqs)*(1-freqs)
        tmp <- 1 - tmp
        tmp <- tmp * nb.sample / (nb.sample - 1)
        final.sumstats[sprintf('mean.het.%s', pop)] <- mean(tmp)
        final.sumstats[sprintf('var.het.%s', pop)] <- var(tmp)
    }
    final.sumstats
}

tryCatch(final.sumstats <- compute.inbreeding.het(), error=function(e){invisible()}, warning=function(w){invisible()})

## compute statistics based on ASD (which was computed by asd)
compute.ASD.stats <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c('s1','s2','adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    for (i in 1:3){
        pop1 <- all.pops[i]
        wanted.ids.1 <- indivs$ids[indivs$pop == pop1]
        for (j in i:3){
            pop2 <- all.pops[j]
            wanted.ids.2 <- indivs$ids[indivs$pop == pop2]
            tmp <- as.matrix(d[wanted.ids.1, wanted.ids.2])
            
            if (i == j){
                tmp <- tmp[lower.tri(tmp)]
                pop <- pop1
                final.sumstats[sprintf('mean.ASD.%s', pop)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s', pop)] <- var(tmp)
            } else {
                tmp <- as.vector(tmp)
                final.sumstats[sprintf('mean.ASD.%s.%s', pop1, pop2)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s.%s', pop1, pop2)] <- var(tmp)
            }
        }
    }
    final.sumstats
}


tryCatch(final.sumstats <- compute.ASD.stats(), error=function(e){invisible()}, warning=function(e){invisible()})

## computes admixture proportions stats using projection in MDS
compute.adm.props <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c('s1','s2','adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    mds <- suppressWarnings(cmdscale(as.matrix(d), k=2))
    ids.s1 <- indivs$ids[indivs$pop == 's1']
    ids.s2 <- indivs$ids[indivs$pop == 's2']
    centroid.s1 <- colMeans(mds[ids.s1,])
    centroid.s2 <- colMeans(mds[ids.s2,])
    centroids.dist <- dist(rbind(centroid.s1,centroid.s2))[1]
    ids.adm <- indivs$ids[indivs$pop == 'adm']
    adm.props <- rep(0, length(ids.adm))
    vec.b <- centroid.s2 - centroid.s1
    for (i in 1:length(adm.props)){
        idx.adm <- which(rownames(mds) == ids.adm[i])
        vec.a <- mds[idx.adm,] - centroid.s1
        proj <- (sum(vec.a*vec.b) / sum(vec.b*vec.b)) * vec.b
        new.point <- centroid.s1 + proj
        adm.props[i] <- 1 - dist(rbind(centroid.s1, new.point))[1] / centroids.dist
    }
    final.sumstats['mean.adm.props'] <- mean(adm.props)
    final.sumstats['var.adm.props'] <- var(adm.props)
    final.sumstats['skew.adm.props'] <- skewness(adm.props)
    final.sumstats['kurt.adm.props'] <- kurtosis(adm.props)
    tmp <- density(adm.props)
    final.sumstats['mode.adm.props'] <- tmp$x[which.max(tmp$y)][1]
    tmp <- sprintf('perc%d.adm.props', seq(0,100,10))
    final.sumstats[tmp] <- quantile(adm.props, seq(0,1,.1))
    final.sumstats
}

tryCatch(final.sumstats <- compute.adm.props(), error=function(e){invisible()}, warning=function(e){invisible()})

compute.adm.angles <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c('s1','s2','adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    mds <- suppressWarnings(cmdscale(as.matrix(d), k=2))
    ids.s1 <- indivs$ids[indivs$pop == 's1']
    ids.s2 <- indivs$ids[indivs$pop == 's2']
    centroid.s1 <- colMeans(mds[ids.s1,])
    centroid.s2 <- colMeans(mds[ids.s2,])
    ids.adm <- indivs$ids[indivs$pop == 'adm']
    adm.angles <- rep(0, length(ids.adm))
    for (i in 1:length(adm.angles)){
        idx.adm <- which(rownames(mds) == ids.adm[i])
        vec.a <- mds[idx.adm,] - centroid.s1
        vec.b <- mds[idx.adm,] - centroid.s2
        cosine <- sum(vec.a * vec.b) / (sqrt(sum(vec.a*vec.a)) * sqrt(sum(vec.b*vec.b)))
        adm.angles[i] <- acos(cosine)
    }
    final.sumstats['mean.adm.angles'] <- mean(adm.angles)
    final.sumstats['var.adm.angles'] <- var(adm.angles)
    final.sumstats['skew.adm.angles'] <- skewness(adm.angles)
    final.sumstats['kurt.adm.angles'] <- kurtosis(adm.angles)
    tmp <- density(adm.angles)
    final.sumstats['mode.adm.angles'] <- tmp$x[which.max(tmp$y)][1]
    tmp <- sprintf('perc%d.adm.angles', seq(0,100,10))
    final.sumstats[tmp] <- quantile(adm.angles, seq(0,1,.1))    
    final.sumstats
}


tryCatch(final.sumstats <- compute.adm.angles(), error=function(e){invisible()}, warning=function(e){invisible()})

write.table(t(final.sumstats), 'final_sumstats.txt', quote=FALSE, row.names=FALSE)
