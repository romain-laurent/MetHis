#!/usr/bin/env Rscript

library(e1071)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx.simu <- args[2]
idx.gen <- args[3]

path <- sprintf('%s/simu_%s/', prefix, idx.simu)

final.sumstats <- c(Fst.s1.adm = NA, Fst.s2.adm = NA, Fst.s1.s2 = NA, mean.ASD.s1.adm = NA, var.ASD.s1.adm = NA, mean.ASD.s2.adm = NA, var.ASD.s2.adm = NA, mean.ASD.s1.s2 = NA, var.ASD.s1.s2 = NA, f3 = NA, mean.ASD.s1 = NA, var.ASD.s1 = NA, mean.ASD.s2 = NA, var.ASD.s2 = NA, mean.ASD.adm = NA, var.ASD.adm = NA, mean.het.s1 = NA, var.het.s1 = NA, mean.het.s2 = NA, var.het.s2 = NA, mean.het.adm = NA, var.het.adm = NA, mean.F.s1 = NA, var.F.s1 = NA, mean.F.s2 = NA, var.F.s2 = NA, mean.F.adm = NA, var.F.adm = NA, mean.adm.props=NA, var.adm.props=NA, skew.adm.props=NA, kurt.adm.props=NA, mode.adm.props=NA)
tmp <- sprintf('perc%d.adm.props', seq(0,100,10))
final.sumstats[tmp] <- NA

## pairwise Fst computed by vcftools
get.Fsts <- function(path){
    d <- read.table(sprintf('%s/result.fst', path))
    final.sumstats['Fst.s1.adm'] <- d[1,]
    final.sumstats['Fst.s2.adm'] <- d[2,]
    final.sumstats['Fst.s1.s2'] <- d[3,]
    final.sumstats
}


## we use 'tryCatch' tricks to prevent problems if files do not exist
tryCatch(final.sumstats <- get.Fsts(path), error=function(e){invisible()}, warning=function(w){invisible()})



## F3 statistic, following Patterson 2012
## does not give exactly the same result as qp3Pop (admixtools) for some reason...
compute.F3 <- function(path, idx.simu, idx.gen){
    ## get allelic frequencies computed by vcftools
    p.adm <- read.table(sprintf('%s/simu_%s_g%s_adm.frq', path, idx.simu, idx.gen), header=TRUE, skip=1)[,6]
    p.s1 <- read.table(sprintf('%s/simu_%s_g%s_s1.frq', path, idx.simu, idx.gen), header=TRUE, skip=1)[,6]
    p.s2 <- read.table(sprintf('%s/simu_%s_g%s_s2.frq', path, idx.simu, idx.gen), header=TRUE, skip=1)[,6]
    d <- data.frame(p.adm, p.s1, p.s2)
    ## compute numerators and denomitors of statistic
    num <- sum((d$p.adm - d$p.s1) * (d$p.adm - d$p.s2))
    denom <- sum(2 * d$p.adm * (1-d$p.adm))
    ## compute statistic
    final.sumstats['f3'] <- num/denom
    final.sumstats
}

tryCatch(final.sumstats <- compute.F3(path, idx.simu, idx.gen), error=function(e){invisible()}, warning=function(w){invisible()})

compute.inbreeding.het <- function(path, idx.simu, idx.gen){
    for (pop in c('s1', 's2', 'adm')){
        ## inbreeding coefficients from vcftools
        d <- read.table(sprintf('%s/%s.het', path, pop), header=TRUE)
        final.sumstats[sprintf('mean.F.%s', pop)] <- mean(d$F)
        final.sumstats[sprintf('var.F.%s', pop)] <- var(d$F)
        ## we compute heterozigosity using allelic frequencies computed by vcftools
        d <- read.table(sprintf('%s/simu_%s_g%s_%s.frq', path, idx.simu, idx.gen, pop), header=TRUE, skip=1)
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

tryCatch(final.sumstats <- compute.inbreeding.het(path, idx.simu, idx.gen), error=function(e){invisible()}, warning=function(w){invisible()})

## compute statistics based on ASD (which was computed by asd)
compute.ASD.stats <- function(path, idx.simu, idx.gen){
    d <- read.table(sprintf('%s/simu_%s_g%s.asd.dist', path, idx.simu, idx.gen), header=TRUE, row.names=1)
    all.pops <- c('s1','s2','adm')
    for (i in 1:3){
        pop1 <- all.pops[i]
        wanted.idxs.1 <- grep(pop1, rownames(d))
        for (j in i:3){
            pop2 <- all.pops[j]
            wanted.idxs.2 <- grep(pop2, rownames(d))
            tmp <- d[wanted.idxs.1, wanted.idxs.2]
            tmp <- tmp[lower.tri(tmp)]
            if (i == j){
                pop <- pop1
                final.sumstats[sprintf('mean.ASD.%s', pop)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s', pop)] <- var(tmp)
            } else {
                final.sumstats[sprintf('mean.ASD.%s.%s', pop1, pop2)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s.%s', pop1, pop2)] <- var(tmp)
            }
        }
    }
    final.sumstats
}

tryCatch(final.sumstats <- compute.ASD.stats(path, idx.simu, idx.gen), error=function(e){invisible()}, warning=function(e){invisible()})

## computes admixture proportions stats using projection in MDS
compute.adm.props <- function(path, idx.simu, idx.gen){
    d <- read.table(sprintf('%s/simu_%s_g%s.asd.dist', path, idx.simu, idx.gen), header=TRUE, row.names=1)
    mds <- suppressWarnings(cmdscale(as.matrix(d), k=2))
    centroid.s1 <- colMeans(mds[grep('s1', rownames(mds)),])
    centroid.s2 <- colMeans(mds[grep('s2', rownames(mds)),])
    centroids.dist <- dist(rbind(centroid.s1,centroid.s2))[1]
    idxs.adm <- grep('adm', rownames(mds))
    adm.props <- rep(0, length(idxs.adm))
    vec.b <- centroid.s2 - centroid.s1
    for (i in 1:length(adm.props)){
        idx.adm <- idxs.adm[i]
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
    final.sumstats['mode.adm.props'] <- tmp$x[which.max(tmp$y)]
    tmp <- sprintf('perc%d.adm.props', seq(0,100,10))
    final.sumstats[tmp] <- quantile(adm.props, seq(0,1,.1))
    final.sumstats
}

tryCatch(final.sumstats <- compute.adm.props(path, idx.simu, idx.gen), error=function(e){invisible()}, warning=function(e){invisible()})

write.table(t(final.sumstats), sprintf('%s/final_sumstats.txt', path), quote=FALSE, row.names=FALSE)
