#' Compare regulons
#'
#' Thins funciton compare regulatory networks based on regulon activity
#'
#' @param regulon1 Regulon to be compare agains the reference regulon \code{regulon2}
#' @param regulon2 Regulon to use as reference to compare \code{regulon1}
#' @param dset Numeric matrix or list of numeric matrixes containing the datasets, with genes in rows and samples in columns
#' @param subreg1 Single regulon or list of two regulons to use as positive conservation control for \code{regulon1}, usually networks obtained from context-matching independent datasets. If a single regulon is given, then \code{regulon1} will be used as second regulon 
#' @param subreg2 Single regulon or list of two regulons to use as positive conservation control for \code{regulon2}, usually networks obtained from context-matching independent datasets. If a single regulon is given, then \code{regulon2} will be used as second regulon
#' @param method Character string indicating the similarity metric, either signatureDistance, pearson, spearman, kendall, euclidean, jaccard, oddsratio or fet
#' @param groups Vector of factors indicating the groups of datasets used during the integration of results
#' @param nn Integer indicating the number of genes used to compute the similarity between regulator's signatures
#' @param minsize Integer indicating the minimum allowed size for the regulons
#' @param ws Number indicating the weight score for \code{signatureDistance}
#' @param rep Integer indicating the number of iterations for the negative controls
#' @return List of six components containing the normalized enrichment scores for the comparison between regulons, for the internal comparion of regulons from network 1, and for the internal comparison of regulons from network 2, the maximum difference in the normalized enrichment score and the associated p-value, a similarity matrix with regulon1 in rows and regulon2 in columns
#' @export
compareRegulon <- function(regulon1, regulon2, dset=NULL, subreg1=NULL, subreg2=NULL, method="pearson", groups=NULL, nn=NULL, minsize=25, ws=2, rep=1) {
  if (is.null(groups)) groups <- 1:length(dset)
  method <- match.arg(method, c("pearson", "spearman", "kendall", "signatureDistance", "euclidean", "jaccard", "oddsratio", "FET"))
  subreg1 <- composeSubreg(regulon1, subreg1)
  subreg2 <- composeSubreg(regulon2, subreg2)
  regul <- list(reg1=regulon1, reg2=regulon2, subreg1=subreg1[[1]], subreg2=subreg1[[2]], subreg3=subreg2[[1]], subreg4=subreg2[[2]])
  rm(regulon1, regulon2, subreg1, subreg2)
  if (method %in% c("jaccard", "oddsratio", "FET")) {
    rr <- lapply(regul[1:2], randomNetwork, rep=rep)
    names(rr) <- c("rr1", "rr2")
    regul <- c(regul, rr)
    rm(rr)
    reg1.reg2 <- discreteSimilarity(regul$reg1, regul$reg2, method=method)
    rr1.rr2 <- discreteSimilarity(regul$rr1, regul$rr2, method=method)
    dnullf <- aecdf(rr1.rr2[is.finite(rr1.rr2)])
    sr1.sr2 <- discreteSimilarity(regul$subreg1, regul$subreg2, method=method)
    sr3.sr4 <- discreteSimilarity(regul$subreg3, regul$subreg4, method=method)
    tmp <- list(regul=reg1.reg2, dnull=rr1.rr2, subreg1=sr1.sr2, subreg2=sr3.sr4, regul.nes=dnullf(reg1.reg2)$nes, subreg1.nes=dnullf(sr1.sr2)$nes, subreg2.nes=dnullf(sr3.sr4)$nes)
    class(tmp) <- "compareRegulon"
    return(tmp)  
  }
  if (!is.list(dset)) dset <- list(dset1=dset)
  dset <- lapply(dset, function(x) t(scale(t(scale(x)))))
  regulators <- table(unlist(lapply(regul, function(regul, dset, minsize) {
    tmp <- apply(sapply(dset, function(dset, regul) sapply(regul, function(x, genes) length(which(names(x$tfmode) %in% genes)), genes=rownames(dset)), regul=regul), 1, min, na.rm=T)
    names(tmp)[tmp>=minsize]
  }, dset=dset, minsize=minsize), use.names=F))
  regulators <- names(regulators)[regulators==6]
  regul <- lapply(regul, function(x, regulators) x[match(regulators, names(x))], regulators=regulators)
  dset <- lapply(dset, function(x, genes) x[rownames(x) %in% genes, ], genes=unique(unlist(lapply(regul, function(x) unique(unlist(lapply(x, function(x) names(x$tfmode)), use.names=F))), use.names=F)))
  rr <- lapply(regul[1:2], randomNetwork, rep=rep)
  names(rr) <- c("rr1", "rr2")
  regul <- c(regul, rr)
  res <- lapply(regul, function(regul, dset) {
    tmp <- lapply(dset, function(dset, regul) ssmarina(dset, regul, method="none", minsize=2, eset.filter=F, nes=F), regul=regul)
    return(tmp)
  }, dset=dset)
  reg1.reg2 <- lapply(1:length(res[[1]]), function(i, res, method, nn, ws) vpProfileDistance(res$reg1[[i]], res$reg2[[i]], method=method, nn=nn, ws=ws), res=res, method=method, nn=nn, ws=ws)
  rr1.rr2 <- lapply(1:length(res[[1]]), function(i, res, method, nn, ws) vpProfileDistance(res$rr1[[i]], res$rr2[[i]], method=method, nn=nn, ws=ws), res=res, method=method, nn=nn, ws=ws)
  dnullf <- lapply(rr1.rr2, aecdf)
  sr1.sr2 <- lapply(1:length(res[[1]]), function(i, res, method, nn, ws) vpProfileDistance(res$subreg1[[i]], res$subreg2[[i]], method=method, nn=nn, ws=ws), res=res, method=method, nn=nn, ws=ws)
  sr3.sr4 <- lapply(1:length(res[[1]]), function(i, res, method, nn, ws) vpProfileDistance(res$subreg3[[i]], res$subreg4[[i]], method=method, nn=nn, ws=ws), res=res, method=method, nn=nn, ws=ws)
  tt <- sapply(1:length(reg1.reg2), function(i, x, dnull) dnull[[i]](x[[i]])$nes, x=reg1.reg2, dnull=dnullf)
  ttc1 <- sapply(1:length(sr1.sr2), function(i, x, dnull) dnull[[i]](x[[i]])$nes, x=sr1.sr2, dnull=dnullf)
  ttc2 <- sapply(1:length(sr3.sr4), function(i, x, dnull) dnull[[i]](x[[i]])$nes, x=sr3.sr4, dnull=dnullf)
  tmp <- c(lapply(list(regul=reg1.reg2, dnull=rr1.rr2, subreg1=sr1.sr2, subreg2=sr3.sr4), function(x, groups) groupMeans(sapply(x, function(x) x), groups), groups=groups), lapply(list(regul.nes=tt, subreg1.nes=ttc1, subreg2.nes=ttc2), function(x, groups) groupStouffer(x, groups), groups=groups))
  class(tmp) <- "compareRegulon"
  return(tmp)
}

groupMeans <- function(x, groups) {
  x1 <- sapply(unique(groups), function(i, x, groups) rowMeans(filterColMatrix(x, groups==i)), x=x, groups=groups)
  if (is.null(ncol(x1))) return(x1)
  return(rowMeans(x1))
}

groupStouffer <- function(x, groups) {
  x1 <- sapply(unique(groups), function(i, x, groups) rowSums(filterColMatrix(x, groups==i))/sqrt(length(which(groups==i))), x=x, groups=groups)
  if (is.null(ncol(x1))) return(x1)
  return(rowSums(x1)/sqrt(ncol(x1)))
}

#' Plot the regulon comparison analysis
#' 
#' This function plots the regulon comparison analysis results
#' 
#' @param x compareRegulon object
#' @param type Character string indicating the type of plot, either 'raw' or 'zscore'
#' @param annot Vector of 2 character strings indicating the name of the networks to compare
#' @param pval Numberr indicating the p-value threshold for Zscore type
#' @param color Vecor of two character strings indicating the color for the positive controls, no positive ontrol is ploted when color is set to NULL
#' @param legend Character string indicating the position of the legend, which is ommited when legend is set to none
#' @param ... Additional parameters to pass to the generic plot function
#' @return Nothing, a plot is generated in the default output device
#' @method plot compareRegulon
#' @S3method plot compareRegulon
plot.compareRegulon <- function(x, type="raw", annot=c("Ntw1", "Ntw2"), pval=.05, color=c("blue", "red"), legend="topleft", ...) {
  switch(match.arg(type, c("raw", "zscore")),
         raw=plotCompareRegulonRaw(x=x, annot=annot, color=color, legend, ...),
         zscore=plotCompareRegulonZscore(x=x, pval=pval, annot=annot, color=color, legend, ...))
}

#' Table of results for regulon comparison analysis
#' 
#' This function generates a table of results for the regulon comparison analysis
#' 
#' @param object compareRegulon object
#' @param pval Number indicating the p-value cutoff for regulon reliability
#' @param annot Vector of 2 character strings indicating the names of the networks
#' @param ... Given for compatibility with the summary generic function
#' @return Data.frame of results
#' @method summary compareRegulon
#' @S3method summary compareRegulon

summary.compareRegulon <- function(object, pval=1, annot=c("Ntw1", "Ntw2"), ...) {
  qval <- qnorm(pval, lower.tail=F)
  pos <- object$subreg1.nes>qval & object$subreg2.nes>qval
  z3 <- object$regul.nes[pos]
  z1=object$subreg1.nes[pos]
  z2=object$subreg2.nes[pos]
  pos <- order(z3, decreasing=T)
  tmp <- data.frame(Gene=names(z3)[pos], round(cbind(z3, z1, z2)[pos, ], 2), signif(pnorm(abs(z3[pos]), lower.tail=F)*2, 3), signif(pnorm(cbind(z1, z2)[pos, ], lower.tail=F), 3))
  colnames(tmp) <- c("Gene", paste(c(paste(annot, collapse="-"), annot), " (z)", sep=""), paste(c(paste(annot, collapse="-"), annot), " (p)", sep=""))
  return(tmp)
}

plotCompareRegulonRaw <- function(x, annot, color=c("blue", "red"), legend, ...) {
  tmp <- lapply(list(regul=x$regul, dnull=x$dnull, subreg1=x$subreg1, subreg2=x$subreg2), function(x) density(x, from=min(x), to=1))
  xl <- ifelse(min(x$dnull)<0, -1, 0)
  yl <- max(unlist(lapply(tmp, function(x) max(x$y))))
  plot(tmp$dnull, col="grey", xlim=c(xl, 1), ylim=c(0, yl*1.05), xlab="Correlation coefficient", main="", ...)
  lines(tmp$subreg1, col=color[1])
  lines(tmp$subreg2, col=color[2])
  lines(tmp$regul, lwd=2, col="black")
  if (legend != "none") legend(legend, c("Random network", paste(annot, collapse=" to "), paste(annot[1], " (+) ctrl", sep=""), paste(annot[2], " (+) ctrl", sep="")), col=c("grey", "black", color), lwd=5, cex=.8)
}

plotCompareRegulonZscore <- function(x, pval, annot, color=c("blue", "red"), legend, ...) {
  qval <- qnorm(pval/2, lower.tail=F)
  pos <- x$subreg1.nes>qval & x$subreg2.nes>qval
  if (min(x$regul.nes) < 0) tmp <- c(lapply(list(x$subreg1.nes[pos], x$subreg2.nes[pos]), function(x) density(x, from=min(x))), list(density(x$regul.nes[pos])))
  else tmp <- lapply(list(x$subreg1.nes[pos], x$subreg2.nes[pos], x$regul.nes[pos]), function(x) density(x, from=range(x)[1]))
  if (is.null(color)) tmp <- tmp[3]
  xlim <- range(sapply(tmp, function(x) x$x))
  ylim <- range(sapply(tmp, function(x) x$y))
  plot(0, 0, type="n", ylim=c(0, ylim[2]*1.05), xlim=xlim, xlab="Z-Score", ylab="Density", main="", ...)
  abline(h=0, col="grey")
  if (min(x$regul.nes, na.rm=T) < 0) lines(c(-qval, -qval), c(0, approx(x=tmp[[length(tmp)]]$x, y=tmp[[length(tmp)]]$y, xout=(-qval))$y), col="grey")
  lines(c(qval, qval), c(0, approx(x=tmp[[length(tmp)]]$x, y=tmp[[length(tmp)]]$y, xout=qval)$y), col="grey")
  if (length(tmp)==3) {
    lines(tmp[[1]], col=color[1])
    lines(tmp[[2]], col=color[2])
  }
  lines(tmp[[length(tmp)]], lwd=2)
  if (length(tmp)==3 & legend != "none") legend(legend, c(paste(annot, collapse=" to "), paste(annot[1], " (+) ctrl", sep=""), paste(annot[2], " (+) ctrl", sep="")), col=c("black", color), lwd=4, cex=.8)
  if (min(x$regul.nes) < 0) text(c(-qval, 0, qval)*2, rep(max(tmp[[length(tmp)]]$y)*1.4, 3), paste(round(c(length(which(x$regul.nes[pos]<(-qval))), length(which(x$regul.nes[pos]>=(-qval) & x$regul.nes[pos] <= qval)), length(which(x$regul.nes[pos]>qval)))*100/length(which(pos))), "%", sep=""))
  else text(c(qval/2, qval*2), rep(max(tmp[[length(tmp)]]$y)*1.4, 2), paste(round(c(length(which(x$regul.nes[pos]>=(-qval) & x$regul.nes[pos] <= qval)), length(which(x$regul.nes[pos]>qval)))*100/length(which(pos))), "%", sep=""))
}

composeSubreg <- function(reg, subreg) {
  if (is.null(subreg)) {
    pos <- lapply(reg, function(x) {tmp <- sample(length(x$tfmode)); nn <- floor(length(tmp)/2); return(list(tmp[1:nn], tmp[(nn+1):(nn*2)]))})
    reg1 <- lapply(1:length(reg), function(i, reg, pos) list(tfmode=reg[[i]]$tfmode[pos[[i]][[1]]], likelihood=reg[[i]]$likelihood[pos[[i]][[1]]]), reg=reg, pos=pos)
    reg2 <- lapply(1:length(reg), function(i, reg, pos) list(tfmode=reg[[i]]$tfmode[pos[[i]][[2]]], likelihood=reg[[i]]$likelihood[pos[[i]][[2]]]), reg=reg, pos=pos)
    names(reg1) <- names(reg2) <- names(reg)
    return(list(reg1, reg2))
  }
  if (length(subreg)==2) return(subreg)
  if (length(subreg)==1) subreg <- subreg[[1]]
  return(list(reg, subreg))
}
  
# signatureDistance, pearson, spearman, kendall, or euclidean
vpProfileDistance <- function(d1, d2, method, nn=nn, ws=ws, diagonal=T) {
  switch(match.arg(method, c("signatureDistance", "pearson", "spearman", "kendall", "euclidean")),
         signatureDistance={
           tmp <- scale(signatureDistance(t(rbind(d1, d2)), nn=nn, ws=ws))
           pos <- which(duplicated(colnames(tmp)))
           tmp <- tmp[1:length(pos), ][, pos]
         },
         pearson={tmp <- cor(t(d1), t(d2))},
         spearman={tmp <- cor(t(d1), t(d2), method="spearman")},
         kendall={tmp <- cor(t(d1), t(d2), method="kendall")},
         euclidean={stop("Euclidean similarity is still not implemented")}
  )
  if (diagonal) return(diag(tmp))
  return(tmp)
}

# networkSimilarity, FET, odds ratio and Jaccard
discreteSimilarity <- function(regulon1, regulon2, method="jaccard") {
  method <- match.arg(method, c("jaccard", "oddsratio", "FET"))
  regulators <- intersect(names(regulon1), names(regulon2))
  regulon1 <- regulon1[match(regulators, names(regulon1))]
  regulon2 <- regulon2[match(regulators, names(regulon2))]
  regulon1 <- lapply(regulon1, function(x) names(x$tfmode))
  regulon2 <- lapply(regulon2, function(x) names(x$tfmode))
  pb <- txtProgressBar(max=length(regulon1), style=3)
  switch(method,
         jaccard={
           dis <- unlist(lapply(1:length(regulon1), function(i, regulon1, regulon2, bkg, pb) {setTxtProgressBar(pb, i); jaccard(bkg %in% regulon1[[i]], bkg %in% regulon2[[i]])}, regulon1=regulon1, regulon2=regulon2, bkg=unique(unlist(c(regulon1, regulon2), use.names=F)), pb=pb), use.names=F)
           names(dis) <- names(regulon1)
         },
         oddsratio={
           dis <- unlist(lapply(1:length(regulon1), function(i, regulon1, regulon2, bkg, pb) {
             setTxtProgressBar(pb, i)
             tmp <- table(bkg %in% regulon1[[i]], bkg %in% regulon2[[i]])
             return(tmp[1, 1]*tmp[2, 2]/tmp[1, 2]/tmp[2, 1])
           }, pb=pb, regulon1=regulon1, regulon2=regulon2, bkg=unique(unlist(c(regulon1, regulon2), use.names=F))), use.names=F)
           names(dis) <- names(regulon1)
           dis <- log2(dis)
         },
         FET={
           dis <- unlist(lapply(1:length(regulon1), function(i, regulon1, regulon2, bkg, pb) {
             setTxtProgressBar(pb, i)
             tmp <- fisher.test(bkg %in% regulon1[[i]], bkg %in% regulon2[[i]], alternative="greater")$p.value
             return(tmp)
           }, pb=pb, regulon1=regulon1, regulon2=regulon2, bkg=unique(unlist(c(regulon1, regulon2), use.names=F))), use.names=F)
           names(dis) <- names(regulon1)
           dis <- qnorm(dis, lower.tail=F)
         })
  return(dis)
}


randomNetwork <- function(regul, rep=1) {
  reg1 <- NULL
  for (i in 1:rep) {
    reg1 <- c(reg1, lapply(regul, function(reg, targets) {
      tfmode <- reg$tfmode
      names(tfmode) <- sample(targets, length(tfmode))
      list(tfmode=tfmode, likelihood=reg$likelihood)
    }, targets=unique(unlist(lapply(regul, function(x) names(x$tfmode)), use.names=F))))
  }
  class(reg1) <- "regulon"
  return(reg1)
}
