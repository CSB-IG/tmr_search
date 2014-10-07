#' MARINA
#'
#' This function performs MAster Regulator INference Analysis
#'
#' @param ges Vector containing the gene expression signature to analyze, or matrix with columns containing bootstraped signatures
#' @param regulon Object of class regulon
#' @param nullmodel Matrix of genes by permutations containing the NULL model signatures. A parametric approach equivalent to shuffle genes will be used if nullmodel is ommitted.
#' @param minsize Number indicating the minimum allowed size for the regulons
#' @param adaptive.size Logical, whether the weight (likelihood) should be used for computing the regulon size
#' @param iterative Logical, whether a two step analysis with adaptive redundancy estimation should be performed
#' @param ges.filter Logical, whether the gene expression signature should be limited to the genes represented in the interactome
#' @param synergy Number indicating the synergy computation mode: (0) for no synergy computation; (0-1) for establishing the p-value cutoff for individual TFs to be included in the synergy analysis; (>1) number of top TFs to be included in the synergy analysis
#' @param level Integer, maximum level of combinatorial regulation
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A MARINA object containing the following components:
#' \describe{
#' \item{signature}{The gene expression signature}
#' \item{regulon}{The final regulon object used}
#' \item{es}{Enrichment analysis results including regulon size, normalized enrichment score and p-value}
#' \item{param}{MARINa parameters, including \code{minsize}, \code{adaptive.size}, \code{iterative}}
#' }
#' @seealso \code{\link{ssmarina}}
#' @export

marina <- function(ges, regulon, nullmodel=NULL, minsize=25, adaptive.size=F, iterative=F, ges.filter=T, synergy=0, level=10, verbose=T) {
# Cleaning the parameters
	if (is.vector(ges)) ges <- matrix(ges, length(ges), 1, dimnames=list(names(ges), NULL))
	regulon <- updateRegulon(regulon)
	if (ges.filter) ges <- filterRowMatrix(ges, rownames(ges) %in% unique(c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names=F))))
	regulon <- lapply(regulon, function(x, genes) {filtro <- names(x$tfmode) %in% genes; x$tfmode <- x$tfmode[filtro]; x$likelihood <- x$likelihood[filtro]; return(x)}, genes=rownames(ges))
	if (!is.null(nullmodel)) nullmodel <- nullmodel[match(rownames(ges), rownames(nullmodel)), ]
# Running MARINa
	if (iterative) {
		tmp <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names=F))
		tw <- rep(1, length(tmp))
		names(tw) <- tmp
		rm(tmp)
		if (adaptive.size) regulon <- regulon[sapply(regulon, function(x, tw) {t1 <- tw[match(names(x$tfmode), names(tw))]; sum(x$likelihood*t1/max(x$likelihood*t1))}, tw=tw)>=minsize]
		else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=minsize]
		if (is.null(nullmodel)) nullf <- ppwea3NULLf(regulon, tw=tw, lw=1)
		else nullf <- pwea3NULLf(pwea3NULLgroups(nullmodel, regulon, tw=tw, lw=1, simsig=T, verbose=verbose), verbose=verbose)
		res <- groupPwea3(ges, regulon, nullf, minsize=1, tw=tw, lw=1, simsig=T, verbose=verbose)
		tmp1 <- cbind(rep(names(regulon), sapply(regulon, function(x) length(x$tfmode))), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names=F))
		tmp1 <- cbind(tmp1, abs(res$nes)[match(tmp1[, 1], names(res$nes))])
		tw <- tapply(tmp1[, 3], tmp1[, 2], function(x) sum(as.numeric(x)))
		tw <- 1-tw/max(tw)
		if (adaptive.size) regulon <- regulon[sapply(regulon, function(x, tw) {t1 <- tw[match(names(x$tfmode), names(tw))]; sum(x$likelihood*t1/max(x$likelihood*t1))}, tw=tw)>=minsize]
		else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=minsize]
		if (is.null(nullmodel)) nullf <- ppwea3NULLf(regulon, tw=tw, lw=1)
		else nullf <- pwea3NULLf(pwea3NULLgroups(nullmodel, regulon, tw=tw, lw=1, simsig=T, verbose=verbose), verbose=verbose)
		res <- groupPwea3(ges, regulon, nullf, minsize=1, tw=tw, lw=1, simsig=T, verbose=verbose)
	}
	else{
		tw <- 1/table(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = F))^.5
		if (adaptive.size) regulon <- regulon[sapply(regulon, function(x, tw) {t1 <- tw[match(names(x$tfmode), names(tw))]; sum(x$likelihood*t1/max(x$likelihood*t1))}, tw=tw)>=minsize]
		else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=minsize]
		if (is.null(nullmodel)) nullf <- ppwea3NULLf(regulon, tw=tw, lw=1)
		else nullf <- pwea3NULLf(pwea3NULLgroups(nullmodel, regulon, tw=tw, lw=1, simsig=T, verbose=verbose), verbose=verbose)
		res <- groupPwea3(ges, regulon, nullf, minsize=1, tw=tw, lw=1, simsig=T, verbose=verbose)
	}
	res <- list(signature=ges, regulon=regulon, es=res, param=list(minsize=minsize, adaptive.size=adaptive.size, iterative=iterative), nullmodel=nullmodel)
  class(res) <- "marina"
	if (synergy>0) res <- marinaSynergy(marinaCombinatorial(mobj=res, regulators=synergy, nullmodel=nullmodel, minsize=minsize, adaptive.size=adaptive.size, iterative=iterative, level=level, verbose=verbose), verbose=verbose)
	return(res)
}

#' MARINA bootstraps integration
#' 
#' This function integrates the bootstrap MARINA results
#' 
#' @param mobj MARINA object
#' @param method Character string indicating the method to use, either mean, median or mode
#' @return MARINA object
#' @seealso \code{\link{marina}}
#' @export
bootstrapMarina <- function(mobj, method="mode") {
  if (ncol(mobj$es$nes.bt)>1) switch(match.arg(method, c("mean", "median", "mode")),
         mean={mobj$es$nes <- rowMeans(mobj$es$nes.bt)},
         median={mobj$es$nes <- apply(mobj$es$nes.bt, 1, median)},
         mode={mobj$es$nes <- apply(mobj$es$nes.bt, 1, distMode)})
  mobj$es$p.value <- pnorm(abs(mobj$es$nes), lower.tail=F)*2
  return(mobj)
}

#' MARINA combinatorial analysis
#'
#' This function performs combinatorial analysis for MARINA objects
#'
#' @param mobj MARINa object generated by \code{marina} function
#' @param regulators Either a number between 0 and 1 indicating the p-value cutoff for individual TFs to be included in the combinations analysis; (>1) indicating the number of top TFs to be included in the combinations analysis; or a vector of character strings indicating the TF IDs to be included in the analysis
#' @param nullmodel Matrix of genes by permutations containing the NULL model signatures. Taken from \code{mobj} by default
#' @param minsize Number indicating the minimum allowed size for the regulons, taken from \code{mobj} by default
#' @param adaptive.size Logical, whether the weight (likelihood) should be used for computing the size, taken from \code{mobj} by default
#' @param iterative Logical, whether a two step analysis with adaptive redundancy estimation should be performed, taken from \code{mobj} by default
#' @param level Integer, maximum level of combinatorial regulation
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A MARINA object
#' @seealso \code{\link{marina}}
#' @export
marinaCombinatorial <- function(mobj, regulators=100, nullmodel=NULL, minsize=NULL, adaptive.size=NULL, iterative=NULL, level=10, verbose=T) {
  synergy <- regulators
  if (is.null(minsize)) minsize <- mobj$param$minsize
	if (is.null(adaptive.size)) adaptive.size <- mobj$param$adaptive.size
	if (is.null(iterative)) iterative <- mobj$param$iterative
	if (is.null(nullmodel)) nullmodel <- mobj$nullmodel
  if (length(synergy)==1 & synergy[1]==0) return(mobj)
	res <- mobj$es
	regulon <- mobj$regulon
	if (iterative) {
		tmp1 <- cbind(rep(names(regulon), sapply(regulon, function(x) length(x$tfmode))), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names=F))
		tmp1 <- cbind(tmp1, abs(res$nes)[match(tmp1[, 1], names(res$nes))])
		tw <- tapply(tmp1[, 3], tmp1[, 2], function(x) sum(as.numeric(x)))
		tw <- 1-tw/max(tw)
	}
	else tw <- 1/table(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = F))^.5
  tfs <- synergy
  if (length(synergy)==1) {
    if (synergy<1) tfs <- names(res$p.value)[res$p.value<synergy]
	  else tfs <- names(res$p.value)[order(res$p.value)[1:synergy]]
  }
  res1 <- list(es=res)
	resp <- res$nes
	res <- list(res)
	regul <- NULL
	n=2
	while(length(tfs)>n & n<=level) {
		if (verbose) cat("\n\n-------------------------------------------\nComputing synergy for combination of ", n, " TFs\n-------------------------------------------\n\n", sep="")
		res1 <- comregulationAnalysis(combn(tfs, n), mobj$signature, regulon, nullmodel, resp, res1$es$p.value, minsize, adaptive.size, tw=tw, verbose=verbose)
		if(length(res1)==0) break
		res <- c(res, list(res1$es))
		regul <- c(regul, res1$regul)
		n <- n+1
		tfs <- unique(unlist(strsplit(as.character(names(res1$es$nes)), "--"), use.names=F))
	}	
	names(res) <- 1:length(res)
	res <- res[sapply(res, function(x) length(x$nes))>0]
	tmp <- names(res[[1]])
	res <- lapply(tmp, function(nom, res) {
		tmp <- unlist(lapply(res, function(x, nom) x[[nom]], nom=nom), use.names=F)
		names(tmp) <- unlist(lapply(res, function(x, nom) names(x[[nom]]), nom=nom), use.names=F)
		return(tmp)
	}, res=res)
	names(res) <- tmp
	regulon <- c(regulon, regul)
	mobj$regulon <- regulon
	mobj$es <- res
  class(mobj) <- "marina"
	return(mobj)
}

#' MARINa synergy analysis
#'
#' This function performs a synergy analysis for combinatorial regulation
#'
#' @param mobj MARINa object containing combinatorial regulation results generated by \code{marinaCombinatorial}
#' @param per Integer indicating the number of permutations
#' @param seed Integer indicating the seed for the permutations, 0 for disable it
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return Updated MARINA object containing the sygergy p-value
#' @seealso \code{\link{marina}}
#' @export

marinaSynergy <- function(mobj, per=1000, seed=1, verbose=T) {
  if (seed>0) set.seed(round(seed))
	pos <- which(sapply(strsplit(names(mobj$regulon), "--"), length)>1)
  pb <- NULL
  if (verbose) {
    cat("\nComputing synergy statistics for ", length(pos), " co-regulons.\n", sep="")
    cat("Process started at ", date(), "\n", sep="")
    pb <- txtProgressBar(max=length(pos), style=3)
  }
	res <- sapply(1:length(pos), function(i, pos, ss, nes, regul, pb, verbose) {
    if (verbose) setTxtProgressBar(pb, i)
    i <- pos[i]
    pos <- match(unlist(strsplit(names(regul)[i], "--")), names(regul))
		tgs <- unique(unlist(lapply(regul[pos], function(x) names(x$tfmode)), use.names=F))
		tfmode <- sapply(regul[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs)
		rownames(tfmode) <- tgs
		tfmode <- t(t(tfmode)*sign(nes[match(names(regul)[pos], names(nes))]))
		tfmode[is.na(tfmode)] <- 0
		ll <- sapply(regul[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs)
		ll[is.na(ll)] <- 0
		ll1 <- ll/rowSums(ll)
		tfmode <- rowSums(tfmode*ll1)
		ll <- apply(ll, 1, max)
		tgs <- names(regul[[i]]$tfmode)
    r2 <- apply(ss, 2, rank)/(nrow(ss)+1)*2-1
    r1 <- abs(r2)*2-1
    r1[r1==(-1)] <- 1-(1/length(r1))
    r1 <- qnorm(r1/2+.5)
    r2 <- qnorm(r2/2+.5)
    dim(r1) <- dim(r2) <- dim(ss)
		rownames(r1) <- rownames(r2) <- rownames(ss)
    pos <- match(names(tfmode), rownames(r1))
    r1 <- filterRowMatrix(r1, pos)
    r2 <- filterRowMatrix(r2, pos)
		dnull <- sapply(1:per, function(i, ll, samp) {
			ll[sample(length(ll), samp)] <- 0
			return(ll)
		}, ll=ll, samp=length(tgs))
		dnull <- cbind(ll, dnull)
		dnull[match(tgs, names(tfmode)), 1] <- 0
		dnull <- t(t(dnull)/colSums(dnull))		
		sum1 <- t(r2) %*% (dnull*tfmode)
    sum2 <- t(r1) %*% (dnull * (1-abs(tfmode)))
    es <- colMeans(abs(sum1)+sum2*(sum2>0))
		x <- es[1]
		es <- es[-1]
    iqr <- quantile(es, c(.5, 5/length(es)))
    pd <- ecdf(es)
    a <- list(x=knots(pd), y=pd(knots(pd)))
    filtro <- a$x<iqr[1] & a$x>=iqr[2] & a$y<1
    spl <- smooth.spline(a$x[filtro], -log(a$y[filtro]), spar=.75)
    p <- exp(-predict(spl, x)$y)
    pos <- which(x>iqr[1])
    if (x>iqr[1]) p <- pd(x)
		return(p)
	}, ss=mobj$signature, nes=mobj$es$nes, regul=mobj$regulon, pos=pos, pb=pb, verbose=verbose)
  if (verbose) cat("\nProcess ended at ", date(), "\n", sep="")
	names(res) <- names(mobj$regulon)[pos]
	mobj$es$synergy <- res
  class(mobj) <- "marina"
	return(mobj)	
}


comregulationAnalysis <- function(tfs, ges, regulon, nullmodel=NULL, ones, resp, minsize=5, adaptive.size=F, tw=.5, verbose=T) {
	reg1 <- apply(tfs, 2, function(x, regulon, res1) {
		pos <- which(names(regulon) %in% x)
		tgs <- table(unlist(lapply(regulon[pos], function(x) names(x$tfmode)), use.names=F))
		tgs <- names(tgs)[tgs==length(x)]
		if (length(tgs)<2) return(list(tfmode=NULL, likelihood=NULL))
		list(tfmode=frmeanna(t(t(sapply(regulon[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs))*sign(res1[pos])))[, 1], likelihood=apply(sapply(regulon[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs), 1, prod))
	}, regulon=regulon, res1=ones)
	names(reg1) <- apply(tfs, 2, paste, collapse="--")
	reg1 <- reg1[sapply(reg1, function(x) length(x$tfmode))>0]
	if (adaptive.size) reg1 <- reg1[sapply(reg1, function(x, tw) {t1 <- tw[match(names(x$tfmode), names(tw))]; sum(x$likelihood*t1/max(x$likelihood*t1))}, tw=tw)>=minsize]
	else reg1 <- reg1[sapply(reg1, function(x) length(x$tfmode))>=minsize]
	if (length(reg1)==0) return(list())
	if (is.null(nullmodel)) nullf <- ppwea3NULLf(reg1, tw=tw, lw=1)
	else nullf <- pwea3NULLf(pwea3NULLgroups(nullmodel, reg1, tw=tw, lw=1, simsig=T, verbose=verbose))
  res1 <- groupPwea3(ges, reg1, nullf, minsize=1, tw=tw, lw=1, simsig=T, verbose=verbose)
	filtro <- sapply(strsplit(names(res1$p.value), "--"), function(x, res1, res) {
		test <- res[sapply(strsplit(names(res), "--"), function(x1, x) all(x1 %in% x), x=x)]
	  if (length(test)>0) return(all(test > res1[paste(x, collapse="--")]))
	  return(F)
	}, res1=res1$p.value, res=resp)
  tmp <- lapply(res1, function(x, filtro) ifelse(is.null(ncol(x)), return(x[filtro]), return(filterRowMatrix(x, filtro))), filtro=filtro)
	return(list(es=tmp, regul=reg1[names(reg1) %in% names(tmp$nes)]))
}

selectTFs <- function(res, n) {
	tfs <- unlist(lapply(strsplit(names(res$nes), "--"), function(x, tfs, n) {
		tfs <- unique(unlist(tfs[sapply(tfs, function(tfs, x) any(tfs %in% x), x=x)], use.names=F))
		if (length(tfs)<n) return(NULL)
		return(combn(tfs, n))
	}, tfs=strsplit(names(res$nes), "--"), n=n), use.names=F)
	tfs <- matrix(tfs, n, length(tfs)/n)
	tfs <- tfs[, !duplicated(apply(tfs, 2, paste, collapse="--"))]
	i <- 1
	while(i<ncol(tfs)) {
		filtro <- apply(tfs, 2, function(x, pat) all(x %in% pat), pat=tfs[, i])
		filtro[i] <- F
		tfs <- filterColMatrix(tfs, !filtro)
		i <- i+1
	}
	return(tfs)
}

#' @method print marina
#' @S3method print marina
print.marina <- function(x, ...) cat("Object of class marina with ", length(x$es$nes), " regulators.\n", sep="")

#' List MARINA results
#' 
#' This function generates a table of MARINA results
#' 
#' @param object MARINA object
#' @param mrs Either number of top MRs to report or vector containing the genes to display
#' @param ... Given for compatibility with the summary generic function
#' @return Data.frame with results
#' @method summary marina
#' @S3method summary marina

summary.marina <- function(object, mrs=10, ...) {
  if (length(mrs)==1) {
    mrs <- names(object$es$nes)[order(object$es$p.value)[1:round(mrs)]]
    mrs <- mrs[order(object$es$nes[match(mrs, names(object$es$nes))], decreasing=T)]
  }
  pos <- match(mrs, names(object$es$nes))
  tmp <- data.frame(Regulon=mrs, Size=object$es$size[pos], NES=round(object$es$nes[pos], 2), p.value=signif(object$es$p.value[pos], 3), FDR=signif(p.adjust(object$es$p.value, "fdr")[pos], 3))
  if (!is.null(object$es$synergy)) {
    synp <- object$es$synergy[match(tmp$Regulon, names(object$es$synergy))]
    tmp <- cbind(tmp, Synergy=signif(synp, 3))
  }
  if (!is.null(object$ledge)) tmp <- cbind(tmp, Ledge=sapply(object$ledge[match(tmp$Regulon, names(object$ledge))], function(x) {
    if (length(x)<6) return(paste(x, collapse=", "))
    return(paste(paste(x[1:4], collapse=", "), ", + ", length(x)-4, " genes", sep=""))
  }))
  if (!is.null(object$shadow)) {
    tmps <- apply(object$shadow, 1, paste, collapse=" -> ")
    tmp <- list("MARINA.results"=tmp, "Shadow.pairs"=tmps)
  }
  return(tmp)
}
