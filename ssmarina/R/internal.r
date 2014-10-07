#Internal functions for marina package
#Check and update regulons to version 2
updateRegulon <- function(regul) {
  if (is.null(names(regul[[1]]))) return(lapply(regul, function(x) {tmp <- rep(0, length(x)); names(tmp) <- x; list(tfmode=tmp, likelihood=rep(1, length(tmp)))}))
  if (names(regul[[1]])[1]=="tfmode") return(regul)
  return(lapply(regul, function(x) list(tfmode=x, likelihood=rep(1, length(x)))))
}

#' Proportionally Weighted Enrichment Analysis for gene-set groups
#'
#' This function performs a Proportionally Weighted Enrichment Analysis on groups of gene-sets
#'
#' @param rlist Named vector containing the scores to rank the expression profile or matrix where columns contains bootstraped signatures
#' @param groups List of gene-sets (regulons), each component is a list of two vectors: \emph{TFmode} containing the TFMoA index (-1; 1) and \emph{likelihood} containing the interaction relative likelihood
#' @param nullpw Numerical matrix representing the null model, with genes as rows (geneID as rownames) and permutations as columns
#' @param alternative Character string indicating the alternative hypothesis, either two.sided, greater or less
#' @param per Integer indicating the number of permutations for the genes in case ???nullpw??? is ommited
#' @param tnorm Logical, whether the ranking vector should be normally transformed after the copula transformation
#' @param minsize Integer indicating the minimum size for the regulons
#' @param tw Number indicating the exponent for the target pleotropy weight
#' @param lw Number indicating the exponent for the edge likelihood weight
#' @param simsig Logical, whether a symetrical distribution of the signature is used
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A list containing four matrices:
#' \describe{
#' \item{es}{Enrichment score}
#' \item{nes}{Normalized Enrichment Score}
#' \item{size}{Regulon size}
#' \item{p.value}{Enrichment p.value}
#' }

groupPwea3 <- function(rlist, groups, nullpw=NULL, alternative=c("two.sided", "less", "greater")[1], per=0, tnorm=T, minsize=5, tw=1, lw=1, simsig=T, verbose=T) {
  if (is.vector(rlist)) rlist <- matrix(rlist, length(rlist), 1, dimnames=list(names(rlist), NULL))
  if (!is.list(groups) | names(groups)[1]=="tfmode") groups <- list(set=groups)
  groups <- updateRegulon(groups)
  groups <- lapply(groups, function(x, lw) {x$likelihood <- x$likelihood^lw; return(x)}, lw=lw)
  if (length(tw)==1) tw <- 1/table(unlist(lapply(groups, function(x) names(x$tfmode)), use.names=F))^tw
  rlist2 <- apply(rlist, 2, rank)/(nrow(rlist)+1)*2-1
  rlist1 <- abs(rlist2)*2-1
  rlist1[rlist1==(-1)] <- 1-(1/length(rlist1))
  if (tnorm) {
    rlist1 <- qnorm(rlist1/2+.5)
    rlist2 <- qnorm(rlist2/2+.5)
  }
  else {
    rlist1 <- rlist1/sqrt(1/3)
    rlist2 <- rlist2/sqrt(1/3)
  }
  groups <- lapply(groups, function(x, genes) {
    tmp <- names(x$tfmode) %in% genes
    x$tfmode <- x$tfmode[tmp]
    if(length(x$likelihood)==length(tmp)) x$likelihood <- x$likelihood[tmp]
    return(x)
  }, genes=rownames(rlist))
  groups <- groups[sapply(groups, function(x) length(x$tfmode))>=minsize]
  if (is.null(nullpw)) {
    if (per>0) {
			nullpw <- list(eset=sapply(1:per, function(i, rlist) sample(rlist), rlist=rlist[, 1]))
			colnames(nullpw$eset) <- 1:per
			rownames(nullpw$eset) <- rownames(rlist)
		}
		else nullpw <- ppwea3NULLf(groups, tw=tw, lw=1)
  }
  if (names(nullpw)[1]=="eset") nullpw <- pwea3NULLgroups(nullpw, groups, tnorm=tnorm, tw=tw, verbose=verbose)
  if (names(nullpw)[1]=="groups") nullpw <- pwea3NULLf(nullpw, verbose=verbose)
  pb <- NULL
  if (verbose) {
    cat("\nComputing enrichment by PWEA3 on ", length(groups), " gene-sets.\n", sep="")
    cat("Process started at ", date(), "\n", sep="")
    pb <- txtProgressBar(max=length(groups), style=3)
  }
  es <- sapply(1:length(groups), function(i, reg, r1, r2, tw, pb, verbose) {
    reg <- reg[[i]]
    pos <- match(names(reg$tfmode), rownames(r1))
    tw <- tw[match(names(reg$tfmode), names(tw))]
		ll <- reg$likelihood * tw / sum(tw * reg$likelihood)
    sum1 <- t(filterRowMatrix(r2, pos)) %*% matrix(reg$tfmode * ll, length(tw), 1)
    ss <- sign(sum1)
    ss[ss==0] <- 1
    if (verbose) setTxtProgressBar(pb, i)
    sum2 <- t(filterRowMatrix(r1, pos)) %*% matrix(ll * (1-abs(reg$tfmode)), length(tw), 1)
    return((abs(sum1)+(sum2*(sum2>0))) * ss)
  }, reg=groups, r1=rlist1, r2=rlist2, tw=tw, pb=pb, verbose=verbose)
  names(es) <- names(groups)
  if(is.vector(es)) es <- matrix(es, 1, length(es), dimnames=list(NULL, names(es)))
	colnames(es) <- names(groups)
  temp <- lapply(colnames(es), function(x, es, pfun, alter) pfun[[x]](es[, x]), es=es, pfun=nullpw)
  names(temp) <- colnames(es)
  nes <- sapply(temp, function(x) qnorm(x$p.value/2, lower.tail=F))*sign(es)
	senes <- sqrt(fcvarna(nes)/nrow(nes))[, 1]
	if (nrow(nes)==1) senes <- NULL
	nes1 <- colMeans(nes)
  pval1 <- pnorm(abs(nes1), lower.tail=F)*2
  if (verbose) cat("\nProcess ended at ", date(), "\n", sep="")
  return(list(es=colMeans(es), nes=nes1, nes.se=senes, size=sapply(groups, function(x) length(x$tfmode)), p.value=pval1, nes.bt=t(nes)))
}

#' Regulon-specific NULL model
#'
#' This function generates the regulon-specific NULL models
#'
#' @param pwnull Numerical matrix representing the null model, with genes as rows (geneID as rownames) and permutations as columns
#' @param groups List containing the regulons
#' @param tnorm Logical, whether the ranking vector should be normally transformed after the copula transformation
#' @param tw Number indicating the exponent for the target pleotropy weight
#' @param lw Number indicating the exponent for the edge likelihood weight
#' @param simsig Logical, whether a symetrical distribution of the signature is used
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A list containing two elements:
#' \describe{
#' \item{groups}{Regulon-specific NULL model containing the enrichment scores}
#' \item{ss}{Direction of the regulon-specific NULL model}
#' }

pwea3NULLgroups <- function(pwnull, groups, tnorm=T, tw=1, lw=1, simsig=T, verbose=T) {
	if (is.null(pwnull)) return(NULL)
  if (is.matrix(pwnull)) pwnull <- list(eset=pwnull)
  groups <- updateRegulon(groups)
  groups <- lapply(groups, function(x, lw) {x$likelihood <- x$likelihood^lw; return(x)}, lw=lw)
  if (length(tw)==1) tw <- 1/table(unlist(lapply(groups, function(x) names(x$tfmode)), use.names=F))^tw
  if (is.list(pwnull$eset)) {
    t <- unlist(pwnull$eset, use.names=F)
    dim(t) <- c(length(pwnull$eset[[1]]), length(pwnull$eset))
    rownames(t) <- names(pwnull$eset[[1]])
    colnames(t) <- 1:length(pwnull$eset)
  }
  else t <- pwnull$eset
	t2 <- apply(t, 2, rank)/(nrow(t)+1)*2-1
	t1 <- abs(t2)*2-1
	t1[t1==(-1)] <- 1-(1/length(t1))
	if (tnorm) {
	  t1 <- qnorm(t1/2+.5)
	  t2 <- qnorm(t2/2+.5)
	}
	else {
	  t1 <- t1/sqrt(1/3)
	  t2 <- t2/sqrt(1/3)
	}
  pb <- NULL
  if (verbose) {
    cat("\nComputing the null distribution for ", length(groups), " gene-sets.\n", sep="")
    cat("Process started at ", date(), "\n", sep="")
    pb <- txtProgressBar(max=length(groups), style=3)
  }
  temp <- lapply(1:length(groups), function(i, groups, t1, t2, tw, pb, verbose=T) {
    x <- groups[[i]]
    pos <- match(names(x$tfmode), rownames(t1))
    tw <- tw[match(names(x$tfmode), names(tw))]
    sum1 <- matrix(x$tfmode * tw * x$likelihood, 1, length(x$tfmode)) %*% t2[pos, ]
    ss <- sign(sum1)
    ss[ss==0] <- 1
    if (verbose) setTxtProgressBar(pb, i)
    sum2 <- matrix((1-abs(x$tfmode)) * tw * x$likelihood, 1, length(x$tfmode)) %*% t1[pos, ]
    return(list(es=as.vector(abs(sum1) + sum2*(sum2>0)) / sum(tw * x$likelihood), ss=ss))
  }, groups=groups, pb=pb, t1=t1, t2=t2, tw=tw, verbose=verbose)
  names(temp) <- names(groups)
  if (verbose) cat("\nProcess ended at ", date(), "\n", sep="")
  es <- t(sapply(temp, function(x) x$es))
  ss <- t(sapply(temp, function(x) x$ss))
  return(list(groups=es, ss=ss))
}

#' Null model function
#'
#' This function generates the NULL model function, which computes the normalized enrichment score and associated p-value
#'
#' @param pwnull Object generated by \code{pwea3NULLgroups} function
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return List of function to compute NES and p-value

pwea3NULLf <- function(pwnull, verbose=T) {
	if (!is.null(pwnull)) {
    pb <- NULL
    if (verbose) {
      cat("\nComputing the null distribution function for ", nrow(pwnull$groups), " TFs.\n", sep="")
		  cat("Process started at ", date(), "\n", sep="")
      pb <- txtProgressBar(max=nrow(pwnull$groups), style=3)
    }
		res <- lapply(1:nrow(pwnull$groups), function(i, pwnull, pb, verbose) {
      if (verbose) setTxtProgressBar(pb, i)
		  x <- pwnull$groups[i, ]
		  iqr <- quantile(x, c(.5, 1-5/length(x)))
		  pd <- ecdf(x)
		  a <- list(x=knots(pd), y=pd(knots(pd)))
		  filtro <- a$x>iqr[1] & a$x<=iqr[2] & a$y<1
		  spl <- smooth.spline(a$x[filtro], -log(1-a$y[filtro]), spar=.75)
		  return(function(x, alternative="greater") {
		    x <- abs(x)
        p <- exp(-predict(spl, x)$y)
		    pos <- which(x<iqr[1])
		    if (length(pos)>0) p[pos] <- 1-pd(x[pos])
		    nes <- qnorm(p/2, lower.tail=F)
		    if (length(p)==1) {p <- as.vector(p); nes <- as.vector(nes)}
		    list(nes=nes*pwnull$ss[i], p.value=p)
		  })
		}, pwnull=pwnull, pb=pb, verbose=verbose)
		if (verbose) cat("\nProcess ended at ", date(), "\n", sep="")
		names(res) <- rownames(pwnull$groups)
		return(res)
	}
}

pwea3NULLf1 <- function(pwnull) {
  apply(pwnull$groups, 1, function(x) {
    xmed <- median(x)
    x <- x-xmed
    iqr <- quantile(x, c(10/length(x), 1-10/length(x)))
    x[x<0] <- x[x<0]/abs(iqr[1])
    x[x>=0] <- x[x>=0]/iqr[2]
    qnx <- abs(qnorm(10/length(x)))
    x <- x * qnx
    pd <- ecdf(x)
    cutoff=abs(qnorm(10/length(x)))
    return(function(x, alternative=c("two.sided", "less", "greater")[3]) {
      x <- x-xmed
      filtro <- x<0
      x[filtro] <- x[filtro]/abs(iqr[1])
      x[!filtro] <- x[!filtro]/iqr[2]
      x <- x * qnx
      res <- sapply(x, function(x, cutoff, alternative) {
        if (x>(-cutoff) & x<cutoff) res <- switch(pmatch(alternative, c("two.sided", "less", "greater")),
          min(pd(x), 1-pd(x))*2,
          pd(x),
          1-pd(x))
        else res <- switch(pmatch(alternative, c("two.sided", "less", "greater")),
          min(pnorm(x), pnorm(x, lower.tail=F))*2,
          pnorm(x),
          pnorm(x, lower.tail=F))
      }, cutoff=cutoff, alternative=alternative)
      list(nes=as.vector(x), p.value=as.vector(res))
    })
  })
}

ppwea3NULLf <- function(regulon, tw, lw) {
  lapply(regulon, function(x, tw, lw) {
    ww <- x$likelihood^lw*tw[match(names(x$tfmode), names(tw))]
    ww <- ww/max(ww)
    ww <- sqrt(sum(ww^2))
    return(function(x, alternative="two.sided") {
      x <- x*ww
      p <- switch(pmatch(alternative, c("two.sided", "less", "greater")),
                  pnorm(abs(x), lower.tail=F)*2,
                  pnorm(x, lower.tail=T),
                  pnorm(x, lower.tail=F))
      list(nes=x, p.value=p)
    })
  }, tw=tw, lw=lw)
}

regulonScore <- function(x, slope=20, inflection=.3) {
  if (is.list(x)) return(lapply(x, regulonScore, slope=slope, inflection=inflection))
  return(1/(1+inflection^(slope*(abs(x)-inflection)))*sign(x))
}

regulonDist <- function (groups, cutoff = 50, method = "FET", copula = T, mode = F) 
{
    if (!mode) 
        groups <- lapply(groups, function(x) {
            res <- abs(sign(x))
            res[res == 0] <- 1
            return(res)
        })
    total <- unique(unlist(lapply(groups, names), use.names = F))
    groups <- groups[sapply(groups, length) >= cutoff]
    group <- groups
    dmatrix <- NULL
    for (i in 1:(length(groups) - 1)) {
        group <- group[-which(names(group) == names(groups)[i])]
        test <- groups[[i]]
        res <- sapply(group, function(x, test, total, method) {
            tmp1 <- tmp2 <- rep(0, length(total))
            tmp1[match(names(test), total)] <- sign(test)
            tmp2[match(names(x), total)] <- sign(x)
            tmp <- abs(tmp1 + tmp2) > 0
            test1 <- abs(tmp1 * tmp) > 0
            test2 <- abs(tmp2 * tmp) > 0
            tmp <- fisher.test(table(test1, test2), alternative = "greater")
            switch(pmatch(method, c("odds", "FET")), res <- tmp$estimate, 
                res <- tmp$p.value)
            res
        }, test = test, total = total, method = method)
        dmatrix <- cbind(dmatrix, c(dmatrix[i, ], 0, res))
    }
    dmatrix <- cbind(dmatrix, c(dmatrix[length(groups), ], 0))
    rownames(dmatrix) <- colnames(dmatrix) <- names(groups)
    if (copula) {
        rdm <- rank(dmatrix)
        dmatrix <- matrix(rdm/max(rdm), nrow(dmatrix), ncol(dmatrix), 
            dimnames = list(names(groups), names(groups)))
    }
    return(dmatrix)
}

#' Filter for rows of a matrix with no loss of col and row names
#'
#' This function filters the rows of a matrix returning always a two dimensional matrix
#'
#' @param x Matrix
#' @param filter Logical or numerical index of rows
#' @return Matrix

filterRowMatrix <- function(x, filter) {
  if (is.logical(filter)) largo <- length(which(filter))
  else largo <- length(filter)
  matrix(x[filter, ], largo, ncol(x), dimnames=list(rownames(x)[filter], colnames(x)))
}

#' Filter for columns of a matrix with no loss of col and row names
#'
#' This function filters the columns of a matrix returning always a two dimensional matrix
#'
#' @param x Matrix
#' @param filter Logical or numerical index of columns
#' @return Matrix
filterColMatrix <- function(x, filter) t(filterRowMatrix(t(x), filter))

#' Mode of continuous distributions
#' 
#' This function computes the mode for continuous distributions
#' 
#' @param x Numeric data vector
#' @param adj Number indicating the adjustment for the kernel bandwidth
#' @return Number
#' @export
distMode <- function (x, adj = 1) 
{
  tmp <- density(x, adjust = adj)
  return(tmp$x[which.max(tmp$y)])
}

jaccard <- function (c1, c2, zerobyzero = NA) 
{
  if (sum(c1) + sum(c2) - sum(c1 & c2) == 0) 
    out <- zerobyzero
  else out <- sum(c1 & c2)/(sum(c1) + sum(c2) - sum(c1 & c2))
  out
}

aecdf <- function(dnull) {
  dnull <- dnull[is.finite(dnull)]
  iqr <- quantile(dnull, c(5/length(dnull), .5, 1-5/length(dnull)))
  pd <- ecdf(dnull)
  a <- list(x=knots(pd), y=pd(knots(pd)))
  filtro <- a$x<iqr[2] & a$x>=iqr[1] & a$y<1
  spl1 <- smooth.spline(abs(a$x[filtro]), -log(a$y[filtro]), spar=.9)
  filtro <- a$x>iqr[2] & a$x<=iqr[3] & a$y<1
  spl2 <- smooth.spline(a$x[filtro], -log(1-a$y[filtro]), spar=.9)
  dnull <- function(x, alternative="two.sided") {
    p1 <- exp(-predict(spl1, abs(x))$y)
    p2 <- exp(-predict(spl2, abs(x))$y)
    p <- p1*(x<iqr[2])+p2*(x>=iqr[2])
    nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
    switch(match.arg(alternative, c("two.sided", "greater", "less")),
           two.sided={p <- p*2},
           greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
           less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]})
    names(nes) <- names(p) <- names(x)
    list(nes=nes, p.value=p)
  }
  return(dnull)
}
