#' Plot MARINA results
#'
#' This function generate a plot for MARINA results showing the enrichment of the target genes for each significant master regulator on the gene expression signature
#'
#' @param x MARINa object produced by \code{marina} function
#' @param mrs Either an integer indicating the number of master regulators to include in the plot, or a character vector containing the names of the master regulators to include in the plot
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param pval Optional matrix of p-values to include in the plot
#' @param bins Number of bins to split the vector of scores in order to compute the density color of the bars
#' @param cex Number indicating the text size scaling, 0 indicates automatic scaling
#' @param density Integrer indicating the number of steps for the kernel density. Zero for not ploting it
#' @param smooth Number indicating the proportion of point for smoothing the density distribution. Zero for not using the smoother
#' @param sep Number indicating the separation from figure and text
#' @param hybrid Logical, whether the 3-tail approach used for computingthe enrichment should be reflected in the plot
#' @param include Vector indicating the information to include as heatmap to the right of the marina plot: expression and activity
#' @param gama Positive number indicating the exponential transformation for the activity and expression color scale
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot marina
#' @S3method plot marina
#' @seealso \code{\link{marina}}

plot.marina <- function(x, mrs=10, color=c("cornflowerblue","salmon"), pval=NULL, bins=500, cex=0, density=0, smooth=0, sep=.2, hybrid=T, include=c("expression", "activity"), gama=2, ...) {
  maobject <- x
  rm(x)
  marg <- par("mai")
  rlist <- maobject$signature
  if (ncol(rlist)>0) rlist <- rowMeans(rlist)
	if (maobject$param$iterative) {
		tmp1 <- cbind(rep(names(maobject$regulon), sapply(maobject$regulon, function(x) length(x$tfmode))), unlist(lapply(maobject$regulon, function(x) names(x$tfmode)), use.names=F))
		tmp1 <- cbind(tmp1, abs(maobject$es$nes)[match(tmp1[, 1], names(maobject$es$nes))])
		tw <- tapply(tmp1[, 3], tmp1[, 2], function(x) sum(as.numeric(x)))
	}
	else tw <- 1/table(unlist(lapply(maobject$regulon, function(x) names(x$tfmode)), use.names = F))^.5
	tw <- 1-tw/max(tw)
  if (length(mrs)==1 & is.numeric(mrs[1])) {
    mrs <- names(maobject$es$nes)[order(maobject$es$p.value)[1:round(mrs)]]
    mrs <- mrs[order(maobject$es$nes[match(mrs, names(maobject$es$nes))], decreasing=T)]
  }
  groups <- maobject$regulon[match(mrs, names(maobject$regulon))]
  if (is.null(pval)) pval <- maobject$es$p.value[match(mrs, names(maobject$es$nes))]
  if (is.data.frame(pval)) pval <- as.matrix(pval)
  if (is.matrix(rlist)) rlist <- rlist[,1]
  if (min(rlist)<0) rlist <- sort(rlist)
  else rlist <- sort(-rlist)
  groups <- groups[length(groups):1]
  color1 <- color
  layout(matrix(1:2, 1, 2), widths=c(10-length(include), length(include)))
  color <- rgb2hsv(col2rgb(color))
  satval <- color[3,]
  color <- color[1,]
  preg <- as.numeric(any(sapply(groups, function(x) any(x$tfmode<0))))+1
  textsize <- 1
  xlimit <- c(0, length(rlist)*(1+.2*max(nchar(names(groups)))/8))
  if (!is.null(pval)) {
    if (is.matrix(pval)) {
      pval <- pval[nrow(pval):1, ]
      pval <- pval[, ncol(pval):1]
      xlimit <- c(-length(rlist)*(sep*ncol(pval)+.02), length(rlist)*(1+.2*max(nchar(names(groups)))/9))
    }
    else {
      pval <- pval[length(pval):1]
      xlimit <- c(-length(rlist)*.12, length(rlist)*(1+.2*max(nchar(names(groups)))/9))
      textsize=.8
    }
  }
  if (length(include)>0 & include[1] != "") par(mai=c(marg[1:3], .05))
  plot(1, type="n", ylim=c(0, length(groups)), xlim=xlimit, axes=F, ylab="", xlab="", yaxs="i")
  if (cex>0) textsize <- (length(groups)<=20)*cex+(length(groups)>20)*(20/length(groups))*cex
  switch(preg,
  {
    for (i in 1:length(groups)) {
      densi <- rep(0, length(rlist))
      x <- which(names(rlist) %in% names(groups[[i]]$tfmode))
      if (length(x)>0) {
	      densi[x] <- 1
	      denStep <- round(length(densi)/bins)
	      x1 <- x[x<denStep]
	      x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
	      x3 <- x[x>(length(rlist)-denStep)]
	      densiRes <- sapply(x2, function(i, densi, denStep)
	      sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
	      densiRes <- densiRes/max(densiRes)
	      temp <- rlist[x]
	      if (satval[1]==0) temp <- hsv((temp<0)*color[1] + (temp>0)*color[2], satval, 1-densiRes)
	      else temp <- hsv((sign(temp)<0)*color[1] + (sign(temp)>0)*color[2], densiRes, satval)
	      for (ii in order(densiRes)) lines(c(x[ii], x[ii]),c(i-1, i), col=temp[ii])
	      if (density>0) {
	        denStep <- round(length(densi)/density)
	        xpos <- seq(denStep, length(rlist)-denStep, length=density)
	        densiRes <- sapply(xpos, function(i, densi, denStep)
	        sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
	        densiRes <- densiRes/max(densiRes)
	        if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
	        lines(xpos, i+densiRes-1)
      	}
      }
    }
    text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
    if (!is.null(pval)) {
      if (is.matrix(pval)) {
	for (i in 1:ncol(pval)) {
	  text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, pval[, i], adj=1, cex=.85*textsize)
	  text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(pval)[i], adj=1, cex=1)
	}
      }
      else {
	text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(pval, 3), adj=1, cex=.85*textsize)
	text(0, length(groups)+.5, ifelse(max(pval)>1, "oddsR", "p-value"), adj=1, cex=1.2)
      }
    }
    text(length(rlist)*1.02, length(groups)+.5, "Set", adj=0, cex=1.2)
  },
  {
    for (i in 1:length(groups)) {
      for (ii in 1:2) {
        tset <- groups[[i]]$tfmode
        tset <- tset[(tset<0 & ii==1)|(!(tset<0) & ii==2)]
        tset1 <- names(tset)[abs(tset)>.5]
        tset2 <- names(tset)[abs(tset)<.5]
        if (length(tset)>1) {
	        densi <- rep(0, length(rlist))
          x <- match(names(tset), names(rlist))
					tw1 <- rep(1, length(x))
	        if (hybrid) {
	          x <- match(tset1, names(rlist))
	          if (ii==1) x <- c(x, match(tset2, names(sort(-abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
	          else x <- c(x, match(tset2, names(sort(abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
						tw1 <- groups[[i]]$likelihood[match(c(tset1, tset2), names(groups[[i]]$tfmode))]
						tw1 <- tw1*tw[match(c(tset1, tset2), names(tw))]
						tw1 <- tw1/max(tw1)
	        }
	        densi[x] <- 1
	        denStep <- round(length(densi)/bins)
	        x1 <- x[x<denStep]
	        x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
	        x3 <- x[x>(length(rlist)-denStep)]
	        densiRes <- sapply(x2, function(i, densi, denStep) sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
					densiRes <- densiRes*(tw1[x>=denStep & x <= (length(rlist)-denStep)])
	        densiRes <- densiRes/max(densiRes)
	        temp <- rlist[x]
	        temp <- hsv(color[ii], densiRes, satval)
	        for (iii in order(densiRes)) lines(c(x[iii], x[iii]),c(i-1+(ii-1)/2, i-1+ii/2), col=temp[iii])
	        if (density>0) {
	          denStep <- round(length(densi)/density)
	          xpos <- seq(denStep, length(rlist)-denStep, length=density)
	          densiRes <- sapply(xpos, function(i, densi, denStep)
	          sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
	          densiRes <- densiRes/max(densiRes)
	          if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
	          lines(xpos, i-1+densiRes/2+(ii-1)/2)
	        }
        	}
      }
    }
    text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
    if (!is.null(pval)) {
      if (is.matrix(pval)) {
	for (i in 1:ncol(pval)) {
	  text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, pval[, i], adj=1, cex=.85*textsize)
	  text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(pval)[i], adj=1, cex=1)
	}
      }
      else {
	text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(pval, 3), adj=1, cex=.85*textsize)
  axis(3, -.05*length(rlist), ifelse(max(pval)>1, "oddsR", "p-value"), adj=1, cex=1.2, tick=F, line=-.5)
      }
    }
    axis(3, length(rlist)*1.05, "Set", adj=0, cex=1.2, line=-.5, tick=F)
  })
  abline(h=0:length(groups))
  lines(c(0, 0), c(0, length(groups)))
  lines(c(length(rlist), length(rlist)), c( 0, length(groups)))
  ss <- maobject$signature
  if (!is.null(dim(ss))) ss <- ss[, 1]
  x <- NULL
  xn <- NULL
  xpos <- NULL
  if("activity" %in% include) {x <- maobject$es$nes[match(mrs, names(maobject$es$nes))]/max(abs(maobject$es$nes)); xn <- "Act"}
  if ("expression" %in% include) {x <- cbind(x, ss[match(mrs, names(ss))]/max(abs(ss))); xn <- c(xn, "Exp"); xpos <- rank(-abs(ss))[match(mrs, names(ss))]}
  if (!is.null(x)) {
    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
    rownames(x) <- colnames(x) <- NULL
    marg[2] <- .05
    par(mai=marg)
    scmax <- max(abs(x), na.rm=T)
    x <- abs(x/scmax)^gama*sign(x)
    x <- filterRowMatrix(x, nrow(x):1)
    x1 <- x
    x1[is.na(x1)] <- 0
    coli <- hsv(ifelse(x1<0, color[1], color[2]), abs(x1), 1)
    coli[is.na(x)] <- hsv(0, 0, .5)
    image(1:ncol(x), 1:nrow(x), t(matrix(1:(ncol(x)*nrow(x)), nrow(x), ncol(x))), col=coli, ylab="", xlab="", axes=F, yaxs="i")
    box()
    grid(ncol(x), nrow(x), col="black", lty=1)
    axis(3, 1:length(xn), xn, las=2, line=-.5, tick=F)
    axis(4, length(xpos):1, xpos, las=1, tick=F, line=-.4, cex.axis=.85*textsize)
  }    
  par(mai=marg)
}
