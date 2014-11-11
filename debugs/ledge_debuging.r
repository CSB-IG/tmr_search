#' Leading-edge analysis
#'
#' This function performs a Leading-Edge analysis on an object of class marina
#'
#' @param mobj marina class object
#' @return marina object updated with a ledge slot
#' @seealso \code{\link{marina}}
#' @export

ledge <- function(mobj) {
  ss <- mobj$signature
  if (!is.null(dim(ss))) ss <- ss[, 1]
  tmp <- lapply(mobj$regulon, function(x, ss) {

    #x <- mobj$regulon

    est <- gsea2ESint(ss, x$tfmode, 1)
    
    #return(est)
    print(class(est$es))
    print(max(est$es))
    print(min(est$es))
    
    res <- max(est$es, na.rm=TRUE) + min(est$es, na.rm=TRUE)
    pos1 <- which.max(abs(est$es))
    #print(res)
    if (res>0) {
      leg <- names(est$rlist)[1:pos1]
      leg <- leg[leg %in% names(x$tfmode)]
    }
    else {
      leg <- names(est$rlist)[pos1:length(est$es)]
      leg <- leg[leg %in% names(x$tfmode)]
    }
    return(leg)
  }, ss=ss)
  mobj$ledge <- tmp
  return(mobj)
}

# Corre todos los análisis de *gsea2ENint*
gsae_result <- lapply(mobj$regulon, function(x, rlist, score) 
  {gsea2ESint(rlist=rlist, x=x$tfmode, score=score)}, rlist=ss, score=1)

# Corre *gsea2ENint* para un gen (p.ej. BARX1)
gsea2ESint(rlist=ss,x=mobj$regulon$BARX1$ftmode,score=1)



gsea2ESint <- function(rlist, x, score) {
  n_inf <- which(is.infinite(rlist))
  rlist <- rlist[-n_inf] # ¡Elimita los infinitos de la lista rlist! estos están presentes en la firma (ss)
  x1 <- names(x)[x>=0]
  x2 <- names(x)[x<0]
  px1 <- match(x1, names(sort(rlist, decreasing=T)))
  px2 <- match(x2, names(sort(rlist, decreasing=F)))
  names(px1) <- x1
  names(px2) <- x2
  rlistr <- rank(rlist)
  rlistr <- rlistr[!(names(rlist) %in% names(x))]
  rlistr <- sort(c(rlistr, px1, px2))
  rlist <- rlist[match(names(rlistr), names(rlist))]
  x <- which(names(rlist) %in% names(x))
  nr <- sum(abs(rlist[x])^score)
  nh <- length(rlist)-length(x)
  es <- rep(-(1/nh),length(rlist))
  es[x] <- abs(rlist[x])^score/nr
  return(list(es=cumsum(es), rlist=rlist))
}


### Debug *gsea2ENint*

y <- mobj$regulon$BARX1$ftmode

gsea2ESint <- function(rlist, y, score) {

  x1 <- names(y)[y>=0]
  x2 <- names(y)[y<0]
  px1 <- match(x1, names(sort(ss, decreasing=T)))
  px2 <- match(x2, names(sort(ss, decreasing=F)))
  names(px1) <- x1
  names(px2) <- x2
  rlistr <- rank(ss)
  rlistr <- rlistr[!(names(ss) %in% names(x))]
  rlistr <- sort(c(rlistr, px1, px2))
  rlist <- ss[match(names(rlistr), names(ss))]
  y <- which(names(rlist) %in% names(y))
  nr <- sum(abs(rlist[y])^1)
  nh <- length(rlist)-length(y)
  es <- rep(-(1/nh),length(rlist))
  es[y] <- abs(rlist[y])^1/nr
  return(list(es=cumsum(es), rlist=rlist))
}