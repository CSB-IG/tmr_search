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
    est <- gsea2ESint(ss, x$tfmode, 1)
    res <- max(est$es) + min(est$es)
    pos1 <- which.max(abs(est$es))
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

gsea2ESint <- function(rlist, x, score) {
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
