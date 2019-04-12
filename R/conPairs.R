# ---- roxygen documentation ----
#
#' @title Identify contact pairs
#'
#' @description
#' Create a dataframe where each row represents a single contact pair.
#'
#' @details
#' This function is used to extract contact pairs following use of the \code{conProcess} or \code{conPhase} function. The returned data frame has two new columns: contact_orig_rowid - the original row id of that particular fix, and contact_pair_id - a unique identifier to show which two fixes are represented by a pair of contacts. The number of unique pairs of contacts is then the highest number in this column, and will be equal to half the number of rows in the data frame.

#' @param ltraj an object of the class \code{ltraj} which is output from the function \code{conProcess} or \code{conPhase}.
#'
#' @return
#' A data frame, where each row represents one of the two fixes in each unique contact pair.
#'
# @references
#'
#' @keywords contacts
#' @seealso conProcess, conPhase
#' 
#' @examples 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doephas <- conPhase(doecons,pc=60*60)
#' prs <- conPairs(doephas)
#' head(prs)
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----
conPairs <- function(ltraj){
  cdf <- ld(ltraj)
  ind <- which(cdf$contacts > 0)
  outcdf <- cdf[rep(rownames(cdf), cdf$contacts), ]
  outcdf$contact_orig_rowid <- NA
  outcdf$contact_id <- as.vector(unlist(sapply(cdf$contact_id[ind],function(x) unlist(strsplit(as.character(x),',')),USE.NAMES = FALSE)))
  outcdf$contact_rowid <- as.vector(unlist(sapply(cdf$contact_rowid[ind],function(x) unlist(strsplit(as.character(x),',')),USE.NAMES = FALSE)))
  outcdf$contact_d <- as.vector(unlist(sapply(cdf$contact_d[ind],function(x) unlist(strsplit(as.character(x),',')),USE.NAMES = FALSE)))
  outcdf$contact_dt <- as.vector(unlist(sapply(cdf$contact_dt[ind],function(x) unlist(strsplit(as.character(x),',')),USE.NAMES = FALSE)))
  outcdf$contact_orig_rowid <- rep(ind,times=cdf$contacts[ind])

  outcdf$contact_pair_id <- NA
  pair_id <- 1
  for (i in 1:(dim(outcdf)[1])){
    if (is.na(outcdf$contact_pair_id[i])){
      ind <- which(outcdf$id == outcdf$contact_id[i] & outcdf$contact_rowid == outcdf$contact_orig_rowid[i])
      outcdf$contact_pair_id[i] <- pair_id
      outcdf$contact_pair_id[ind] <- pair_id
      pair_id <- pair_id + 1
    }
  }
  outcdf <- outcdf[order(outcdf$contact_pair_id),]
  outcdf$contact_d <- as.numeric(outcdf$contact_d)
  outcdf$contact_dt <- as.numeric(outcdf$contact_dt)
  return(outcdf)
}
