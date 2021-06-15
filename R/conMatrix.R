# ---- roxygen documentation ----
#
#' @title conMatrix
#'
#' @description
#' Create a matrix that can be used for social network analysis.
#'
#' @details
#' This function is used to calculate the contact rates between individuals and output them in the form of a matrix.
#' NOTE: This function is designed to be used only when a single ltraj object is input into the conProcess function. I.e., single-species contact networks.

#' @param mtraj an object of the class \code{ltraj} which is output from the function \code{conProcess}.
#' @param idcol column id associated with IDs of the individuals, default is the 'burst'
#' @param output ('count' or 'rate') whether to compute the counts of contacts in the contact matrix (default) or the contact-rate.
#'
#' @return
#' A matrix, with the contact rates between individuals.
#'
# @references
#'
#' @keywords contacts
#' @seealso conProcess, conPairs
#' 
#' @examples 
#' \dontrun{
#' data(does)
#' doecons <- conProcess(does,tc=15*60,dc=50)
#' doemat <- conMatrix(doecons)
#' doemat_rate <- conMatrix(doecons,output='rate')
#' }
#' 
#' @export
#
# ---- End of roxygen documentation ----

conMatrix <- function(mtraj,idcol='burst',output='count'){
  
  conp = conPairs(mtraj)
  
  dfr <- ld(mtraj)
  col1 <- which(names(dfr)==idcol)
  ids <- as.character(unique(dfr[,col1]))
  #Get all the unique combinations between ids
  mat <- matrix(nrow=length(ids),ncol=length(ids),
               dimnames=list(ids,ids))
  
  for (i in 1:length(ids)){
    id1 = ids[i]
    traj_sub = dfr[dfr[,idcol] == id1,]
    n_fix = nrow(traj_sub)
    
    df_sub = conp[conp[,idcol] == id1,]
    for (j in 1:length(ids)){
      id2 = ids[j]
      df_con = df_sub[df_sub[,'contact_id'] == id2,]
      n_con = nrow(df_con)
      
      if (output == 'rate'){
        mat[i,j] = n_con / n_fix
      } else {
        mat[i,j] = n_con
      }
    }
  }
  return(mat)
}
