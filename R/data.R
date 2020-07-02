
pA.toy2 <- read.table("C:/Science/Projects/PolyA site selection/APAlog/data/pAsite.toy2.txt", header = TRUE)
pA_design <- read.table("C:/Science/Projects/PolyA site selection/APAlog/data/sgHNRNPC_sgCTRL_design.txt", header = TRUE)

usethis::use_data(pA.toy2, pA_design, overwrite = TRUE)

#' Examples of read count dataset and design matrix for alternative poly A site usage analysis

#' pA.toy2
#' @format A data frame typically with a minimum of four columns
#' \describe{
#'  \item{transcript}{transcript IDs}
#'  \item{pA.site}{pA site names or codes}
#'  \item{sample}{sample IDs}
#'  \item{count}{read counts mapping to the the pA site - transcript - sample specified in previous columns}
#'  }
#'  Note: The column containing read counts MUST be called "count".
#'
#'  Sample name can be used directly as the predictor variable in the regression model. Alternatively,
#'  the user can provide any number of sample attributes in a design matrix to be used as regressors. Below is
#'  an example of a simple design matrix. There must be a key variable with the same name in the count matrix
#'  and the design matrix to connect the two (similar to relational databases). In our example data,
#'  it is called "sample".
#'
#'  pA_design
#'  @format A data frame describing sample attributes typically with a minimum of two columns
#'  \describe{
#'   \item{sample}{sample IDs}
#'   \item{cell_line}{Cell line names}
#'   \item{rep}{replicate number}
#'  }
