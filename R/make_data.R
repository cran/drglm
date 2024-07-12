#' Reading Data File Larger than Memory for Fitting GLMs Using \code{big.drglm} Function
#'
#' @param filename Path to the data set on disk.
#' @param chunksize Size of the chunk or subset to be read from the large file for fitting GLMs.
#' @param ... Additional arguments to be passed to \code{read.csv}.
#' @importFrom utils read.csv
#' @export make.data
#' @return A function that reads chunks of the data set.
#' @examples
#' # Create a toy dataset
#' set.seed(123)
#' # Number of rows to be generated
#' n <- 10000
#'
#' # Creating dataset
#' dataset <- data.frame(
#'   Var_1 = round(rnorm(n, mean = 50, sd = 10)),
#'   Var_2 = round(rnorm(n, mean = 7.5, sd = 2.1)),
#'   Var_3 = as.factor(sample(c("0", "1"), n, replace = TRUE)),
#'   Var_4 = as.factor(sample(c("0", "1", "2"), n, replace = TRUE)),
#'   Var_5 = as.factor(sample(0:15, n, replace = TRUE)),
#'   Var_6 = round(rnorm(n, mean = 60, sd = 5))
#' )
#'
#' # Save the dataset to a temporary file
#' temp_file <- tempfile(fileext = ".csv")
#' write.csv(dataset, file = temp_file, row.names = FALSE)
#'
#' # Path to the temporary file
#' dataset_path <- temp_file
#' dataset_path  # Display the path to the temporary file
#'
#' # Initialize the data reading function with the data set path and chunk size
#' da <- drglm::make.data(dataset_path, chunksize = 1000)
#'
#' # Fitting MLR Models
#' nmodel <- drglm::big.drglm(da,
#' formula = Var_1 ~ Var_2 + factor(Var_3) + factor(Var_4) + factor(Var_5) + Var_6,
#' 10, family = "gaussian")
#' # View the results table
#' print(nmodel)
#'
#' # Fitting logistic Regression Model
#' bmodel <- drglm::big.drglm(da,
#' formula = factor(Var_3) ~ Var_1 + Var_2 + factor(Var_4) + factor(Var_5) + Var_6,
#' 10, family = "binomial")
#' # View the results table
#' print(bmodel)
#'
#' # Fitting Poisson Regression Model
#' pmodel <- drglm::big.drglm(da,
#' formula = Var_5 ~ Var_1 + Var_2 + factor(Var_3) + factor(Var_4) + Var_6,
#' 10, family = "poisson")
#' # View the results table
#' print(pmodel)


make.data <- function(filename, chunksize, ...) {
  conn <- NULL
  col_names <- NULL

  function(reset = FALSE) {
    if (reset) {
      if (!is.null(conn)) close(conn)
      conn <<- file(filename, open = "r")
      chunk <- read.csv(conn, nrows = chunksize, header = TRUE, ...)
      col_names <<- names(chunk)
      return(chunk)
    } else {
      rval <- read.csv(conn, nrows = chunksize, header = FALSE, ...)
      if (nrow(rval) == 0) {
        close(conn)
        conn <<- NULL
        rval <<- NULL
      } else {
        names(rval) <- col_names
      }
      return(rval)
    }
  }
}
