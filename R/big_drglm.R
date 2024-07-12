#' Fitting Linear and Generalized Linear Models to out of the memory data sets in "Divide and Recombine" approach
#' @description
#' Function \code{big.drglm} aimed to fit GLMs to datasets larger in size that can not be  stored in memory. It uses popular divide and recombine technique to handle large data sets efficiently.
#' @param data.generator Using the function \code{\link{make.data}} to initialize the data reading function with the data set path and chunk size, then the data.generate is used directly as data source for the \code{\link{big.drglm}} function.
#' @param formula An entity belonging to the "formula" class (or one that can be transformed into that class) represents a symbolic representation of the model that needs to be adjusted. Specifics about how the model is defined can be found in the 'Details' section.
#' @param chunks Number of subsets to be divided.
#' @param family An explanation of the error distribution that will be implemented in the model.
#'
#' @return A Generalized Linear Model is fitted in "Divide & Recombine" approach using preferred number of chunks to data set. A list of model coefficients is estimated using divide and recombine method with the respective standard error of estimates.
#' @references
#' - Xi, R., Lin, N., & Chen, Y. (2009). Compression and aggregation for logistic regression analysis in data cubes. IEEE Transactions on Knowledge and Data Engineering, 21(4).
#' - Chen, Y., Dong, G., Han, J., Pei, J., Wah, B. W., & Wang, J. (2006). Regression cubes with losseless compression and aggregation. IEEE Transactions on Knowledge and Data Engineering, 18(12).
#' - Zuo, W., & Li, Y. (2018). A New Stochastic Restricted Liu Estimator for the Logistic Regression Model. Open Journal of Statistics, 08(01).
#' - Karim, M. R., & Islam, M. A. (2019). Reliability and Survival Analysis. In Reliability and Survival Analysis.
#' - Enea, M. (2009) Fitting Linear Models and Generalized Linear Models with large data sets in R.
#' - Bates, D. (2009) Technical Report on Least Square Calculations.
#' - Lumley, T. (2009) \emph{biglm} package documentation.
#' @author MH Nayem
#' @seealso \code{\link{drglm}}, \code{\link{drglm.multinom}}
#' @import nnet
#' @import stats
#' @import speedglm speedglm
#' @importFrom stats gaussian model.matrix model.response model.frame qnorm pnorm binomial coef qt poisson
#' @importFrom nnet multinom
#' @importFrom speedglm speedglm
#' @export big.drglm
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
#'# Initialize the data reading function with the data set path and chunk size
#'da <- drglm::make.data(dataset_path, chunksize = 1000)
#'# Fitting MLR Models
#'nmodel <- drglm::big.drglm(da,
#'formula = Var_1 ~ Var_2+ factor(Var_3)+factor(Var_4)+ factor(Var_5)+ Var_6,
#'10, family="gaussian")

#'# View the results table
#'print(nmodel)
#'# Fitting logistic Regression Model
#'bmodel <- drglm::big.drglm(da,
#'formula = factor(Var_3) ~ Var_1+ Var_2+ factor(Var_4)+ factor(Var_5)+ Var_6,
#'10, family="binomial")
#'# View the results table
#'print(bmodel)

#'# Fitting Poisson Regression Model
#'pmodel <- drglm::big.drglm(da,
#'formula = Var_5 ~ Var_1+ Var_2+ factor(Var_3)+ factor(Var_4)+ Var_6,
#'10, family="poisson")
#'# View the results table
#'print(pmodel)


big.drglm <- function(data.generator, formula, chunks, family)
{
  if(family=="gaussian" )
  {
    XtX_combined <- NULL
    XtY_combined <- NULL
    vcov_combined <- NULL

    for (i in 1:chunks) {
      chunk <- data.generator(reset = i == 1)

      model <- speedglm::speedglm(formula, family = gaussian(), data = chunk)

      X <- model.matrix(formula, data = chunk)

      Xt <- t(X)

      Y <- as.matrix(model.response(model.frame(formula, data = chunk)))

      if (is.null(XtX_combined)) {
        XtX_combined <- Xt %*% X
        XtY_combined <- Xt %*% Y
        vcov_combined <- vcov(model) / chunks
      } else {
        XtX_combined <- XtX_combined + (Xt %*% X)
        XtY_combined <- XtY_combined + (Xt %*% Y)
        vcov_combined <- vcov_combined + (vcov(model) / chunks)
      }

      rm(chunk, model, X, Xt, Y)
      gc()
    }

    B <- solve(XtX_combined) %*% XtY_combined

    v_com <- vcov_combined / chunks
    se_com <- sqrt(diag(v_com))

    alpha <- 0.05
    z <- qnorm(1 - alpha / 2)
    l_normal <- B - z * se_com
    u_normal <- B + z * se_com
    Z <- B / se_com
    p_value <- 2 * (1 - pnorm(abs(Z)))

    # Create a data frame with the results
    table <- data.frame(
      "Estimate" = B,
      "standard.error" = se_com,
      "t_value" = Z,
      "Pr(>|z|)" = p_value,
      "normal.CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]")
    )

    # Return the results table
    return(table)
  }
  else if(family=="binomial")
  {

    vcov <- list()
    v <- list()
    H <- list()
    b <- list()

    for (i in 1:chunks) {
      chunk <- data.generator(reset = i == 1)

      model <- speedglm::speedglm(formula, family = binomial(), data = chunk)

      vcov[[i]] <- vcov(model)

      v[[i]] <- as.matrix(diag(vcov[[i]]))/chunks

      H[[i]] <- as.matrix(solve(vcov[[i]]))

      b[[i]] <- as.matrix(coef(model))

      rm(chunk, model)
      gc()
    }

    v_com <- Reduce("+", v)/chunks
    se_com <- sqrt(v_com)

    HE <- Reduce("+", H)
    IH=solve(HE)

    HB <- lapply(1:chunks, function(i) H[[i]]%*%b[[i]])

    HB <- Reduce("+", HB)

    B=(IH%*%HB)
    OR=exp(B)

    alpha <- 0.05
    z<-qnorm(1-alpha/2)
    t<-qt(1-alpha/2,df = length(B))

    l_normal= B-z*se_com
    u_normal=B+z*se_com

    lower_t=B-t*se_com
    upper_t=B+t*se_com

    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    table <- data.frame("Estimate"=B,"Odds Ratio"=OR,"standard.error"=se_com , "t value"= Z, "Pr(>|z|)"=p_value,"normal.CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"))

    # Return the results table
    return(table)
  }
  else if(family=="poisson")
  {

    model <- list()
    b <- list()
    vcov <- list()
    v <- list()

    for (i in 1:chunks) {
      chunk <- data.generator(reset = i == 1)

      model[[i]] <- speedglm::speedglm(formula, chunk, family = poisson())

      b[[i]] <- as.matrix(coef(model[[i]]))

      vcov[[i]] <- vcov(model[[i]])

      v[[i]] <- as.matrix(diag(vcov[[i]]))/chunks

      rm(chunk)
      gc()
    }

    B = Reduce("+", b) / chunks
    OR = exp(B)

    v_com = Reduce("+", chunks) / chunks
    se_com = sqrt(v_com)

    alpha <- 0.05
    z <- qnorm(1 - alpha / 2)

    l_normal = B - z * se_com
    u_normal = B + z * se_com

    Z = B / se_com
    p_value = 2 * (1 - pnorm(abs(Z)))


    table <- data.frame(
      "Estimate" = B,
      "Odds Ratio" = OR,
      "standard.error" = se_com,
      "Z_value" = Z,
      "Pr(>|z|)" = p_value,
      "normal.CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]")
    )

    # Return the results table
    return(table)

  } else
  {
    stop("Unsupported family")
  }
}



make.data <- function(filename, chunksize, ...)
{
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
