#' Fitting Linear and Generalized Linear Model in "Divide and Recombine" approach to Large Data Sets
#' @description
#' Function \code{drglm} aimed to fit GLMs to datasets larger in size that can be stored in memory. It uses popular divide and recombine technique to handle large data sets efficiently.Function \code{drglm} optimizes performance when linked with optimized BLAS libraries like ATLAS.The function \code{drglm} requires defining the number of chunks K and the fitfunction.The rest of the arguments are almost identical with the speedglm or biglm package.
#'
#' @param formula An entity belonging to the "formula" class (or one that can be transformed into that class) represents a symbolic representation of the model that needs to be adjusted. Specifics about how the model is defined can be found in the 'Details' section.
#' @param family An explanation of the error distribution that will be implemented in the model.
#' @param data A data frame, list, or environment that is not required but can be provided if available.
#' @param fitfunction The function to be utilized for model fitting. \code{glm} or \code{speedglm} should be used.For Multinomial models, \code{multinom} function is preferred.
#' @param k Number of subsets to be used.
#'
#' @return A Generalized Linear Model is fitted in "Divide & Recombine" approach using "k" chunks to data set. A list of model coefficients is estimated using divide and recombine method with the respective standard error of estimates.
#' @references
#' - Xi, R., Lin, N., & Chen, Y. (2009). Compression and aggregation for logistic regression analysis in data cubes. IEEE Transactions on Knowledge and Data Engineering, 21(4).
#' - Chen, Y., Dong, G., Han, J., Pei, J., Wah, B. W., & Wang, J. (2006). Regression cubes with lossless compression and aggregation. IEEE Transactions on Knowledge and Data Engineering, 18(12).
#' - Zuo, W., & Li, Y. (2018). A New Stochastic Restricted Liu Estimator for the Logistic Regression Model. Open Journal of Statistics, 08(01).
#' - Karim, M. R., & Islam, M. A. (2019). Reliability and Survival Analysis. In Reliability and Survival Analysis.
#' - Enea, M. (2009) Fitting Linear Models and Generalized Linear Models with large data sets in R.
#' - Bates, D. (2009) Technical Report on Least Square Calculations.
#' - Lumley, T. (2009) \emph{biglm} package documentation.
#' @author MH Nayem
#' @seealso \code{\link{big.drglm}}, \code{\link{drglm.multinom}}
#' @import nnet
#' @import speedglm
#' @import stats
#' @importFrom stats glm model.matrix model.response model.frame qnorm pnorm coef qt gaussian binomial poisson
#' @importFrom speedglm speedglm
#' @importFrom nnet multinom
#' @export drglm
#' @examples
#' set.seed(123)
#' #Number of rows to be generated
#' n <- 10000
#' #creating dataset
#' dataset <- data.frame( pred_1 = round(rnorm(n, mean = 50, sd = 10)),
#' pred_2 = round(rnorm(n, mean = 7.5, sd = 2.1)),
#' pred_3 = as.factor(sample(c("0", "1"), n, replace = TRUE)),
#' pred_4 = as.factor(sample(c("0", "1", "2"), n, replace = TRUE)),
#' pred_5 = sample(0:15, n, replace = TRUE),
#' pred_6 = round(rnorm(n, mean = 60, sd = 5)))
#' #fitting MLRM
#' nmodel= drglm::drglm(pred_1 ~ pred_2+ pred_3+ pred_4+ pred_5+ pred_6,
#' data=dataset, family="gaussian", fitfunction="speedglm", k=10)
#' #Output
#' nmodel
#' #fitting simple logistic regression model
#' bmodel=drglm::drglm(pred_3~ pred_1+ pred_2+ pred_4+ pred_5+ pred_6,
#' data=dataset, family="binomial", fitfunction="speedglm", k=10)
#' #Output
#' bmodel
#' #fitting poisson regression model
#' pmodel=drglm::drglm(pred_5~ pred_1+ pred_2+ pred_3+ pred_4+ pred_6,
#' data=dataset, family="poisson", fitfunction="speedglm", k=10)
#' #Output
#' pmodel
#' #fitting multinomial logistic regression model
#' mmodel=drglm::drglm(pred_4~ pred_1+ pred_2+ pred_3+ pred_5+ pred_6,
#' data=dataset, family="multinomial", fitfunction="multinom", k=10)
#' #Output
#' mmodel


drglm<-function(formula,family,data,k,fitfunction)
{
  if(family=="gaussian" & fitfunction=="glm")
  {

    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {
      model[[i]] <- glm(formula ,split_data[[i]], family = gaussian)
    }

    X <- list()
    Xt <- list()

    for (i in 1:k) {
      X[[i]] <- model.matrix(model[[i]])
      Xt[[i]] <- t(X[[i]])
    }

    Y <- list()

    for (i in 1:k) {
      Y[[i]] <- as.matrix(model.response(model.frame(model[[i]])))
    }

    XtX=lapply(1:k, function(i) Xt[[i]]%*%X[[i]])
    XtX=Reduce("+",XtX)
    XtY=lapply(1:k, function(i) Xt[[i]]%*%Y[[i]])
    XtY=Reduce("+",XtY)


    IXtX=solve(XtX)
    B=IXtX%*%XtY

    vcov <- list()

    for (i in 1:length(model)) {

      vcov[[i]] <- vcov(model[[i]])
    }


    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com)

    alpha <- 0.05

    z<-qnorm(1-alpha/2)


    #lower and upper bounds using gaussian distribution
    l_normal= B-z*se_com
    u_normal=B+z*se_com


    #calculating z values
    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    # creating a data frame with four columns
    table <- data.frame("Estimate"=B,
                        "standard error"=se_com ,
                        "t value"= Z,
                        "Pr(>|t|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    # viewtable
    table
  }
  else if(family=="gaussian" & fitfunction=="speedglm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      # fit the glm() function for each element
      model[[i]] <- speedglm::speedglm(formula ,data=split_data[[i]], family = gaussian())
    }

    X <- list()
    Xt <- list()

    for (i in 1:k) {
      X[[i]] <- model.matrix(formula, data=split_data[[i]])
      attr(X[[i]], "contrasts") <- NULL
      attr(X[[i]], "assign") <- NULL
      Xt[[i]] <- t(X[[i]])
    }

    Y <- list()

    for (i in 1:k) {
      Y[[i]]<- as.matrix(model.response(model.frame(formula, data=split_data[[i]])))
    }


    XtX=lapply(1:k, function(i) Xt[[i]]%*%X[[i]])
    XtX=Reduce("+",XtX)
    XtY=lapply(1:k, function(i) Xt[[i]]%*%Y[[i]])
    XtY=Reduce("+",XtY)

    IXtX=solve(XtX)
    B=IXtX%*%XtY

    vcov <- list()

    for (i in 1:length(model)) {

      vcov[[i]] <- vcov(model[[i]])
    }


    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com)

    alpha <- 0.05

    z<-qnorm(1-alpha/2)


    #lower and upper bounds using gaussian distribution
    l_normal= B-z*se_com
    u_normal=B+z*se_com


    #calculating z values
    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    # creating a data frame with four columns
    table <- data.frame("Estimate"=B,
                        "standard error"=se_com ,
                        "t value"= Z,
                        "Pr(>|t|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    # viewtable
    table
  }
  else if(family=="binomial" & fitfunction=="glm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      # fit the glm() function for each element
      model[[i]] <- glm(formula,split_data[[i]], family = binomial)
    }

    vcov <- list()

    for (i in 1:length(model)) {
      # calculate the vcov() function for each element
      vcov[[i]] <- vcov(model[[i]])
    }

    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com) # Take the square root of v_com

    H <- list()

    for (i in 1:length(vcov))
    {
      H[[i]] <- as.matrix(solve(vcov[[i]]))
    }

    # Beta's of models
    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }

    HE <- Reduce("+", H)

    IH=solve(HE)

    HB <- lapply(1:k, function(i) H[[i]]%*%b[[i]])
    # Sum the elements of the list
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

    # creating a data frame with four columns
    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    table

  }
  else if (family=="binomial" &  fitfunction=="speedglm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      model[[i]] <- speedglm::speedglm(formula,split_data[[i]], family = binomial())
    }

    vcov <- list()

    for (i in 1:length(model)) {
      # calculate the vcov() function for each element
      vcov[[i]] <- vcov(model[[i]])
    }

    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com) # Take the square root of v_com

    H <- list()

    for (i in 1:length(vcov))
    {
      H[[i]] <- as.matrix(solve(vcov[[i]]))
    }

    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }

    HE <- Reduce("+", H)
    IH=solve(HE)

    HB <- lapply(1:k, function(i) H[[i]]%*%b[[i]])
    # Sum the elements of the list
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

    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    table
  }
  else if(family=="poisson" & fitfunction=="glm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      model[[i]] <- glm(formula,split_data[[i]], family = poisson)
    }

    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }

    B={Reduce("+",b)}/k
    OR=exp(B)

    vcov <- list()

    for (i in 1:length(model)) {

      vcov[[i]] <- vcov(model[[i]])
    }


    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com)

    alpha <- 0.05

    z<-qnorm(1-alpha/2)

    l_normal= B-z*se_com
    u_normal=B+z*se_com

    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    table
  }
  else if(family=="poisson" & fitfunction=="speedglm")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      model[[i]] <- speedglm::speedglm(formula,split_data[[i]], family = poisson())
    }

    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }

    B={Reduce("+",b)}/k
    OR=exp(B)

    vcov <- list()

    for (i in 1:length(model)) {

      vcov[[i]] <- vcov(model[[i]])
    }


    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- as.matrix(diag(vcov[[i]]))/k
    }

    v_com <- Reduce("+", v)/k

    se_com <- sqrt(v_com)

    alpha <- 0.05

    z<-qnorm(1-alpha/2)

    l_normal= B-z*se_com
    u_normal=B+z*se_com

    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    table <- data.frame("Estimate"=B,
                        "Odds Ratio"=OR,
                        "standard error"=se_com ,
                        "z value"= Z,
                        "Pr(>|z|)"=p_value,
                        "95% CI" = paste("[", round(l_normal, 2), ",", round(u_normal, 2), "]"),
                        check.names = FALSE)

    table
  }
  else if(family=="multinomial" & fitfunction=="multinom")
  {
    n = nrow(data)
    rows_per_chunk = ceiling(n / k)
    split_data = lapply(1:k, function(i) {
      start_row = (i-1) * rows_per_chunk + 1
      end_row = min(i * rows_per_chunk, n)
      return(data[start_row:end_row, ])
    })

    model <- list()

    for (i in 1:length(split_data)) {

      model[[i]] <- nnet::multinom(formula=formula,data=split_data[[i]],Hess = TRUE)
    }
    response=as.character(formula(model[[1]])[[2]])
    b <- list()
    for (i in 1:length(model))
    {
      b[[i]] <- as.matrix(coef(model[[i]]))
    }

    B={Reduce("+",b)}/k
    row_names_B <- rownames(B)
    OR=exp(B)


    H<-list()
    for (i in 1:length(model))
    {
      H[[i]]<-(model[[i]])$Hessian
    }

    vcov <- list()

    for (i in 1:length(model)) {

      vcov[[i]] <- as.matrix(solve(H[[i]]))
    }

    v <- list()
    for (i in 1:length(vcov))
    {
      v[[i]] <- matrix(diag(vcov[[i]]),nrow=(nlevels(data[[response]]))-1,byrow=TRUE)/k
    }

    v_com <- Reduce("+", v)/k



    se_com <- sqrt(v_com)
    rownames(se_com) <- row_names_B

    alpha <- 0.05

    z<-qnorm(1-alpha/2)

    l_normal= B-z*se_com
    u_normal=B+z*se_com

    Z=B/se_com

    p_value=2*(1- pnorm(abs(Z)))

    table <- data.frame("Estimate"=t(B),
                        "Odds Ratio"=t(OR),
                        "standard error"=t(se_com) ,
                        "z value"= t(Z),
                        "Pr(>|z|)"=t(p_value),
                        "95% lower CI"=t(l_normal),
                        "95% upper CI" = t(u_normal ),
                        check.names = FALSE)

    # viewtable
    table
  }
  else
  {
    stop("Unsupported family")
  }
}
