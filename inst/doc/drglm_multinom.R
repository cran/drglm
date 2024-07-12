## ----setup, include=FALSE-----------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
set.seed(123)
#Number of rows to be generated
n <- 1000000
#creating dataset
dataset <- data.frame( 
Var_1 = round(rnorm(n, mean = 50, sd = 10)), 
Var_2 = round(rnorm(n, mean = 7.5, sd = 2.1)), 
Var_3 = as.factor(sample(c("0", "1"), n, replace = TRUE)), 
Var_4 = as.factor(sample(c("0", "1", "2"), n, replace = TRUE)), 
Var_5 = as.factor(sample(0:15, n, replace = TRUE)), 
Var_6 = round(rnorm(n, mean = 60, sd = 5))
)


## -----------------------------------------------------------------------------

mmodel=drglm::drglm.multinom(Var_4~ Var_1+ Var_2+ Var_3+ Var_5+ Var_6, 
                             data=dataset, k=10)
#Output
print(mmodel)

