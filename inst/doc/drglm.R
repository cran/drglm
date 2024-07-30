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
Var_5 = sample(0:6, n, replace = TRUE), 
Var_6 = round(rnorm(n, mean = 60, sd = 5))
)


## -----------------------------------------------------------------------------

nmodel= drglm::drglm(Var_1 ~ Var_2+ Var_3+ Var_4+ Var_5+ Var_6,  
                     data=dataset, family="gaussian", 
                     fitfunction="speedglm", k=10)
#Output
print(nmodel)

## -----------------------------------------------------------------------------

bmodel=drglm::drglm(Var_3~ Var_1+ Var_2+ Var_4+ Var_5+ Var_6, 
                    data=dataset, family="binomial",
                    fitfunction="speedglm", k=10)
#Output

print(bmodel)

## -----------------------------------------------------------------------------

pmodel=drglm::drglm(Var_5~ Var_1+ Var_2+ Var_3+ Var_4+ Var_6, 
                    data=dataset, family="poisson", 
                    fitfunction="speedglm", k=10)

#Output
print(pmodel)

## -----------------------------------------------------------------------------

mmodel=drglm::drglm(Var_4~ Var_1+ Var_2+ Var_3+ Var_5+ Var_6, 
              data=dataset,family="multinomial",
              fitfunction="multinom", k=10)
#Output
print(mmodel)

