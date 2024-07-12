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

# Save the dataset to a temporary file
temp_file <- tempfile(fileext = ".csv")
write.csv(dataset, file = temp_file, row.names = FALSE)

# Path to the temporary file
dataset_path <- temp_file
dataset_path  # Display the path to the temporary file


## -----------------------------------------------------------------------------
# Path to the temporary file
dataset_path <- temp_file
dataset_path  # Display the path to the temporary file


# Initialize the data reading function with the data set path and chunk size
da <- drglm::make.data(dataset_path, chunksize = 100000)

## -----------------------------------------------------------------------------
# Fitting MLR Model
nmodel <- 
  drglm::big.drglm(da, formula = Var_1 ~ Var_2+ factor(Var_3)+
                     factor(Var_4)+ factor(Var_5)+ Var_6, 10, family="gaussian")

# View the results table
print(nmodel)

## -----------------------------------------------------------------------------
# Fitting Logistic Model
bmodel <- drglm::big.drglm(da,formula = factor(Var_3) ~ Var_1+ Var_2+
                             factor(Var_4)+ factor(Var_5)+ Var_6, 
                           10, family="binomial")
# View the results table
print(bmodel)

## -----------------------------------------------------------------------------
# Fitting Poisson Regression Model
pmodel <- drglm::big.drglm(da,
                           formula = Var_5 ~ Var_1+ Var_2+ factor(Var_3)+ 
                             factor(Var_4)+ Var_6, 10, family="poisson")
# View the results table
print(pmodel)

