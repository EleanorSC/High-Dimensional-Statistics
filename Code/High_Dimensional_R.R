## ----------------------------------------------------------------------
##
## Script name: High_Dimensional_R
##
## Purpose of script: 
##      # Full code for High Dimensional Statistics course
##
## Author: Dr. Eleanor Conole
##
## Date Created: 2024-11-05 
##
##
##
## ----------------------------------------------------------------------
## Setup
## ----------------------------------------------------------------------
##
## Ensure you have completed any necessary setup instructions before proceeding.
## ----------------------------------------------------------------------

install.packages("BiocManager", quiet = TRUE)

# Download and read dependencies file
download.file(
  "https://raw.githubusercontent.com/EleanorSC/High-Dimensional-Statistics/main/dependencies.csv",
  destfile = 'dependencies.csv'
)
table <- read.table('dependencies.csv')

# Install dependencies using BiocManager
BiocManager::install(table[[1]])

# Create a directory for data files
dir.create("data", showWarnings = FALSE)

# List of data files to download
data_files <- c(
  "prostate.rds"
  "coefHorvath.rds",
  "methylation.rds",
)

# Download data files into the "data" directory
for (file in data_files) {
  download.file(
    url = file.path(
      "https://raw.githubusercontent.com/EleanorSC/High-Dimensional-Statistics/main/Data",
      file
    ),
    destfile = file.path("data", file)
  )
}

# ----------------------------------------------------------------------
# High-Dimensional Data Analysis: Challenge 2
# 
# Challenge introduces how we can work with a dataset that contains multiple features
# and explores the common challenges of high-dimensional data analysis.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Load the Prostate Dataset
# ----------------------------------------------------------------------

# Load the dataset using the `here` package. 
# Ensure the 'prostate.rds' file is located in the 'data' folder of your project directory.

prostate <- readRDS(here("data/prostate.rds"))

# ----------------------------------------------------------------------
# Examine the Dataset
# ----------------------------------------------------------------------

# The dataset represents clinical data where each row corresponds to a single patient.
# Let's first examine the dimensions of the dataset to understand its structure.

# Determine the number of observations (n) and features (p)
dimensions <- dim(prostate)  # Returns a vector: [number of rows, number of columns]
print(dimensions)  # Print the dimensions

# Examine the variables measured in the dataset
variable_names <- names(prostate)  # Returns the column names (variables)
print(variable_names)

# View the first few rows of the dataset to get a sense of its contents
head(prostate)

# ----------------------------------------------------------------------
# Visualize Relationships Between Variables
# ----------------------------------------------------------------------

# Plot pairwise relationships between variables using the `pairs()` function.
# This function creates a scatterplot matrix to visualize the relationships
# between all pairs of numeric variables.

pairs(prostate)

# ----------------------------------------------------------------------
# Reflection on High-Dimensional Data Challenges
# ----------------------------------------------------------------------

# Observations:
# - Pairwise scatterplots reveal correlations between variables, which can
#   complicate analyses due to multicollinearity.
# - High-dimensional datasets often have a greater number of features (p)
#   than observations (n), which may lead to overfitting or difficulties in
#   model interpretation.


# ----------------------------------------------------------------------
# High-Dimensional Data Analysis: Challenge 3
# 
# 1). Use the cor() function to examine correlations between all variables 
#  in the prostate dataset. Are some pairs of variables highly correlated using a threshold of 0.75 for the correlation coefficients?
# 2). Use the lm() function to fit univariate regression models to predict patient age using two variables that are highly correlated as predictors. Which of these variables are statistically significant predictors of age? Hint: the summary() function can help here.
# 3). Fit a multiple linear regression model predicting patient age using both variables. What happened?
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# High-Dimensional Data Analysis: Challenge 3
#
# This script explores correlations, univariate regression, and multivariate
# regression in the prostate dataset.
# ----------------------------------------------------------------------
# 1. Examine Correlations Between Variables
# ----------------------------------------------------------------------
# Use the cor() function to compute pairwise correlations between all variables 
# in the prostate dataset. Identify pairs of variables with a correlation
# coefficient greater than 0.75 (absolute value).

# Compute the correlation matrix for the dataset
# NOTE: save your matrix if you want to come back to it: cor_matrix <- cor(prostate, use = "complete.obs")

cor(prostate)

# View the correlation matrix
print(cor_matrix)
round(cor(prostate), 2) # rounding helps to visualise the correlations

# As seen above, some variables are highly correlated.
# In particular, the correlation between gleason and pgg45 is equal to 0.75
# NOTE: writing code to to this is easy; e.g. to find variable pairs with high correlations (> 0.75 or < -0.75):
# high_correlations <- which(abs(cor_matrix) > 0.75 & abs(cor_matrix) < 1, arr.ind = TRUE)

# Extract and display the highly correlated variable pairs
high_cor_pairs <- tibble::tibble(
  Var1 = rownames(cor_matrix)[high_correlations[, 1]],
  Var2 = colnames(cor_matrix)[high_correlations[, 2]],
  Correlation = cor_matrix[high_correlations]
)
print(high_cor_pairs)


# ----------------------------------------------------------------------
# 2. Fit Univariate Regression Models
# ----------------------------------------------------------------------
# Fitting univariate regression models to predict age using `gleason` and `pgg45` as predictors.
# Use lm() to fit univariate regression models predicting patient age 
# using two highly correlated variables. Evaluate statistical significance
# of each predictor using summary().

model_gleason <- lm(age ~ gleason, data = prostate)
model_pgg45 <- lm(age ~ pgg45, data = prostate)


# Summarize results for each model
summary(model_gleason)  # Check significance of gleason
summary(model_pgg45)  # Check significance of pgg45

# ----------------------------------------------------------------------
# 3. Fit a Multiple Regression Model
# ----------------------------------------------------------------------
# Fit a multiple regression model predicting patient age using both variables.
# Examine what happens when both correlated variables are included.

# Fit a multiple regression model
model_multivar <- lm(age ~ gleason + pgg45, data = prostate)

# Summarize the multiple regression model
summary(model_multivar)

# ----------------------------------------------------------------------
# Using Bioconductor
# ----------------------------------------------------------------------

# This workshop focuses on statistical methods for visualizing and analyzing 
# high-dimensional biological data using Bioconductor packages.

# Bioconductor is an open-source platform designed for analyzing high-throughput genomic data.
# It provides a variety of useful packages and example datasets.
# More details and resources are available at: https://www.bioconductor.org/

# Bioconductor packages can be installed and managed using the BiocManager package.

# Let's load the "minfi" package, a Bioconductor package specifically designed
# for analyzing Illumina Infinium DNA methylation arrays.
# BiocManager::install("minfi")

library("minfi")
browseVignettes("minfi")

#methylation <- readRDS("methylation.rds")
methylation <- readRDS(here("data/methylation.rds"))
head(colData(methylation))

methyl_mat <- t(assay(methylation))
## calculate correlations between cells in matrix
cor_mat <- cor(methyl_mat)

cor_mat[1:10, 1:10] # print the top-left corner of the correlation matrix

# ----------------------------------------------------------------------
# Observations
# ----------------------------------------------------------------------
# - Highly correlated variables may cause issues in multiple regression due
#   to multicollinearity.
# - Look for inflated standard errors, changes in significance, or opposite
#   coefficient signs when both variables are included.
# - Dimensionality reduction techniques (e.g., PCA) or regularization methods
#   (e.g., ridge regression, LASSO) can help address multicollinearity.


# ----------------------------------------------------------------------
# Part 2: Regression with many outcomes using DNAm data
# ----------------------------------------------------------------------

library("here")
library("minfi")
methylation <- readRDS(here("data/methylation.rds"))

# Check the dimensions of the dataset using dim().
# Note: In computational biology data structures in R, observations are stored as columns 
# and features (e.g., genomic sites) are stored as rows.
# This contrasts with typical tabular data where features are columns and observations are rows.

# Access assay data (e.g., normalised methylation levels) using assay().

# Retrieve sample-level information using colData().
dim(methylation)

# The output shows that the object has dimensions of 5000 × 37,
# indicating 5000 features and 37 observations.
# To extract the matrix of methylation M-values, use the assay() function.

methyl_mat <- assay(methylation)

hist(methyl_mat, xlab = "M-value")

head(colData(methylation))

# Association between age and DNA methylation:
# The following heatmap summarises age and methylation levels available in the methylation dataset:

library("ComplexHeatmap")

age <- methylation$Age

# sort methylation values by age 
order <- order(age)
age_ord <- age[order]
methyl_mat_ord <- methyl_mat[, order]

# plot heatmap
Heatmap(methyl_mat_ord,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "M-value",
        row_title = "Feature", 
        column_title =  "Sample", 
        top_annotation = columnAnnotation(age = age_ord))


# ----------------------------------------------------------------------
# CHALLENGE 1
# ----------------------------------------------------------------------
# Why can't we simply fit many linear regression models for every combination of features (colData and assays)
# and draw conclusions based on significant p-values?

# Solution:
# There are several problems with this approach:

# 1. Multiple Testing Problem:
#    - If we perform 5000 tests for each of 14 variables, even with no true associations,
#      random noise would produce some significant (spurious) results.

# 2. Small Sample Size:
#    - Some covariates may have very small sample sizes (e.g., certain ethnicities),
#      leading to unreliable or spurious findings.

# 3. Lack of Research Focus:
#    - Without a clear research question, interpreting each model becomes ambiguous.
#    - Rationalising findings after the fact can lead to creating unsupported "stories"
#      based solely on significant results.

# Arbitrarily select the first CpG in the methyl_mat matrix (the one on its first row):

age <- methylation$Age

# methyl_mat[1, ] indicates that the 1st CpG will be used as outcome variable
lm_age_methyl1 <- lm(methyl_mat[1, ] ~ age)
lm_age_methyl1

# We now have estimates for the expected methylation level when age equals 0 (the intercept) 
# and the change in methylation level for a unit change in age (the slope).
# We could plot this linear model:

plot(age, methyl_mat[1, ], xlab = "Age", ylab = "Methylation level", pch = 16)
abline(lm_age_methyl1)

# For this linear model, we can use tidy() from the broom package to extract detailed
# information about the coefficients and the associated hypothesis tests in this model:

library("broom")
tidy(lm_age_methyl1)

# ----------------------------------------------------------------------
# CHALLENGE 2
# ----------------------------------------------------------------------

# QUESTION: In the model we fitted, the estimate for the intercept is 0.902 
# and its associated p-value is 0.0129. What does this mean?

# SOLUTION:
# The first coefficient, the intercept, represents the mean methylation value 
# for the first CpG when age is zero (estimated as 0.902).
# However, this is not meaningful since there are no observations with age zero 
# or below 20 in the dataset.
# The p-value tests whether the intercept (β₀) is zero, 
# but this is not relevant since methylation levels at age zero are not of interest.
# The key focus is on the regression coefficient for age, 
# as it shows whether there is a linear relationship between age and methylation.

# ----------------------------------------------------------------------
# What is a model matrix?

design_age <- model.matrix(lm_age_methyl1) # model matrix
head(design_age)
dim(design_age)

# The model matrix has the same number of rows as the methylation data has samples.
# It has two columns: one for the intercept (like in the linear model we fit earlier)
# and one for age.

# When using lm(), this step happens automatically, but here we specify the model matrix directly.
# For more complex experimental designs, the limma user manual provides guidance on creating model matrices.
# In this case, we use a simple two-variable model.

# We use the lmFit() function to fit the model, passing the methylation data and model matrix.
# Internally, lmFit() efficiently runs lm() for each row of the data.

# The eBayes() function, applied to the lmFit() output, performs pooled standard error estimation,
# resulting in moderated t-statistics and p-values.

library("limma")
design_age <- model.matrix(lm_age_methyl1) # model matrix
fit_age <- lmFit(methyl_mat, design = design_age)
fit_age <- eBayes(fit_age)

# Use the topTable() function to obtain the results of the linear models.
# By default, topTable() returns results for the first coefficient in the model.
# The first coefficient corresponds to the intercept term, which is not of interest here.
# Specify coef = 2 to focus on the second coefficient (age in this case).

# By default, topTable() returns only the top 10 results.
# To see all results, set number = nrow(fit_age), ensuring a row for every input row.

toptab_age <- topTable(fit_age, coef = 2, number = nrow(fit_age))
head(toptab_age)

# The output of topTable includes several columns:
# - logFC: The coefficient estimate (log fold change), representing effect size.
# - aveExpr: The average expression level.
# - t: The t-statistic for the coefficient.
# - P.Value: The p-value for the test.
# - adj.P.Val: The adjusted p-value (we'll discuss adjusted p-values shortly).
# - B: The log-odds that a feature is significantly different (a transformation of the p-value, not covered here).
# The term logFC is used for historical reasons from microarray experiments.

# These results provide effect sizes and p-values for the association between methylation levels
# at each locus and age across the 37 samples.

# To visualise effect sizes (coefficients) and statistical significance (p-values),
# we can plot effect sizes against p-values for all linear models.
# Such plots are called "volcano plots" because their shape resembles an eruption.

plot(toptab_age$logFC, -log10(toptab_age$P.Value),
     xlab = "Effect size", ylab = bquote(-log[10](p-value)),
     pch = 19
)

# In this figure:
# - Each point represents a feature of interest.
# - The x-axis shows the effect size from a linear model.
# - The y-axis shows −log10(p-value), where higher values indicate stronger statistical evidence
#   of a non-zero effect size.

# - Positive effect sizes indicate increasing methylation with age.
# - Negative effect sizes indicate decreasing methylation with age.
# - Points higher on the y-axis represent features with results unlikely under the null hypothesis.

# The goal is to identify features with different methylation levels across age groups.
# Ideally, there would be a clear separation between “null” (no effect) and “non-null” (effect exists) features.
# However, in practice, we often see a continuum of effect sizes and p-values without clear separation.

# Statistical methods can provide insights from these continuous measures,
# but it is often useful to generate a list of features with confident non-zero effect sizes.
# This is challenging due to the large number of tests performed.

# ----------------------------------------------------------------------
# CHALLENGE 3
# ----------------------------------------------------------------------

# 1. Try fitting a linear model using smoking status as a covariate instead of age,
# and create a volcano plot.
# Note: Smoking status is stored as methylation$smoker.

# 2. In the lecture example, we saw that information sharing can result in larger p-values.
# Why might this be preferable?

# SOLUTION:
# 1. 
design_smoke <- model.matrix(~methylation$smoker)
fit_smoke <- lmFit(methyl_mat, design = design_smoke)
fit_smoke <- eBayes(fit_smoke)
toptab_smoke <- topTable(fit_smoke, coef = 2, number = nrow(fit_smoke))
plot(toptab_smoke$logFC, -log10(toptab_smoke$P.Value),
     xlab = "Effect size", ylab = bquote(-log[10](p)),
     pch = 19
)

# 2. 
# Being more conservative when identifying features can help reduce false discoveries.
# It is also important to be cautious when rejecting the null hypothesis based on 
# small standard errors caused by abnormally low variability for certain features.

# ----------------------------------------------------------------------
# The problem of multiple tests:
# With a large number of features, it is important to determine which features are "interesting" 
# or "significant" for further study.

# Using a standard significance threshold of 0.05 may lead to many false positives.
# A p-value threshold of 0.05 means there is a 1 in 20 chance of observing results 
# as extreme or more extreme under the null hypothesis (no association between age and methylation).

# If we perform many more than 20 tests, we are likely to observe significant p-values 
# purely due to random chance, even when the null hypothesis is true.

# To illustrate this, we can permute (scramble) the age values and rerun the test 
# to see how random chance affects the results.

set.seed(123) 
age_perm <- age[sample(ncol(methyl_mat), ncol(methyl_mat))]
design_age_perm <- model.matrix(~age_perm)

fit_age_perm <- lmFit(methyl_mat, design = design_age_perm)
fit_age_perm <- eBayes(fit_age_perm)
toptab_age_perm <- topTable(fit_age_perm, coef = 2, number = nrow(fit_age_perm))

plot(toptab_age_perm$logFC, -log10(toptab_age_perm$P.Value),
     xlab = "Effect size", ylab = bquote(-log[10](p)),
     pch = 19
)
abline(h = -log10(0.05), lty = "dashed", col = "red")

# A random sequence of ages was generated, so there is no reason to expect 
# a true association between methylation levels and this random sequence.

# Despite this, many features still have p-values lower than the traditional significance threshold (p = 0.05).

# In this example, 226 features are significant at p < 0.05.
# Using this fixed threshold in a real experiment could lead to identifying many features 
# as associated with age, even though the results are simply due to random chance.

# When performing multiple tests, features are classified as "significant" or "non-significant".
# However, with many tests, some classifications will inevitably be incorrect.

# Results can fall into four categories:
# - True Positive: Truly different and labeled as different.
# - False Negative: Truly different but labeled as not different.
# - False Positive (False Discovery): Not truly different but labeled as different.
# - True Negative: Not truly different and labeled as not different.

# At a 5% significance level, 5% of results will be false positives (false discoveries) by chance,
# as p-values are uniformly distributed under the null hypothesis.

# To control false discoveries, one method is the Bonferroni correction:
# - Adjust the significance threshold by dividing it by the number of tests (n).
# - Alternatively, multiply p-values by the number of tests.
# - Bonferroni is highly conservative, especially with many features.

p_raw <- toptab_age$P.Value
p_fwer <- p.adjust(p_raw, method = "bonferroni")
plot(p_raw, p_fwer, pch = 16, log="xy")
abline(0:1, lty = "dashed")
abline(v = 0.05, lty = "dashed", col = "red")
abline(h = 0.05, lty = "dashed", col = "red")

# ----------------------------------------------------------------------
# CHALLENGE 4
# ----------------------------------------------------------------------
# 1. At a significance level of 0.05, with 100 tests, what is the Bonferroni significance threshold?

# A = threshold = 0.05 / 100.

# 2. In a gene expression experiment, after FDR correction with an adjusted p-value threshold of 0.05,
# 500 significant genes are observed. What proportion of these genes are truly different?

# A = We can’t say what proportion of these genes are truly different. However, if we repeated this
# experiment and statistical test over and over, on average 5% of the results from each run
# would be false discoveries.

# Try applying FDR correction to the p_raw vector.
# Hint: Use the p.adjust() function and check help("p.adjust") for details on the method.
# A = The following code runs FDR correction and compares it to non-corrected values and to Bonferroni:

p_fdr <- p.adjust(p_raw, method = "BH")
plot(p_raw, p_fdr, pch = 16, log="xy")
abline(0:1, lty = "dashed")
abline(v = 0.05, lty = "dashed", col = "red")
abline(h = 0.05, lty = "dashed", col = "red")

# plot of chunk plot-fdr-fwer
plot(p_fwer, p_fdr, pch = 16, log="xy")
abline(0:1, lty = "dashed")
abline(v = 0.05, lty = "dashed", col = "red")
abline(h = 0.05, lty = "dashed", col = "red")

# ----------------------------------------------------------------------
# SUMMARY
# ----------------------------------------------------------------------

# Key Points

# 1. Performing linear regression in high-dimensional data requires hypothesis testing,
#    which is less critical in low-dimensional regression.

# 2. Sharing information between features can improve statistical power and reduce false positives.

# 3. Multiple testing correction is essential when performing many hypothesis tests
#    on high-dimensional data, helping to retain power while avoiding costly false discoveries.

# 4. The choice of multiple testing method depends on goals:
#    - Conservative methods (e.g., Bonferroni) minimize false positives.
#    - Liberal methods (e.g., FDR) balance power and error rates.

# ----------------------------------------------------------------------
# REGULARISATION
# ----------------------------------------------------------------------

library("here")
library("minfi")
methylation <- readRDS(here("data/methylation.rds"))

## here, we transpose the matrix to have features as rows and samples as columns
methyl_mat <- t(assay(methylation))
age <- methylation$Age

# by using methyl_mat in the formula below, R will run a multivariate regression
# model in which each of the columns in methyl_mat is used as a predictor. 
fit <- lm(age ~ methyl_mat)
summary(fit)

# Singularities
xtx <- t(methyl_mat) %*% methyl_mat
det(xtx)

# Correlated features

library("ComplexHeatmap")
small <- methyl_mat[, 1:500]
cor_mat <- cor(small)
Heatmap(cor_mat,
        column_title = "Feature-feature correlation in methylation data",
        name = "Pearson correlation",
        show_row_dend = FALSE, show_column_dend = FALSE,
        show_row_names = FALSE, show_column_names = FALSE
)

### Model selection using training and test sets

methylation <- readRDS(here::here("data/methylation.rds"))

library("SummarizedExperiment")
age <- methylation$Age
methyl_mat <- t(assay(methylation))

# Subset CpGs

cpg_markers <- c("cg16241714", "cg14424579", "cg22736354", "cg02479575", "cg00864867", 
                 "cg25505610", "cg06493994", "cg04528819", "cg26297688", "cg20692569", 
                 "cg04084157", "cg22920873", "cg10281002", "cg21378206", "cg26005082", 
                 "cg12946225", "cg25771195", "cg26845300", "cg06144905", "cg27377450"
)

horvath_mat <- methyl_mat[, cpg_markers]

## Generate an index to split the data
set.seed(42)
train_ind <- sample(nrow(methyl_mat), 25)

## Split the data 
train_mat <- horvath_mat[train_ind, ]
train_age <- age[train_ind]
test_mat <- horvath_mat[-train_ind, ]
test_age <- age[-train_ind]

## Fit a linear model
# as.data.frame() converts train_mat into a data.frame
# Using the `.` syntax above together with a `data` argument will lead to
# the same result as using `train_age ~ train_mat`: R will fit a multivariate 
# regression model in which each of the columns in `train_mat` is used as 
# a predictor. We opted to use the `.` syntax because it will help us to 
# obtain model predictions using the `predict()` function. 

fit_horvath <- lm(train_age ~ ., data = as.data.frame(train_mat))

## Function to calculate the (mean squared) error
mse <- function(true, prediction) { 
  mean((true - prediction)^2) 
} 

## Calculate the training error 
err_lm_train <- mse(train_age, fitted(fit_horvath)) 
err_lm_train

# ----------------------------------------------------------------------
# CHALLENGE 3
# ----------------------------------------------------------------------

pred_lm <- predict(fit_horvath, newdata = as.data.frame(test_mat)) 

err_lm <- mse(test_age, pred_lm)
err_lm

par(mfrow = c(1, 1))
plot(test_age, pred_lm, pch = 19)
abline(coef = 0:1, lty = "dashed")


# Ridge regression
library("glmnet")

## glmnet() performs scaling by default, supply un-scaled data:
horvath_mat <- methyl_mat[, cpg_markers] # select the same 20 sites as before
train_mat <- horvath_mat[train_ind, ] # use the same individuals as selected before
test_mat <- horvath_mat[-train_ind, ]

ridge_fit <- glmnet(x = train_mat, y = train_age, alpha = 0)
plot(ridge_fit, xvar = "lambda")
abline(h = 0, lty = "dashed")

# Obtain a matrix of predictions from the ridge model,
# where each column corresponds to a different lambda value
pred_ridge <- predict(ridge_fit, newx = test_mat)

# Calculate MSE for every column of the prediction matrix against the vector of true ages
err_ridge <- apply(pred_ridge, 2, function(col) mse(test_age, col)) 
min_err_ridge <- min(err_ridge)

# Identify the lambda value that results in the lowest MSE (ie, the "best" lambda value)
which_min_err <- which.min(err_ridge)
pred_min_ridge <- pred_ridge[, which_min_err]

## Return errors
min_err_ridge

err_lm 

chosen_lambda <- ridge_fit$lambda[which.min(err_ridge)]
plot(ridge_fit, xvar = "lambda")
abline(v = log(chosen_lambda), lty = "dashed")

# ----------------------------------------------------------------------
# CHALLENGE 4
# ----------------------------------------------------------------------

# 1. Which performs better, ridge or OLS?
# 2. Plot predicted ages for each method against the true ages.
# How do the predictions look for both methods? Why might ridge be performing better?

#   1. Ridge regression performs significantly better on unseen data, 
# despite being “worse” on the training data.
min_err_ridge
err_lm


# 2. The ridge ones are much less spread out with far fewer extreme predictions.

all <- c(pred_lm, test_age, pred_min_ridge)
lims <- range(all)
par(mfrow = 1:2)
plot(test_age, pred_lm,
     xlim = lims, ylim = lims,
     pch = 19
)
abline(coef = 0:1, lty = "dashed")
plot(test_age, pred_min_ridge,
     xlim = lims, ylim = lims,
     pch = 19
)
abline(coef = 0:1, lty = "dashed")


### LASSO Regularization

# LASSO regularization uses the **L1 norm** (sum of absolute coefficient values) to shrink coefficients:
  
# To balance model complexity and information retention, we choose an optimal λ. 
# A common method is cross-validation, which splits the data into K chunks. 
# K−1 chunks are used for training, and the remaining chunk is for testing.
# This process rotates through all chunks, providing a reliable estimate of how 
# well each λ value generalizes to new data.
#
# We can use this new idea to choose a lambda value by finding the lambda that minimises 
# the error across each of the test and training splits. In R:

# fit lasso model with cross-validation across a range of lambda values
lasso <- cv.glmnet(methyl_mat, age, alpha = 1)
plot(lasso)

# Extract the coefficients from the model with the lowest mean squared error from cross-validation
coefl <- coef(lasso, lasso$lambda.min)
# select only non-zero coefficients
selection <- which(coefl != 0)
# and convert to a normal matrix
selected_coefs <- as.matrix(coefl)[selection, 1]
selected_coefs

# We can see that cross-validation has selected a value of λ
# resulting in 44 features and the intercept.



# ----------------------------------------------------------------------
# ELASTIC NET REGRESSION
# ----------------------------------------------------------------------

# Elastic net regression combines the properties of ridge (alpha = 0) 
# and LASSO (alpha = 1) regression, blending their advantages. 
# It drops uninformative variables like LASSO while maintaining 
# conservative coefficient estimates like ridge, leading to improved predictions.

# ----------------------------------------------------------------------
# CHALLENGE 5
# ----------------------------------------------------------------------

# 1. Fit an elastic net model (hint: alpha = 0.5) without cross-validation and plot the model object.
# 2. Fit an elastic net model with cross-validation and plot the error. Compare with LASSO.

elastic <- glmnet(methyl_mat[, -1], age, alpha = 0.5)
plot(elastic)

#The process of model selection is similar for elastic net models as for LASSO models.

elastic_cv <- cv.glmnet(methyl_mat[, -1], age, alpha = 0.5)
plot(elastic_cv)

