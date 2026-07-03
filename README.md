## README Document for stabel

This readme file explains three R scripts, two for simulation studies and one for case study analysis for stabel (Stability Selection with Ensemble Learning) package paper. The method can be run by sourcing the scripts and setting the required input parameters in the R environment. The original methodology paper can be downloaded from https://link.springer.com/article/10.1007/s12561-026-09516-w .

## 1.	Script 1 (Simulation Gaussian additive model.R)

stabel_wrapper(n, p, rho, cutoff, iter, train.pcntg, SL.library, seed, p1, corr.type = c("independent","ar1"))

## Inputs for stable_wrapper

n	Number of observations (sample size) generated in each simulation replicate.

p	Total number of predictor variables generated in the simulated dataset.

rho	Correlation coefficient controlling the dependence among predictors.

cutoff	Stability selection cutoff used by stabel.vs() to determine the final selected variables.

iter	Number of simulation replications to perform.

train.pcntg	Proportion of observations randomly allocated to the training dataset.

SL.library	Vector of Super Learner algorithms used for the prediction stage in stabel.pred().

seed	Random seed used to ensure reproducibility of the STABEL procedures.

p1	Number of true signal variables used when computing TPR and FDR (here, the first 8 predictors).

corr.type	Correlation structure for the simulated predictors. Options are "independent" and "ar1".

## Outputs for stable_wrapper

TPR.avg	True Positive Rate for stabel using the average aggregation rule.

FDR.avg	False Discovery Rate for stabel using the average aggregation rule.

TPR.union	True Positive Rate for stabel using the union aggregation rule.

FDR.union	False Discovery Rate for stabel using the union aggregation rule.

TPR.SS.lasso	True Positive Rate for stability selection with LASSO.

FDR.SS.lasso	False Discovery Rate for stability selection with LASSO.

TPR.t.lasso	True Positive Rate for the traditional LASSO model.

FDR.t.lasso	False Discovery Rate for the traditional LASSO model.

rmse.avg	Root Mean Squared Error (RMSE) of the prediction model built using variables selected by stabel (average aggregation).

mae.avg	Mean Absolute Error (MAE) of the prediction model built using variables selected by stabel (average aggregation).

rmse.union	RMSE of the prediction model built using variables selected by stabel (union aggregation).

mae.union	MAE of the prediction model built using variables selected by STABEL (union aggregation).

rmse.SS.lasso	RMSE of the prediction model built using variables selected by stability selection with LASSO.

mae.SS.lasso	MAE of the prediction model built using variables selected by stability selection with LASSO.

rmse.t.lasso	RMSE of the prediction model built using variables selected by the traditional LASSO model.

mae.t.lasso	MAE of the prediction model built using variables selected by the traditional LASSO model.

## Running time for stable_wrapper

It takes approximately 44 seconds to complete a single iteration (iter=1) when n=200, p=1000, rho=0, and p1=8.

## 2.	Script 2 (Simulation Gaussian additive model custom correlation.R)

All of the arguments described for Script 1 are also applicable to Script 2, with the exception of the “corr.type” argument. In Script 2, no predefined correlation structure is used to generate the simulated data. Instead, the correlation structure is estimated empirically from the training set of the PTCL-NOS dataset and then used to generate the synthetic predictors. The output variables produced by stabel_wrapper() are identical to those in Script 1.

## 3.	Script 3 (PTCL-NOS analysis.R)

This script reproduces the real-data analysis presented in the manuscript using the PTCL-NOS gene expression dataset. The objective is to identify a parsimonious set of genes for distinguishing the GATA3 and TBX21 molecular subtypes of PTCL-NOS and to compare the predictive performance of stabel with several competing variable selection methods.

## Analysis Workflow:

Step	Description

Load gene expression data	Import the normalized gene expression matrix and associated clinical annotations.

Data preprocessing	Replace missing gene names with unique identifiers, sort genes alphabetically, and remove duplicate gene entries while retaining the first occurrence.

Data transformation	Transpose the expression matrix so that rows represent patients and columns represent genes.

Outcome construction	Extract PTCL-NOS samples, remove unclassified cases, and encode the response variable as binary (GATA3 = 0, TBX21 = 1).

Merge datasets	Merge the gene expression matrix with the clinical outcome information.

Training/testing split	Randomly divide the data into 60% training and 40% testing sets using a fixed random seed.

Variance filtering	Calculate gene variances using only the training data and retain the 10,000 most variable genes. The testing set is restricted to the same genes.

Sure Independence Screening (SIS)	Apply SIS to the training data and retain the 30 most informative genes. The same subset is extracted from the testing data.

Variable selection	Apply stabel using LASSO, Random Forest, and sparseSVM base learners. Comparison methods include traditional LASSO, modified LASSO, and stability selection with LASSO (stabs).

Prediction	Build prediction models using the selected variables and evaluate performance on the independent testing dataset using Super Learner.

Performance evaluation	Compare methods using prediction accuracy, AUC, sensitivity, specificity, and sensitivity at the target specificity.

## Variable Selection Methods:

Method	Description

stabel	Stability-based ensemble variable selection using LASSO, Random Forest, and sparseSVM with average aggregation.

Random Forest Model	Prediction using the stabel-selected variables with Random Forest as the only Super Learner algorithm.

Traditional LASSO	Penalized logistic regression with the tuning parameter selected by cross-validation.

Modified LASSO	LASSO with a manually adjusted penalty parameter to select the same number of variables as STABEL.

Stability Selection with LASSO	Stability selection implemented using the stabs package.

## Performance Metrics:

Output	Description

prediction.accuracy	Overall classification accuracy on the testing dataset.

auc	Area under the ROC curve.

sensitivity	True positive rate.

specificity	True negative rate.

sensitivity.at.specificity	Sensitivity evaluated at the pre-specified target specificity (0.985).

