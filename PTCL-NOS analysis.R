library(dplyr)
library(tidyr)
library(tibble)
library(caret)
library(SIS)
library(stabel)
library(SuperLearner)
library(pROC)
library(glmnet)
library(stabs)

### ---------------------------
### 1. Load raw data
### ---------------------------
df <- read.table(
  "C:/Users/.../Normalized GEP Expression Data.txt",
  header = TRUE,
  sep = "\t"
)

## Replace empty gene names with unique identifiers
df1 <- df %>%
  mutate(Gene = ifelse(Gene == "" | is.na(Gene),
                       paste0("XX", row_number()),
                       Gene))

## Sort by gene name
df2 <- df1[order(df1$Gene), ]

## Remove duplicate genes (keep first occurrence)
df_cleaned <- df2 %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

### ---------------------------
### 2. Transpose data
### ---------------------------
df_transposed <- df_cleaned %>%
  select(-ProbeSet) %>%
  column_to_rownames(var = "Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")


### ---------------------------
### 3. Outcome data
### ---------------------------
df_res <- read.csv(
  "C:/Users/.../Gene Expression Descritpions with Cohorts.csv"
)

molecular.ptcl_nos <- subset(df_res, Final.Molecular == "PTCL-NOS")

molecular.Tbet_GATA <- subset(
  molecular.ptcl_nos,
  PTCL.NOS..GATA3.TBX21 != "Unclass"
)[, c("Experiment.Name", "PTCL.NOS..GATA3.TBX21")]

colnames(molecular.Tbet_GATA) <- c("Sample", "Response")
molecular.Tbet_GATA$Response <- ifelse(molecular.Tbet_GATA$Response == "GATA", 0, 1)

## Merge outcome after transposing
final.data <- merge(df_transposed, molecular.Tbet_GATA,
                    by = "Sample",
                    all.y = TRUE)


### ---------------------------
### 4. Train-test split
### ---------------------------
set.seed(199055)

n <- nrow(final.data)
train_idx <- sample(seq_len(n), floor(0.6 * n))

train_data_full <- final.data[train_idx, ]
test_data_full  <- final.data[-train_idx, ]


### ---------------------------
### 5. Variance filtering on TRAINING set only
### ---------------------------
X_train_full <- train_data_full[, !(names(train_data_full) %in% c("Sample", "Response"))]
Y_train <- train_data_full$Response

## Calculate variance on training data only
gene_var <- apply(X_train_full, 2, var, na.rm = TRUE)

## Keep top 10,000 most variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:10000]

X_train <- X_train_full[, top_genes, drop = FALSE]

## Apply same genes to testing set
X_test <- test_data_full[, top_genes, drop = FALSE]
Y_test <- test_data_full$Response

### ---------------------------
### 6. SIS on training set only
### ---------------------------
sis_mod <- SIS(
  as.matrix(X_train),
  Y_train,
  family = "binomial",
  penalty = "lasso",
  tune = "cv",
  nsis = 30,
  iter.max = 100,
  iter = FALSE
)

sig_var <- sis_mod$sis.ix0

## Final training and testing datasets with same 30 genes
X_train_sis <- as.matrix(X_train[, sig_var, drop = FALSE])
X_test_sis  <- X_test[, sig_var, drop = FALSE]

cutoff=0.70
########################################################
##############Data analysis#############################

#cv=cv.stabel(X=X_train_sis,Y=Y_train,family = "binomial",fit.fun = c("Lasso","sparseSVM"),eval.metric.Lasso = "mse",
             #eval.metric.sparseSVM = "me", nfolds = 4, seed=500)

vs=stabel.vs(X=X_train_sis,Y=Y_train,cutoff = cutoff,B = 50, bestlam.Lasso = 0.0003, maxit.Lasso = 1e+05,
             bestlam.sparseSVM = 0.63, maxit.sparseSVM=1000, gamma=0.1,
             family = "binomial", dfmax = 10, ntree=300, mcAdj=TRUE, maxRuns = 300, pValue = 0.05,
             fit.fun = c("Lasso", "sparseSVM", "RF"),comb.method = "average", seed = 10)

vs.avg=colnames(X_train_sis)[vs$final.selected.set] #STABEL with average

pred.avg=stabel.pred(X=X_train_sis[, vs.avg, drop = FALSE],Y=Y_train, newX = as.matrix(X_test_sis[, vs.avg, drop = FALSE]),
                     newY = Y_test, method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), 
                     family = "binomial",thr.prob=NULL, use.youden = TRUE, nfolds = 4, target.specificity=0.985, seed = 10)

pred.RF=stabel.pred(X=X_train_sis[, vs.avg, drop = FALSE],Y=Y_train, newX = as.matrix(X_test_sis[, vs.avg, drop = FALSE]),
                    newY = Y_test, method = "method.NNLS", SL.library="SL.randomForest", 
                    family = "binomial",thr.prob=NULL, use.youden = TRUE, nfolds = 4, target.specificity=0.985, seed = 10)

#######Lasso
set.seed(100)
cvfit.lasso <- cv.glmnet(x=X_train_sis, y= Y_train, family = "binomial", alpha=1, type.measure="mse", nfolds=4) #finding the best value of lambda
bestlam.t.lasso <- cvfit.lasso$lambda.min
t.lasso=glmnet(x=X_train_sis,y=Y_train, family="binomial", alpha=1, lambda=bestlam.t.lasso, maxit=100000, intercept=TRUE)
vs.t.lasso=colnames(X_train_sis)[which(coef(t.lasso)[-1] !=0)] #non-zero coefficients

pred.lasso=stabel.pred(X=X_train_sis[, vs.t.lasso, drop = FALSE],Y=Y_train,newX = as.matrix(X_test_sis[, vs.t.lasso, drop = FALSE]),newY = Y_test,
                       method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), family = "binomial", 
                       thr.prob=NULL, use.youden = TRUE, nfolds = 4,target.specificity=0.985, seed = 10)

####Modified Lasso
set.seed(100)
lambda=930*bestlam.t.lasso #selects exactly 5 genes same as STABEL
m.lasso=glmnet(x=X_train_sis,y=Y_train, family="binomial", alpha=1, lambda=lambda, maxit=100000, intercept=TRUE)
vs.m.lasso=colnames(X_train_sis)[which(coef(m.lasso)[-1] !=0)] #non-zero coefficients

pred.m.lasso=stabel.pred(X=X_train_sis[, vs.m.lasso, drop = FALSE],Y=Y_train,newX = as.matrix(X_test_sis[, vs.m.lasso, drop = FALSE]),newY = Y_test,
                         method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), family = "binomial", 
                         thr.prob=NULL, use.youden = TRUE, nfolds = 4,target.specificity=0.985, seed = 10)

###########Stability selection with LASSO (stabs)################
set.seed(100)
SS.lasso <- stabsel(x = X_train_sis, y = Y_train, B=50, fitfun = glmnet.lasso, cutoff = 0.70,
                    PFER = 3, sampling.type = "SS")
vs.stabs <- colnames(X_train_sis)[SS.lasso$selected]

pred.SS.lasso=stabel.pred(X=X_train_sis[, vs.stabs, drop = FALSE],Y=Y_train, newX = as.matrix(X_test_sis[, vs.stabs, drop = FALSE]),newY = Y_test,
                          method = "method.NNLS", SL.library="SL.glm", family = "binomial", 
                          thr.prob=NULL, use.youden = TRUE, nfolds = 4,target.specificity=0.985, seed = 10)

pred.avg$predition.accuracy; pred.RF$predition.accuracy; pred.lasso$predition.accuracy; pred.m.lasso$predition.accuracy; pred.SS.lasso$predition.accuracy
pred.avg$auc; pred.RF$auc; pred.lasso$auc; pred.m.lasso$auc; pred.SS.lasso$auc
pred.avg$sensitivity; pred.RF$sensitivity; pred.lasso$sensitivity; pred.m.lasso$sensitivity; pred.SS.lasso$sensitivity
pred.avg$specificity; pred.RF$specificity; pred.lasso$specificity; pred.m.lasso$specificity; pred.SS.lasso$specificity
pred.avg$sensitivity.at.specificity; pred.RF$sensitivity.at.specificity; pred.lasso$sensitivity.at.specificity;  pred.m.lasso$sensitivity.at.specificity; pred.SS.lasso$sensitivity.at.specificity

