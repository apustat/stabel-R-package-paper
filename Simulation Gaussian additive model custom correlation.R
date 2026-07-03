library(stabel)
library(SuperLearner)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)
library(SIS)
library(MASS)


stabel_wrapper=function(n, p, cutoff, iter, train.pcntg, SL.library, seed, p1){
  cl <- parallel::makeCluster(min(4, detectCores()-1), setup_strategy = "sequential")
  # Activate cluster for foreach library
  registerDoParallel(cl)
  
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  #####This part is just for calculating correlation matrix#####
  df=read.csv("C:/Users/.../train_10000_genes.csv") #training data with 10000 genes
  n0 <- nrow(df)
  p0 <- ncol(df) #including Y
  X=as.matrix(df[, -p0])
  Y=df[, p0]
  sis_mod<- SIS(X, Y, family='binomial', penalty = "lasso", nsis=p, iter.max =100, iter = FALSE)
  sig_var=sis_mod$sis.ix0 #significant variables by SIS
  
  sel_X=X[, sig_var, drop = FALSE] # data with the variables from SIS]
  corr.str=cor(sel_X)
  ######################################
  
  stabel_iter=foreach(i = 1:iter, .packages = c("MASS", "sparseSVM", "Boruta", "SuperLearner", "stabel", "glmnet", "stabs", "SIS")) %dopar%{
    set.seed(1000+i)
    
    X <- mvrnorm(
      n = n,
      mu = rep(0, p),
      Sigma = corr.str
    )
    
    # Define the functions
    f1 <- function(x) sin(2 * pi * x) + 0.5 * x^3 # f1: Oscillating function
    f2 <- function(x) x^2 # f2: Quadratic function
    f3 <- function(x) sin(2 * x) + x^3 - x # f3: Arbitrary smooth function
    f4 <- function(x) 1 / (1 + exp(-x)) - 0.5 # f4: Smooth increasing function with plateau (scaled sigmoid function)
    f5 <- function(x) cos(2 * pi * x) # f5: Cosine function
    f6 <- function(x) ifelse(x < -1, -1, ifelse(x > 1, 1, x)) # f6: Piecewise linear function
    f7 <- function(x) x # f7: Linear function with positive slop
    f8 <- function(x) -x # f8: Linear function with negative slope
    
    # Apply the functions to selected covariates
    eta <- f1(X[, 1]) + f2(X[, 2]) + f3(X[, 3]) + f4(X[, 4]) + f5(X[, 5]) + f6(X[, 6]) + f7(X[, 7]) + f8(X[, 8])
    
    colnames(X)=paste0("X", 1:p)
    # Generate Y
    sigma2 <- 1  # Variance of Y
    Y <- rnorm(n, mean = eta, sd = sqrt(sigma2))
    
    train=sample(1:n, n*train.pcntg, replace=FALSE)
    test=(-train)
    
    pre_train_X=as.matrix(X[train, , drop = FALSE]) #training data with all variables
    train_Y=Y[train]
    
    sis_mod=SIS(pre_train_X, train_Y, family='gaussian', iter = FALSE, penalty = "lasso", nsis=n/log(n)) 
    sig_var=sort(sis_mod$sis.ix0) #significant variables by SIS
    
    train_X=as.matrix(X[train, sig_var, drop = FALSE]) #training data with the variables from SIS
    
    test_X=as.matrix(X[test, sig_var, drop = FALSE])
    test_Y=Y[test]
    
    ########Cross-validation
    #cv=cv.stabel(X=train_X,Y=train_Y,family = "gaussian",fit.fun = "Lasso",eval.metric.Lasso = "mse",nfolds = 5,seed = seed)
    
    ##########Variable selection
    vs=stabel.vs(X=train_X,Y=train_Y,cutoff = cutoff,B = 50,bestlam.Lasso = 0.14, maxit.Lasso = 1e+05,
                 family = "gaussian",dfmax = 10, ntree=200, mcAdj=TRUE, maxRuns = 100, pValue = 0.05,
                 fit.fun = c("Lasso", "RF"),comb.method = "average",seed = seed)
    
    vs.avg <- colnames(train_X)[
      if (length(vs$final.selected.set) == 0) 1 else vs$final.selected.set
    ] #STABEL with average
    
    
    vs.union <- colnames(train_X)[
      if (length(sort(union(vs$selected.variables$Lasso,
                            vs$selected.variables$RF))) == 0) {
        1
      } else {
        sort(union(vs$selected.variables$Lasso,
                   vs$selected.variables$RF))
      }
    ] #STABEL with union
    
    vs.SS.lasso <- colnames(train_X)[
      if (length(sort(vs$selected.variables$Lasso)) == 0) {
        1
      } else {
        sort(vs$selected.variables$Lasso)
      }
    ] #stability selection with Lasso
    # if no variable is selected, we will use the first variable to avoid error in the prediction step.
    ##################################
    true.b=paste0("X", 1:8)
    
    tp.avg <- length(intersect(vs.avg, true.b)) #True Positives
    fp.avg <- length(setdiff(vs.avg, true.b)) #False Positives
    fn.avg <- length(setdiff(true.b, vs.avg)) #False Negatives
    tn.avg <- length(setdiff(1:p1, union(true.b, vs.avg))) #True negatives
    TPR.avg <- tp.avg/(tp.avg+fn.avg) #proportion of true covariates selected in the final model.
    FDR.avg <- fp.avg/(fp.avg+tp.avg) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    
    tp.union <- length(intersect(vs.union, true.b)) #True Positives
    fp.union <- length(setdiff(vs.union, true.b)) #False Positives 
    fn.union <- length(setdiff(true.b, vs.union)) #False Negatives
    tn.union <- length(setdiff(1:p1, union(true.b, vs.union))) #True negatives
    TPR.union <- tp.union/(tp.union+fn.union) #proportion of true covariates selected in the final model.
    FDR.union <- fp.union/(fp.union+tp.union) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    
    tp.SS.lasso <- length(intersect(vs.SS.lasso, true.b)) #True Positives
    fp.SS.lasso <- length(setdiff(vs.SS.lasso, true.b)) #False Positives 
    fn.SS.lasso <- length(setdiff(true.b, vs.SS.lasso)) #False Negatives
    tn.SS.lasso <- length(setdiff(1:p1, union(true.b, vs.SS.lasso))) #True negatives
    TPR.SS.lasso <- tp.SS.lasso/(tp.SS.lasso+fn.SS.lasso) #proportion of true covariates selected in the final model.
    FDR.SS.lasso <- fp.SS.lasso/(fp.SS.lasso+tp.SS.lasso) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    
    ###################Prediction
    pred.avg=stabel.pred(X=train_X[, vs.avg, drop = FALSE],Y=train_Y, newX = test_X[, vs.avg, drop = FALSE],newY = test_Y,method = "method.NNLS",
                         SL.library=SL.library,
                         family = "gaussian", nfolds = 5, seed = seed)
    
    pred.union=stabel.pred(X=train_X[, vs.union, drop = FALSE],Y=train_Y,newX = test_X[, vs.union, drop = FALSE],newY = test_Y,method = "method.NNLS",
                           SL.library=SL.library,
                           family = "gaussian", nfolds = 5, seed = seed)
    
    pred.SS.lasso=stabel.pred(X=train_X[, vs.SS.lasso, drop = FALSE],Y=train_Y,newX = test_X[, vs.SS.lasso, drop = FALSE],newY = test_Y,method = "method.NNLS",
                              SL.library="SL.lm",
                              family = "gaussian", nfolds = 5, seed = seed)
    
    ######################Traditional Lasso
    set.seed(200+i)
    t.lasso=glmnet(x=train_X,y=train_Y, family="gaussian", alpha=1, lambda=0.14, maxit=100000, intercept=TRUE)
    vs.t.lasso=colnames(train_X)[which(coef(t.lasso)[-1] !=0)] #non-zero coefficients
    
    tp.t.lasso <- length(intersect(vs.t.lasso, true.b)) #True Positives
    fp.t.lasso <- length(setdiff(vs.t.lasso, true.b)) #False Positives
    fn.t.lasso <- length(setdiff(true.b, vs.t.lasso)) #False Negatives
    tn.t.lasso <- length(setdiff(1:p1, union(true.b, vs.t.lasso))) #True negatives
    TPR.t.lasso= tp.t.lasso/(tp.t.lasso+fn.t.lasso) #proportion of true covariates selected in the final model.
    FDR.t.lasso= fp.t.lasso/(fp.t.lasso+tp.t.lasso) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    
    pred.t.lasso=stabel.pred(X=train_X[, vs.t.lasso, drop = FALSE],Y=train_Y, newX = test_X[, vs.t.lasso, drop = FALSE],newY = test_Y,method = "method.NNLS",
                             SL.library="SL.lm",
                             family = "gaussian", nfolds = 5, seed = seed)
    
    c(TPR.avg, FDR.avg, TPR.union, FDR.union, TPR.SS.lasso, FDR.SS.lasso, TPR.t.lasso, FDR.t.lasso,
      pred.avg$rmse, pred.avg$mae,
      pred.union$rmse, pred.union$mae,
      pred.SS.lasso$rmse, pred.SS.lasso$mae,
      pred.t.lasso$rmse, pred.t.lasso$mae)
    
  }
  pred.sum=as.data.frame(do.call(rbind, stabel_iter))
}

result=stabel_wrapper(n=400, p=1000, cutoff=0.7, iter=10, train.pcntg=0.6, SL.library=c("SL.randomForest", "SL.lm", "SL.ksvm"), seed=10, p1=8)

colnames(result) <- c("TPR.avg", "FDR.avg", "TPR.union", "FDR.union", "TPR.SS.lasso", "FDR.SS.lasso", "TPR.t.lasso", "FDR.t.lasso",
                      "rmse.avg", "mae.avg",
                      "rmse.union", "mae.union",
                      "rmse.SS.lasso", "mae.SS.lasso",
                      "rmse.t.lasso", "mae.t.lasso")