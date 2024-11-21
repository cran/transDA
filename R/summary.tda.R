summary.tda <-
function(object, ...){
  # object = MDA
  class <- object$true_ID
  Predicted <- object$predict_ID
  
  
  tab <- try(table(class, Predicted))
  
  names(dimnames(tab)) <- c("Class", "Predicted") 
  
  
  
  loglik = max(object$loglik)
  
  bic = min(na.omit(as.vector(object$BIC)))
  classes <- levels(class)
  nclass <- length(classes)
  n <- as.vector(table(class))
  
  if( !is.null(object$lambda)== TRUE){
    type <- "Transformation discriminant analysis"
    
  }else{
    
    type <- "Discriminant analysis"
  }
  
  if( !is.null(object$lambda)== TRUE){
    par <- list(object$sub_prior,object$mu, object$sigma,object$lambda)
    names(par) <- c("subgroup_prior", "Mean", "Variance","lambda")
  }else{
    par <- list(object$sub_prior,object$mu, object$sigma)
    names(par) <- c("subgroup_prior", "Mean", "Variance")
  }
  
  
  
  
  if(is.null(object$test_true_ID) ){ 
    
    obj <- list(type = type, n = n, 
                loglik = loglik,  bic = bic,
                nclass = nclass, classes = classes,
                prior = object$prior, 
                parameters = par, 
                tab = tab,
                Misclassification_Rate = object$misclassification_rate,
                ARI = object$ARI)
    
  } else{
    
    
    testset_class <- object$test_true_ID
    testset_Predicted <- object$test_pred_ID
    testset_tab <- try(table(testset_class,testset_Predicted))
    names(dimnames(testset_tab)) <- c("Class", "Predicted") 
    
    obj <- list(type = type, n = n, 
                loglik = loglik,  bic = bic,
                nclass = nclass, classes = classes,
                prior = object$prior, 
                parameters = par, 
                tab = tab,
                Misclassification_Rate = object$misclassification_rate,
                ARI = object$ARI,
                Testset_tab = testset_tab, 
                Testset_Misclassification_Rate = object$test_misrate,
                Testset_ARI_test = object$test_ARI)
    
    
  }
  
  
  
  
  class(obj) <- "summary.tda"
  return(obj)
  
  
  
}
