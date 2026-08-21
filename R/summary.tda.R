summary.tda <-
function(object, ...){
  # object = MDA
  class <- object$true_ID
  Predicted <- object$predict_ID
  
  
  tab <- try(table(class, Predicted))
  
  names(dimnames(tab)) <- c("Class", "Predicted") 
  
  
  
  loglik = max(object$loglik)
  
  criterion_name <- attr(object, "criterion_name")
  if (is.null(criterion_name)) {
    criterion_name <- attr(object, "criterion")
  }
  if (is.null(criterion_name)) {
    available_criteria <- intersect(c("BIC", "AIC", "ICL"), names(object))
    criterion_name <- if (length(available_criteria) > 0) available_criteria[1] else "BIC"
  }
  criterion_values <- object[["criterion", exact=TRUE]]
  if (is.null(criterion_values)) {
    criterion_values <- object[[criterion_name]]
  }
  if (is.null(criterion_values)) {
    stop("The fitted object does not contain values for ", criterion_name)
  }
  criterion <- min(na.omit(as.vector(criterion_values)))
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
                loglik = loglik,
                criterion = criterion,
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
                loglik = loglik,
                criterion = criterion,
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
  
  attr(obj, "criterion_name") <- criterion_name
  
  
  
  class(obj) <- "summary.tda"
  return(obj)
  
  
  
}
