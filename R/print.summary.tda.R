print.summary.tda <-
function(x,...){
  # x = obj
    cat("----------------------------------------------\n")
  if(is.null(x$parameters$lambda)) { 
    cat("         Discriminant analysis\n")
  } else{
    
    cat("  Transformation discriminant analysis\n ")
    
  }
  
   cat("-----------------------------------------------\n")
  
  cat("\n")
  cat("Model summary:")
  cat("\n")
  criterion_name <- attr(x, "criterion_name")
  if (is.null(criterion_name)) {
    criterion_name <- attr(x, "criterion")
  }
  if (is.null(criterion_name)) {
    stored_criterion <- x[["criterion", exact=TRUE]]
    criterion_name <- if (is.character(stored_criterion)) stored_criterion else "BIC"
  }
  criterion_value <- x[["criterion", exact=TRUE]]
  
  # Support summary objects created by earlier versions of these methods.
  if (is.character(criterion_value)) {
    criterion_value <- x[["criterion_value", exact=TRUE]]
  }
  if (is.null(criterion_value)) {
    criterion_value <- x[["criterion_value", exact=TRUE]]
  }
  if (is.null(criterion_value)) {
    old_criterion <- attr(x, "criterion_name")
    if (is.null(old_criterion)) {
      old_criterion <- attr(x, "criterion")
    }
    if (is.null(old_criterion)) {
      old_criterion <- "BIC"
    }
    criterion_value <- x[[tolower(old_criterion)]]
  }
  if (is.null(criterion_value)) {
    stop("The summary object does not contain a criterion value")
  }
  
  tab1 <- data.frame("log-likelihood" = x$loglik,
                     "n" = sum(x$n),  
                     "criterion" = criterion_value,
                     row.names = "", check.names = FALSE)
  names(tab1)[3] <- criterion_name
  print(tab1)
  
  
  
  tab2 <-  data.frame("n" = x$n, "%" = round(x$n/sum(x$n)*100,2), 
                      "sub_G" = as.vector(sapply(x$parameters$subgroup_prior, length)),
                      check.names = FALSE,
                      row.names = x$classes)
  
  
  tab2 <- as.matrix(tab2)
  names(dimnames(tab2)) <- c("Classes", "")
  
  cat("\nTraining set confusion matrix:\n")
  print(x$tab)
  
  cat("Training set Misclassification rate =", round(x$Misclassification_Rate,4),"\n")
  cat( "Training set ARI = ", round(x$ARI,4),"\n")

    
  
  if (!is.null(x$Testset_Misclassification_Rate)){ 
    cat("-----------------------------------------------\n")
    
    cat("\nTesting set confusion matrix:\n")
    print(x$Testset_tab )
    cat("\n")
    cat("Testing set Misclassification rate =", round(x$Testset_Misclassification_Rate,4),"\n")
    cat( "Testing set ARI = ", round(x$Testset_ARI_test,4),"\n")
    
  }
  cat("------------------------------------------------\n")
  invisible(x)
  
}
