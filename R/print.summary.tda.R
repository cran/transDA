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
  tab1 <- data.frame("log-likelihood" = x$loglik,
                     "n" = sum(x$n),  
                     "BIC" = x$bic, 
                     row.names = "", check.names = FALSE)
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
