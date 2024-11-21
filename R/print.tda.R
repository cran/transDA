print.tda <-
function(x,...) {
  
  cat("-----------------------------------\n")
  cat("Discriminant Analysis", "\n")
  # cat(" ", "\n")
  cat("-----------------------------------\n")
  
  if(length(c(x$BIC))==1 ) {
    
    cat("BIC:",x$BIC,"\n")
    cat(" ", "\n")
    cat("prior:",  "\n")
    
    print(x$prior)
    cat(" ", "\n")
    cat("sub_prior:", "\n")
    
    print(x$sub_prior)
    
    
    cat("ARI:", x$ARI, "\n")
    cat(" ", "\n")
    cat("Training Misclassification Error:", x$misclassification_rate, "\n")
    cat(" ", "\n")
    cat("loglik:", x$loglik[length(x$loglik)], "\n")
    cat(" ", "\n")
    cat("iteration:", length(x$loglik), "\n")
    
  } else {
    
    cat("BIC:" , "\n")
    print(x$BIC)
    cat(" ", "\n")
    cat("minimum BIC:", min(na.omit(as.vector(x$BIC))), "\n")
    
    cat(" ", "\n")
    
    cat("prior:",  "\n")
    
    print(x$prior)
    cat(" ", "\n")
    cat("sub_prior:", "\n")
    
    print(x$sub_prior)
    
    cat("ARI:", x$ARI, "\n")
    cat(" ", "\n")
    cat("Training Misclassification Error:", x$misclassification_rate, "\n")
    cat(" ", "\n")
    cat("loglik:", x$loglik[length(x$loglik)], "\n")
    cat(" ", "\n")
    cat("iteration:", length(x$loglik), "\n")
    cat(" ", "\n")
    cat("the subgroup combination of best model:",
        sapply(x$sub_prior, length),"\n")
    
  }
  invisible()
}
