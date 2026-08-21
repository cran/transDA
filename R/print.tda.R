print.tda <-
function(x,...) {
  
  cat("-----------------------------------\n")
  cat("Discriminant Analysis", "\n")
  # cat(" ", "\n")
  cat("-----------------------------------\n")
  
  criterion_name <- attr(x, "criterion_name")
  if (is.null(criterion_name)) {
    criterion_name <- attr(x, "criterion")
  }
  if (is.null(criterion_name)) {
    available_criteria <- intersect(c("BIC", "AIC", "ICL"), names(x))
    criterion_name <- if (length(available_criteria) > 0) available_criteria[1] else "BIC"
  }
  criterion_values <- x[["criterion", exact=TRUE]]
  if (is.null(criterion_values)) {
    criterion_values <- x[[criterion_name]]
  }
  if (is.null(criterion_values)) {
    stop("The fitted object does not contain criterion values")
  }
  
  if(length(c(criterion_values))==1 ) {
    
    cat(paste0(criterion_name, ":"),criterion_values,"\n")
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
    
    cat(paste0(criterion_name, ":") , "\n")
    print(criterion_values)
    cat(" ", "\n")
    cat(paste0("minimum ", criterion_name, ":"),
        min(na.omit(as.vector(criterion_values))), "\n")
    
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
