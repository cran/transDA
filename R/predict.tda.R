predict.tda <-
function(object, newdata, ...)
  
{
  
  if(!inherits(object, "tda")) 
    stop("object not of class 'tda'")
  
  if (missing(newdata)) {
    stop("newdata must be provided")
  }
  
  
  
  pred<-function(x,p,mu,sigma,prior,lambda){
    
    N=nrow(x)
    lambda1 <- lambda
    
    
    data_t<-vector("list", length = length(mu))
    
    
    for (g in 1:length(mu)) {
      data_k<-vector("list", length = nrow(mu[[g]]))
      for (k in 1:nrow(mu[[g]])) {
        
        
        data_p<-matrix(NA,nrow = nrow(x), ncol = ncol(x))
        
        for (j in 1:ncol(mu[[g]]) ) {
          
          if (lambda1[[g]][k,j] == 0) {
            data_p[,j] <- x[,j]
          } else {
            data_p[,j]<- (exp(lambda1[[g]][k,j]*x[,j] )-
                            1)/(lambda1[[g]][k,j])
          }
          
        }
        data_k[[k]]<-data_p
        
      }
      data_t[[g]]<-data_k
    }
    
    
    pred_matrix <- matrix(NA, nrow = N, ncol = length(mu))
    
    for (g in 1:length(mu)) {
      
      prob<-matrix(NA,nrow = N,ncol = nrow(mu[[g]]))
      
      for (k in 1:nrow(mu[[g]]) ) {
        
        prob[,k]<-as.vector(p[[g]][k]*dmvnorm(data_t[[g]][[k]],mean = mu[[g]][k,],
                                              sigma = sigma[[g]][[k]]))*
          as.vector(exp(t(as.matrix(lambda1[[g]][k,]))%*%t(x)))
      }
      
      pred_matrix[,g] <- rowSums(prob)*prior[g]
      
    }
    
    pred_matrix <- pred_matrix/rowSums(pred_matrix)
    pred_ID <- max.col(pred_matrix)
    
    est <- list(pred_ID,pred_matrix)
    
    return(est)
    
  }
  
  
  p_list <- object$sub_prior
  
  mu_list <- object$mu
  
  sigma_list <- object$sigma
  
  prior <- object$prior
  
  
  if (is.null(object$lambda)==TRUE){
    
    lambda <- list()
    mu <- mu_list
    start <- 1
    lambda0 <- rep(0,length(unlist(mu)))
    for (dims in 1:length(mu)) {
      size <- prod(dim(mu[[dims]]))
      end <- start + size - 1
      
      lambda[[length(lambda) + 1]] <- matrix(lambda0[start:end],
                                             nrow=dim(mu[[dims]])[[1]],
                                             ncol = dim(mu[[dims]])[[2]],
                                             byrow = T)
      start <- end + 1
    }
    
    
  }else{
    
    lambda <- object$lambda
  }
  
  
  pred_ID<-pred(x= newdata,p = p_list,
                mu= mu_list,sigma =  sigma_list, 
                lambda = lambda,
                prior=prior)
  
  ### it should the original gourp ID 
  classNames <-  levels(as.factor(object$true_ID))
  
  
  Z <- pred_ID[[2]]
  cl <- apply(Z, 1, which.max)
  class <- factor(classNames[cl], levels = classNames)
  colnames(Z ) <- classNames
  
  out <- list(classification = class, Z = Z )
  return(out) 
  
  
  
}
