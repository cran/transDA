





tda<- function(x, max_k, ID, trans = TRUE, common_lambda= FALSE, 
               common_sigma = FALSE, iter = 50, 
               subgroup = NULL, tol= 0.001, lambda0 = 0.015){ 
  
  
  # E_step0
  
  
  E_step0<-function(x,k){
    
    
    kmean<-list()
    for (g in 1:length(x)) {
      
      kmean[[g]]<-kmeans(x[[g]],k[g])
      
    }
    
    p_mix_list<-list()
    
    
    for (g in 1:length(x)) {
      
      p_mix<-rep(NA,length(kmean[[g]]$size))
      
      for (j in 1: length(kmean[[g]]$size)) {
        
        p_mix[j]<-kmean[[g]]$size[j]/nrow(x[[g]])
      }
      
      p_mix_list[[g]]<-p_mix
    }
    
    
    mu_list<-list()
    for (j in 1:length(x)) {
      
      mu_list[[j]]<-kmean[[j]]$centers
      
    }
    
    gamma_list<-list()
    for (g in 1:length(x)) {
      
      gamma <- matrix(NA, nrow = nrow(x[[g]]), ncol = nrow(mu_list[[g]]))
      
      for (k in 1:nrow(mu_list[[g]])) {
        
        gamma[,k] <- p_mix_list[[g]][k] * dmvnorm(x[[g]],mean = mu_list[[g]][k,], sigma = cov(x[[g]]))
        
        
      }
      
      gamma <- gamma / rowSums(gamma)
      gamma_list[[g]]<-gamma
      
    }
    return(gamma_list)
  }
  
  
  ########################
  #                      #
  #  Estimate lambda     #
  #.                     #
  ########################
  
  fun_opt_lambda <- function(lambda,gamma,x) {
    
    lambda1 <- list()
    start <- 1
    
    for (dims in 1:length(x)) {
      size <- ncol(gamma[[dims]])*ncol(x[[dims]])
      end <- start + size - 1
      
      lambda1[[length(lambda1) + 1]] <- matrix(lambda[start:end],
                                               nrow=ncol(gamma[[dims]]),
                                               ncol = ncol(x[[dims]]),
                                               byrow = T)
      start <- end + 1
    }
    
    
    
    
    
    subsets_t<-vector("list", length = length(x))
    
    for (g in 1:length(x)) {
      subsets_k<-vector("list", length = ncol(gamma[[g]]) )
      for (k in 1:ncol(gamma[[g]])  ) {
        
        
        subsets_p<-matrix(NA,nrow = nrow(x[[g]]), ncol = ncol(x[[g]]))
        
        for (p in 1:ncol(x[[g]]) ) {
          
          
          if (lambda1[[g]][k,p] == 0) {
            subsets_p[,p] <- x[[g]][,p]
          } else {
            subsets_p[,p]<- (exp(lambda1[[g]][k,p]*x[[g]][,p] )-
                               1)/(lambda1[[g]][k,p])
          }
          
        }
        subsets_k[[k]]<-subsets_p
        
      }
      subsets_t[[g]]<-subsets_k
    }
    
    
    #  mu_t_list
    
    mu_t_list<-vector("list", length = length(x))
    
    for (g in 1:length(x)) {
      
      mu_k<-matrix(NA,nrow = ncol(gamma[[g]]), ncol = ncol(x[[g]]))
      
      for (k in 1:ncol(gamma[[g]])  ) {
        
        mu_k[k,]<-colSums(gamma[[g]][,k]*subsets_t[[g]][[k]])/sum(gamma[[g]][,k])
        
      }
      
      mu_t_list[[g]]<-mu_k
      
    }
    
    
    
    
    sigma_t_list<-vector("list", length = length(x))
    
    if( common_sigma == FALSE){
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                             t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))/sum(gamma[[g]][,k])
          
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
    } else{
      
      sigma_t_list0 <- vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                         t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))
          
        }
        
        sigma<-Reduce("+", sigma)
        
        sigma_t_list0[[g]] <- sigma
        
      }
      
      
      sigma_pooled <- Reduce("+", sigma_t_list0)/N
      
      
      
      sigma_t_list<-vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-sigma_pooled
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
      
      
      
    }
    
    
    
    loglik_g<-vector("list", length = length(mu_t_list))
    for (g in 1:length(subsets_t) ) {
      loglik_k<-vector("list", length = nrow(mu_t_list[[g]]))
      for (k in 1:length(subsets_t[[g]])) {
        
        
        loglik_k[[k]]<-gamma[[g]][,k]*(
          log(dmvnorm(subsets_t[[g]][[k]],
                      mean = mu_t_list[[g]][k,],
                      sigma =sigma_t_list[[g]][[k]] ))+
            as.vector(t(as.matrix(lambda1[[g]][k,]))%*%t(x[[g]]))
        )
        
      }
      loglik_g[[g]]<- unlist(loglik_k)
    }
    
    
    
    loglik_t <- -sum(unlist(loglik_g))
    
    return(loglik_t)
    
  }
  
  
  
  
  fun_opt_lambda1 <- function(lambda,gamma,x) {
    
    
    lambda1 = lambda
    
    
    subsets_t<-vector("list", length = length(x))
    
    for (g in 1:length(x)) {
      subsets_k<-vector("list", length = ncol(gamma[[g]]) )
      for (k in 1:ncol(gamma[[g]])  ) {
        
        
        subsets_p<-matrix(NA,nrow = nrow(x[[g]]), ncol = ncol(x[[g]]))
        
        for (p in 1:ncol(x[[g]]) ) {
          
          
          if (lambda1[p] == 0) {
            subsets_p[,p] <- x[[g]][,p]
          } else {
            subsets_p[,p]<- (exp(lambda1[p]*x[[g]][,p] )-
                               1)/(lambda1[p])
          }
          
        }
        subsets_k[[k]]<-subsets_p
        
      }
      subsets_t[[g]]<-subsets_k
    }
    
    
    #  mu_t_list
    
    mu_t_list<-vector("list", length = length(x))
    
    for (g in 1:length(x)) {
      
      mu_k<-matrix(NA,nrow = ncol(gamma[[g]]), ncol = ncol(x[[g]]))
      
      for (k in 1:ncol(gamma[[g]])  ) {
        
        mu_k[k,]<-colSums(gamma[[g]][,k]*subsets_t[[g]][[k]])/sum(gamma[[g]][,k])
        
      }
      
      mu_t_list[[g]]<-mu_k
      
    }
    
    
    
    
    sigma_t_list<-vector("list", length = length(x))
    
    if( common_sigma == FALSE){
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                             t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))/sum(gamma[[g]][,k])
          
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
    } else{
      
      sigma_t_list0 <- vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                         t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))
          
        }
        
        sigma<-Reduce("+", sigma)
        
        sigma_t_list0[[g]] <- sigma
        
      }
      
      
      sigma_pooled <- Reduce("+", sigma_t_list0)/N
      
      
      
      sigma_t_list<-vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-sigma_pooled
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
      
      
      
    }
    
    
    
    loglik_g<-vector("list", length = length(mu_t_list))
    for (g in 1:length(subsets_t) ) {
      loglik_k<-vector("list", length = nrow(mu_t_list[[g]]))
      for (k in 1:length(subsets_t[[g]])) {
        
        
        loglik_k[[k]]<-gamma[[g]][,k]*(
          log(dmvnorm(subsets_t[[g]][[k]],
                      mean = mu_t_list[[g]][k,],
                      sigma =sigma_t_list[[g]][[k]] ))+
            as.vector(t(as.matrix(lambda1))%*%t(x[[g]]))
        )
        
      }
      loglik_g[[g]]<- unlist(loglik_k)
    }
    
    
    
    loglik_t <- -sum(unlist(loglik_g))
    
    return(loglik_t)
    
  }
  
  
  
  ##########################
  #                        #
  #.      M_step           #
  #                        #
  ##########################
  
  
  M_step<-function(x,gamma,lambda){
    
    Mix_p<-list()
    
    for (i in 1:length(gamma)) {
      
      Mix_p[[i]]<-colSums(gamma[[i]])/nrow(gamma[[i]])
      
    }
    
    
    
    
    
    
    lambda1 <- list()
    start <- 1
    
    for (dims in 1:length(gamma)) {
      size <- ncol(gamma[[dims]])*ncol(x[[1]])
      end <- start + size - 1
      
      lambda1[[length(lambda1) + 1]] <- matrix(lambda[start:end],
                                               nrow=dim(gamma[[dims]])[[2]],
                                               ncol = dim(x[[dims]])[[2]],
                                               byrow = T)
      start <- end + 1
    }
    
    
    subsets_t<-vector("list", length = length(x))
    
    
    for (g in 1:length(x)) {
      subsets_k<-vector("list", length = ncol(gamma[[g]]))
      for (k in 1:ncol(gamma[[g]])) {
        
        
        subsets_p<-matrix(NA,nrow = nrow(x[[g]]), ncol = ncol(x[[g]]))
        
        for (j in 1:ncol(x[[g]]) ) {
          
          if (lambda1[[g]][k,j] == 0) {
            subsets_p[,j] <- x[[g]][,j]
          } else {
            subsets_p[,j]<- (exp(lambda1[[g]][k,j]*x[[g]][,j] )-
                               1)/(lambda1[[g]][k,j])
          }
          
        }
        subsets_k[[k]]<-subsets_p
        
      }
      subsets_t[[g]]<-subsets_k
    }
    
    
    
    mu_t_list<-vector("list", length = length(x))
    
    for (g in 1:length(x)) {
      
      mu_k<-matrix(NA,nrow = ncol(gamma[[g]]), ncol = ncol(x[[g]]))
      
      for (k in 1:ncol(gamma[[g]])) {
        
        mu_k[k,]<-colSums(gamma[[g]][,k]*subsets_t[[g]][[k]])/sum(gamma[[g]][,k])
        
      }
      
      mu_t_list[[g]]<-mu_k
      
    }
    
    
    
    
    sigma_t_list<-vector("list", length = length(x))
    
    if( common_sigma == FALSE){
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                             t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))/sum(gamma[[g]][,k])
          
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
    } else {
      
      sigma_t_list0 <- vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma[[k]]<-((t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,])%*%diag(gamma[[g]][,k])%*%
                         t(t(subsets_t[[g]][[k]])-mu_t_list[[g]][k,]))
          
        }
        
        sigma<-Reduce("+", sigma)
        
        sigma_t_list0[[g]] <- sigma
        
      }
      
      
      sigma_pooled <- Reduce("+", sigma_t_list0)/N
      
      
      
      sigma_t_list<-vector("list", length = length(x))
      
      for (g in 1:length(x)) {
        
        sigma_t_k<-vector("list", length = nrow(mu_t_list[[g]]))
        
        for (k in 1:nrow(mu_t_list[[g]])) {
          
          sigma_t_k[[k]]<-sigma_pooled
        }
        
        sigma_t_list[[g]]<-sigma_t_k
        
      }
      
      
      
    }
    
    
    
    return(list(p=Mix_p,
                subsets_t=subsets_t,
                mu=mu_t_list,
                sigma=sigma_t_list))
    
  }
  
  
  #######################
  #                     #
  #    E_step.          #
  #                     #
  #######################
  
  E_step<-function(p, mu,subsets_t,
                   x,sigma,lambda
                   
  ){
    
    
    lambda1 <- list()
    start <- 1
    
    for (dims in 1:length(mu)) {
      size <- prod(dim(mu[[dims]]))
      end <- start + size - 1
      
      lambda1[[length(lambda1) + 1]] <- matrix(lambda[start:end],
                                               nrow=dim(mu[[dims]])[[1]],
                                               ncol = dim(mu[[dims]])[[2]],
                                               byrow = T)
      start <- end + 1
    }
    
    
    gamma_list<-list()
    
    for (g in 1:length(x)) {
      
      gamma<-matrix(NA,nrow = nrow(subsets_t[[g]][[1]]),ncol = nrow( mu[[g]]))
      
      for (k in 1:nrow(mu[[g]]) ) {
        
        
        gamma[,k]<-as.vector(p[[g]][[k]]*dmvnorm(subsets_t[[g]][[k]],mean = mu[[g]][k,],
                                                 sigma = sigma[[g]][[k]]))*
          as.vector(exp(t(as.matrix(lambda1[[g]][k,]))%*%t(x[[g]])))
        
      }
      
      gamma_list[[g]]<-gamma/rowSums(gamma)
    }
    
    return(gamma=gamma_list)
    
  }
  
  
  ##########################
  #                        #
  #.    loglihood.         #
  #                        #
  ##########################
  
  
  
  log_likelihood <- function(gamma,
                             mu,
                             sigma,
                             x,
                             subsets_t,
                             lambda,
                             p,
                             N
                             
  ){
    
    lambda1 <- list()
    start <- 1
    
    for (dims in 1:length(mu)) {
      size <- prod(dim(mu[[dims]]))
      end <- start + size - 1
      
      lambda1[[length(lambda1) + 1]] <- matrix(lambda[start:end],
                                               nrow=dim(mu[[dims]])[[1]],
                                               ncol = dim(mu[[dims]])[[2]],
                                               byrow = T)
      start <- end + 1
    }
    
    
    
    
    prior<-rep(NA,length(x))
    
    for (i in 1:length(x)) {
      
      prior[i]<-nrow(x[[i]])/N
      
    }
    
    
    
    loglik_list<-vector("list",length = length(sigma))
    
    for (g in 1:length(x)) {
      loglik <- vector("list",length = length(sigma[[g]]))
      for (k in 1:length(sigma[[g]])){
        
        
        loglik[[k]]<-(p[[g]][k]*dmvnorm(subsets_t[[g]][[k]],mean=mu[[g]][k,],
                                        sigma=sigma[[g]][[k]]))*
          exp(as.vector(t(as.matrix(lambda1[[g]][k,]))%*%t(x[[g]])))
      }
      
      loglik_list[[g]]<-sum(log(Reduce("+",loglik)))+nrow(x[[g]])*log(prior[g])
    }
    
    Loglik<- sum(unlist(loglik_list))
    
    
    return(Loglik)
    
  }
  
  
  
  ### calculate BIC   #####
  
  BIC_est<-function(loglik,mu,x){
    
    
    N_k=nrow(do.call(rbind,mu))
    N_p=ncol(do.call(rbind,mu))
    
    if(common_sigma==FALSE && trans==TRUE && common_lambda==FALSE){
      M=N_k-1+2*N_k*N_p +N_k*N_p*(N_p+1)/2 }
    else if( common_sigma==TRUE && trans==TRUE && common_lambda==FALSE){
      
      M=N_k-1+2*N_k*N_p +N_p*(N_p+1)/2
      
    } else  if(common_sigma==FALSE && trans==TRUE && common_lambda==TRUE){
      
      M=N_k-1+(N_k+1)*N_p +N_k*N_p*(N_p+1)/2 }
    
    else if( common_sigma==TRUE && trans==TRUE && common_lambda==TRUE){
      
      M=N_k-1+(N_k+1)*N_p +N_p*(N_p+1)/2
    }
    
    
    else if( common_sigma==TRUE && trans==FALSE){
      
      M=N_k-1+N_k*N_p +N_p*(N_p+1)/2
      
    }else if( common_sigma==FALSE && trans==FALSE){
      
      M=N_k-1+N_k*N_p +N_k*N_p*(N_p+1)/2
    }
    
    BIC_est=-2*loglik+log(nrow(do.call(rbind, x)))*M
    
    return(BIC_est)
  }
  
  
  
  ######################################
  #                                    #
  #               predict function     #
  #                                    #
  ######################################
  
  
  pred<-function(x,p,mu,sigma,prior,lambda){
    
    N=nrow(x)
    lambda1 <- list()
    start <- 1
    
    for (dims in 1:length(mu)) {
      size <- nrow(mu[[dims]])*ncol(mu[[dims]])  # calculates xi * p
      end <- start + size - 1
      
      lambda1[[length(lambda1) + 1]] <- matrix(lambda[start:end],
                                               nrow=dim(mu[[dims]])[[1]],
                                               ncol = dim(mu[[dims]])[[2]],
                                               byrow = T)
      start <- end + 1
    }
    
    
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
    
    # pred_ID <- max.col(pred_matrix)
    
    return(pred_matrix)
    
    
  }
  
  
  # convert a character column to numeric ID
  
  # data = x column_name = ID  id_column_name = ID
  
  convert_to_numeric_id <- function(data,  id_column_name) {
    data= data.frame(data,id_column_name)
    #  data[[column_name]] <- as.factor(data[[column_name]])
    
    data$id_column_name <- as.factor(data$id_column_name)
    
    original_levels <- levels(data$id_column_name)
    
    
    #   convert_back <- function(numeric_ids) {
    #      sapply(numeric_ids, function(id) original_levels[id])}
    
    
    data$id_column_name <- as.numeric(data$id_column_name )
    #    if (id_column_name != column_name) {
    #   data[[column_name]] <- NULL  }
    
    #   return(list("data" = data, "convert_back" = convert_back))
    return(list("data" = data, classnames = original_levels))
  }
  
  
  
  
  result <- convert_to_numeric_id(x, ID)
  
  
  
  convert_back <- result$classnames
  
  
  x <- data.frame(x,ID)
  
  
  
  
  traindata <- x 
  subsets <- split(traindata, traindata$ID)
  
  
  
  N=nrow(traindata)
  
  
  ##############################
  #                            #
  #.  permutations function    #
  #                            #
  ##############################
  
  permutations<- function(n, r) {
    if (r == 0) {
      return(matrix(numeric(0), nrow = 1, ncol = 0))
    }
    
    perms <- matrix(nrow = n^r, ncol = r)
    for (i in 1:r) {
      
      times <- n^(r-i)
      perms[, i] <- rep(rep(1:n, each=times), n^(i-1))
    }
    
    return(perms)
  }
  
  ####################################################################
  
  subsets <- lapply(subsets, function(sub_df) {
    sub_df$ID <- NULL
    return(sub_df)
  })
  
  if( length(subgroup) != length(subsets) && missing(max_k)){
    
    stop( "The length of subgroup combination should equal to length of unique ID " )
    
  }
  
  is_positive_integer_vector <- function(vector) {
    is_positive <- all(vector > 0)
    is_whole_number <- all(vector == as.integer(vector))
    
    return(is_positive && is_whole_number)
  }
  
  
  
  if (is_positive_integer_vector(subgroup)!= TRUE ){
    
    stop( "subgroup should be all positive and integer " )
    
  }
  
  
  
  
  if (missing(max_k) && is.null(subgroup)) {
    stop("Both max_k and subgroup are missing. Please provide one.")
  } else if (!missing(max_k) && !is.null(subgroup)) {
    stop("Both max_k and subgroup  are provided. Please provide only one.")
  }
  
  if (!is.null(subgroup)) {
    subgroup_combination <- matrix(subgroup, nrow=1 , byrow=T)
  } else {
    # subgroup_combination<-gtools::permutations(n = max_k,r = length(subsets), repeats.allowed=TRUE)
    subgroup_combination <-  as.matrix(permutations(n = max_k,r = length(subsets)))
    
    subgroup_combination <- as.matrix(subgroup_combination[,ncol(subgroup_combination):1,drop = FALSE])
  }
  
  
  
  est<-list()
  
  #  lambda_0 <- lambda0
  # j= 512   lambda0 = 0.015
  for (j in 1:nrow(subgroup_combination)) {
    
    K<-subgroup_combination[j,]
    
    gamma_list<-E_step0(x=subsets,k=K)
    
    count=0
    loglik_old <- 0
    
    
    if(trans== TRUE && common_lambda==TRUE){
      
      lambda_int<-rep(lambda0, ncol(subsets[[1]]))
    }else{
      lambda_int<-rep(lambda0,sum(K)*ncol(subsets[[1]]))
    }
    
    
    
    params<-mu<-p<-sigma_list<-subsets_t<-list()
    loglik<-rep(NA,iter+1)
    
    
    repeat{
      
      if( trans== TRUE && common_lambda==FALSE){
        lambda00<-try(suppressWarnings(optim(par=lambda_int,fn=fun_opt_lambda,x=subsets,
                                             gamma=gamma_list,method = "Nelder-Mead"))$par,
                      silent = T)
        
        lambda_int<-lambda00
        
      } else if(trans== TRUE && common_lambda==TRUE){
        
        lambda000<-try(suppressWarnings(optim(par=lambda_int,fn=fun_opt_lambda1,x=subsets,
                                              gamma=gamma_list,method = "Nelder-Mead"))$par,
                       silent = T)
        lambda00 <- rep(lambda000,sum(K))
        lambda_int<-lambda000
      } else{
        
        lambda00 <- rep(0,sum(K)*ncol(subsets[[1]]))
      }
      
      # (!inherits(lambda000, "try-error"))
      # (!inherits(lambda00, "try-error"))
      # if (class(lambda0) != "try-error") {
      # if (!inherits(lambda00, "try-error")&&!inherits(lambda000, "try-error")) {
      # if (!inherits(lambda00, "try-error") || !inherits(lambda000, "try-error")) {
      if (!inherits(lambda_int, "try-error")) {
        count<-count+1
        ### lambda_int<-lambda0
        params<- M_step(x=subsets,gamma=gamma_list,lambda = lambda00)
        if( any(is.nan(unlist(params)))){
          break
        }
        
        
        mu<-params$mu
        p<-params$p
        sigma_list<-params$sigma
        subsets_t<-params$subsets_t
        
        gamma_list<-E_step(p=p, mu=mu, subsets_t=subsets_t,
                           x=subsets,sigma =sigma_list,lambda=lambda00)
        
        
        
        
        loglik[count]<- loglik_new <- log_likelihood (gamma=gamma_list,mu=mu,
                                                      sigma=sigma_list,
                                                      x=subsets,
                                                      subsets_t = subsets_t,
                                                      lambda = lambda00,
                                                      p=p,N=N )
        
        
        
      } else {
        count<-iter+1
      }
      
      if (count>iter || (abs(loglik_new - loglik_old) < tol) || any(is.nan(unlist(params))) ){
        loglik1 <- loglik[length(na.omit(loglik))-1]
        break
      }
      loglik_old <- loglik_new
      
    }
    
    BIC<-BIC_est(loglik = loglik1,mu=mu,x=subsets)
    
    est[[j]]<-list(params,gamma_list,loglik,BIC,lambda00)
    
  }
  
  if (any(is.nan(unlist(params))) && !is.null(subgroup)){
    stop("Sample size is too small, can't estimate, try to use smaller size of subgroup ")
  }
  
  
  for ( i in 1: length(est)){
    if (inherits(est[[i]][[5]], "try-error")  ) {
  
      est[[i]][[4]] = Inf
    }
    
  }
  
  for ( i in 1: length(est)){
 
    loglik_tt<- na.omit(est[[i]][[3]]) 
 
    if (length(loglik_tt) >=2 && length(loglik_tt) < iter && abs(loglik_tt[length(loglik_tt)] - loglik_tt[length(loglik_tt) - 1]) > tol) {
      est[[i]][[4]] <- Inf
    }else if(length(loglik_tt)==1){
      est[[i]][[4]] <- Inf
    }
    
  }
  
  
  
  

  est_bic<-rep(NA,length = length(est))
  
  for ( i in 1:length(est)){
    
    est_bic[i]<-est[[i]][[4]]
  }
  
  if(!is.null(subgroup)){
    
    est_bic <- est_bic
    
  }else if (length(subsets)==2 && is.null(subgroup) ) {
    est_bic<-matrix(est_bic,max_k,max_k, byrow = F)
    
    colnames(est_bic) <- paste0("G2_k =", 1:max_k)
    rownames(est_bic) <- paste0("G1_k =", 1:max_k)
    
  } else if ( length(subsets) > 2 && is.null(subgroup))
  {
    est_bic<-array(est_bic, dim=c(rep(max_k,length(subsets))) ) # names for array
    
    dim_names <- function(est_bic) {
      lapply(seq_along(dim(est_bic)), function(dim_index) {
        labels <- seq_len(dim(est_bic)[dim_index])
        paste0("G", dim_index, "_k=", labels)
      })
    }
    
    
    dim_names <- dim_names(est_bic)
    
    
    dimnames(est_bic) <- dim_names
    
  }
  
  
  
  est_bic[is.infinite(est_bic)] <- NA
  
  
  
  min_value <- Inf
  min_BIC_est <- NULL
  
  
  for (i in seq_along(est)) {
    current_BIC_value <- est[[i]][[4]]
    if (current_BIC_value < min_value) {
      min_value <- current_BIC_value
      min_BIC_est<- est[[i]]
    }
  }
  
  
  ## prediction function
  
  
  prior<-c(NA,rep(length(subsets)))
  
  for ( i in 1:length(subsets)) {
    
    prior[i]<-nrow(subsets[[i]])/N
    
  }
  
  col_name <- "ID"
  
  ID_col_number <- which(names(x) == col_name)
  
  pred_data <- traindata[,-ID_col_number]
  
  p_list <- min_BIC_est[[1]]$p
  
  
  if (is.null(p_list) &&inherits(lambda00, "try-error")) {
    stop("Can't estimate transformation parameter --lambda, try to change the innial lambda or increase the sample size")
  }else{ names(p_list)<-c(paste0("G", 1:length(subsets)))}
  
  mu_list <- min_BIC_est[[1]]$mu
  
  names(mu_list) <- c(paste0("G", 1:length(subsets)))
  
  sigma_list<-min_BIC_est[[1]]$sigma
  
  names(sigma_list) <- c(paste0("G", 1:length(subsets)))
  
  lambda00<-min_BIC_est[[5]]
  
  
  pred_matrix <- pred(x=pred_data,p = p_list,
                      mu=mu_list,sigma = sigma_list,
                      lambda=lambda00,prior=prior)
  
  z <- pred_matrix
  
  pred_ID <- max.col(pred_matrix)
  
  ID <- as.numeric(factor(traindata$ID))
  

  original_pred_ID <- factor(convert_back[as.vector(pred_ID)], levels = convert_back)
  original_ID <- factor(convert_back[as.vector(ID)], levels = convert_back)
  
  classNames <-  levels(as.factor(original_ID))
  colnames(z ) <- classNames
  
  
  ### alignment for groups
  
  calculate_accuracy <- function(cluster_assignments, true_labels) {
    
    
    true_labels <- as.factor(true_labels)
    
    
    
    unique_clusters <- as.character( unique(cluster_assignments))
    
    cluster_assignments <- as.character(cluster_assignments)
    
    cluster_to_label_map <- numeric(length(unique_clusters))
    
    
    
    for (cluster in unique_clusters)  {
      labels_in_cluster <- true_labels[cluster_assignments == cluster]
      most_common_label <- as.numeric(names(sort(table(labels_in_cluster), decreasing = TRUE)[1]))
      cluster_to_label_map[cluster] <- most_common_label
    }
    
    
    predicted_labels <- cluster_to_label_map[cluster_assignments]
    
 
    
    calculate_ARI <- function(true_labels, pred_labels) {
      
      true_labels <- factor(true_labels)
      pred_labels <- factor(pred_labels)
      
      contingency_table <- table(true_labels, pred_labels)
      
      sum_comb_nij <- sum(choose(contingency_table, 2))
      
      sum_comb_ai <- sum(choose(rowSums(contingency_table), 2))
      sum_comb_bj <- sum(choose(colSums(contingency_table), 2))
      
      
      N <- choose(sum(contingency_table), 2)
      
      expected_index <- sum_comb_ai * sum_comb_bj / N
      
      
      max_index <- 0.5 * (sum_comb_ai + sum_comb_bj)
      
      if (max_index == expected_index) {
        ARI <- 1
      } else {
        ARI <- (sum_comb_nij - expected_index) / (max_index - expected_index)
      }
      
      return(ARI)
    }
    

    ARI <- calculate_ARI(true_labels, predicted_labels)
    accuracy <- sum(predicted_labels == as.numeric(true_labels)) / length(true_labels)

    return(list(accuracy,ARI))
    
  }
  
  
  
  accuracy <- calculate_accuracy(cluster_assignments=pred_ID, true_labels=ID)
  
  

  
  ARI <-  accuracy[[2]]
  

  
  misclassification_rate <- 1 - accuracy[[1]]
  
  
  
  
  lambda <- list()
  mu <- mu_list
  start <- 1
  
  for (dims in 1:length(mu)) {
    size <- prod(dim(mu[[dims]]))
    end <- start + size - 1
    
    lambda[[length(lambda) + 1]] <- matrix(lambda00[start:end],
                                           nrow=dim(mu[[dims]])[[1]],
                                           ncol = dim(mu[[dims]])[[2]],
                                           byrow = T)
    start <- end + 1
  }
  
  
  
  names(lambda) <- c(paste0("G", 1:length(subsets)))
  
  loglik0 <- min_BIC_est[[3]]
  
  loglik  <- loglik0[!is.na(loglik0)]
  

  if (trans == FALSE)  {
    fit<-list(est_bic,p_list,mu_list,sigma_list,
              loglik,original_pred_ID,original_ID,prior,misclassification_rate,ARI,z)
    
    
    names(fit)<-c("BIC","sub_prior","mu","sigma",
                  "loglik","predict_ID", "true_ID", "prior",
                  "misclassification_rate","ARI","Z")
    
    
  } else if (trans == TRUE)    {
    
    fit<-list(est_bic,p_list,mu_list,sigma_list,lambda,
              loglik,original_pred_ID,original_ID,prior,misclassification_rate,ARI,z)
    
    
    names(fit)<-c("BIC","sub_prior","mu","sigma", "lambda",
                  "loglik","predict_ID","true_ID","prior",
                  "misclassification_rate","ARI","Z")
    
    
  }   
  
  
  class(fit) <- "tda"
  
  return(fit) # get the estimated parameter only with min BIC
  
  
  
}

#' Print Method for TMDA Objects
#'
#' This function is a method for the generic `print` function for objects of class `'TMDA'`.
#' It formats the printing of TMDA objects for better readability.
#'
#' @param x An object of class `tda`.
#' @param ... Further arguments passed to or from other methods.
#' @export

print.tda <- function(x,...) {
  
  cat("Discriminant Analysis", "\n")
  cat(" ", "\n")
  
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
    cat("Training Misclassification rate:", x$misclassification_rate, "\n")
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




