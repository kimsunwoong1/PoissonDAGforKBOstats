# Copyright (c) 2018 - 2019  Gunwoong Park  [gw.park23@gmail.com]
# All rights reserved.  See the file COPYING for license terms. 

##### ODS algorithm

ODS_Algorithm = function(data, N_min = 1, sparsity_level = 2){

  library(glmnet) 
  library(dplyr)
  library(plyr)
  
  data = as.matrix(data)
  p = ncol(data)
  n = nrow(data)
  
  #### Step 1: Lasso ####
  temp_CP = lapply( 1:ncol(data), candidate.parents.function, data = data, lambda.lasso = NULL, sparsity_level = sparsity_level )
  CandidateParents = matrix(0, p ,p)
  for(i in 1:p){ 
    if(is.numeric(temp_CP[[i]]) ){
      CandidateParents[i, temp_CP[[i]]] <- CandidateParents[temp_CP[[i]], i] <-1
    } 
  }
  
  #### Step 2: given CandidateParents
  pi_ODS = matrix(NA, ncol = p, nrow = 1)
  RemNode = 1:p
  condition_set = setdiff(1:p, RemNode)
  
  ## Overdispersion Scoring
  ODSscore1 = ODSscore_fun_S(data, condition_set = condition_set, RemNode = RemNode, minimum_N = 1)
  pi_ODS[1,1] = which.min( ODSscore1$scores )
  for( i in 2:(p-1)){
    RemNode = setdiff( 1:p,  pi_ODS[1,1:(i-1)] )
    ODSscore1 = sapply(RemNode, function(j) 
      ODSscore_fun_S(data,
                     condition_set = intersect(pi_ODS[1, 1:(i-1)], which(CandidateParents[,j] == 1) ), 
                     RemNode = j, 
                     minimum_N = N_min )[[1]] )
    pi_ODS[1,i] = RemNode[which.min(ODSscore1)] 
  }
  pi_ODS[1,p] = setdiff(1:p, pi_ODS[1,1:(p-1)])
  
  
  #### Step 3) ####
  ODS_DAG = matrix(0, p, p)
  ODS_lasso = glmnet( as.matrix(data[,pi_ODS[1,rep(1,2)]]), data[,pi_ODS[1,2]] , family = "poisson", alpha=1, lambda= log(p)/40, thresh = 1e-04 )$beta[,1]
  pos001 = which( ODS_lasso != 0 )
  if( !is.null(pos001) ){
    ODS_DAG[pi_ODS[1,1], pi_ODS[1,2] ] = 1 ;
  }
  for( i in 3:p){
    lambda_lasso = cv.glmnet( as.matrix(data[,pi_ODS[1,1:(i-1)]]), data[,pi_ODS[1,i]] , family = "poisson", alpha=1, nfolds = n, thresh = 1e-04, pmax = 500 )
    lambda_optimal = lambda_lasso$lambda.min + sparsity_level*( lambda_lasso$lambda.1se - lambda_lasso$lambda.min )
    ODS_lasso = glmnet( as.matrix(data[,pi_ODS[1,1:(i-1)]]), data[,pi_ODS[1,i]] , family = "poisson", alpha=1, thresh = 1e-04, lambda = lambda_optimal, pmax = p )$beta
     pos001 = which( ODS_lasso[,1] != 0 )
    if( length(pos001) > 0 ){
      ODS_DAG[ pi_ODS[1, pos001 ]  , pi_ODS[1,i] ] =1 ;
    }
  }
  ODS_DAG = t(ODS_DAG)
  
  return( list( DAG = ODS_DAG, 
                Ordering = pi_ODS) 
          )
}
################  
#### candidate parents function via lasso for Step 1) ####
candidate.parents.function = function(data, index, lambda.lasso = NULL, sparsity_level = 2.0 ) {
  p = ncol(data)
  n = nrow(data)
  class = setdiff(1:p, index)
  if( is.null(lambda.lasso) ){
    cv.lambda.lasso = cv.glmnet( as.matrix(data[,class]), data[,index] , family = "poisson", alpha=1, nfolds = n, type.measure="mse", thresh = 1e-05, pmax = 100)
    lambda.lasso = cv.lambda.lasso$lambda.min + sparsity_level*(cv.lambda.lasso$lambda.1se -cv.lambda.lasso$lambda.min )
    result_beta = cv.lambda.lasso$glmnet.fit$beta[,which.min( cv.lambda.lasso$glmnet.fit$lambda > lambda.lasso )]
  }else{
    result = glmnet( as.matrix( data[, class] ), data[, index] , family = "poisson" , alpha=1, lambda= lambda.lasso , thresh = 1e-04)
    result_beta = result$beta[,1]
  }
  return( class[which(result_beta!= 0 )] )
}

#########Score Functions: ODSscore_fun_S for Step 2) ######################
ODS_fun = function(x){
  ODS = var(x)-mean(x)
  return(ODS)
} 

ODSscore_fun_S = function(data, condition_set = NULL, RemNode, minimum_N = 2){
  #score_fun calculates the ratio of moments scores for given nodes (RemNode) provided parents set (conditions_set)
  #condition_set is a candidate parents set. 
  #RemNodes is a set of nodes to be considered as the the element of the ordering.   #moment.power is the power parameter for f_fun
  small_number = 10^(-20)
  data = data.frame(data)
  p = ncol(data)
  n = nrow(data)
  #  scores = scores_pos = rep(NA, p )

  if( length(condition_set)== 0 ){
    scores = sapply(RemNode, function(j) ODS_fun(data[,j]) )
  }else{
    
    #############
    id_result = id_fun_S(data, condition_set, minimum_N)
    if( is.null(id_result$id) ){ print("Warning!! Not enough sample size.")}
    if( !is.null(id_result$id) ){
      N = id_result$N
      res = matrix(NA, nrow = length(id_result$id), ncol = p)
      for( m in 1:id_result$M  ){
        if( sum( data[id_result$id[[m]],RemNode] >=  1 ) > minimum_N ){
          res[m,RemNode] = sapply(RemNode, function(j) ODS_fun( data[id_result$id[[m]],j] ) )
        }else{
          res[m,RemNode] = 0 #undefined score it depends on the power of the momoment. Hence it should be zero. 
        }
      }
      scores = sapply(RemNode, function(j) 
        sum( N[which(res[,j]!=0)]/sum(N[which(res[,j]!=0)]) * res[which(res[,j]!=0), j] ) )
    }
    #    colnames(scores[[i]]) =  colnames(scores_pos[[i]]) = 1:p
  }
  return( list( scores = scores) )
}



id_fun_S = function(data, condition_set, minimum_N){ 
  #id_fun finds indexes of samples such that conditioning set is the same.  
  #data is a sample_size by p matrix.
  #j is a condition node index.
  
  data = data.frame(data)
  p = ncol(data)
  n = nrow(data)
  rownames(data) = 1:nrow(data)
  #colnames(data) = 1:ncol(data)
  id = NULL
  set_S = plyr::count(data, vars= colnames(data)[condition_set]) #it depends on the order of loading libaries plyr must be later than dplyr.data must be dataframe instead of matrix.
  
  #print(max(set_S[,ncol(set_S)]))
  
  if( max(set_S[,ncol(set_S)]) > minimum_N ){
    set_S = set_S[set_S[,ncol(set_S)] > minimum_N, ]
    N_S = set_S[,ncol(set_S)]
    M = nrow(set_S)
    #set_S = as.data.frame( count(data, data[,j]) )
    sample_id = set_S
    if( length(condition_set) == 1){
      id = lapply(1:M, function(m) which( data[,condition_set] ==  set_S[m,-ncol(set_S)]) )
    }
    if( length(condition_set) > 1){
      id = lapply(1:M, function(m) as.numeric( rownames( match_df(data[,condition_set], set_S[m,-ncol(set_S)], on = NULL) ) ) )
    }
  }
  return( list(M = nrow(set_S), id = id, N = N_S) )
}


# match_df is from dplyr
match_df = function (x, y, on = NULL) 
{
  if (is.null(on)) {
    on <- intersect(names(x), names(y))
  }
  keys <- join.keys(x, y, on)
  x[keys$x %in% keys$y, , drop = FALSE]
}

