library(pracma)
library(latex2exp)
library(dplyr)
### Kendall' tau
library(Kendall)
library(igraph)
library(parallel)
num_cores=7

library(glmnet)
#  Graph functions
min_vertex_set <- function(edges_matrix) {
  edges_matrix=as.matrix(as.matrix(edges_matrix ))
  #  edges_matrix=as.matrix(cleaned.data[,c(2,3)])
  V0 <- sample(edges_matrix, 1)
  # Initialize set S
  S <- V0
  # Initialize edge_index vector with NA values
  edge_index <- c()
  # Loop until no more edges can be found
  for(j in 1:(length(unique(edges_matrix))-1)) {
    iinndd=which( (edges_matrix[,1] %in% S & !( edges_matrix[,2] %in% S))|(edges_matrix[,2] %in% S & !( edges_matrix[,1] %in% S))  ) 
    if(length(iinndd) == 0) {
      # If no such neighbors exist, exit loop
      break
    }
    edge_index[j]=iinndd[sample(length(iinndd),1)]
    V1=setdiff(edges_matrix[edge_index[j],],S)
    # Add this vertex to S
    S <- c(S, V1)
  }
  # Return set S and edge_index vector
  return(edge_index)
}

I_pi_truncate <- function(design_X, p0, n_d, epsilon = 1e-11) {
  # Adjust p0 values based on epsilon
  p0[p0 < epsilon|p0 >1- epsilon] <- epsilon
  # Compute the weights
  w <- n_d * p0 * (1 - p0)
  # Scale the rows of design_X by the weights
  X_scaled <- design_X * w
  # Compute the matrix product
  result <- t(design_X) %*% X_scaled
  return(result)
}

find_closest_indices <- function(theta,dataset) {
  n <- length(theta)
  if (n < 2) {
    stop("Vector must contain at least two values.")
  }
  dis=abs(theta[unlist(dataset[,1])]-theta[unlist(dataset[,2])])
  dis[(dataset[ ,3]+dataset[ ,4])==0]=max(dis)
  ii=which(dis == min(dis, na.rm = TRUE))
  if(length(ii)>1){ii=sample(ii,1)}
  #  ii=which.min(dis)
  return(ii)
}


########################################################


### 100sushi
#We asked the respondent to sort the objects in the item set "B" 
#according to the respondent's preference in an ascending order.
{
  raw.data100=read.table("sushi3b.5000.10.order.txt", header = F)
  raw.data_100b=raw.data100[,3:12]
  names(raw.data_100b)=names(raw.data100)[1:10]
  first_elements  = c()
  second_elements = c()
  for(i in 1:9){
    first_elements  = c(first_elements , as.vector(as.matrix(raw.data_100b[, 1:i])       ) )
    second_elements = c(second_elements, as.vector(as.matrix(raw.data_100b[, (11-i):10] )) )
  }
  result_df100b <- data.frame(Object1 = first_elements, Object2 = second_elements)
  result_df100b=result_df100b+1
  result_df100b$obj1_won = result_df100b$Object1>result_df100b$Object2
  result_df100b[, c(1, 2)] <- t(apply(result_df100b[, c(1, 2)], 1, function(x) sort(as.numeric(x))) )
  
  # Aggregating the data
  aggregated_data_100b <- result_df100b %>%
    group_by(Object1, Object2) %>%
    summarise(
      count_A_won_TRUE = sum(obj1_won, na.rm = TRUE),
      count_A_won_FALSE = sum(!obj1_won, na.rm = TRUE)
    )
}

### estimate full model as true model
{
  p_plus=100
  n_full=nrow(result_df100b)
  design_matrix <-  matrix(0, nrow = n_full, ncol = p_plus)
  for (row in 1:n_full ) {
    design_matrix[row, as.numeric(result_df100b[ row,2]) ]  <-  -1
    design_matrix[row, as.numeric(result_df100b[ row,1]) ] <- 1
  }
  design_X=design_matrix[,-1]
  
  ttttt_1=Sys.time()
  #  names(data_full)[1:99] = c(names(data_full)[2:99], "V100")
  model_full_100_net <- glmnet(design_X,result_df100b$obj1_won ,  
                               family = "binomial", alpha = 0, 
                               standardize=FALSE,intercept=FALSE,lambda=0)
  #                               lower.limits=-5,upper.limits=5,
  theta_full_100=c(0,as.numeric(model_full_100_net$beta))
  
  ttttt_2=Sys.time()
  print(  ttttt_2-ttttt_1)
  
  
}

###Initial sample 0
{ 
  p_plus=100
  p=99
  data_old=result_df100b
  # full design matrix X
  nn=nrow(data_old)
  design_matrix <-  matrix(0, nrow = nn, ncol = p_plus)
  for (row in 1:nn ) {
    design_matrix[row, data_old[row,2]]  <-  -1
    design_matrix[row, data_old[row,1]] <- 1
  }
  design_X_full=design_matrix[,-1]
  design_matrix_full=design_matrix
  data_full=as.data.frame(design_X_full)
  data_full$Y=data_old[  ,3]
  ###  data_full <--> result_df100b
  
  #full design matrix X (compact version)
  n_summary=nrow(aggregated_data_100b)
  design_matrix <-  matrix(0, nrow = n_summary, ncol = p_plus)
  for (row in 1:n_summary ) {
    design_matrix[row, as.numeric(aggregated_data_100b[row,2])]  <-  -1
    design_matrix[row, as.numeric(aggregated_data_100b[row,1])]  <- 1
  }
  design_matrix_full_compact=design_matrix
  design_X=design_matrix[,-1]
}

###Initial sample 1
{
  epsilon = 1e-11
  set.seed(365)
  n_max  = 4000-1
  n_init = 100-1
  l_max=100
  Kendall_rank_Gi1_full=matrix(0, nrow = l_max, ncol = n_max-n_init)
  Kendall_rank_Gi0_full=matrix(0, nrow = l_max, ncol = n_max-n_init)
  Kendall_rank_unif_full_project=matrix(0, nrow = l_max, ncol = n_max-n_init)
  Kendall_rank_UC_full=matrix(0, nrow = l_max, ncol = n_max-n_init)
  
  Kendall_rank_Gi0=rep(0,n_max-n_init)
  Kendall_rank_Gi1=rep(0,n_max-n_init)
  Kendall_rank_unif_project=rep(0,n_max-n_init)  
  Kendall_rank_UC=rep(0,n_max-n_init)
  #  cal_time=c(0,0)
}


kk=1
tt_1=Sys.time() 
# initial data inside loop
for(ii in 1:l_max)
{  
  {  
    init_index<-min_vertex_set(as.matrix(data_old[,c(1,2)]))
    init_index=c(init_index,sample(setdiff(1:nn,init_index), n_init-length(init_index)) )
    init_index=sort(init_index)
    init_data=data_full[init_index,]
    aggregated_data_adp_init=aggregated_data_100b
    k = nrow( design_X  )
    n_d_init=rep(0,k)
    for (i in init_index){
      row=which(aggregated_data_100b$Object1==data_old$Object1[i]&aggregated_data_100b$Object2==data_old$Object2[i])
      #    row=sample(row,1)
      n_d_init[row]=n_d_init[row]+1
      if(data_full$Y[i]==1){
        aggregated_data_adp_init$count_A_won_TRUE[row]=aggregated_data_adp_init$count_A_won_TRUE[row]-1
      }
      else{
        aggregated_data_adp_init$count_A_won_FALSE[row]=aggregated_data_adp_init$count_A_won_FALSE[row]-1
      }
    }
  }
  
  ### GI1
  {
    method.order = 1
    {
      aggregated_data_adp=aggregated_data_adp_init
      adp_index=init_index
      k = nrow( design_X )
      k0 = length(init_index)
      n = n_max
      dt_adp=init_data 
      n_d=n_d_init
      { 
        theta_loop=rep(0,p)
        j=0
        for (i in (k0+1):n){
          zero_first_order_index=rep(0,k)
          model_loop <- glmnet(design_X_full[adp_index,],data_full$Y[adp_index],  
                               family = "binomial", alpha = 0, 
                               standardize=FALSE,intercept=FALSE,
                               lower.limits=-3,upper.limits=3,lambda=0)
          theta_loop=as.numeric(model_loop$beta)
          #c(0,as.numeric(model_full_100_net$beta))
          Kendall_rank_Gi1[i-k0]= cor(theta_full_100, c(0,theta_loop), method = "kendall") 
          t0 <- -abs(as.numeric(design_X %*% theta_loop ) )
          #expit function
          p0 <- plogis(t0)
          return_val=I_pi_truncate(design_X, p0, n_d, epsilon = 0)
          #        if (method.order ==1)
          {
            inv_X=solve(return_val , t(design_X) ) 
            zero_first_order_index= -colSums(inv_X^2)*p0*(1-p0)
          }
          zero_first_order_index[(aggregated_data_adp[,3]+aggregated_data_adp[,4])==0]=max(zero_first_order_index)+1
          action=which.min( zero_first_order_index )
          n_d[action]=n_d[action] +1
          {
            row_col_ind =as.numeric( aggregated_data_adp[action, c(1,2)] )
            action_set=setdiff(which(data_old$Object1==row_col_ind[1]&data_old$Object2==row_col_ind[2]),init_index)
            if (length(action_set)==1){ind_new = action_set}else{ind_new=sample(action_set ,1)}
            adp_index=c(adp_index, ind_new)
            if(data_full$Y[ind_new]==1){aggregated_data_adp$count_A_won_TRUE[action]=aggregated_data_adp$count_A_won_TRUE[action]-1}
            else{aggregated_data_adp$count_A_won_FALSE[action]=aggregated_data_adp$count_A_won_FALSE[action]-1}
          }
          # Generate the response variable Y based on the probabilities
          dt_adp = data_full[adp_index,]
        }
      }
    }
    Kendall_rank_Gi1_full[kk,]=Kendall_rank_Gi1
  }
  # Existing txt file (create it if it doesn't exist)
  gi1_file <- "kendall_full_GI1.txt"
  # Check if the file exists; if not, create it with headers
  if (!file.exists(gi1_file)) {
    write.table(Kendall_rank_Gi1, file = gi1_file, sep = ",", row.names = FALSE)
  } else {
    # Append data to the existing file without writing headers again
    write.table(Kendall_rank_Gi1, file = gi1_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }
  #  write.csv(Kendall_rank_Gi1_full,"kendall_full_GI1.csv")

  ### GI0 par
  {
    method.order = 0
    {
      aggregated_data_adp=aggregated_data_adp_init
      adp_index=init_index
      k = nrow( design_X )
      k0 = length(init_index)
      n = n_max
      dt_adp=init_data 
      n_d=n_d_init
      cls <- makeCluster(num_cores)
      { 
        theta_loop=rep(0,p)
        j=0
        for (i in (k0+1):n){
          zero_first_order_index=rep(0,k)
          model_loop <- glmnet(design_X_full[adp_index,],data_full$Y[adp_index],  
                               family = "binomial", alpha = 0, 
                               standardize=FALSE,intercept=FALSE,
                               lower.limits=-3,upper.limits=3,lambda=0)
          theta_loop=as.numeric(model_loop$beta)
          Kendall_rank_Gi0[i-k0]= cor(theta_full_100, c(0,theta_loop) , method = "kendall") 
          # Obtain the vector of fitted probabilities
          t0 <- -abs(as.numeric(design_X %*% theta_loop ) )
          #expit function
          p0 <- plogis(t0)
          return_val= I_pi_truncate(design_X, p0, n_d, epsilon = 0) 
          #design_X_0=design_matrix_full_compact%*%VV
          {
            clusterExport(cls,c('return_val','p0','design_X'))
            zero_first_order_index=parSapply(cls, 1:k, function(a) {
              return(
                sum(diag(solve(return_val +design_X[a,] %*% t(design_X[a,]) * p0[a] * (1 - p0[a]))  ))
              )
            })
          }
          zero_first_order_index[(aggregated_data_adp[,3]+aggregated_data_adp[,4])==0]=max(zero_first_order_index)+1
          action=which.min( zero_first_order_index )
          n_d[action]=n_d[action] +1
          {
            row_col_ind =as.numeric( aggregated_data_adp[action, c(1,2)] )
            action_set=setdiff(which(data_old$Object1==row_col_ind[1]&data_old$Object2==row_col_ind[2]),init_index)
            if (length(action_set)==1){ind_new = action_set}else{ind_new=sample(action_set ,1)}
            adp_index=c(adp_index, ind_new)
            if(data_full$Y[ind_new]==1){
              aggregated_data_adp$count_A_won_TRUE[action]=aggregated_data_adp$count_A_won_TRUE[action]-1
            }else{aggregated_data_adp$count_A_won_FALSE[action]=aggregated_data_adp$count_A_won_FALSE[action]-1}
          }
          # Generate the response variable Y based on the probabilities
          dt_adp = data_full[adp_index,]
        }
      }
      stopCluster(cls)
    }
    Kendall_rank_Gi0_full[kk,]=Kendall_rank_Gi0
  }
  # Existing txt file (create it if it doesn't exist)
  gi0_file <- "kendall_full_GI0.txt"
  # Check if the file exists; if not, create it with headers
  if (!file.exists(gi0_file)) {
    write.table(Kendall_rank_Gi0, file = gi0_file, sep = ",", row.names = FALSE)
  } else {
    # Append data to the existing file without writing headers again
    write.table(Kendall_rank_Gi0, file = gi0_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }
  #  write.csv(Kendall_rank_Gi0_full,"kendall_full_GI0.csv")
  
  ### UC
  {
    {
      aggregated_data_adp=aggregated_data_adp_init
      adp_index=init_index
      k = nrow( design_X )
      k0 = length(init_index)
      n = n_max
      dt_adp=init_data 
      n_d=n_d_init
      { 
        for (i in (k0+1):n){
          model_loop <- glmnet(design_X_full[adp_index,],data_full$Y[adp_index],  
                               family = "binomial", alpha = 0, 
                               standardize=FALSE,intercept=FALSE,
                               lower.limits=-3,upper.limits=3,lambda=0)
          #model_loop <- glmnet(design_matrix_full[adp_index,],data_full$Y[adp_index],  
          #                     family = "binomial", alpha = 0, 
          #                     standardize=FALSE,intercept=FALSE,
          #                     lambda=0.1)
          theta_loop=as.numeric(model_loop$beta)
          Kendall_rank_UC[i-k0]= cor(theta_full_100, c(0,theta_loop) , method = "kendall") 
          action=find_closest_indices(theta=c(0,theta_loop), aggregated_data_adp)
          n_d[action]=n_d[action] +1
          {
            row_col_ind =as.numeric( aggregated_data_adp[action, c(1,2)] )
            action_set=setdiff(which(data_old$Object1==row_col_ind[1]&data_old$Object2==row_col_ind[2]),init_index)
            if (length(action_set)==1){ind_new = action_set}else{ind_new=sample(action_set ,1)}
            adp_index=c(adp_index, ind_new)
            if(data_full$Y[ind_new]==1){aggregated_data_adp$count_A_won_TRUE[action]=aggregated_data_adp$count_A_won_TRUE[action]-1}
            else{aggregated_data_adp$count_A_won_FALSE[action]=aggregated_data_adp$count_A_won_FALSE[action]-1}
          }
          # Generate the response variable Y based on the probabilities
          dt_adp = data_full[adp_index,]
        }
      }
    }
    Kendall_rank_UC_full[kk,]=Kendall_rank_UC
  }
  # Existing txt file (create it if it doesn't exist)
  US_file <- "kendall_full_US.txt"
  # Check if the file exists; if not, create it with headers
  if (!file.exists(US_file)) {
    write.table(Kendall_rank_UC, file = US_file, sep = ",", row.names = FALSE)
  } else {
    # Append data to the existing file without writing headers again
    write.table(Kendall_rank_UC, file = US_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }
  #  write.csv(Kendall_rank_UC_full,"Kendall_rank_UC_full.csv")

  ### uniform
  {
    ### uniform
    init_index_adp=c(init_index, sample(setdiff(1:nrow(data_full),init_index)) )
    for (i in (k0+1):n){
      model_loop_project <- glmnet(design_X_full[init_index_adp[1:i],],data_full$Y[init_index_adp[1:i]],  
                                   family = "binomial", alpha = 0, 
                                   standardize=FALSE,intercept=FALSE,
                                   lower.limits=-3,upper.limits=3,lambda=0)
      theta_loop_project=as.numeric(model_loop_project$beta)
      Kendall_rank_unif_project[i-k0]= cor(theta_full_100, c(0,theta_loop_project) , method = "kendall") 
      #      Kendall_rank_unif[i-k0]=as.numeric(Kendall(as.numeric( theta_full_100),c(0, theta_loop) )$tau)
    }
    #    Kendall_rank_unif_full[kk,]=Kendall_rank_unif
    Kendall_rank_unif_full_project[kk,]=Kendall_rank_unif_project
  }
  # Existing txt file (create it if it doesn't exist)
  unif_file <- "kendall_full_unif_project.txt"
  # Check if the file exists; if not, create it with headers
  if (!file.exists(unif_file)) {
    write.table(Kendall_rank_unif_project, file = unif_file, sep = ",", row.names = FALSE)
  } else {
    # Append data to the existing file without writing headers again
    write.table(Kendall_rank_unif_project, file = unif_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  }
  #  write.csv(Kendall_rank_unif_full_project,"kendall_full_unif_project.csv")
  kk=kk+1
}
tt_2=Sys.time() 
print(tt_2-tt_1)

