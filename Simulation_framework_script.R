# Final experimentation script
# This is designed to generate a csv file with the results from combinatorial  
# manipulation of the parameters of this experiment
# We will change method, n, p, structure, heritability and number of causal variants,
# and measure allele frequency of causal SNP and causal rank/p

################################Preliminaries###################################
#install.packages('tidyverse')
#install.packages('glmnet')

library(tidyverse)
library(glmnet)


#################################Parameters#####################################
heritabilities <- c(0.1,0.2,0.3,0.4,0.5) # to determine phenotype vec
dependency_weights <- c(seq(0,0.7,by = 0.05),seq(0.7,0.98,by = 0.02)) # to create structure
num_causal_list <- c(1,2,5,10) # how many causal variants
num_causal_large_p <- c(1,2,5,10,100) # how many causal for p = 1000  
reps <- 2 # how many reps to try of each combination of parameters


###############################(n,p)=(20,100)###################################
n <- 20
p <- 100

# Ewens allele frequency distribution
maf <- 1/(1:(n/2))
maf <- maf/sum(maf)

# Let k be the number of subpopulations of size m, where n = m*k
# We will use this to create X in a structured way
k <- 5
m <- n/k

# We will use this to generate the noise component of the Y
sigma <- 1
noise <- rnorm(n,0,sigma)

# Fix the heritability of the causal SNP. This will impact the generation of Y
for(h2 in heritabilities){
  
  # Generating populations with different levels of structure. 
  # d = 0 is no structure, d -> 1 is more structure
  for(d in dependency_weights){
    
    # we repeat the process reps many times, generating a new X and Y each time
    for(iter in 1:reps){
      
      
      # Generate the matrix X to represent genotype
      X <- matrix(0,nrow=n,ncol=p)
      for (snp in 1:p)
      {
        maf_count <- sample(1:(n/2),1,prob=maf)
        n1 <- maf_count # number of ones
        n0 <- n-maf_count # number of zeroes
        W <- matrix(0,nrow=m,ncol=k) # matrix that will become a column in X
        
        for (subp in 1:k)
        {
          next_subpopulation <- rep(0,m)
          
          #Assign the first 0/1 to the subpopulation based off the overall probability of selecting it
          firstval <- rbinom(1,1,(n1/(n0+n1)))
          W[1,subp] <- firstval
          next_subpopulation[1] <- firstval
          ifelse(firstval == 1, n1 <- n1-1, n0 <- n0-1)
          
          for (i in 2:m)
          {
            ones_so_far <- sum(next_subpopulation[1:i-1]) 
            total_so_far <- i-1
            
            #calculate probability of assigning 1 or 0 to the subpopulation depending on what is already in the population
            if(n1 == 0){
              cur_prob <- 0
            } else if(n0 == 0){
              cur_prob <- 1
            } else {
              cur_prob <- (n1/(n0+n1))*(1-d)+(ones_so_far/total_so_far)*d
            }
            
            nextval <- rbinom(1,1,cur_prob)
            W[i,subp] <- nextval
            next_subpopulation[i] <- nextval
            ifelse(nextval==1,n1<-n1-1,n0<-n0-1)
            
          }
        }
        
        # This subtle line of code limits the extent of linkage disequilibrium. Toggle the comment and observe what happens.
        #W <- W[,sample(1:k,k)]
        
        X_vector <- c(W)
        X[,snp] <- X_vector
      }
      
      # Measure the amount of structure as the average LD
      LDmat <- cor(X)^2
      average_LD <- sum(LDmat)/(p^2)
      
      #Calculate r2 to use later in Evalue
      mean_r2_perSNP <- (rowSums(LDmat)-1)/(p-1)
      
      
      #Vectors to store results - rank of first causal variant to be discovered
      # and maf of that variant
      GWAS_ranks <- rep(0,length(num_causal_list))
      GWAS_mafs <- rep(0,length(num_causal_list))
      Eval_ranks <- rep(0,length(num_causal_list))
      Eval_mafs <- rep(0,length(num_causal_list))
      LASSO_ranks <- rep(0,length(num_causal_list))
      LASSO_mafs <- rep(0,length(num_causal_list))
      
      
      # Now we have filled in the X matrix, we use it to generate Y
      # This depends on heritability and the allele frequency of the causal SNP(s)
      # So first, we choose how many casual SNPs we want
      # Once we have generated Y we apply our methods to find causal SNP(s)
      for(t in 1:length(num_causal_list)){
        
        num_causal <- num_causal_list[t]
        
        # We need to divide the hertiability between the causal variants
        h2 <- h2/num_causal
        
        signal <- 0
        for(i in 1:num_causal){
          # determine the coefficients of the causal vector(s) in the model by their maf
          f <- sum(X[,i])/n
          beta <- (((h2/(1-h2))*(1/(f*(1-f))))^.5) * sigma
          # signal = sum(beta_i*X_i) 
          signal <- signal + beta*X[,i]
        }
        
        h2 <- h2*num_causal
        
        # Generate the Y vector
        Y <- signal + noise
        
        # We will use this correlation in finding Evalue
        rXY_2 <- cor(X,Y)^2
        
        # The next step is to successively apply each method and see where they 
        # rank the causal SNPs, recording the lowest (best) rank
        
        ###########################GWAS and Evalue##############################
        # Compute the p-values and E-values
        gwaseval_tibble <- tibble(SNP_index = 1:p,GWAS_p = 0, E_val = 0)
        for (focalSNP in 1:p)
        {
          
          F_val <- (n-2)*rXY_2[focalSNP]/(1-rXY_2[focalSNP])
          gwaseval_tibble$GWAS_p[focalSNP] <- pf(F_val,1,n-2,0,lower.tail = FALSE)
          e_val <- 0
          for (otherSNP in 1:p)
          {
            if (LDmat[focalSNP,otherSNP] == 1)
            {
              e_val <- e_val + .5
            }
            else 
            {
              statistic <- F_val * (1-sqrt(LDmat[focalSNP,otherSNP]))^2/(1 - LDmat[focalSNP,otherSNP])
              e_val <- e_val + pf(statistic,1,n-2,0,lower.tail = FALSE)
            }
          }
          gwaseval_tibble$E_val[focalSNP] <- e_val/p
        }
        
        # Extract the columns of pvals and evals now they've been calculated
        GWAS_pvals <- gwaseval_tibble$GWAS_p
        E_vals <- gwaseval_tibble$E_val
        
        # Record the ranking of the causal SNP in GWAS
        sortedp <- sort(GWAS_pvals,decreasing=FALSE,index.return=TRUE)
        # Find the highest ranked causal variant and it's maf 
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        GWAS_ranks[t] <- causal_rank
        GWAS_mafs[t] <- causal_maf
        
        # Record the ranking of the causal SNP in Evalue
        sortedE <- sort(abs(E_vals),decreasing=FALSE,index.return=TRUE)
        # Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        Eval_ranks[t] <- causal_rank
        Eval_mafs[t] <- causal_maf
        
        
        ##################################LASSO#######################################
        #apply the glm and find the coefficients of the columns of the matrix
        fit <- cv.glmnet(x=as.matrix(X), y=Y, alpha=1)
        cvec<-as.matrix(coef(fit,s='lambda.1se'))
        num_skip <- length(cvec)-p
        new_vec <- cvec[(1+num_skip):(p+num_skip),]
        
        # Record the ranking of the causal SNP
        sortedp <- sort(abs(new_vec),decreasing=TRUE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        LASSO_ranks[t] <- causal_rank
        LASSO_mafs[t] <- causal_maf
        
        
      }
      
      
      #Store the results in a tibble! Including the measured structure for X
      if(iter == 1){
        fixed_dh2_tibble <- tibble(rep = rep(iter,length.out= length(num_causal_list)),
                                   num_causal = num_causal_list, GWAS_ranks = GWAS_ranks, 
                                   GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                                   Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                                   LASSO_mafs = LASSO_mafs, average_LD 
                                   = rep(average_LD,length.out = length(num_causal_list)))
      } else {
        next_block <- tibble(rep = rep(iter,length.out= length(num_causal_list)),
                             num_causal = num_causal_list, GWAS_ranks = GWAS_ranks, 
                             GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                             Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                             LASSO_mafs = LASSO_mafs, average_LD 
                             = rep(average_LD,length.out = length(num_causal_list)))
        fixed_dh2_tibble <- bind_rows(fixed_dh2_tibble,next_block)
      }
      
      
      
    }
    
    # We record and compile our results for each value of d
    # Note the important assumption that the first value in dependency_weights is 0
    if(d == 0){
      fixed_h2_tibble <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
    } else {
      next_block <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
      fixed_h2_tibble <- bind_rows(fixed_h2_tibble,next_block)
    }
    
    
  }
  
  # Finally, compile the results for heritability
  if(h2 == heritabilities[1]){
    final_results <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
  } else {
    next_block <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
    final_results <- bind_rows(final_results,next_block)
  }
  

}


# Now that we have produced our results we want to manipulate them so they are in our desired form.
results_n20_p100 <- final_results  %>%  
  pivot_longer(., GWAS_ranks:LASSO_mafs,names_to = c("method", ".value"),names_pattern = "([^_]+)_([^_]+)") %>% 
  arrange(., heritability,dependency_weight,num_causal,rep) %>% 
  mutate(., rank_over_p = ranks/p) %>% 
  select(., -c(ranks))

# And save them to a csv file!
file_name <- paste("Results_n",as.character(n),"_p",as.character(p),".csv",sep = "")

write.csv(results_n20_p100,file_name,row.names = FALSE)




###############################(n,p)=(20,1000)###################################
n <- 20
p <- 1000

# Ewens allele frequency distribution
maf <- 1/(1:(n/2))
maf <- maf/sum(maf)

# Let k be the number of subpopulations of size m, where n = m*k
# We will use this to create X in a structured way
k <- 5
m <- n/k

# We will use this to generate the noise component of the Y
sigma <- 1
noise <- rnorm(n,0,sigma)

# Fix the heritability of the causal SNP. This will impact the generation of Y
for(h2 in heritabilities){
  
  # Generating populations with different levels of structure. 
  # d = 0 is no structure, d -> 1 is more structure
  for(d in dependency_weights){
    
    # we repeat the process reps many times, generating a new X and Y each time
    for(iter in 1:reps){
      
      
      # Generate the matrix X to represent genotype
      X <- matrix(0,nrow=n,ncol=p)
      for (snp in 1:p)
      {
        maf_count <- sample(1:(n/2),1,prob=maf)
        n1 <- maf_count # number of ones
        n0 <- n-maf_count # number of zeroes
        W <- matrix(0,nrow=m,ncol=k) # matrix that will become a column in X
        
        for (subp in 1:k)
        {
          next_subpopulation <- rep(0,m)
          
          #Assign the first 0/1 to the subpopulation based off the overall probability of selecting it
          firstval <- rbinom(1,1,(n1/(n0+n1)))
          W[1,subp] <- firstval
          next_subpopulation[1] <- firstval
          ifelse(firstval == 1, n1 <- n1-1, n0 <- n0-1)
          
          for (i in 2:m)
          {
            ones_so_far <- sum(next_subpopulation[1:i-1]) 
            total_so_far <- i-1
            
            #calculate probability of assigning 1 or 0 to the subpopulation depending on what is already in the population
            if(n1 == 0){
              cur_prob <- 0
            } else if(n0 == 0){
              cur_prob <- 1
            } else {
              cur_prob <- (n1/(n0+n1))*(1-d)+(ones_so_far/total_so_far)*d
            }
            
            nextval <- rbinom(1,1,cur_prob)
            W[i,subp] <- nextval
            next_subpopulation[i] <- nextval
            ifelse(nextval==1,n1<-n1-1,n0<-n0-1)
            
          }
        }
        
        # This subtle line of code limits the extent of linkage disequilibrium. Toggle the comment and observe what happens.
        #W <- W[,sample(1:k,k)]
        
        X_vector <- c(W)
        X[,snp] <- X_vector
      }
      
      # Measure the amount of structure as the average LD
      average_LD <- sum(cor(X)^2)/(p^2)
      
      
      #Vectors to store results - rank of first causal variant to be discovered
      # and maf of that variant
      GWAS_ranks <- rep(0,length(num_causal_large_p))
      GWAS_mafs <- rep(0,length(num_causal_large_p))
      Eval_ranks <- rep(0,length(num_causal_large_p))
      Eval_mafs <- rep(0,length(num_causal_large_p))
      LASSO_ranks <- rep(0,length(num_causal_large_p))
      LASSO_mafs <- rep(0,length(num_causal_large_p))
      
      
      # Now we have filled in the X matrix, we use it to generate Y
      # This depends on heritability and the allele frequency of the causal SNP(s)
      # So first, we choose how many casual SNPs we want
      # Once we have generated Y we apply our methods to find causal SNP(s)
      for(t in 1:length(num_causal_large_p)){
        
        num_causal <- num_causal_large_p[t]
        
        # We need to divide the hertiability between the causal variants
        h2 <- h2/num_causal
        
        signal <- 0
        for(i in 1:num_causal){
          # determine the coefficients of the causal vector(s) in the model by their maf
          f <- sum(X[,i])/n
          beta <- (((h2/(1-h2))*(1/(f*(1-f))))^.5) * sigma
          # signal = sum(beta_i*X_i) 
          signal <- signal + beta*X[,i]
        }
        
        h2 <- h2*num_causal
        
        # Generate the Y vector
        Y <- signal + noise
        
        
        # The next step is to successively apply each method and see where they 
        # rank the causal SNPs, recording the lowest (best) rank
        
        ###################################GWAS###################################
        # Compute the p-values
        GWAS_pvals <- rep(0,p)
        for (j in 1:p)
        {
          GWAS_model <- lm(Y~X[,j])
          GWAS_pvals[j] <- summary(GWAS_model)$coefficients[2,4]
        }
        
        # Record the ranking of the causal SNP
        sortedp <- sort(GWAS_pvals,decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and it's maf 
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        GWAS_ranks[t] <- causal_rank
        GWAS_mafs[t] <- causal_maf
        
        
        
        
        ###################################Evalue#####################################
        #Find the list of E values
        E_vals <- rep(0,p)
        for (j in 1:p)
        {
          cond_pvals <- rep(0,p)
          for (s in 1:p)
          {
            if((sum(X[,j]!=X[,s]) == 0)|(sum(X[,j]==X[,s]) == 0))
            {
              cond_pvals[s] <- NA
            }
            else
            {
              two_variable_model <- lm(Y~X[,j]+X[,s])
              cond_pvals[s] <- summary(two_variable_model)$coefficients[2,4]
            }
            E_vals[j] <- (sum(cond_pvals[!is.na(cond_pvals)]) + (1/2)*sum(is.na(cond_pvals)))/p
          }
        }
        
        # Record the ranking of the causal SNP
        sortedE <- sort(abs(E_vals),decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        Eval_ranks[t] <- causal_rank
        Eval_mafs[t] <- causal_maf
        
        
        ##################################LASSO#######################################
        #apply the glm and find the coefficients of the columns of the matrix
        fit <- cv.glmnet(x=as.matrix(X), y=Y, alpha=1)
        cvec<-as.matrix(coef(fit,s='lambda.1se'))
        num_skip <- length(cvec)-p
        new_vec <- cvec[(1+num_skip):(p+num_skip),]
        
        # Record the ranking of the causal SNP
        sortedp <- sort(abs(new_vec),decreasing=TRUE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        LASSO_ranks[t] <- causal_rank
        LASSO_mafs[t] <- causal_maf
        
        
      }
      
      
      #Store the results in a tibble! Including the measured structure for X
      if(iter == 1){
        fixed_dh2_tibble <- tibble(rep = rep(iter,length.out= length(num_causal_large_p)),
                                   num_causal = num_causal_large_p, GWAS_ranks = GWAS_ranks, 
                                   GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                                   Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                                   LASSO_mafs = LASSO_mafs, average_LD 
                                   = rep(average_LD,length.out = length(num_causal_large_p)))
      } else {
        next_block <- tibble(rep = rep(iter,length.out= length(num_causal_large_p)),
                             num_causal = num_causal_large_p, GWAS_ranks = GWAS_ranks, 
                             GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                             Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                             LASSO_mafs = LASSO_mafs, average_LD 
                             = rep(average_LD,length.out = length(num_causal_large_p)))
        fixed_dh2_tibble <- bind_rows(fixed_dh2_tibble,next_block)
      }
      
      
      
    }
    
    # We record and compile our results for each value of d
    # Note the important assumption that the first value in dependency_weights is 0
    if(d == 0){
      fixed_h2_tibble <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
    } else {
      next_block <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
      fixed_h2_tibble <- bind_rows(fixed_h2_tibble,next_block)
    }
    
    
  }
  
  # Finally, compile the results for heritability
  if(h2 == heritabilities[1]){
    final_results <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
  } else {
    next_block <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
    final_results <- bind_rows(final_results,next_block)
  }
  
  
}


# Now that we have produced our results we want to manipulate them so they are in our desired form.
results_n20_p1000 <- final_results  %>%  
  pivot_longer(., GWAS_ranks:LASSO_mafs,names_to = c("method", ".value"),names_pattern = "([^_]+)_([^_]+)") %>% 
  arrange(., heritability,dependency_weight,num_causal,rep) %>% 
  mutate(., rank_over_p = ranks/p) %>% 
  select(., -c(ranks))

# And save them to a csv file!
file_name <- paste("Results_n",as.character(n),"_p",as.character(p),".csv",sep = "")

write.csv(results_n20_p1000,file_name,row.names = FALSE)



###############################(n,p)=(100,100)###################################
n <- 100
p <- 100

# Ewens allele frequency distribution
maf <- 1/(1:(n/2))
maf <- maf/sum(maf)

# Let k be the number of subpopulations of size m, where n = m*k
# We will use this to create X in a structured way
k <- 10
m <- n/k

# We will use this to generate the noise component of the Y
sigma <- 1
noise <- rnorm(n,0,sigma)

# Fix the heritability of the causal SNP. This will impact the generation of Y
for(h2 in heritabilities){
  
  # Generating populations with different levels of structure. 
  # d = 0 is no structure, d -> 1 is more structure
  for(d in dependency_weights){
    
    # we repeat the process reps many times, generating a new X and Y each time
    for(iter in 1:reps){
      
      
      # Generate the matrix X to represent genotype
      X <- matrix(0,nrow=n,ncol=p)
      for (snp in 1:p)
      {
        maf_count <- sample(1:(n/2),1,prob=maf)
        n1 <- maf_count # number of ones
        n0 <- n-maf_count # number of zeroes
        W <- matrix(0,nrow=m,ncol=k) # matrix that will become a column in X
        
        for (subp in 1:k)
        {
          next_subpopulation <- rep(0,m)
          
          #Assign the first 0/1 to the subpopulation based off the overall probability of selecting it
          firstval <- rbinom(1,1,(n1/(n0+n1)))
          W[1,subp] <- firstval
          next_subpopulation[1] <- firstval
          ifelse(firstval == 1, n1 <- n1-1, n0 <- n0-1)
          
          for (i in 2:m)
          {
            ones_so_far <- sum(next_subpopulation[1:i-1]) 
            total_so_far <- i-1
            
            #calculate probability of assigning 1 or 0 to the subpopulation depending on what is already in the population
            if(n1 == 0){
              cur_prob <- 0
            } else if(n0 == 0){
              cur_prob <- 1
            } else {
              cur_prob <- (n1/(n0+n1))*(1-d)+(ones_so_far/total_so_far)*d
            }
            
            nextval <- rbinom(1,1,cur_prob)
            W[i,subp] <- nextval
            next_subpopulation[i] <- nextval
            ifelse(nextval==1,n1<-n1-1,n0<-n0-1)
            
          }
        }
        
        # This subtle line of code limits the extent of linkage disequilibrium. Toggle the comment and observe what happens.
        #W <- W[,sample(1:k,k)]
        
        X_vector <- c(W)
        X[,snp] <- X_vector
      }
      
      # Measure the amount of structure as the average LD
      average_LD <- sum(cor(X)^2)/(p^2)
      
      
      #Vectors to store results - rank of first causal variant to be discovered
      # and maf of that variant
      GWAS_ranks <- rep(0,length(num_causal_list))
      GWAS_mafs <- rep(0,length(num_causal_list))
      Eval_ranks <- rep(0,length(num_causal_list))
      Eval_mafs <- rep(0,length(num_causal_list))
      LASSO_ranks <- rep(0,length(num_causal_list))
      LASSO_mafs <- rep(0,length(num_causal_list))
      
      
      # Now we have filled in the X matrix, we use it to generate Y
      # This depends on heritability and the allele frequency of the causal SNP(s)
      # So first, we choose how many casual SNPs we want
      # Once we have generated Y we apply our methods to find causal SNP(s)
      for(t in 1:length(num_causal_list)){
        
        num_causal <- num_causal_list[t]
        
        # We need to divide the hertiability between the causal variants
        h2 <- h2/num_causal
        
        signal <- 0
        for(i in 1:num_causal){
          # determine the coefficients of the causal vector(s) in the model by their maf
          f <- sum(X[,i])/n
          beta <- (((h2/(1-h2))*(1/(f*(1-f))))^.5) * sigma
          # signal = sum(beta_i*X_i) 
          signal <- signal + beta*X[,i]
        }
        
        h2 <- h2*num_causal
        
        # Generate the Y vector
        Y <- signal + noise
        
        
        # The next step is to successively apply each method and see where they 
        # rank the causal SNPs, recording the lowest (best) rank
        
        ###################################GWAS###################################
        # Compute the p-values
        GWAS_pvals <- rep(0,p)
        for (j in 1:p)
        {
          GWAS_model <- lm(Y~X[,j])
          GWAS_pvals[j] <- summary(GWAS_model)$coefficients[2,4]
        }
        
        # Record the ranking of the causal SNP
        sortedp <- sort(GWAS_pvals,decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and it's maf 
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        GWAS_ranks[t] <- causal_rank
        GWAS_mafs[t] <- causal_maf
        
        
        
        
        ###################################Evalue#####################################
        #Find the list of E values
        E_vals <- rep(0,p)
        for (j in 1:p)
        {
          cond_pvals <- rep(0,p)
          for (s in 1:p)
          {
            if((sum(X[,j]!=X[,s]) == 0)|(sum(X[,j]==X[,s]) == 0))
            {
              cond_pvals[s] <- NA
            }
            else
            {
              two_variable_model <- lm(Y~X[,j]+X[,s])
              cond_pvals[s] <- summary(two_variable_model)$coefficients[2,4]
            }
            E_vals[j] <- (sum(cond_pvals[!is.na(cond_pvals)]) + (1/2)*sum(is.na(cond_pvals)))/p
          }
        }
        
        # Record the ranking of the causal SNP
        sortedE <- sort(abs(E_vals),decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        Eval_ranks[t] <- causal_rank
        Eval_mafs[t] <- causal_maf
        
        
        ##################################LASSO#######################################
        #apply the glm and find the coefficients of the columns of the matrix
        fit <- cv.glmnet(x=as.matrix(X), y=Y, alpha=1)
        cvec<-as.matrix(coef(fit,s='lambda.1se'))
        num_skip <- length(cvec)-p
        new_vec <- cvec[(1+num_skip):(p+num_skip),]
        
        # Record the ranking of the causal SNP
        sortedp <- sort(abs(new_vec),decreasing=TRUE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        LASSO_ranks[t] <- causal_rank
        LASSO_mafs[t] <- causal_maf
        
        
      }
      
      
      #Store the results in a tibble! Including the measured structure for X
      if(iter == 1){
        fixed_dh2_tibble <- tibble(rep = rep(iter,length.out= length(num_causal_list)),
                                   num_causal = num_causal_list, GWAS_ranks = GWAS_ranks, 
                                   GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                                   Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                                   LASSO_mafs = LASSO_mafs, average_LD 
                                   = rep(average_LD,length.out = length(num_causal_list)))
      } else {
        next_block <- tibble(rep = rep(iter,length.out= length(num_causal_list)),
                             num_causal = num_causal_list, GWAS_ranks = GWAS_ranks, 
                             GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                             Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                             LASSO_mafs = LASSO_mafs, average_LD 
                             = rep(average_LD,length.out = length(num_causal_list)))
        fixed_dh2_tibble <- bind_rows(fixed_dh2_tibble,next_block)
      }
      
      
      
    }
    
    # We record and compile our results for each value of d
    # Note the important assumption that the first value in dependency_weights is 0
    if(d == 0){
      fixed_h2_tibble <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
    } else {
      next_block <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
      fixed_h2_tibble <- bind_rows(fixed_h2_tibble,next_block)
    }
    
    
  }
  
  # Finally, compile the results for heritability
  if(h2 == heritabilities[1]){
    final_results <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
  } else {
    next_block <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
    final_results <- bind_rows(final_results,next_block)
  }
  
  
}


# Now that we have produced our results we want to manipulate them so they are in our desired form.
results_n100_p100 <- final_results  %>%  
  pivot_longer(., GWAS_ranks:LASSO_mafs,names_to = c("method", ".value"),names_pattern = "([^_]+)_([^_]+)") %>% 
  arrange(., heritability,dependency_weight,num_causal,rep) %>% 
  mutate(., rank_over_p = ranks/p) %>% 
  select(., -c(ranks))

# And save them to a csv file!
file_name <- paste("Results_n",as.character(n),"_p",as.character(p),".csv",sep = "")

write.csv(results_n100_p100,file_name,row.names = FALSE)




###############################(n,p)=(100,1000)###################################
n <- 100
p <- 1000

# Ewens allele frequency distribution
maf <- 1/(1:(n/2))
maf <- maf/sum(maf)

# Let k be the number of subpopulations of size m, where n = m*k
# We will use this to create X in a structured way
k <- 10
m <- n/k

# We will use this to generate the noise component of the Y
sigma <- 1
noise <- rnorm(n,0,sigma)

# Fix the heritability of the causal SNP. This will impact the generation of Y
for(h2 in heritabilities){
  
  # Generating populations with different levels of structure. 
  # d = 0 is no structure, d -> 1 is more structure
  for(d in dependency_weights){
    
    # we repeat the process reps many times, generating a new X and Y each time
    for(iter in 1:reps){
      
      
      # Generate the matrix X to represent genotype
      X <- matrix(0,nrow=n,ncol=p)
      for (snp in 1:p)
      {
        maf_count <- sample(1:(n/2),1,prob=maf)
        n1 <- maf_count # number of ones
        n0 <- n-maf_count # number of zeroes
        W <- matrix(0,nrow=m,ncol=k) # matrix that will become a column in X
        
        for (subp in 1:k)
        {
          next_subpopulation <- rep(0,m)
          
          #Assign the first 0/1 to the subpopulation based off the overall probability of selecting it
          firstval <- rbinom(1,1,(n1/(n0+n1)))
          W[1,subp] <- firstval
          next_subpopulation[1] <- firstval
          ifelse(firstval == 1, n1 <- n1-1, n0 <- n0-1)
          
          for (i in 2:m)
          {
            ones_so_far <- sum(next_subpopulation[1:i-1]) 
            total_so_far <- i-1
            
            #calculate probability of assigning 1 or 0 to the subpopulation depending on what is already in the population
            if(n1 == 0){
              cur_prob <- 0
            } else if(n0 == 0){
              cur_prob <- 1
            } else {
              cur_prob <- (n1/(n0+n1))*(1-d)+(ones_so_far/total_so_far)*d
            }
            
            nextval <- rbinom(1,1,cur_prob)
            W[i,subp] <- nextval
            next_subpopulation[i] <- nextval
            ifelse(nextval==1,n1<-n1-1,n0<-n0-1)
            
          }
        }
        
        # This subtle line of code limits the extent of linkage disequilibrium. Toggle the comment and observe what happens.
        #W <- W[,sample(1:k,k)]
        
        X_vector <- c(W)
        X[,snp] <- X_vector
      }
      
      # Measure the amount of structure as the average LD
      average_LD <- sum(cor(X)^2)/(p^2)
      
      
      #Vectors to store results - rank of first causal variant to be discovered
      # and maf of that variant
      GWAS_ranks <- rep(0,length(num_causal_large_p))
      GWAS_mafs <- rep(0,length(num_causal_large_p))
      Eval_ranks <- rep(0,length(num_causal_large_p))
      Eval_mafs <- rep(0,length(num_causal_large_p))
      LASSO_ranks <- rep(0,length(num_causal_large_p))
      LASSO_mafs <- rep(0,length(num_causal_large_p))
      
      
      # Now we have filled in the X matrix, we use it to generate Y
      # This depends on heritability and the allele frequency of the causal SNP(s)
      # So first, we choose how many casual SNPs we want
      # Once we have generated Y we apply our methods to find causal SNP(s)
      for(t in 1:length(num_causal_large_p)){
        
        num_causal <- num_causal_large_p[t]
        
        # We need to divide the hertiability between the causal variants
        h2 <- h2/num_causal
        
        signal <- 0
        for(i in 1:num_causal){
          # determine the coefficients of the causal vector(s) in the model by their maf
          f <- sum(X[,i])/n
          beta <- (((h2/(1-h2))*(1/(f*(1-f))))^.5) * sigma
          # signal = sum(beta_i*X_i) 
          signal <- signal + beta*X[,i]
        }
        
        h2 <- h2*num_causal
        
        # Generate the Y vector
        Y <- signal + noise
        
        
        # The next step is to successively apply each method and see where they 
        # rank the causal SNPs, recording the lowest (best) rank
        
        ###################################GWAS###################################
        # Compute the p-values
        GWAS_pvals <- rep(0,p)
        for (j in 1:p)
        {
          GWAS_model <- lm(Y~X[,j])
          GWAS_pvals[j] <- summary(GWAS_model)$coefficients[2,4]
        }
        
        # Record the ranking of the causal SNP
        sortedp <- sort(GWAS_pvals,decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and it's maf 
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        GWAS_ranks[t] <- causal_rank
        GWAS_mafs[t] <- causal_maf
        
        
        
        
        ###################################Evalue#####################################
        #Find the list of E values
        E_vals <- rep(0,p)
        for (j in 1:p)
        {
          cond_pvals <- rep(0,p)
          for (s in 1:p)
          {
            if((sum(X[,j]!=X[,s]) == 0)|(sum(X[,j]==X[,s]) == 0))
            {
              cond_pvals[s] <- NA
            }
            else
            {
              two_variable_model <- lm(Y~X[,j]+X[,s])
              cond_pvals[s] <- summary(two_variable_model)$coefficients[2,4]
            }
            E_vals[j] <- (sum(cond_pvals[!is.na(cond_pvals)]) + (1/2)*sum(is.na(cond_pvals)))/p
          }
        }
        
        # Record the ranking of the causal SNP
        sortedE <- sort(abs(E_vals),decreasing=FALSE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        Eval_ranks[t] <- causal_rank
        Eval_mafs[t] <- causal_maf
        
        
        ##################################LASSO#######################################
        #apply the glm and find the coefficients of the columns of the matrix
        fit <- cv.glmnet(x=as.matrix(X), y=Y, alpha=1)
        cvec<-as.matrix(coef(fit,s='lambda.1se'))
        num_skip <- length(cvec)-p
        new_vec <- cvec[(1+num_skip):(p+num_skip),]
        
        # Record the ranking of the causal SNP
        sortedp <- sort(abs(new_vec),decreasing=TRUE,index.return=TRUE)
        #Find the highest ranked causal variant and its maf
        causal_rank <- p
        causal_maf <- 1
        for(index in 1:num_causal){
          causal_rank <- ifelse(causal_rank < which(sortedp$ix == index),causal_rank,which(sortedp$ix == index))
          causal_maf <- ifelse(causal_rank < which(sortedp$ix == index),causal_maf,sum(X[,index])/n)
        }
        LASSO_ranks[t] <- causal_rank
        LASSO_mafs[t] <- causal_maf
        
        
      }
      
      
      #Store the results in a tibble! Including the measured structure for X
      if(iter == 1){
        fixed_dh2_tibble <- tibble(rep = rep(iter,length.out= length(num_causal_large_p)),
                                   num_causal = num_causal_large_p, GWAS_ranks = GWAS_ranks, 
                                   GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                                   Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                                   LASSO_mafs = LASSO_mafs, average_LD 
                                   = rep(average_LD,length.out = length(num_causal_large_p)))
      } else {
        next_block <- tibble(rep = rep(iter,length.out= length(num_causal_large_p)),
                             num_causal = num_causal_large_p, GWAS_ranks = GWAS_ranks, 
                             GWAS_mafs = GWAS_mafs, Eval_ranks = Eval_ranks,
                             Eval_mafs = Eval_mafs, LASSO_ranks = LASSO_ranks,
                             LASSO_mafs = LASSO_mafs, average_LD 
                             = rep(average_LD,length.out = length(num_causal_large_p)))
        fixed_dh2_tibble <- bind_rows(fixed_dh2_tibble,next_block)
      }
      
      
      
    }
    
    # We record and compile our results for each value of d
    # Note the important assumption that the first value in dependency_weights is 0
    if(d == 0){
      fixed_h2_tibble <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
    } else {
      next_block <- mutate(fixed_dh2_tibble,dependency_weight = rep(d,nrow(fixed_dh2_tibble)))
      fixed_h2_tibble <- bind_rows(fixed_h2_tibble,next_block)
    }
    
    
  }
  
  # Finally, compile the results for heritability
  if(h2 == heritabilities[1]){
    final_results <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
  } else {
    next_block <- mutate(fixed_h2_tibble, heritability = rep(h2,nrow(fixed_h2_tibble)))
    final_results <- bind_rows(final_results,next_block)
  }
  
  
}


# Now that we have produced our results we want to manipulate them so they are in our desired form.
results_n100_p1000 <- final_results  %>%  
  pivot_longer(., GWAS_ranks:LASSO_mafs,names_to = c("method", ".value"),names_pattern = "([^_]+)_([^_]+)") %>% 
  arrange(., heritability,dependency_weight,num_causal,rep) %>% 
  mutate(., rank_over_p = ranks/p) %>% 
  select(., -c(ranks))

# And save them to a csv file!
file_name <- paste("Results_n",as.character(n),"_p",as.character(p),".csv",sep = "")

write.csv(results_n100_p1000,file_name,row.names = FALSE)




