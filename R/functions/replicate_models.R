

replicate_phylo_glmm <- function(tree_list, formula, data,family ="binomial"){
        out_list<-list()
       
        for(i in 1:length(tree_list)){
          tree <- read.tree(tree_list[i])

          out_list[[i]] <- 
          
          pglmm(formula = formula ,
                data = data,
                family = family,
                cov_ranef = list(species_i = tree))    
         
          
          out_list[[i]]$partial_R2 <- suppressMessages(rr2::R2(mod = out_list[[i]]))  
           
        }
        
        
    return(out_list)
  
}


standard_error <- function(x){sd(x)/sqrt(length(x))}


summarize_replicates <- function( phylo_reps_output){
  
  out <- array(dim = c(length(phylo_reps_output[[1]]$B),
                        3,
                        length(phylo_reps_output)))

  colnames(out)<- c("estimate","se","p_val")
  rownames(out)<- rownames(phylo_reps_output[[1]]$B)
  aic_list <- NULL
  R2_list <- NULL
  
  for(i in 1:length(phylo_reps_output)){
    
    out[,,i] <- cbind(phylo_reps_output[[i]]$B, #estimate
                      phylo_reps_output[[i]]$B.se, #standard error 
                      phylo_reps_output[[i]]$B.pvalue #pvale
                      )
    aic_list <- c(aic_list,phylo_reps_output[[i]]$AIC)
    R2_list <-  c(R2_list,phylo_reps_output[[i]]$partial_R2)
    
  } #i loop 
  

  
  out_mean <- apply(X = out,MARGIN = c(1,2),FUN = mean)
  aic_mean <- mean(aic_list)
  R2 <- cbind(mean(R2_list),standard_error(R2_list))
  colnames(R2)<- c("mean","SE")
  out_list<-list(means=out_mean,aic_mean=aic_mean,R2 = R2)    
  return(out_list)  
}



