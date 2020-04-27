#code to summarize the phylo models
#note that code is very specific to this project and is not generally useful


summarize_phyl0_models <- function(digits_to_round = 3){
  
    phylo_means  <- ls(pattern = "\\.mean",name = globalenv())
    phylo_means <- phylo_means[grep(pattern = "phylo\\.full",x = phylo_means)]

      #First, get all model parameters to build an output matrix to hold things  
    all_parms <- NULL
    for(i in phylo_means){
      
      model_i <- get(i)
      all_parms <-unique(c(all_parms,row.names(model_i$means)))
    
    }
    
  rm(i,model_i)
  
  #make output file
  output <- matrix(nrow = length(phylo_means),ncol = length(all_parms)+2)
  output <- as.data.frame(output)
  colnames(output) <- c("model","aic",all_parms)
  
  #populate output file
    for(i in 1:length(phylo_means)){
      
      model_i <- get(phylo_means[i])
      output$model[i] <- phylo_means[i]
      output$aic[i] <- round(x = model_i$aic_mean,digits = digits_to_round)      
        #Append bits to designate p-val
          for( j in 1: nrow(model_i$means)){
           
            as.numeric(model_i$means[j,"estimate"])
            if(as.numeric(model_i$means[j,"p_val"]) > 0.1){output[i,rownames(model_i$means)[j]] <-   round(x = as.numeric(model_i$means[j,"estimate"]),digits = digits_to_round);next}
            as.numeric(model_i$means[j,"estimate"])
            if(as.numeric(model_i$means[j,"p_val"]) < 0.1 &  model_i$means[j,"p_val"] > 0.05){output[i,rownames(model_i$means)[j]] <-   paste(round(x = as.numeric(model_i$means[j,"estimate"]),digits = digits_to_round),".");next  }
            as.numeric(model_i$means[j,"estimate"])
            if(as.numeric(model_i$means[j,"p_val"]) < 0.05 &  model_i$means[j,"p_val"] > 0.01 ){  output[i,rownames(model_i$means)[j]] <-   paste(round(x = as.numeric(model_i$means[j,"estimate"]),digits = digits_to_round),"*") ;next            }
            as.numeric(model_i$means[j,"estimate"])
            if(as.numeric(model_i$means[j,"p_val"]) < 0.01 ){   output[i,rownames(model_i$means)[j]] <-   paste(round(x = as.numeric(model_i$means[j,"estimate"]),digits = digits_to_round),"**") ;next            }
            as.numeric(model_i$means[j,"estimate"])
            
            #output[i,rownames(model_i$means)[j]] <- model_i$means[j,"estimate"]
            
            
          }
        
    
      
      
      
      
      
      
    }# i loop

  
  
  colnames(output)<- qdap::multigsub(pattern = c("s_range_size","s_richness","r_e_kpd","r_n_kpd","r_ne_kpd","s_kpd","r_e_mpd","r_n_mpd","r_ne_mpd","s_mpd","r_e_nnd","r_n_nnd","r_ne_nnd","s_nnd","r_e_spd",
                              "r_n_spd","r_ne_spd","s_spd","r_e_vpd","r_n_vpd","r_ne_vpd","s_vpd","s_pd" ),
                  replacement = c("Range size s","Richness s",
                                  "KPD r,e","KPD r,n","KPD r,ne","KPD s",
                                  "MPD r,e","MPD r,n","MPD r,ne","MPD s",
                                  "NND r,e","NND r,n","NND r,ne","NND s",
                                  "SPD r,e","SPD r,n","SPD r,ne","SPD s",
                                  "VPD r,e","VPD r,n","VPD r,ne","VPD s",
                                  "PD s" ),text.var = colnames(output)
                  )
  

    return(output)

}


##########################################################################



#code to summarize the phylo models
#note that code is very specific to this project and is not generally useful


summarize_nonphylo_models <- function(digits_to_round = 3){
  
  model_objects  <- ls(pattern = "full_",name = globalenv())
  model_objects <- model_objects[grep(pattern = ".",fixed = T,x = model_objects,invert = T)]
  model_objects <- model_objects[grep(pattern = "data",fixed = T,x = model_objects,invert = T)]
  model_objects <- model_objects[grep(pattern = "model",fixed = T,x = model_objects,invert = T)]
  
  #First, get all model parameters to build an output matrix to hold things  
  all_parms <- NULL
  for(i in model_objects){
    
    model_i <- get(i)
    x<-summary(model_i)
    all_parms <-unique(c(all_parms,rownames(x$coefficients)))
    rm(x)
  }
  
  rm(i,model_i)
  
  #make output file
  output <- matrix(nrow = length(model_objects),ncol = length(all_parms)+2)
  output <- as.data.frame(output)
  colnames(output) <- c("model","aic",all_parms)
  
  #populate output file
  for(i in 1:length(model_objects)){
    
    model_i <- get(model_objects[i])
    output$model[i] <- model_objects[i]
    
    sum_i <- summary(model_i)
    output$aic[i] <- sum_i$AICtab[1]
    sum_i$coefficients
    
    
    #Append bits to designate p-val
    for( j in 1: nrow(sum_i$coefficients)){
      
      as.numeric(sum_i$coefficients[j,"Estimate"])
      if(as.numeric(sum_i$coefficients[j,"Pr(>|z|)"]) > 0.1){output[i,rownames(sum_i$coefficients)[j]] <-   round(x = as.numeric(sum_i$coefficients[j,"Estimate"]),digits = digits_to_round);next}
      as.numeric(sum_i$coefficients[j,"Estimate"])
      if(as.numeric(sum_i$coefficients[j,"Pr(>|z|)"]) < 0.1 &  sum_i$coefficients[j,"Pr(>|z|)"] > 0.05){output[i,rownames(sum_i$coefficients)[j]] <-   paste(round(x = as.numeric(sum_i$coefficients[j,"Estimate"]),digits = digits_to_round),".");next  }
      as.numeric(sum_i$coefficients[j,"Estimate"])
      if(as.numeric(sum_i$coefficients[j,"Pr(>|z|)"]) < 0.05 &  sum_i$coefficients[j,"Pr(>|z|)"] > 0.01 ){  output[i,rownames(sum_i$coefficients)[j]] <-   paste(round(x = as.numeric(sum_i$coefficients[j,"Estimate"]),digits = digits_to_round),"*") ;next            }
      as.numeric(sum_i$coefficients[j,"Estimate"])
      if(as.numeric(sum_i$coefficients[j,"Pr(>|z|)"]) < 0.01 ){   output[i,rownames(sum_i$coefficients)[j]] <-   paste(round(x = as.numeric(sum_i$coefficients[j,"Estimate"]),digits = digits_to_round),"**") ;next            }
      as.numeric(sum_i$coefficients[j,"Estimate"])
      
      #output[i,rownames(model_i$means)[j]] <- model_i$means[j,"Estimate"]
      
      
    }
    
    
    
    
    
    
    
    
  }# i loop
  
  
  
  colnames(output)<- qdap::multigsub(pattern = c("s_range_size","s_richness","r_e_kpd","r_n_kpd","r_ne_kpd","s_kpd","r_e_mpd","r_n_mpd","r_ne_mpd","s_mpd","r_e_nnd","r_n_nnd","r_ne_nnd","s_nnd","r_e_spd",
                                                 "r_n_spd","r_ne_spd","s_spd","r_e_vpd","r_n_vpd","r_ne_vpd","s_vpd","s_pd" ),
                                     replacement = c("Range size s","Richness s",
                                                     "KPD r,e","KPD r,n","KPD r,ne","KPD s",
                                                     "MPD r,e","MPD r,n","MPD r,ne","MPD s",
                                                     "NND r,e","NND r,n","NND r,ne","NND s",
                                                     "SPD r,e","SPD r,n","SPD r,ne","SPD s",
                                                     "VPD r,e","VPD r,n","VPD r,ne","VPD s",
                                                     "PD s" ),text.var = colnames(output)
  )
  
  
  return(output)
  
}
