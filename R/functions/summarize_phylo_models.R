#code to summarize the phylo models
#note that code is very specific to this project and is not generally useful


summarize_phyl0_models <- function(digits_to_round = 3){
  
    phylo_means  <- ls(pattern = "\\.mean",name = globalenv())
    
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


