library(phyr)
library(ape)

#Load Sol data

sol<-read.csv("Data/1221523Database1.csv",stringsAsFactors = F)
sol_nz<-subset(sol,subset = sol$Location_of_introduction == "New_Zealand")
us_nz<-subset(intros,intros$site_i=="N")

sol_hi<-subset(sol,subset = sol$Location_of_introduction == "Hawaiian_Islands")
us_hi<-subset(intros,intros$site_i=="H")


#The data in NZ is broken down by different introductions, need to convert to #events, # individuals

summarized_nz_output<-NULL
for(i in 1:length(unique(sol_nz$Species))){
  species_i<-as.character(unique(sol_nz$Species)[i] ) 
  data_i<-subset(sol_nz,sol_nz$Species==species_i)  
  
  n_events<-nrow(data_i)
  n_individuals<-sum(data_i$Propagule_size,na.rm = T)
  
  if(sum(data_i$Outcome_binary)>0){success<-1}else{success<-0}  
  
  species_data_i<-unique(data_i[,c(9:32,34:35)])  
  climate_match_mean<-mean(data_i$Climate_matching,na.rm = T)  
  
  output_i<-cbind(species_i,success,n_individuals,n_events,species_data_i,climate_match_mean)
  summarized_nz_output<-rbind(summarized_nz_output,output_i)
  
  
}

rm(species_i,data_i,n_events,n_individuals,success,species_data_i,climate_match_mean,output_i,i)



#The data in HI is broken down by different introductions, need to convert to #events, # individuals

summarized_hi_output<-NULL
for(i in 1:length(unique(sol_hi$Species))){
  species_i<-as.character(unique(sol_hi$Species)[i] ) 
  data_i<-subset(sol_hi,sol_hi$Species==species_i)  
  
  n_events<-nrow(data_i)
  n_individuals<-sum(data_i$Propagule_size,na.rm = T)
  
  if(sum(data_i$Outcome_binary)>0){success<-1}else{success<-0}  
  
  species_data_i<-unique(data_i[,c(9:32,34:35)])  
  climate_match_mean<-mean(data_i$Climate_matching,na.rm = T)  
  
  output_i<-cbind(species_i,success,n_individuals,n_events,species_data_i,climate_match_mean)
  summarized_hi_output<-rbind(summarized_hi_output,output_i)
  
  
}

rm(species_i,data_i,n_events,n_individuals,success,species_data_i,climate_match_mean,output_i,i)



#Combine with our data:
full_data_set_nz <- merge(x = us_nz,y = summarized_nz_output,by = "species_i")
full_data_set_nz <- full_data_set_nz[which(!is.na(full_data_set_nz$site_i)),]
full_data_set_nz <- full_data_set_nz[which(!is.na(full_data_set_nz$success)),]

full_data_set_hi <- merge(x = us_hi,y = summarized_hi_output,by = "species_i")
full_data_set_hi <- full_data_set_hi[which(!is.na(full_data_set_hi$site_i)),]
full_data_set_hi <- full_data_set_hi[which(!is.na(full_data_set_hi$success)),]

full_data_set<- rbind(full_data_set_hi,full_data_set_nz)
rm(full_data_set_hi,full_data_set_nz,us_nz,us_hi,sol_hi,sol_nz,sol)


# Traits included in model selection were those identified by Sol et al. (2012) in their best-fit model and included 
  # body mass, the residuals of brain mass against body mass, the value of a brood relative to expected lifetime reproductive output, and habitat generalism (the number of habitat types used by a species). 
# We used stepwise model selection (R Core Team 2016) to test hypotheses based on 
  # 1) Propagule pressure only; 
  # 2) Propagule pressure and MPD metrics; and 
  # 3) Propagule pressure and species’ traits. 
# Stepwise model selection was based on AIC and included both forward and backward selection.  
# Propagule pressure metrics included both the number of individuals and number of introduction events. 
# MPD metrics included both source and recipient community MPD.



full_model = r_success ~ s_range_size + s_richness + s_pd + s_nnd + s_mpd + s_vpd + s_spd + s_kpd + r_n_nnd + r_n_mpd + r_n_vpd + r_n_spd + r_n_kpd + 
r_e_nnd + r_e_mpd + r_e_vpd + r_e_spd + r_e_kpd + r_ne_nnd + r_ne_mpd + r_ne_vpd + r_ne_spd + r_ne_kpd + 
n_individuals + n_events + Body_mass + Brain_residual + Brood_value + Habitat_generalism + (1|site_i) + (1|species_i__)


#phylogenetic model selection

#tree_list is a character vector of phylogenies to use
phylo_model_selection <- function(tree = tree_list[1], model=full_model,data=full_data_set,family ="binomial",species_column = "species_i"){
      
      tree <- read.tree(tree)
      
      
      stop("add code to prune tree to relevant species here")
  
      #Get model components
      rhs <- as.character(full_model)[grep(pattern = "+",x = as.character(full_model),fixed = T)]    
      lhs <- as.character(full_model)[grep(pattern = "+",x = as.character(full_model),fixed = T,invert = T)]    
      lhs <- lhs[which(lhs != "~")]
      
      #Get model variables
      rhs_model_terms <- strsplit(x = rhs, split = " + ",fixed = T)[[1]]      
      
      model_variables <- sapply(X = rhs_model_terms,FUN = function(x){
              x <- strsplit(x = x,split = "(1 | ",fixed = T)[[1]]
              x <- x[which(x != "")]
              x <- strsplit(x = x,split = ")",fixed = T)[[1]]
              x <- gsub(pattern = "__",replacement = "",x = x)        
          
          })
      
      model_variables <- c(lhs,model_variables)
      
      #Drop species that are missing variables    
  
      data <- data[model_variables]    
      data <- na.omit(data)
      
      
      
      #Prepare output structure
      output <- matrix(ncol = length(grep(rhs_model_terms,pattern = "(",fixed = T,invert = T))+2,nrow = 1 )
      output <- as.data.frame(output)
      colnames(output) <- c("model","AIC",rhs_model_terms[grep(rhs_model_terms,pattern = "(",fixed = T,invert = T)])
      
      
      #Make initial model
          rhs_permanent_terms <- rhs_model_terms[grep(rhs_model_terms,pattern = "(",fixed = T,invert = F)]
          
          current_formula <- paste(lhs," ~ ",paste(rhs_model_terms,collapse = " + "))
          #as.formula(current_formula)
          
          model_out_i <- pglmm(formula = as.formula(current_formula) ,
                                data = data,
                                family = family,
                                cov_ranef = list(species_i = tree))    
          
          
          #populate output
          
          output$AIC[1] <- model_out_i$AIC
          output$model[1] <- current_formula
      
        #Compare model to alternatives
      
          
          for(i in 3:ncol(output)){
            
            formula_i <- current_formula
            formula_i2<- gsub(pattern = paste(colnames(output)[i]," \\+",sep = ""),replacement = "",x = formula_i)
            
            if(formula_i == formula_i2){stop("problem updating formulas")}else{formula_i <- formula_i2; rm(formula_i2)}
            model_out_i <- NULL
            
            try(model_out_i <- pglmm(formula = as.formula(formula_i) ,
                                 data = data,
                                 family = family,
                                 cov_ranef = list(species_i = tree)) )
            
            
            if(is.null(model_out_i)){output[1,i] <- NA}else{output[1,i] <- model_out_i$AIC}
            
            
            
              
          }
          
      
          #While loop
            #note, to save computational time, compare formulas to previously tested formulas prior to testing.  no need to re-fit the same model
            #to decide whether term needs to be added or removed from model, grep model term against current formula
            #which.min aics to decide which to use next.  if the lowest is equal to or higher than the current model, stop (can't do any better)
      
        
              current_aic <- output$AIC[1]  
              min_aic <- min(output[1,3:ncol(output)],na.rm = T)
              iteration <- 1
              
              
              
              
                      
              while(current_aic > min_aic){
                
                #keep track of the row we're populating.  Also good to print to make sure things are still working
                iteration <- iteration + 1 
                cat("model iteration ", iteration," \n current model: ",output$model[iteration-1], " \n with an AIC of",output$AIC[iteration-1],"\n",sep = "")
                
                #add to data.frame
                new_row <- output[1,]
                new_row[1:ncol(new_row)] <- NA
                output <- rbind(output,new_row)
                
                #update overall model based on the best available option (adding or removing)
                
                current_formula <- output$model[iteration-1]
                focal_variable  <- colnames(output)[3:ncol(output)][which.min(output[(iteration-1),3:ncol(output)])]
                
                #if the variable is present, remove it
                if(length(grep(pattern = focal_variable,x = current_formula))>0){new_current_formula <- gsub(pattern = paste(focal_variable," \\+",sep = ""),replacement = "",x = current_formula) }
                
                #if the variable is lacking, add it
                if(length(grep(pattern = focal_variable,x = current_formula))==0){new_current_formula <- gsub(pattern = " ~ ",replacement = paste(" ~ ",focal_variable," \\+",sep = ""),x = current_formula)}
                
                #update formula
                current_formula <- new_current_formula
                
                #fit current model
                
                model_out_i <- pglmm(formula = as.formula(current_formula) ,
                                     data = data,
                                     family = family,
                                     cov_ranef = list(species_i = tree))    
                
                #populate output
                
                output$AIC[iteration] <- model_out_i$AIC
                output$model[iteration] <- current_formula
                
                
                  #Compare model to alternatives
                  
                  
                  for(i in 3:ncol(output)){
                    
                    formula_i <- current_formula
                    focal_variable <- colnames(output)[i]
                    
                    
                    #if the variable is present, remove it
                    if(length(grep(pattern = focal_variable,x = current_formula))>0){formula_i2 <- gsub(pattern = paste(focal_variable," \\+",sep = ""),replacement = "",x = current_formula) }
                    
                    #if the variable is lacking, add it
                    if(length(grep(pattern = focal_variable,x = current_formula))==0){formula_i2 <- gsub(pattern = " ~ ",replacement = paste(" ~ ",focal_variable," \\+",sep = ""),x = current_formula)}
                    
                    
                    
                    if(formula_i == formula_i2){stop("problem updating formulas")}else{formula_i <- formula_i2; rm(formula_i2)}
                    model_out_i <- NULL
                    
                    try(model_out_i <- pglmm(formula = as.formula(formula_i) ,
                                             data = data,
                                             family = family,
                                             cov_ranef = list(species_i = tree)) )
                    
                    #add AIC to output file
                    if(is.null(model_out_i)){output[iteration,i] <- NA }else{ output[iteration,i] <- model_out_i$AIC }
                    
                    
                    
                    
                  }#i loop
                
                #check whether you can do better in another iteration, if not,  leave the loop
                min_aic <- min(output[iteration,3:ncol(output)],na.rm = T)  
                
                if(output$AIC[iteration] <= min_aic ){break()}  
                
                
                
              }#while loop  
          
          
              
  
  
  
  
}


###########################################################################################################################



#phylogenetic model selection

#tree_list is a character vector of phylogenies to use
phylo_model_selection_backward <- function(tree = tree_list[1], model=full_model,data=full_data_set,family ="binomial",species_column = "species_i"){
  
  tree <- read.tree(tree)
  
  #Get model components
  rhs <- as.character(full_model)[grep(pattern = "+",x = as.character(full_model),fixed = T)]    
  lhs <- as.character(full_model)[grep(pattern = "+",x = as.character(full_model),fixed = T,invert = T)]    
  lhs <- lhs[which(lhs != "~")]
  
  #Get model variables
  rhs_model_terms <- strsplit(x = rhs, split = " + ",fixed = T)[[1]]      
  
  model_variables <- sapply(X = rhs_model_terms,FUN = function(x){
    x <- strsplit(x = x,split = "(1 | ",fixed = T)[[1]]
    x <- x[which(x != "")]
    x <- strsplit(x = x,split = ")",fixed = T)[[1]]
    x <- gsub(pattern = "__",replacement = "",x = x)        
    
  })
  
  model_variables <- c(lhs,model_variables)
  
  #Drop species that are missing variables    
  
  data <- data[model_variables]    
  data <- na.omit(data)
  
  
  
  #Prepare output structure
  output <- matrix(ncol = length(grep(rhs_model_terms,pattern = "(",fixed = T,invert = T))+2,nrow = 1 )
  output <- as.data.frame(output)
  colnames(output) <- c("model","AIC",rhs_model_terms[grep(rhs_model_terms,pattern = "(",fixed = T,invert = T)])
  
  
  #Make initial model
  rhs_permanent_terms <- rhs_model_terms[grep(rhs_model_terms,pattern = "(",fixed = T,invert = F)]
  
  current_formula <- paste(lhs," ~ ",paste(rhs_model_terms,collapse = " + "))
  #as.formula(current_formula)
  
  model_out_i <- pglmm(formula = as.formula(current_formula) ,
                       data = data,
                       family = family,
                       cov_ranef = list(species_i = tree))    
  
  
  #populate output
  
  output$AIC[1] <- model_out_i$AIC
  output$model[1] <- current_formula
  
  #Compare model to alternatives
  
  for(i in 3:ncol(output)){
    
    formula_i <- current_formula
    formula_i2<- gsub(pattern = paste(colnames(output)[i]," \\+",sep = ""),replacement = "",x = formula_i)
    
    if(formula_i == formula_i2){stop("problem updating formulas")}else{formula_i <- formula_i2; rm(formula_i2)}
    model_out_i <- NULL
    
    try(model_out_i <- pglmm(formula = as.formula(formula_i) ,
                             data = data,
                             family = family,
                             cov_ranef = list(species_i = tree)) )
    
    
    if(is.null(model_out_i)){output[1,i] <- NA}else{output[1,i] <- model_out_i$AIC}
    
    
    
    
  }
  
  
  #While loop
      #note, to save computational time, compare formulas to previously tested formulas prior to testing.  no need to re-fit the same model
  current_aic <- output$AIC[1]  
  min_aic <- min(output[1,3:ncol(output)],na.rm = T)
  iteration <- 1
  
  
  while(current_aic > min_aic){
    
    #keep track of the row we're populating.  Also good to print to make sure things are still working
    iteration <- iteration + 1 
    cat("model iteration ", iteration," \n current model: ",output$model[iteration-1], " \n with an AIC of",output$AIC[iteration-1],"\n",sep = "")
    
    #add to data.frame
    new_row <- output[1,]
    new_row[1:ncol(new_row)] <- NA
    output <- rbind(output,new_row)
    
    #update overall model based on the best available option (adding or removing)
    
    current_formula <- output$model[iteration-1]
    focal_variable  <- colnames(output)[3:ncol(output)][which.min(output[(iteration-1),3:ncol(output)])]
    
    #if the variable is present, remove it
    if(length(grep(pattern = focal_variable,x = current_formula))>0){new_current_formula <- gsub(pattern = paste(focal_variable," \\+",sep = ""),replacement = "",x = current_formula) }
    
    #if the variable is lacking, add it
    if(length(grep(pattern = focal_variable,x = current_formula))==0){stop("BRIAN FIX CODE")}
    
    #update formula
    current_formula <- new_current_formula
    
    #fit current model
    
    model_out_i <- pglmm(formula = as.formula(current_formula) ,
                         data = data,
                         family = family,
                         cov_ranef = list(species_i = tree))    
    
    #populate output
    
    output$AIC[iteration] <- model_out_i$AIC
    output$model[iteration] <- current_formula
    
    
    #Compare model to alternatives
    
    
    for(i in 3:ncol(output)){
      
      formula_i <- current_formula
      focal_variable <- colnames(output)[i]
      
      
      #if the variable is present, remove it
      if(length(grep(pattern = focal_variable,x = current_formula))>0){formula_i2 <- gsub(pattern = paste(focal_variable," \\+",sep = ""),replacement = "",x = current_formula) }
      
      #if the variable is lacking, skip to the next variable
      #if(length(grep(pattern = focal_variable,x = current_formula))==0){formula_i2 <- gsub(pattern = " ~ ",replacement = paste(" ~ ",focal_variable," \\+",sep = ""),x = current_formula)}
      if(length(grep(pattern = focal_variable,x = current_formula))==0){next}
      
      
      if(formula_i == formula_i2){stop("problem updating formulas")}else{formula_i <- formula_i2; rm(formula_i2)}
      model_out_i <- NULL
      
      try(model_out_i <- pglmm(formula = as.formula(formula_i) ,
                               data = data,
                               family = family,
                               cov_ranef = list(species_i = tree)) )
      
      #add AIC to output file
      if(is.null(model_out_i)){output[iteration,i] <- NA }else{ output[iteration,i] <- model_out_i$AIC }
      
      
      
      
    }#i loop
    
    #check whether you can do better in another iteration, if not,  leave the loop
    min_aic <- min(output[iteration,3:ncol(output)],na.rm = T)  
    
    if(output$AIC[iteration] <= min_aic ){break()}  
    
    
    
  }#while loop  
  
  
  
  
  
  
  
}














