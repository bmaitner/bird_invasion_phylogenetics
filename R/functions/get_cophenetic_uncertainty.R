#effects of uncertainty

  trees <- list.files(path = "data/bird_phylogeny_updated_names/",full.names = T)
  species_list <- as.character(unique(intros$species_i))
  
get_cohpenetic_uncertainty <- function(trees,species_list){

  #make output structure
  tree_i <- read.tree(file = trees[1])
  tree_i <- drop.tip(phy = tree_i,tip = setdiff(x = tree_i$tip.label,y = species_list))
  coph_i <- cophenetic(tree_i)
  coph_i<-coph_i[order(colnames(coph_i)),order(colnames(coph_i))]
  names_coph_i <- colnames(coph_i)
  
  out_array <- array(dim = c(length(tree_i$tip.label),length(tree_i$tip.label),length(trees)))

        for(i in 1:length(trees)){
          tree_i <- read.tree(file = trees[i])
          tree_i <- drop.tip(phy = tree_i,tip = setdiff(x = tree_i$tip.label,y = species_list))
          coph_i <- cophenetic(tree_i)
          coph_i <- coph_i[order(colnames(coph_i)),order(colnames(coph_i))]
          
          if(!identical(names_coph_i,colnames(coph_i)   )  ){stop("Names wrong")}else{names_coph_i <- colnames(coph_i)}
          out_array[,,i]<- coph_i
            
          
        }#for loop  
        
    
          #generate mean and variance
          coph_means <- apply(X = out_array,MARGIN = c(1,2),FUN = mean)
          coph_se <- apply(X = out_array,MARGIN = c(1,2),FUN = function(x){ sd(x)/ sqrt(length(x) ) })

      coph_out_list <- list()
      coph_out_list[[1]] <- out_array
      coph_out_list[[2]] <- coph_means
      coph_out_list[[3]] <- coph_se
      names(coph_out_list) <- c("coph_dist_array","coph_means","coph_standard_error")      
      rm(coph_i,names_coph_i,tree_i,out_array,coph_means,coph_se,i)
      return(coph_out_list) 
}


#write.csv(x = coph_out_list$coph_means,file = "figures_and_tables/mean_cophenetic_distance_100_phylogenies.csv")
#write.csv(x = coph_out_list$coph_standard_error,file = "figures_and_tables/se_cophenetic_distance_100_phylogenies.csv")

