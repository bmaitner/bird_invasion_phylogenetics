#Code to ensure phylogenetic names match.
#To save time, I'm only going to correct the relavent species names, under the assumption Birdlife is correct

#So, I'll correct:
#Introduced species names
#Species that occur in focal regions
#Species that co-occur with introduced species

#load in birdlife taxonomy
bl_taxonomy <- read.csv("data/BirdLife_Checklist_Version_80/Checklist v8 Oct15/BirdLife_Checklist_Version_8.csv",stringsAsFactors = F)

#rename species to match format of phylo, introduction data
bl_taxonomy$Scientific.name <- gsub(pattern = " ",replacement = "_",x = bl_taxonomy$Scientific.name)

#load in introduction data
introductions <- read.csv("data/invasion_data.csv",stringsAsFactors = F)

#cut australia data (too large an area)
introductions<-introductions[which(introductions$Site!="A"),]

#load occs to check names against
occs <- readRDS("data/cea_occurrences.rds")
occs$species<-as.character(occs$species)
occs$species <- gsub(pattern = " ",replacement = "_",x = occs$species)

#####################################

#Check names of introduced species
  #manually corrected introduced names to birdlife checklist

intro_names_to_check <- unique(introductions$Species)
intro_names_to_check <- intro_names_to_check[which(!intro_names_to_check %in% occs$species)]

#one bad name ( hybrid, drop from analysis)
introductions<-introductions[grep(pattern = "_x_",x = introductions$Species,invert = T),]

#no bad names remaining
intro_names_to_check <- unique(introductions$Species)
intro_names_to_check <- intro_names_to_check[which(!intro_names_to_check %in% bl_taxonomy$Scientific.name)]
rm(intro_names_to_check)

#all intro species names are now alright (relative to BL)

#######################################

#Now, check that phylo names match bl names
#first, we need all relevant species from source or recipient regions

occs <- readRDS("data/cea_occurrences.rds")
occs$species <- as.character(occs$species)
occs$species <- gsub(pattern = " ",replacement = "_",x = occs$species)


#get a list of the cells that the introduced species occur in
native_cells <- occs$cell[which(occs$species %in% unique(introductions$Species))]

#use the native cell list to get a list of all relevant species
native_species_to_check <- unique(occs$species[which(occs$cell %in% native_cells)])
rm(native_cells)

#introduced species to check
intro_species_to_check <- unique(introductions$Species)

# recipient region names to check
library(maptools)
library(sf)

gadm <- sf::read_sf("C:/Users/Brian/Desktop/Range maps/gadm28/gadm28.shp") #read in gadm
gadm <- gadm[which((gadm$NAME_0 =="United States" & gadm$NAME_1 %in% c("Hawaii", "Florida")) |
                 gadm$NAME_0 == "New Zealand"),]#cut the gadm down to just relevant spots to save space
gadm <- sf:::as_Spatial(gadm)#convert to spatial

#reproject to crs of template raster
template <- raster()
template <- projectRaster(from = template,crs = crs("+proj=cea +units=km"),res = 110)
gadm <- spTransform(x = gadm,CRSobj = template@crs)
plot(gadm)
gadm@data$name01 <- paste(gadm@data$NAME_0,gadm@data$NAME_1)
gadm@data$name01

focal_cells <- NULL
for(i in 1:length(unique(gadm@data$name01))){

  gadm_raster_i <- rasterize(x = gadm[which(gadm@data$name01==unique(gadm@data$name01)[i]),],y = template,getCover=T)  
  plot(gadm_raster_i)
  focal_cells <- rbind(focal_cells,cbind(which(getValues(gadm_raster_i>0)),gadm_raster_i[which(getValues(gadm_raster_i>0))],unique(gadm@data$name01)[i]))

}
focal_cells <- as.data.frame(focal_cells)
colnames(focal_cells) <- c("cell","pct_cover","political_unit")
focal_cells$cell<-as.numeric(as.character(focal_cells$cell))

#cleanup
rm(gadm_raster_i,i,gadm,template)

recipient_species_to_check <- unique(occs$species[which(occs$cell%in%focal_cells$cell)])

combined_species_to_check <- unique(c(recipient_species_to_check, native_species_to_check))

rm(intro_species_to_check,native_species_to_check,recipient_species_to_check)
##################################################################################################

#Now, we need to check that all of the species in the "combined_species_to_check" file are present in the phylos
#if they are not, need to correct their names

trees <- list.files("data/Bird_Phylogeny/",pattern = "tree_",full.names = T)#these are the birdtree.org phylogenies
tree1 <- read.tree(trees[1])
tree <- read.tree(trees[1])

#try to correct names in an automated fashion first, then manually if need be
#inputs needed:
#tree
#combined_species_to_check
#bl_taxonomy

to_correct <- setdiff(x = combined_species_to_check,y =   tree$tip.label)

bl_taxonomy$Synonyms <- gsub(pattern = " ",replacement = "_",x = bl_taxonomy$Synonyms)  


#first, we'll use the bl synonyms

for(i in 1:length(to_correct)){
  
  name_i<-to_correct[i]  
  print(name_i)
  
  #check if synonym in phylo
  syn_i <- bl_taxonomy$Synonyms[which(bl_taxonomy$Scientific.name==name_i)]
  syn_i <- syn_i[which(syn_i!="")]
  
  
  if(length(syn_i)>1){stop()}
  
  if(length(syn_i)>0){
    #stop()
    syn_i <- unlist(strsplit(x = syn_i,split = ";_"))
    
    
    if(length(syn_i)>0){
      
      
      for(s in 1:length(syn_i)){
        if(length(grep(pattern = syn_i[s],x = tree$tip.label))>0){
          tree$tip.label[grep(pattern = syn_i[s],x = tree$tip.label)] <- name_i  
        }}}}
  
  
  
}# i loop
rm(i,name_i,s,syn_i,to_correct)

# see how many we got

to_correct <- setdiff(x = combined_species_to_check,y =   tree$tip.label)# well, that got ~ half

####################################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#try taxize

library(taxize)

#itis first
itis_syns <- synonyms(to_correct,db="itis")

itis_syns_df <- synonyms_df(x = itis_syns)


for(i in 1:length(unique(itis_syns_df$.id))){
  
  #id field is the good name
  #need to find a match to one of the synonyms in the phylo
  good_name_i <- unique(itis_syns_df$.id)[i]
  bad_names_i <- unique(na.omit(c(itis_syns_df$acc_name[which(itis_syns_df$.id==good_name_i)],itis_syns_df$syn_name[which(itis_syns_df$.id==good_name_i)])))
  
  for(s in 1:length(bad_names_i)){
          
        
      
      if(length(which(tree$tip.label == gsub(pattern = " ",replacement = "_",x = bad_names_i[s],)))>0){
      
      tree$tip.label[which(tree$tip.label == gsub(pattern = " ",replacement = "_",x = bad_names_i[s]))] <- good_name_i
      
    }
    
    
  }#s in bad names 
  
  
  
}#i in names returned by itis

# see how many we got

to_correct <- setdiff(x = combined_species_to_check,y =   tree$tip.label)# well, that got a few
#skipping col because they limit the number of queries...which means they're basically useless here

to_correct  
out<-as.data.frame(to_correct)
out$new_name_phylo <- NA
out$new_name_db <- NA  
colnames(out)[1]<-"bl_name"

#output here and manually correct others
#write.csv(x = out,file = "data/manual_corrections.csv",row.names = F  ) #commented this out so I don't accidentally overwrite

#re-import data
name_corrections <- read.csv("data/manual_corrections.csv",stringsAsFactors = F)

for(i in 1: nrow(name_corrections)){
  
  name_to_replace <- gsub(pattern = " ",replacement = "_",x = name_corrections$new_name_phylo[i])
  new_name <- gsub(pattern = " ",replacement = "_",x = name_corrections$bl_name[i])
  
  if(!is.na(name_to_replace) & name_to_replace !=""){
    
    tree$tip.label[which(tree$tip.label==name_to_replace)] <- new_name
  }
  
  
  
}

to_correct <- setdiff(x = combined_species_to_check,y =   tree$tip.label)# well, that got ~ 60

name_corrections[which(name_corrections$bl_name %in% to_correct),] #all remaining species are either extinct or represent recent taxonomic revisions or new descriptions.


name_translator <- cbind(tree1$tip.label,tree$tip.label)
name_translator <- as.data.frame(name_translator)
name_translator$corrected <- as.character(name_translator$corrected)
colnames(name_translator) <- c("original","corrected")
saveRDS(name_translator,file = "data/name_translator.rds")

#######################################

#Check name of introduced species against phylo

name_translator <- readRDS("data/name_translator.rds")
name_translator$corrected <- as.character(name_translator$corrected)

introductions$Species[which(!introductions$Species %in% name_translator$corrected)]

name_translator$corrected[which(name_translator$original=="Anas_poecilorhyncha")] <- "Anas_poecilorhyncha"
name_translator$corrected[which(name_translator$original=="Paroaria_gularis")] <- "Paroaria_gularis"
name_translator$corrected[which(name_translator$original=="Pitta_guajana")] <- "Pitta_guajana"
name_translator$corrected[which(name_translator$original=="Ramphastos_vitellinus")] <- "Ramphastos_citrolaemus"

saveRDS(name_translator,file = "data/name_translator.rds")


######################################################################################################
######################################################################################################
