#3 Phylodiversity measurements

#for each introduced species x recipient region, calculate:

    #pd,  nnd, mpd, vpd, spd, kpd
      # recipient region metrics relative to : natives only, introduced only, native introduced, 
      # native region metrics

      #review empirical paper first to make sure things are consistent

#inputs: 
  #list of 100 random phylos
  #skinny occs
  #intro data
  #gadm cells of intro regions
  #name translator

#output: data.frame with rows:
  #species
  #region
  #recipient community 
    #n,e,ne
  #native community
  #metric (pd, nnd, mpd, vpd, spd, kpd)
  #range size
  #species richness native
  #species richness intro

#load in data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #load in introduction data
    introductions <- read.csv("data/invasion_data.csv",stringsAsFactors = F)
    
    #cut australia data (too large an area)
    introductions<-introductions[which(introductions$Site!="A"),]
    #cut individual hawaiian islands
    introductions<-introductions[which(introductions$Site!="I"),]
    
    #one bad name ( hybrid, drop from analysis)
    introductions<-introductions[grep(pattern = "_x_",x = introductions$Species,invert = T),]
    
  #load in translator
    translator <- readRDS("data/name_translator.rds")
    translator$original<-as.character(translator$original)
    translator$corrected<-as.character(translator$corrected)
    translator$corrected[which(translator$original=="Eupodotis_afraoides")]<-"Afrotis_afraoides"
    
  #load occurrences
    occs <- readRDS("data/cea_occurrences.rds")
    occs$species<-as.character(occs$species)
    occs$species <- gsub(pattern = " ",replacement = "_",x = occs$species)
  
  #load in recipient region polygons
    library(maptools)
    library(sf)
    library(moments)
    library(PhyloMeasures)
    library(ape)
    
    gadm <- sf::read_sf("C:/Users/Brian/Desktop/Range maps/gadm28/gadm28.shp") #read in gadm
    gadm <- gadm[which((gadm$NAME_0 =="United States" & gadm$NAME_1 %in% c("Hawaii", "Florida")) |
                         gadm$NAME_0 == "New Zealand"),]#cut the gadm down to just relevant spots to save space
    gadm <- sf:::as_Spatial(gadm)#convert to spatial
    
    #reproject to crs of template raster
    template<-raster()
    template<-projectRaster(from = template,crs = crs("+proj=cea +units=km"),res = 110)
    gadm <- spTransform(x = gadm,CRSobj = template@crs)
    gadm@data$name01<-paste(gadm@data$NAME_0,gadm@data$NAME_1)
    
    focal_cells<-NULL
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
  
  #load in list of phylos
    trees <- list.files("data/Bird_Phylogeny/",pattern = "tree_",full.names = T)
    set.seed(2005)
    trees <- trees[sample(x = 1:length(trees),size = 100,replace = F)] #sample 100 phylos to use
  
#calculate phylo metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phy_out <- NULL        

for(p in 1:length(trees)){
  tree_p <- read.tree(trees[p])
    tree_p_taxa<- as.data.frame(tree_p$tip.label)
    tree_p_corrected <- merge(x = tree_p_taxa,y = translator,by.x="tree_p$tip.label",by.y="original",sort=F)
    tree_p$tip.label <- as.character(tree_p_corrected$corrected)
    rm(tree_p_taxa,tree_p_corrected)
    coph_p <- cophenetic(x = tree_p)
    
for(i in 1:nrow(unique(introductions[,c("Species","Site")]))){
  species_i <- unique(introductions[i,c("Species","Site")])$Species
  site_i <- unique(introductions[i,c("Species","Site")])$Site
  
  
  #source region phylo stuff
    
    #get a list of all raster cells where the species is found, where the species is native, reintroduced, or the origin is uncertain
    source_cells <- occs$cell[which(occs$species==species_i & occs$origin %in% c(1,2,5))]
    
    #get a list of all species that occur in the source_cells, again limiting it to those that are native or possible native
    source_species <- unique(occs$species[which(occs$cell %in% source_cells & occs$origin %in% c(1,2,5))])
          
    source_phydist <- coph_p[which(colnames(coph_p) %in% source_species),which(colnames(coph_p) %in% source_species)]
    source_phydist <- source_phydist[,which(row.names(source_phydist)==species_i)]
    source_phydist <- source_phydist[which(source_phydist != 0)]
    
    #prune tree
    #tree_i <- drop.tip(phy = tree_p,tip = which(!tree_p$tip.label %in% source_species))
    comm_s <- matrix(nrow = 1,ncol = length(tree_p$tip.label),data = 0)
    colnames(comm_s) <- tree_p$tip.label
    comm_s[1,][which(colnames(comm_s) %in% source_species)] <- 1
    
    #source output
    s_pd <- pd.query(tree = tree_p,matrix = comm_s)
    s_nnd <- min(source_phydist)
    s_mpd <- mean(source_phydist)  
    s_vpd <- var(source_phydist)
    s_spd <- skewness(source_phydist)
    s_kpd <- kurtosis(source_phydist)
    s_richness <- length(source_species)
    s_range_size <- length(source_cells)
  
  #recipient region stuff
    
    #get a list of all raster cells in recipient region
    recipient_cells <- NULL
    if(site_i=="N"){recipient_cells <- focal_cells$cell[grep(pattern = "New Zealand",x = focal_cells$political_unit)]}
    if(site_i=="F"){recipient_cells <- focal_cells$cell[grep(pattern = "Florida",x = focal_cells$political_unit)]}
    if(site_i=="H"){recipient_cells <- focal_cells$cell[grep(pattern = "Hawaii",x = focal_cells$political_unit)]}
  
    #get a list of all species that occur in the recipient_cells, again limiting it to those that are native or possible native
    recipient_n_species <- c(unique(occs$species[which(occs$cell %in% recipient_cells & occs$origin %in% c(1,2,5))]),species_i)
    
    #get a list of all species that occur in the re_cells, again limiting it to those that are established
    recipient_e_species <- unique(occs$species[which(occs$cell %in% recipient_cells & occs$origin %in% c(3))])
    recipient_e_species <- unique(c(recipient_e_species,introductions$Species[which(introductions$Site==site_i & introductions$Failure==0)],species_i))
    
    #get a list of all species that occur in the source_c
    recipient_ne_species <- unique(c(recipient_n_species,recipient_e_species))
    
    #recipient n phydist
    recipient_n_phydist <- coph_p[which(colnames(coph_p) %in% recipient_n_species),which(colnames(coph_p) %in% recipient_n_species)]
    recipient_n_phydist <- recipient_n_phydist[,which(row.names(recipient_n_phydist)==species_i)]
    recipient_n_phydist <- recipient_n_phydist[which(recipient_n_phydist != 0)]
    
    #recipient e phydist
    recipient_e_phydist <- coph_p[which(colnames(coph_p) %in% recipient_e_species),which(colnames(coph_p) %in% recipient_e_species)]
    recipient_e_phydist <- recipient_e_phydist[,which(row.names(recipient_e_phydist)==species_i)]
    recipient_e_phydist <- recipient_e_phydist[which(recipient_e_phydist != 0)]
    
    #recipient ne phydist
    recipient_ne_phydist <- coph_p[which(colnames(coph_p) %in% recipient_ne_species),which(colnames(coph_p) %in% recipient_ne_species)]
    recipient_ne_phydist <- recipient_ne_phydist[,which(row.names(recipient_ne_phydist)==species_i)]
    recipient_ne_phydist <- recipient_ne_phydist[which(recipient_ne_phydist != 0)]
    
    #recipient communities for PD
    comm_r <- matrix(nrow = 3,ncol = length(tree_p$tip.label),data = 0)
    colnames(comm_r) <- tree_p$tip.label
    rownames(comm_r)<-c("n","e","ne")
    comm_r["n",][which(colnames(comm_r) %in% recipient_n_species)] <- 1
    comm_r["e",][which(colnames(comm_r) %in% recipient_e_species)] <- 1
    comm_r["ne",][which(colnames(comm_r) %in% recipient_ne_species)] <- 1
    
    #recipient output general
    r_range_size <- length(recipient_cells)
    r_failure <- introductions$Failure[which(introductions$Site==site_i & introductions$Species==species_i)]
    
    #recipient output n
    r_n_pd <- pd.query(tree = tree_p,matrix = comm_r)[1]
    r_n_nnd <- min(recipient_n_phydist)
    r_n_mpd <- mean(recipient_n_phydist)  
    r_n_vpd <- var(recipient_n_phydist)
    r_n_spd <- skewness(recipient_n_phydist)
    r_n_kpd <- kurtosis(recipient_n_phydist)
    r_n_richness <- length(recipient_n_species)
    
    #recipient output e
    r_e_pd <- pd.query(tree = tree_p,matrix = comm_r)[2]
    r_e_nnd <- min(recipient_e_phydist)
    r_e_mpd <- mean(recipient_e_phydist)  
    r_e_vpd <- var(recipient_e_phydist)
    r_e_spd <- skewness(recipient_e_phydist)
    r_e_kpd <- kurtosis(recipient_e_phydist)
    r_e_richness <- length(recipient_e_species)
    
    #recipient output ne
    r_ne_pd <- pd.query(tree = tree_p,matrix = comm_r)[3]
    r_ne_nnd <- min(recipient_ne_phydist)
    r_ne_mpd <- mean(recipient_ne_phydist)  
    r_ne_vpd <- var(recipient_ne_phydist)
    r_ne_spd <- skewness(recipient_ne_phydist)
    r_ne_kpd <- kurtosis(recipient_ne_phydist)
    r_ne_richness <- length(recipient_ne_species)
    
        
    phy_out <- rbind(phy_out,
          cbind(species_i,site_i,
          s_range_size,
          s_pd,s_nnd,s_mpd,s_vpd,s_spd,s_kpd,s_richness,
          r_failure,r_range_size,
          r_n_pd,r_n_nnd,r_n_mpd,r_n_vpd,r_n_spd,r_n_kpd,r_n_richness,
          r_e_pd,r_e_nnd,r_e_mpd,r_e_vpd,r_e_spd,r_e_kpd,r_e_richness,
          r_ne_pd,r_ne_nnd,r_ne_mpd,r_ne_vpd,r_ne_spd,r_ne_kpd,r_ne_richness
          ))
  
  #cleanup
  rm(species_i,site_i,
     s_range_size,
     s_pd,s_nnd,s_mpd,s_vpd,s_spd,s_kpd,s_richness,
     r_failure,r_range_size,
     r_n_pd,r_n_nnd,r_n_mpd,r_n_vpd,r_n_spd,r_n_kpd,r_n_richness,
     r_e_pd,r_e_nnd,r_e_mpd,r_e_vpd,r_e_spd,r_e_kpd,r_e_richness,
     r_ne_pd,r_ne_nnd,r_ne_mpd,r_ne_vpd,r_ne_spd,r_ne_kpd,r_ne_richness,comm_r,comm_s,
     source_cells,recipient_cells,
     recipient_e_phydist,recipient_e_species,recipient_n_phydist,recipient_n_species,recipient_ne_phydist,recipient_ne_species,
     source_phydist,source_species)
    
  #Print output of completeness to watch for frozen process
  cat(round(i/nrow(unique(introductions[,c("Species","Site")]))*100),"percent of species complete for", p, "of",length(trees)," of trees\n"  )
    
  
  
}#i loop
}#p loop

#saveRDS(object = phy_out,file = "data/phy_metrics_output.rds")
#phy_out <- readRDS("data/phy_metrics_output.rds")


#Now, compile code by taking means across replicates
sp_x_site <- unique(phy_out[,c('species_i',"site_i")])

mean_phy_out <- NULL

for( i in 1:nrow(sp_x_site)){
species <- sp_x_site[i,]['species_i']  
site <- sp_x_site[i,]['site_i']  
data_i <- phy_out[which(phy_out[,"species_i"]==species & phy_out[,"site_i"]==site),]  
mean_phy_out<-rbind(mean_phy_out,c(species,site,apply(X = data_i[,3:33],MARGIN = 2,FUN = function(x){mean(as.numeric(x))})))
}

rm(sp_x_site,species,site,data_i,i)
saveRDS(object = mean_phy_out,file = "data/mean_phy_metrics_output.rds")

#V temp code V##################

    

meanout<-readRDS("data/mean_phy_metrics_output.rds")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    
    
