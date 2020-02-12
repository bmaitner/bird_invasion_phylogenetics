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
  #recipient community (n, e, ne)
  #metric (pd, nnd, mpd, vpd, spd, kpd)
  #range size
  #species richness native
  #species richness intro

#load in data

  #load in introduction data
    introductions <- read.csv("data/invasion_data.csv",stringsAsFactors = F)
    
    #cut australia data (too large an area)
    introductions<-introductions[which(introductions$Site!="A"),]
    
    #one bad name ( hybrid, drop from analysis)
    introductions<-introductions[grep(pattern = "_x_",x = introductions$Species,invert = T),]
    
  #load in translator
    translator <- readRDS("data/name_translator.rds")

  #load occurrences
    occs <- readRDS("data/cea_occurrences.rds")
    occs$species<-as.character(occs$species)
    occs$species <- gsub(pattern = " ",replacement = "_",x = occs$species)
  
  #load in recipient region polygons
    library(maptools)
    library(sf)
    
    gadm <- sf::read_sf("C:/Users/Brian/Desktop/Range maps/gadm28/gadm28.shp") #read in gadm
    gadm <- gadm[which((gadm$NAME_0 =="United States" & gadm$NAME_1 %in% c("Hawaii", "Florida")) |
                         gadm$NAME_0 == "New Zealand"),]#cut the gadm down to just relevant spots to save space
    gadm <- sf:::as_Spatial(gadm)#convert to spatial
    
    #reproject to crs of template raster
    template<-raster()
    template<-projectRaster(from = template,crs = crs("+proj=cea +units=km"),res = 110)
    gadm <- spTransform(x = gadm,CRSobj = template@crs)
    plot(gadm)
    gadm@data$name01<-paste(gadm@data$NAME_0,gadm@data$NAME_1)
    gadm@data$name01
    
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
    
