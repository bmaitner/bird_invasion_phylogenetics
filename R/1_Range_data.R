# Code to swap out ranges with a "skinny format" dataframe
# Input: Birdlife Maps (individual, since smaller size, easier bug fixing)
# Output: Dataframe with 4 columns: Species, Cell,percent of cell occupied, origin


##Origin
#1	Native
#The species is/was native to the area
#2	Reintroduced
#The species is/was reintroduced through either direct or indirect human activity
#3	Introduced
#The species is/was introduced outside of its historical distribution range through either direct or indirect human activity.
#4	Vagrant
#The species is/was recorded once or sporadically, but it is known not to be native to the area.
#5	Origin Uncertain
#The species’ provenance in an area is not known (it may be native, reintroduced or introduced).
#6	Assisted Colonisation
#Species subject to intentional movement and release outside its native range to reduce the extinction risk of the taxon.


library(sp)
library(maptools)
library(sf)
library(raster)

cont<-read_sf("C:/Users/Brian/Desktop/Range maps/landmass/continents_equal_area.shp")
cont<-as_Spatial(cont)


shapefiles <- list.files(path = "data/All_Shapefiles_2_5_2016/All_Shapefiles",pattern = ".shp",full.names = T)
template<-raster()
template<-projectRaster(from = template,crs = crs("+proj=cea +units=km"),res = 110)

#Prep empty output file
skinny_occurrences<-NULL


for(i in 1:length(shapefiles)){
cat(round(i/length(shapefiles),digits = 2)*100," percent done")

sp_i<-NULL
try(sp_i <- sf::st_read(dsn = shapefiles[i]))
if(is.null(sp_i)){next}

#cut out empty geometries(these prevent conversion to spatial)
sp_i<-sp_i[which(sapply(X = sp_i$geometry,FUN = length)>0),]

#convert to spatial(needed for rasterization, which is needed to get pct cover of cells)
sp_i<-as_Spatial(sp_i)

#if there is no data associated with the geometries, skip the file
if(ncol(sp_i@data)==0){next}

#Excluding origin = 4 since this is unlikely to be a major part of the range or community
sp_i <- sp_i[which(sp_i@data$ORIGIN != 4),]


#rasterize by ORIGIN(this contains native/introduced data)
raster_i<-rasterize(x = spTransform(x = sp_i,CRSobj = template@crs),y = template,field="ORIGIN")
raster_i_cover<-rasterize(x = spTransform(x = sp_i,CRSobj = template@crs),y = template,field="ORIGIN",getCover=T)

#plot(raster_i)
#plot(spTransform(x = cont,CRSobj = template@crs),add=T)

#only try to add data if there is any to add
if(length(which(getValues(raster_i_cover>0)))>0){

skinny_occurrences <- rbind(skinny_occurrences,cbind(as.character(sp_i@data$SCINAME[1]),which(getValues(raster_i_cover>0)),raster_i[which(getValues(raster_i_cover>0))] ,raster_i_cover[which(getValues(raster_i_cover>0))]))

}


}
rm(raster_i,raster_i_cover,sp_i,shapefiles,i)

skinny_occurrences<-as.data.frame(skinny_occurrences)
colnames(skinny_occurrences)<-c("species","cell","origin","percent_cover")

write.csv(x = skinny_occurrences,file = "data/cea_occurrences.csv",row.names = F)
