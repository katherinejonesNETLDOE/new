#Create a function for writing out polyline shapefiles 

#define the function to be called elsewhere
export.spatial.points<-function (dataframe.withxy,layer.name){
 
  projection.info<-"+proj=utm +zone=15+datum=WGS84"
  
  spdf <- SpatialPointsDataFrame(dataframe.withxy[,c("x", "y")], dataframe.withxy[,1:8],proj4string = CRS(projection.info))  
   
  writeOGR(spdf,dsn=paste0(substr(layer.name,1,2),"output_clusters"),
           layer=layer.name,driver="ESRI Shapefile",morphToESRI = T)
  
}
  