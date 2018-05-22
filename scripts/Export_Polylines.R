#Create a function for writing out polyline shapefiles 

#define the function to be called elsewhere
export.spatial.lines<-function (begin.coords,end.coords,layer.name){

#create an empty list that is the length of the # lines to be exported (i.e. - read length of beginning file, every
#every beginning has an end)
lns <- vector("list", nrow(begin.coords))

#for the elements in the list representing the number of lines
for (i in seq_along(lns)) {
  #
  lns[[i]] <- Lines(list(Line(rbind(begin.coords[i, ], end.coords[i,]))), ID=as.character(i))
}

#didn't right in projection because having trouble getting rgdal to recognize utm zone 15
sp_lns<-SpatialLines(lns, proj4string = CRS(as.character(NA)))



df <- data.frame(len = sapply(1:length(sp_lns), function(i) gLength(sp_lns[i, ])))
rownames(df) <- sapply(1:length(sp_lns), function(i) sp_lns@lines[[i]]@ID)

sp_lns_df <- SpatialLinesDataFrame(sp_lns, data = df)

writeOGR(sp_lns_df,dsn=paste0(substr(layer.name,1,2),"output_shapefiles"),
         layer=layer.name,driver="ESRI Shapefile",morphToESRI = T)

}