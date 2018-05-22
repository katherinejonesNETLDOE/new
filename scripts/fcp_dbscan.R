install.packages("fpc")
install.packages("dbscan")
install.packages("factoextra")
install.packages("mapproj", dependencies=TRUE)

library(fpc)
library(dbscan)
library(factoextra)
library(mapproj)

setwd("R:\\D_Users\\68_Jones\\Projects\\SIMPA\\SIMPA_analysis")

##example x y data
# data("multishapes", package = "factoextra")
# df <- multishapes[, 1:2]

##read in my own wellbore (x,y) data
df<-read.csv("DBSCAN\\okwells_GT.csv")

#get xy columns
lat<-data.frame(df[ ,which(colnames(df)=="Surface_Latitude")])
long<-data.frame(df[ ,which(colnames(df)=="Surface_Longitude")])

#project x,y data to dataframe

#get first column of lat and long df's
df.xy<-as.data.frame(cbind(x=lat[,1],y=long[,1]))
#project both to project of choice - sometimes last parameter doesn't matter
df.xy.proj <- mapproject(df.xy$x, df.xy$y, "rectangular",0)
x.proj<-df.xy.proj[["x"]]
y.proj<-df.xy.proj[["y"]]
df.xy.proj<-as.data.frame (cbind(x.proj,y.proj))

#check to see if it's properly reading all points and projection
plot(df.xy.proj)

# Compute DBSCAN using fpc package
library("fpc")
set.seed(123)
db <- fpc::dbscan(df.xy.proj, eps = 0.00001, MinPts = 5)
# Plot DBSCAN results
library("factoextra")

fviz_cluster(db, data = df.xy.proj, stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE,
             geom = "point", ggtheme = theme_classic())

#generates plot to determine most appropriate eps value (i.e. look for the 'knee' in this plot)
dbscan::kNNdistplot(df.xy.proj, k =  100)
abline(h = 0.0002, lty = 10, col="red")
