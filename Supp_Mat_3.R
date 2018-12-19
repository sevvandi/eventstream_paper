#  THIS FILE CONTAINS EXPERIMENTS USING USING REAL APPLICATION 2 DATA (NO2 DATA)
# -----------------------------------------------------------------------------------
# Output 3.1 - NO2 data for 2008 and associated map
# Output 3.2 - NO2 event features and events at March 2008
# Output 3.3 - NO2 event features and events at March 2018
# Output 3.4 - Analysis of the biggest cluster
# Output 3.5 - Match clusters in 2008 and 2018 and analyse differences
# -----------------------------------------------------------------------------------

library("raster")
library("maps")
library("eventstream")
library("ggplot2")
library("reshape2")

# -----------------------------------------------------------------------------------
# Output 3.1 - NO2 data for 2008 and associated map
# -----------------------------------------------------------------------------------
# Load 2008 NO2 data
data("NO2_2008")
tt <- NO2_2008[1, , ]
r <- raster(tt,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r)
# SAVE AS NO2_March_2008_No_Bndry.pdf
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS NO2_March_2008_With_Bndry.pdf

# Load 2018 NO2 data
data("NO2_2018")
tt <- NO2_2018[1, , ]
r <- raster(tt,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r)
# SAVE AS NO2_March_2008_No_Bndry.pdf
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS NO2_March_2018_With_Bndry.pdf

# -------------------------------------------------------------------------
# Output 3.2 - NO2 event features and events at March 2008
# -------------------------------------------------------------------------
output <- get_clusters_3d(NO2_2008, flag="N", filename="Nothing",thres=0.97, epsilon = 2, miniPts = 20, vis=FALSE)
cluster.all <- output$clusters
xyz.high <- output$data
all_no2_clusters_march <- xyz.high[xyz.high[,1]==1,-1]
all_cluster_ids_march <- cluster.all$cluster[xyz.high[,1]==1]
cluster_ids_march <- all_cluster_ids_march[all_cluster_ids_march!=0]
no2_clusters_march <- all_no2_clusters_march[all_cluster_ids_march!=0,]

march_map <- matrix(0, nrow=180, ncol=360)
set.seed(123)
new_ids <- sample( length(unique(cluster_ids_march)),length(unique(cluster_ids_march)) )
new_cluster_ids <- cluster_ids_march
for(i in 1:length(unique(cluster_ids_march)) ) {
  new_cluster_ids[ cluster_ids_march== unique(cluster_ids_march)[i]] <- new_ids[i]
}

march_map[no2_clusters_march[,1:2]] <- new_cluster_ids
r <- raster(march_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r, legend=F)
# SAVE AS Clusters_NO2_March_2008_No_Bndry.pdf

map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS Clusters_NO2_March_2008_With_Bndry.pdf
length(unique(cluster_ids_march))

cluster.all <- output$clusters
xyz.high <- output$data
ll <- ceiling(dim(NO2_2008)[1]/5)
normal.stats <- stats_3d(NO2_2008[1:ll,,])
output_2008 <- output

# Compute features
all.basic.features.2008 <- get_features_3d(xyz.high, cluster.all$cluster , normal.stats,dim(arr)[1], 1)


# -------------------------------------------------------------------------
# Output 3.3 - NO2 event features and events at March 2018
# -------------------------------------------------------------------------
# Load 2018 NO2 data
data("NO2_2018")

output <- get_clusters_3d(NO2_2018, flag="N", filename="Nothing",thres=0.97, epsilon = 2, miniPts = 20, vis=FALSE)
cluster.all <- output$clusters
xyz.high <- output$data
all_no2_clusters_march <- xyz.high[xyz.high[,1]==1,-1]
all_cluster_ids_march <- cluster.all$cluster[xyz.high[,1]==1]
cluster_ids_march <- all_cluster_ids_march[all_cluster_ids_march!=0]
no2_clusters_march <- all_no2_clusters_march[all_cluster_ids_march!=0,]

march_map <- matrix(0, nrow=180, ncol=360)
set.seed(123)
new_ids <- sample( length(unique(cluster_ids_march)),length(unique(cluster_ids_march)) )
new_cluster_ids <- cluster_ids_march
for(i in 1:length(unique(cluster_ids_march)) ) {
  new_cluster_ids[ cluster_ids_march== unique(cluster_ids_march)[i]] <- new_ids[i]
}

march_map[no2_clusters_march[,1:2]] <- new_cluster_ids
r <- raster(march_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r, legend=F)
# SAVE AS Clusters_NO2_March_2018_No_Bndry.pdf

map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS Clusters_NO2_March_2018_With_Bndry.pdf
length(unique(cluster_ids_march))

cluster.all <- output$clusters
xyz.high <- output$data
ll <- ceiling(dim(NO2_2018)[1]/5)
normal.stats <- stats_3d(NO2_2018[1:ll,,])
output_2018 <- output


# Compute features
all.basic.features.2018 <- get_features_3d(xyz.high, cluster.all$cluster , normal.stats,dim(arr)[1], 1)


# -------------------------------------------------------------------------
# Output 3.4 - Analysis of the biggest cluster
# -------------------------------------------------------------------------
mm1 <- which.max(all.basic.features.2008[,2,1])
mm2 <- which.max(all.basic.features.2018[,2,1])

all.basic.features.2008[mm1,,]
all.basic.features.2018[mm2,,]

all.basic.features.2008[mm1,9:10,]
all.basic.features.2018[mm2,9:10,]


## Plot this event for 2008
xyz.high <- output_2008$data
high_event <- which( (output_2008$clusters$cluster==all.basic.features.2008[mm1,1,1]) & (xyz.high[,1]==1 ) )
high_event_data<-xyz.high[high_event,]


high_map <- matrix(0, nrow=180, ncol=360)
high_map[high_event_data[,2:3]] <- high_event_data[,4]
r <- raster(high_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r)
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS High_Event_2008.pdf



## Plot this event for 2018
xyz.high <- output_2018$data
high_event <- which( (output_2018$clusters$cluster==all.basic.features.2018[mm2,1,1]) & (xyz.high[,1]==1 ) )
high_event_data<-xyz.high[high_event,]


high_map <- matrix(0, nrow=180, ncol=360)
high_map[high_event_data[,2:3]] <- high_event_data[,4]
r <- raster(high_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r)
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS High_Event_2018.pdf


month <- c("Mar","Apr", "May", "Jun")
month <- factor(month, levels = month.abb)
NO2_2008 <- all.basic.features.2008[mm1,11,]
NO2_2018 <- all.basic.features.2018[mm2,11,]
df <- cbind.data.frame(month,NO2_2008,NO2_2018)
mdf <- melt(df)
ggplot(mdf, aes(month,value)) + geom_point(aes(color =variable) )+ geom_line(aes(group =variable, color = variable))  + ylab("Average NO2") + theme_bw() 
## SAVE AS Highest_Event_NO2_2008_2018.pdf

# -------------------------------------------------------------------------
# Output 3.5 - Match clusters in 2008 and 2018 and analyse differences
# -------------------------------------------------------------------------

matches <- rep(0,dim(all.basic.features.2008)[1])
for(i in 1:dim(all.basic.features.2008)[1]){
  coords_y <- all.basic.features.2008[i,9,1]
  coords_z <- all.basic.features.2008[i,10,1]
  rec <- which( abs(coords_y - all.basic.features.2018[,9,1]) + abs(coords_z - all.basic.features.2018[,10,1]) < 5 )
  if(length(rec)>0){
    matches[i] <- rec
  }
}
matches
clus1 <- which(matches>0)
clus2 <- matches[clus1]
clus1


mean_vals_2008 <- all.basic.features.2008[clus1,11,]
mean_vals_2018 <- all.basic.features.2018[clus2,11,]
yy <- mean_vals_2018 - mean_vals_2008
yy <- cbind.data.frame(as.factor(paste(clus2)), yy)
colnames(yy) <- c("Event","Mar","Apr", "May", "Jun")
meltdf <- melt(yy, id="Event")
ggplot(meltdf,aes(x=variable,y=value,color=Event,group=Event)) + geom_line() + theme_bw() + ylab("NO2 level difference 2018 - 2008") + xlab("Month")
# SAVE AS NO2_Level_Difference_2018_2008.pdf

meltdf[meltdf$Event==15,]
positive_yy <- yy[which(apply(yy[,2:5],1,function(x)sum(x>0))==4),]
evs <- which(meltdf$Event %in% positive_yy$Event)
meltdf[evs,]
ggplot(meltdf[evs,],aes(x=variable,y=value,color=Event,group=Event)) + geom_line() + theme_bw() + ylab("NO2 level increase 2018 - 2008 ") + xlab("Month")
# SAVE AS NO2_Level_Increase_2018_2008.pdf


out2 <- clus2[which(yy[,2]>200)]
all.basic.features.2018[out2,,]

out1 <- clus1[which(yy[,2]>200)]
all.basic.features.2008[out1,,]

## Plot this event
odd_event <- which( (output$clusters$cluster==out2) & (xyz.high[,1]==1 ) )
odd_event_data<-xyz.high[odd_event,]


diff_map <- matrix(0, nrow=180, ncol=360)
diff_map[odd_event_data[,2:3]] <- odd_event_data[,1]
r <- raster(diff_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r, legend=FALSE)
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS Different_Event_2018_2008.pdf


positive <- clus2[which(apply(yy[,2:5],1,function(x)sum(x>0))==4)]
## Plot this event
pos_event <- which( (output$clusters$cluster%in%positive) & (xyz.high[,1]==1 ) )
pos_event_data<-xyz.high[pos_event,]


pos_map <- matrix(0, nrow=180, ncol=360)
pos_map[pos_event_data[,2:3]] <- output$clusters$cluster[pos_event] 
r <- raster(pos_map,xmn=-179.5,xmx=179.5,ymn=-89.5,ymx=89.5,crs="+proj=longlat +datum=WGS84")
plot(r)


break_points <- sort( unique(c (unique(output$clusters$cluster[pos_event])-1, unique(output$clusters$cluster[pos_event])+1)) )
colors <- c("red", "darkgrey", "blue", "red", "darkgreen", "brown")
plot(r, breaks=break_points, col=colors, legend="F")
map("world",add=T, fill=FALSE, col="darkgrey")
# SAVE AS Positive_NO2_Events_2018_2008.pdf
