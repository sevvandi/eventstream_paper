#  THIS FILE CONTAINS OTHER FIGURES IN THE PAPER
# -----------------------------------------------------------------------------
# TASK 1 - Real world dataset - Figure 2 in paper
# TASK 2 - evenstream events and clusters - Figure 3 in paper
# TASK 3 - Figure with splines - Figure 4 in paper
# -----------------------------------------------------------------------------

library("eventstream")
library("ggplot2")
library("reshape2")
library("AtmRay")

# -----------------------------------------------------------------------------
# TASK 1 - Real world dataset - Figure 2 in paper
# -----------------------------------------------------------------------------
data("real_stream")
par(mfrow=c(1,1))
dat1 <- real_stream[1:100,]
dat <- as.data.frame(t(dat1))
dat.x <- 1:dim(dat)[2]
dat.y <- 1:dim(dat)[1]
mesh.xy <- AtmRay::meshgrid(dat.x,dat.y)
xyz.dat <- cbind(as.vector(mesh.xy$x), as.vector(mesh.xy$y), as.vector(as.matrix(dat)) )
xyz.dat <- as.data.frame(xyz.dat)
colnames(xyz.dat) <- c("Time", "Location", "Value")
ggplot(xyz.dat, aes(Time, Location)) + geom_raster(aes(fill=Value)) +   scale_fill_gradientn(colours=topo.colors(12)) + theme_bw()


ftrs <- extract_event_ftrs(dat1,win_size=100, step_size = 20)
df <- ftrs[c(1,3),11,]
colnames(df) <- c("T1", "T2", "T3", "T4")
df_m <- melt(df)
colnames(df_m) <- c("Event", "Age", "Slope")
df_m$Event <- as.factor(df_m$Event)
levels(df_m$Event) <- c("A", "B")
ggplot(df_m, aes(Age, Slope)) + geom_point(aes(color=Event)) + geom_path(aes(group=Event, color=Event)) + xlab("Event Age") + ylab("Slope of fitted line")  + theme_bw()


# -----------------------------------------------------------------------------
# TASK 2 - evenstream events and clusters - Figure 3 in paper
# -----------------------------------------------------------------------------
out <- gen_stream(n=1,sd=15)
dat <- out$data
image(1:nrow(dat), 1:ncol(dat),  as.matrix(dat),col=topo.colors(12), xlab="Time", ylab="Location", main="Window")

events <- get_clusters(dat, flag="N", thres=0.98, vis=FALSE)
xyz.high.xy <- events$data[,c(1,2)]
res <- events$clusters
plot(xyz.high.xy[res$cluster!=0,c(2,1)], col=res$cluster[res$cluster!=0]+1L, pch=19, xlim=c(0,dim(dat)[1]), ylim=c(0,dim(dat)[2]), xlab="Time", ylab="Location", main="Events"  )


# -----------------------------------------------------------------------------
# TASK 3 - Figure with splines - Figure 4 in paper
# -----------------------------------------------------------------------------
data("real_stream")
dat <- real_stream[1:20,]
normal.stats.splines <- spline_stats(dat)

par(mfrow=c(1,1))
image(1:ncol(dat), 1:nrow(dat), as.matrix(t(dat)),col=topo.colors(12), ylab="Time", xlab="Location")
plot(normal.stats.splines[[1]], type="l", xlab="Location", ylab="Median pixel value")
plot(normal.stats.splines[[2]], type="l", xlab="Location", ylab="IQR pixel value")


