# ----------------------------------------------------------------
# TASK 1 - Kulldorff Scan Statistic Comparison - Fibre optic data
# TASK 2 - Kulldorff Scan Statistic Comparison - Synthetic data
# TASK 3 - Kulldorf Scan Statistic Comparison - NO2 dataset
# ----------------------------------------------------------------

library("eventstream")
library("changepoint")
library("ggplot2")
library("SpatialEpi")
library("dplyr")
library("tibble")


# IMPORTANT NOTE
# THE FUNCTION SpatialEpi::kulldorff TAKES A REALLY LONG TIME TO RUN, AND OFTEN RStudio HAD TO BE RESTARTED TO PERFORM THE NEXT ITERATION. 
# SO, THE OUTPUT OF SpatialEpi::kulldorff ARE PROVIDED AS SEPARATE FILES, WHICH WILL BE USED IN THIS SECTION.
# PLEASE DOWNLOAD THESE FILES FROM https://monash.figshare.com/articles/RData_files/10007978
# AND CHANGE FOLDER TO THE FOLDER HAVING ALL kuldorff_*.RData FILES
folder <- ""  # SET FOLDER TO APPROPRIATE LOCAL FOLDER CONTAINING ALL kuldorff_*.RData FILES

# ----------------------------------------------------------------
# TASK 1 - Kulldorff Scan Statistic Comparison - Fibre optic data
# ----------------------------------------------------------------

st <- 1
en <- 40
for(i in 1:10){
  if(en==379){
    en1 <- 380
  }else{
    en1 <- en
  }
  fname <- paste("kuldorff_real_stream_", st, "_", en1, ".RData", sep="")
  load(paste(folder, fname, sep=""))
  pos <- regexec(".RData",fname)
  kuldorff_scan_name <- substring(fname, 1, (pos[[1]][1]-1))
  kuldorff_scan <- eval(parse(text = kuldorff_scan_name))
  
  dat <- real_stream[st:en,]
  
  # CPDBEE EVENTS
  out <- get_clusters(dat, vis=FALSE, rolling=FALSE)
  dat <- t(dat)
  dat.x <- 1:dim(dat)[2]
  dat.y <- 1:dim(dat)[1]
  mesh.xy <- AtmRay::meshgrid(dat.x,dat.y)
  xyz.dat <- cbind(as.vector(mesh.xy$x), as.vector(mesh.xy$y), as.vector(as.matrix(dat)) )
  xyz.dat <- as.data.frame(xyz.dat)
  colnames(xyz.dat) <- c("Time", "Location", "Value")
  
  events <- as_tibble( cbind.data.frame( out$data[ ,1:2], as.vector(rep("Event", dim(out$data)[1] )) ) )
  colnames(events)[3] <- "State"
  tt <- full_join(as_tibble(xyz.dat),events)
  tt <- tt %>% mutate(State = factor(State, levels = c("Background", "Event")))
  tt[is.na(tt[ ,4]), 4] <- "Background"

  
  # KULDORFF'S SCAN STATISTIC CLUSTERS
  cluster <- kuldorff_scan$most.likely.cluster$location.IDs.included
  cluster_sec <- kuldorff_scan$secondary.clusters[[1]]$location.IDs.included
  
  State <- rep("Background", dim(xyz.dat)[1])
  State[cluster] <- "Event"
  xyz.dat2 <- cbind.data.frame(xyz.dat[ ,1:2], State)
  xyz.ori <- cbind.data.frame(xyz.dat, as.vector( rep("Original", dim(xyz.dat)[1]) ))
  colnames(xyz.ori) <- c(colnames(xyz.dat), "Method")
  xyz.clus <- cbind.data.frame(tt[ ,c(1,2,4)],  as.vector( rep("CPDBEE", dim(tt)[1]) ) )
  colnames(xyz.clus) <- colnames(xyz.ori)
  xyz.scanstat <- cbind.data.frame(xyz.dat2,  as.vector( rep("Scan Statistic", dim(tt)[1]) ) )
  colnames(xyz.scanstat) <- colnames(xyz.ori)
  
  if(i==10){
    rat <- 0.15
  }else{
    rat <- 0.3
  }
  xyz.ori2 <- xyz.ori
  xyz.ori2[xyz.ori$Value >40000, 3] <- "Event"
  xyz.ori2[xyz.ori$Value <= 40000, 3] <- "Background"
  df1 <- rbind.data.frame(xyz.ori2 , xyz.clus, xyz.scanstat) 
  
  filename <- paste("Event_Comparison_", st, "_", en, sep="")
  print(filename)
  print( ggplot(df1, aes(Time, Location)) + geom_raster(aes(fill=Value)) +  facet_grid(.~Method)  +   scale_fill_manual(values=c("grey90", "black")) + coord_fixed(ratio=rat) + theme_bw() )
 Sys.sleep(1)
  
  st <- en + 1
  en <- min(en+40, 379)
  
}



# --------------------------------------------------------------
# TASK 2 - Kulldorff Scan Statistic Comparison - Synthetic data
# --------------------------------------------------------------

out <- gen_stream(n=1,sd=15)
dat_ori <- out$data
dat_ori <- dat_ori - min(dat_ori)

st <- 1
en <- 50
for(i in 1:6){

  fname <- paste("kuldorff_synthetic_", st, "_", en, ".RData", sep="")
  load(paste(folder, fname, sep=""))
  pos <- regexec(".RData",fname)
  kuldorff_scan_name <- substring(fname, 1, (pos[[1]][1]-1))
  kuldorff_scan <- eval(parse(text = kuldorff_scan_name))
  
  dat <- dat_ori[st:en,]
  
  # CPDBEE EVENTS
  out <- get_clusters(dat, vis=FALSE, rolling=FALSE)
  dat <- t(dat)
  dat.x <- 1:dim(dat)[2]
  dat.y <- 1:dim(dat)[1]
  mesh.xy <- AtmRay::meshgrid(dat.x,dat.y)
  xyz.dat <- cbind(as.vector(mesh.xy$x), as.vector(mesh.xy$y), as.vector(as.matrix(dat)) )
  xyz.dat <- as.data.frame(xyz.dat)
  colnames(xyz.dat) <- c("Time", "Location", "Value")
  
  events <- as_tibble( cbind.data.frame( out$data[ ,1:2], as.vector(rep("Event", dim(out$data)[1] )) ) )
  colnames(events)[3] <- "State"
  tt <- full_join(as_tibble(xyz.dat),events)
  tt <- tt %>% mutate(State = factor(State, levels = c("Background", "Event")))
  tt[is.na(tt[ ,4]), 4] <- "Background"
  
  
  # KULDORFF'S SCAN STATISTIC CLUSTERS
  cluster <- kuldorff_scan$most.likely.cluster$location.IDs.included
  cluster_sec <- kuldorff_scan$secondary.clusters[[1]]$location.IDs.included
  
  State <- rep("Background", dim(xyz.dat)[1])
  State[cluster] <- "Event"
  xyz.dat2 <- cbind.data.frame(xyz.dat[ ,1:2], State)
  xyz.ori <- cbind.data.frame(xyz.dat, as.vector( rep("Original", dim(xyz.dat)[1]) ))
  colnames(xyz.ori) <- c(colnames(xyz.dat), "Method")
  xyz.clus <- cbind.data.frame(tt[ ,c(1,2,4)],  as.vector( rep("CPDBEE", dim(tt)[1]) ) )
  colnames(xyz.clus) <- colnames(xyz.ori)
  xyz.scanstat <- cbind.data.frame(xyz.dat2,  as.vector( rep("Scan Statistic", dim(tt)[1]) ) )
  colnames(xyz.scanstat) <- colnames(xyz.ori)
  
   xyz.ori2 <- xyz.ori
  xyz.ori2[xyz.ori$Value >10, 3] <- "Event"
  xyz.ori2[xyz.ori$Value <= 10, 3] <- "Background"
  df1 <- rbind.data.frame(xyz.ori2 , xyz.clus, xyz.scanstat) 
  
  filename <- paste("Event_Comparison_synththetic", st, "_", en, sep="")
  print(filename)
  print( ggplot(df1, aes(Time, Location)) + geom_raster(aes(fill=Value)) +  facet_grid(.~Method)  +   scale_fill_manual(values=c("grey90", "black")) + coord_fixed(ratio=0.5) + theme_bw() )
  Sys.sleep(1)
  
  st <- en + 1
  en <- en + 50
  
}




# --------------------------------------------------------------
# TASK 3 - NO2 dataset - Kulldorf Scan Statistic Comparison
# --------------------------------------------------------------

dat_ori <- NO2_2018[1, , ]


st <- 1
en <- 90
for(i in 1:4){
  
  fname <- paste("kuldorff_no2_", st, "_", en, ".RData", sep="")
  load(paste(folder, fname, sep=""))
  pos <- regexec(".RData",fname)
  kuldorff_scan_name <- substring(fname, 1, (pos[[1]][1]-1))
  kuldorff_scan <- eval(parse(text = kuldorff_scan_name))
  
  dat <- dat_ori[ ,st:en]
  
  # CPDBEE EVENTS
  out <- get_clusters(dat, vis=FALSE, rolling=FALSE)
  dat <- t(dat)
  dat.x <- 1:dim(dat)[2]
  dat.y <- 1:dim(dat)[1]
  mesh.xy <- AtmRay::meshgrid(dat.x,dat.y)
  xyz.dat <- cbind(as.vector(mesh.xy$x), as.vector(mesh.xy$y), as.vector(as.matrix(dat)) )
  xyz.dat <- as.data.frame(xyz.dat)
  colnames(xyz.dat) <- c("Time", "Location", "Value")
  
  events <- as_tibble( cbind.data.frame( out$data[ ,1:2], as.vector(rep("Event", dim(out$data)[1] )) ) )
  colnames(events)[3] <- "State"
  tt <- full_join(as_tibble(xyz.dat),events)
  tt <- tt %>% mutate(State = factor(State, levels = c("Background", "Event")))
  tt[is.na(tt[ ,4]), 4] <- "Background"
  
  
  # KULDORFF'S SCAN STATISTIC CLUSTERS
  cluster <- kuldorff_scan$most.likely.cluster$location.IDs.included
  cluster_sec <- kuldorff_scan$secondary.clusters[[1]]$location.IDs.included
  
  State <- rep("Background", dim(xyz.dat)[1])
  State[cluster] <- "Event"
  xyz.dat2 <- cbind.data.frame(xyz.dat[ ,1:2], State)
  xyz.ori <- cbind.data.frame(xyz.dat, as.vector( rep("Original", dim(xyz.dat)[1]) ))
  colnames(xyz.ori) <- c(colnames(xyz.dat), "Method")
  xyz.clus <- cbind.data.frame(tt[ ,c(1,2,4)],  as.vector( rep("CPDBEE", dim(tt)[1]) ) )
  colnames(xyz.clus) <- colnames(xyz.ori)
  xyz.scanstat <- cbind.data.frame(xyz.dat2,  as.vector( rep("Scan Statistic", dim(tt)[1]) ) )
  colnames(xyz.scanstat) <- colnames(xyz.ori)
  
  xyz.ori2 <- xyz.ori
  xyz.ori2[xyz.ori$Value >100, 3] <- "Event"
  xyz.ori2[xyz.ori$Value <= 100, 3] <- "Background"
  df1 <- rbind.data.frame(xyz.ori2 , xyz.clus, xyz.scanstat) 
  df1[ ,1] <- 180 -df1[ ,1] 
  colnames(df1) <- c("Y", "X", "Value", "Method")
  
  filename <- paste("Event_Comparison_NO2_", st, "_", en, sep="")
  print(filename)
  print( ggplot(df1, aes(X, Y)) + geom_raster(aes(fill=Value)) +  facet_grid(.~Method)  +   scale_fill_manual(values=c("grey90", "black")) + coord_fixed(ratio=1) + theme_bw() )
  Sys.sleep(1)
  
  st <- en + 1
  en <- en + 90
  
}

