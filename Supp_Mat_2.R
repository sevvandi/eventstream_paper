#  THIS FILE CONTAINS EXPERIMENTS USING REAL APPLICATION 1 DATA
# -----------------------------------------------------------------------------------
# Experiment 2.1 - App 1 data - SAVEC and Logistic Regression
# Experiment 2.2 - App 1 data - DAVEC
# Results in Table 3 
# -----------------------------------------------------------------------------------

library("eventstream")
library("dma")
library("abind")
library("AtmRay")
library("pROC")
library("reshape2")
library("ggplot2")

# -----------------------------------------------------------------------------------
# Experiment 2.1 - App 1 data - SAVEC and Logistic Regression
# -----------------------------------------------------------------------------------

# Load Application 1 data and class labels
data("real_stream")
data("real_details")

# Extracting features
# This is done in the supervised setting. Class labels of events are stored in real_details.
# Making vis=TRUE plots the window data and the extracted events side by side.

feat_all <- extract_event_ftrs(real_stream, supervised=TRUE, details=real_details, win_size=40, step_size=10, thres=0.95, folder="None", vis=FALSE, tt=10)

# Extract time, other features and class labels from feat_all
class_col <- dim(feat_all)[2]
tt <- feat_all[,"width",]
Y <- feat_all[, class_col, ]

# Remove unnecessary variables, cluster_id, centroid_x, centroid_y,  and class
X0 <- feat_all[, -c(1, 7, 8, class_col), ]  
X0[is.na(X0)] <- 0
rm_inds <- which(apply(tt,1,function(x) sum(is.na(x))>0 ))
if(length(rm_inds) >0 ){
  X0 <- X0[-rm_inds, , ]
  tt <- tt[-rm_inds, ]
  Y <- Y[-rm_inds, ]
}

all_data <- rbind.data.frame( X0[,,1], X0[,,2], X0[,,3], X0[,,4])
all_data_sc <- all_data 

ll <- dim(all_data_sc)[1]/4
X <- abind::abind( all_data_sc[1:ll, ], all_data_sc[(ll+1):(2*ll), ], all_data_sc[(2*ll+1):(3*ll), ], all_data_sc[(3*ll+1):(4*ll), ], along=3 )


# 4-fold Cross validation
nn <- 4
folds <- cut(seq(1,dim(X)[1]),breaks=nn,labels=FALSE)

# Initialize different accuracy-metric matrices
ppv_td_stages <- matrix(0,nrow=nn, ncol=4)
npv_td_stages <- matrix(0,nrow=nn, ncol=4)
roc_td_stages <- matrix(0,nrow=nn, ncol=4)

ppv_glm_stages <- matrix(0,nrow=nn, ncol=4)
npv_glm_stages <- matrix(0,nrow=nn, ncol=4)
roc_glm_stages <- matrix(0,nrow=nn, ncol=4)

for(i in 1:nn){
  
  # test data
  testinds <- which(folds==i,arr.ind=TRUE)
  test_data_1 <- X[testinds,,]
  test_data <- rbind.data.frame( test_data_1[,,1],test_data_1[,,2],test_data_1[,,3], test_data_1[,,4] )
  test_labels <- c(Y[testinds,1], Y[testinds,2], Y[testinds,3], Y[testinds,4] )
  test_t <- c(tt[testinds,1], tt[testinds,2], tt[testinds,3], tt[testinds,4] )
  
  # train data
  train_data_1 <- X[-testinds,,]
  train_data <- rbind.data.frame( train_data_1[,,1],train_data_1[,,2],train_data_1[,,3], train_data_1[,,4] )
  train_labels <- c(Y[-testinds,1], Y[-testinds,2], Y[-testinds,3], Y[-testinds,4] )
  train_t <- c(tt[-testinds,1], tt[-testinds,2], tt[-testinds,3], tt[-testinds,4] )
  
  # train model
  model_td <- eventstream::td_logistic(train_t, train_data, train_labels, quad = TRUE, interact=TRUE,  lambda=0.005)
  
  # test model
  pr <- predict_tdl(model_td ,test_t, test_data)
  test_Y <- test_labels
  
  # ppv at each stage
  q_1 <- length(pr)/4
  ppv_td_stages[i,1] <- sum( (pr[1:q_1]==1) & (test_Y[1:q_1]==1) )/sum(test_Y[1:q_1], na.rm=TRUE)
  ppv_td_stages[i,2] <- sum( (pr[(q_1+1):(2*q_1)]==1) & (test_Y[(q_1+1):(2*q_1)]==1) )/sum(test_Y[(q_1+1):(2*q_1)], na.rm=TRUE)
  
  ppv_td_stages[i,3] <- sum( (pr[(2*q_1+1):(3*q_1)]==1) & (test_Y[(2*q_1+1):(3*q_1)]==1) )/sum(test_Y[(2*q_1+1):(3*q_1)], na.rm=TRUE)
  
  ppv_td_stages[i,4] <- sum( (pr[(3*q_1+1):(4*q_1)]==1) & (test_Y[(3*q_1+1):(4*q_1)]==1) )/sum(test_Y[(3*q_1+1):(4*q_1)], na.rm=TRUE)
  
  
  # npv at each stage
  npv_td_stages[i,1] <- sum( (pr[1:q_1]==0) & (test_Y[1:q_1]==0) )/( q_1 - sum(test_Y[1:q_1], na.rm=TRUE) )
  npv_td_stages[i,2] <- sum( (pr[(q_1+1):(2*q_1)]==0) & (test_Y[(q_1+1):(2*q_1)]==0) )/(q_1 - sum(test_Y[(q_1+1):(2*q_1)], na.rm=TRUE) )
  
  npv_td_stages[i,3] <- sum( (pr[(2*q_1+1):(3*q_1)]==0) & (test_Y[(2*q_1+1):(3*q_1)]==0) )/(q_1 - sum(test_Y[(2*q_1+1):(3*q_1)], na.rm=TRUE) )
  
  npv_td_stages[i,4] <- sum( (pr[(3*q_1+1):(4*q_1)]==0) & (test_Y[(3*q_1+1):(4*q_1)]==0) )/(q_1 - sum(test_Y[(3*q_1+1):(4*q_1)], na.rm=TRUE) )
  
  # ROC AUC at each stage
  roc_obj <- roc(test_Y[1:q_1], pr[1:q_1])
  roc_td_stages[i,1] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(q_1+1):(2*q_1)], pr[(q_1+1):(2*q_1)])
  roc_td_stages[i,2] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(2*q_1+1):(3*q_1)], pr[(2*q_1+1):(3*q_1)])
  roc_td_stages[i,3] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(3*q_1+1):(4*q_1)], pr[(3*q_1+1):(4*q_1)])
  roc_td_stages[i,4] <- auc(roc_obj)
  
  
  train_data <- as.matrix(train_data)
  rownames(train_data) <- NULL
  test_data <- as.matrix(test_data)
  rownames(test_data) <- NULL
  
  # Train and test GLM
  xglm <- rbind(train_data, test_data)
  yglm <- c(train_labels, test_labels)
  xyglm <- cbind.data.frame(xglm, yglm)
  num <- dim(train_data)[1]
  en <- dim(xyglm)[1]
  ccol <- dim(xyglm)[2]
  
  
  model_glm <- glm(as.factor(yglm)~ ., data=xyglm[1:num,], family="binomial")
  pr_glm <- predict(model_glm , newdata = xyglm[(num+1):en,-ccol], type="response" )
  pred_glm <- rep(0,length(pr_glm))
  pred_glm[pr_glm > 0.5] <- 1
  pred_glm[pr_glm <= 0.5] <- 0
  
  pr <- pred_glm
  q_1 <- length(pr)/4
  test_Y <- test_labels
  
  # ppv at each stage
  q_1 <- length(pr)/4
  ppv_glm_stages[i,1] <- sum( (pr[1:q_1]==1) & (test_Y[1:q_1]==1) )/sum(test_Y[1:q_1], na.rm=TRUE)
  ppv_glm_stages[i,2] <- sum( (pr[(q_1+1):(2*q_1)]==1) & (test_Y[(q_1+1):(2*q_1)]==1) )/sum(test_Y[(q_1+1):(2*q_1)], na.rm=TRUE)
  
  ppv_glm_stages[i,3] <- sum( (pr[(2*q_1+1):(3*q_1)]==1) & (test_Y[(2*q_1+1):(3*q_1)]==1) )/sum(test_Y[(2*q_1+1):(3*q_1)], na.rm=TRUE)
  
  ppv_glm_stages[i,4] <- sum( (pr[(3*q_1+1):(4*q_1)]==1) & (test_Y[(3*q_1+1):(4*q_1)]==1) )/sum(test_Y[(3*q_1+1):(4*q_1)], na.rm=TRUE)
  
  
  # npv at each stage
  q_1 <- length(pr)/4
  npv_glm_stages[i,1] <- sum( (pr[1:q_1]==0) & (test_Y[1:q_1]==0) )/( q_1 - sum(test_Y[1:q_1], na.rm=TRUE) )
  npv_glm_stages[i,2] <- sum( (pr[(q_1+1):(2*q_1)]==0) & (test_Y[(q_1+1):(2*q_1)]==0) )/(q_1 - sum(test_Y[(q_1+1):(2*q_1)], na.rm=TRUE) )
  
  npv_glm_stages[i,3] <- sum( (pr[(2*q_1+1):(3*q_1)]==0) & (test_Y[(2*q_1+1):(3*q_1)]==0) )/(q_1 - sum(test_Y[(2*q_1+1):(3*q_1)], na.rm=TRUE) )
  
  npv_glm_stages[i,4] <- sum( (pr[(3*q_1+1):(4*q_1)]==0) & (test_Y[(3*q_1+1):(4*q_1)]==0) )/(q_1 - sum(test_Y[(3*q_1+1):(4*q_1)], na.rm=TRUE) )
  
  # roc at each stage
  roc_obj <- roc(test_Y[1:q_1], pr[1:q_1])
  roc_glm_stages[i,1] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(q_1+1):(2*q_1)], pr[(q_1+1):(2*q_1)])
  roc_glm_stages[i,2] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(2*q_1+1):(3*q_1)], pr[(2*q_1+1):(3*q_1)])
  roc_glm_stages[i,3] <- auc(roc_obj)
  roc_obj <- roc(test_Y[(3*q_1+1):(4*q_1)], pr[(3*q_1+1):(4*q_1)])
  roc_glm_stages[i,4] <- auc(roc_obj)
  
}




# -----------------------------------------------------------------------------------
# Experiment 2.2 - App 1 data - DAVEC
# -----------------------------------------------------------------------------------

# Load Application 1 data and class labels
data("real_stream")
data("real_details")

# For 4 - fold cross validation
nn <- 4

feat_all_0 <- extract_event_ftrs(real_stream, supervised=TRUE, details=real_details, win_size=40, step_size=10, thres=0.95, folder="None", vis=FALSE, tt=10)
class.col <- dim(feat_all)[2]


ppv_dma_stages <- matrix(0,ncol=4, nrow=nn)
npv_dma_stages <- matrix(0,ncol=4, nrow=nn)
roc_dma_stages <- matrix(0,ncol=4, nrow=nn)

colnames(roc_dma_stages) <- colnames(npv_dma_stages) <- colnames(ppv_dma_stages) <- c("dma_T1", "dma_T2", "dma_T3", "dma_T4")

for(kk in 1:nn){
  
  # Getting a quarter of data to the end of the array for testing
  mm <- floor(dim(feat_all_0)[1]/nn)
  blk <- ((kk-1)*mm+1):(kk*mm)
  feat_all <- abind(feat_all_0[-blk,,], feat_all_0[blk,,], along=1)
  
  
  Y <- feat_all[, class.col, ]
  X0 <- feat_all[, -c(1, 7, 8, class.col), ]  
  
  X1.1 <- scale(X0[,,1], center=TRUE, scale=TRUE)
  X1.2 <- scale(X0[,,2], center=TRUE, scale=TRUE)
  X1.3 <- scale(X0[,,3], center=TRUE, scale=TRUE)
  X1.4 <- scale(X0[,,4], center=TRUE, scale=TRUE)
  
  X1  <- abind(X1.1,X1.2,X1.3,X1.4, along=3)
  
  # mmat1 is a model matrix needed for dma package
  mmat1 <- matrix(  rep(1, (dim(X1)[2])),  nrow=1, byrow=TRUE)
  tt <- floor(3*dim(X1)[1]/4) 
  en <- dim(X1)[1]
  
  mod1 <- logdma.init(X1[1:tt,,1],Y[1:tt,1], mmat1)

  B <- X1[1:tt, , 2]
  mod2 <- logdma.init(B,Y[1:tt,2], mmat1)

  B <- X1[1:tt, , 3]
  mod3 <- logdma.init(B,Y[1:tt,3], mmat1)

  B <- X1[1:tt, , 4]
  mod4 <- logdma.init(B,Y[1:tt,4], mmat1)
  
  ## Prediction
  preds <- array(0, dim =c(en,4))
  hhh <- (tt+1):en
  preds[hhh,1] <- logdma.predict(mod1,X1[hhh, , 1])

  X2 <- X1[hhh, , 2]
  preds[hhh,2] <- logdma.predict(mod2,X2)

  X3 <- X1[hhh, , 3]
  preds[hhh,3] <- logdma.predict(mod3,X3)

  X4 <- X1[hhh, , 4]
  preds[hhh,4] <- logdma.predict(mod4,X4)

  preds.onezero.1 <- array(0, dim=c(en-tt, 4))
  preds.onezero.1 <-   ifelse(preds[(tt+1):en,]>0.5,1,0)
  
 
  
  # ppv values for each stage
  ppv_dma_stages[kk,1] <- sum( ( preds.onezero.1[,1]==1 ) & ( Y[(tt+1):en,1] ==1 ))/sum( preds.onezero.1[,1]==1 )
  ppv_dma_stages[kk,2] <- sum( ( preds.onezero.1[,2]==1 ) & ( Y[(tt+1):en,2] ==1 ))/sum( preds.onezero.1[,2]==1 )
  ppv_dma_stages[kk,3] <- sum( ( preds.onezero.1[,3]==1 ) & ( Y[(tt+1):en,3] ==1 ))/sum( preds.onezero.1[,3]==1 )
  ppv_dma_stages[kk,4] <- sum( ( preds.onezero.1[,4]==1 ) & ( Y[(tt+1):en,4] ==1 ))/sum( preds.onezero.1[,4]==1 )
  
  
  
  # npv values for each stage
  npv_dma_stages[kk,1] <- sum( ( preds.onezero.1[,1]==0 ) & ( Y[(tt+1):en,1] ==0 ))/sum( preds.onezero.1[,1]==0 )
  npv_dma_stages[kk,2] <- sum( ( preds.onezero.1[,2]==0 ) & ( Y[(tt+1):en,2] ==0 ))/sum( preds.onezero.1[,2]==0 )
  npv_dma_stages[kk,3] <- sum( ( preds.onezero.1[,3]==0 ) & ( Y[(tt+1):en,3] ==0 ))/sum( preds.onezero.1[,3]==0 )
  npv_dma_stages[kk,4] <- sum( ( preds.onezero.1[,4]==0 ) & ( Y[(tt+1):en,4] ==0 ))/sum( preds.onezero.1[,4]==0 )
  
  
  
  # roc values for each stage
  roc_dma_stages[kk,1] <- auc(roc(Y[(tt+1):en,1], preds.onezero.1[,1] ))
  roc_dma_stages[kk,2] <- auc(roc(Y[(tt+1):en,2], preds.onezero.1[,2] ))
  roc_dma_stages[kk,3] <- auc(roc(Y[(tt+1):en,3], preds.onezero.1[,3] ))
  roc_dma_stages[kk,4] <- auc(roc(Y[(tt+1):en,4], preds.onezero.1[,4] ))
  
  
  
}

# -----------------------------------------------------------------------------------
# Results in Table 3 
# -----------------------------------------------------------------------------------


# ---- Averages 
apply(ppv_td_stages,2,mean)
apply(ppv_dma_stages,2,mean)
apply(ppv_glm_stages,2,mean)

apply(npv_td_stages,2,mean)
apply(npv_dma_stages,2,mean)
apply(npv_glm_stages,2,mean)

apply(roc_td_stages,2,mean)
apply(roc_dma_stages,2,mean)
apply(roc_glm_stages,2,mean)


# ---- Standard deviations
apply(ppv_td_stages,2,sd)
apply(ppv_dma_stages,2,sd)
apply(ppv_glm_stages,2,sd)

apply(npv_td_stages,2,sd)
apply(npv_dma_stages,2,sd)
apply(npv_glm_stages,2,sd)

apply(roc_td_stages,2,sd)
apply(roc_dma_stages,2,sd)
apply(roc_glm_stages,2,sd)

