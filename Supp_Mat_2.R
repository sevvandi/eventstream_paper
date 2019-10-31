#  THIS FILE CONTAINS EXPERIMENTS USING REAL APPLICATION 1 DATA
# -----------------------------------------------------------------------------------
# Experiment 2 - App 1 data - Age-varying events Classifier & glm
# -----------------------------------------------------------------------------------

library("eventstream")
library("abind")
library("AtmRay")
library("pROC")
library("ggplot2")

# -----------------------------------------------------------------------------------
# Experiment 2 - App 1 data - Age-varying events Classifier & glm
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

ppv_four_glms <- matrix(0,nrow=nn, ncol=4)
npv_four_glms <- matrix(0,nrow=nn, ncol=4)
roc_four_glms <- matrix(0,nrow=nn, ncol=4)

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
  model_td <- eventstream::td_logistic(train_t, train_data, train_labels, quad = FALSE,   lambda=0.05)
  
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
  xyglm_tr <- xyglm[1:num, ]
  
  
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
  
  
  # TRAIN AND TEST FOUR SEPARATE GLMS
  # USE THE SAME TRAIN AND TEST DATA AS BEFORE
  q1_tr <- dim(train_data)[1]/4
  age_1_tr <- 1:q1_tr
  age_2_tr <- (q1_tr+1):(2*q1_tr)
  age_3_tr <- (2*q1_tr+1):(3*q1_tr)
  age_4_tr <- (3*q1_tr+1):(4*q1_tr)
  
  model_glm_t1 <- glm(yglm ~ ., data=xyglm_tr, subset = age_1_tr, family="binomial") 
  model_glm_t2 <- glm(yglm ~ ., data=xyglm_tr, subset = age_2_tr, family="binomial") 
  model_glm_t3 <- glm(yglm ~ ., data=xyglm_tr, subset = age_3_tr, family="binomial") 
  model_glm_t4 <- glm(yglm ~ ., data=xyglm_tr, subset = age_4_tr,  family="binomial")
  
  age_1_ts <- 1:q_1
  age_2_ts <- (q_1+1):(2*q_1)
  age_3_ts <- (2*q_1+1):(3*q_1)
  age_4_ts <- (3*q_1+1):(4*q_1)
  
  pr_glm_t1 <- predict(model_glm_t1, newdata = xyglm[(num+age_1_ts), -ccol], type="response")
  pr_glm_t2 <- predict(model_glm_t2, newdata = xyglm[(num+age_2_ts), -ccol], type="response" )
  pr_glm_t3 <- predict(model_glm_t3, newdata = xyglm[(num+age_3_ts), -ccol], type="response" )
  pr_glm_t4 <- predict(model_glm_t4, newdata = xyglm[(num+age_4_ts), -ccol], type="response" )
  
  pr_glm_t1 <- ifelse(pr_glm_t1>0.5,1,0)
  pr_glm_t2 <- ifelse(pr_glm_t2>0.5,1,0)
  pr_glm_t3 <- ifelse(pr_glm_t3>0.5,1,0)
  pr_glm_t4 <- ifelse(pr_glm_t4>0.5,1,0)
  
  # ppv at each stage
  ppv_four_glms[i,1] <- sum( (pr_glm_t1==1) & (test_Y[age_1_ts]==1) )/sum(test_Y[age_1_ts], na.rm=TRUE)
  ppv_four_glms[i,2] <- sum( (pr_glm_t2==1) & (test_Y[age_2_ts]==1) )/sum(test_Y[age_2_ts], na.rm=TRUE)
  ppv_four_glms[i,3] <- sum( (pr_glm_t3==1) & (test_Y[age_3_ts]==1) )/sum(test_Y[age_3_ts], na.rm=TRUE)
  ppv_four_glms[i,4] <- sum( (pr_glm_t4==1) & (test_Y[age_4_ts]==1) )/sum(test_Y[age_4_ts], na.rm=TRUE)
  
  
  # npv at each stage
  npv_four_glms[i,1] <- sum( (pr_glm_t1==0) & (test_Y[age_1_ts]==0) )/( q_1 - sum(test_Y[age_1_ts], na.rm=TRUE) )
  npv_four_glms[i,2] <- sum( (pr_glm_t2==0) & (test_Y[age_2_ts]==0) )/( q_1 - sum(test_Y[age_2_ts], na.rm=TRUE) )
  npv_four_glms[i,3] <- sum( (pr_glm_t3==0) & (test_Y[age_3_ts]==0) )/( q_1 - sum(test_Y[age_3_ts], na.rm=TRUE) )
  npv_four_glms[i,4] <- sum( (pr_glm_t4==0) & (test_Y[age_4_ts]==0) )/( q_1 - sum(test_Y[age_4_ts], na.rm=TRUE) )
  
  
  # ROC AUC at each stage
  roc_obj <- roc(test_Y[age_1_ts], pr_glm_t1)
  roc_four_glms[i,1] <- auc(roc_obj)
  roc_obj <- roc(test_Y[age_2_ts], pr_glm_t2)
  roc_four_glms[i,2] <- auc(roc_obj)
  roc_obj <- roc(test_Y[age_3_ts], pr_glm_t3)
  roc_four_glms[i,3] <- auc(roc_obj)
  roc_obj <- roc(test_Y[age_4_ts], pr_glm_t4)
  roc_four_glms[i,4] <- auc(roc_obj)
  
  
  
}

apply(ppv_td_stages,2,mean)
apply(ppv_td_stages,2,sd)
apply(ppv_glm_stages,2,mean)
apply(ppv_glm_stages,2,sd)
apply(ppv_four_glms,2,mean)
apply(ppv_four_glms,2,sd)


apply(npv_td_stages,2,mean)
apply(npv_td_stages,2,sd)
apply(npv_glm_stages,2,mean)
apply(npv_glm_stages,2,sd)
apply(npv_four_glms,2,mean)
apply(npv_four_glms,2,sd)


apply(roc_td_stages,2,mean)
apply(roc_td_stages,2,sd)
apply(roc_glm_stages,2,mean)
apply(roc_glm_stages,2,sd)
apply(roc_four_glms,2,mean)
apply(roc_four_glms,2,sd)



all_accs_2 <- cbind.data.frame(as.vector(roc_td_stages), as.vector(roc_glm_stages), as.vector(roc_four_glms) )
colnames(all_accs_2) <- c("CC-Log", "1-Log", "n-Log")

accs <- -1*as.matrix(all_accs_2)
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.9)
order_methods <- names(out1$means)
accs <- accs[ ,order_methods]
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.95, main=paste("Nemenyi test ranks using AUC "))
out1$fpval




all_accs_2 <- cbind.data.frame(as.vector(ppv_td_stages), as.vector(ppv_glm_stages), as.vector(ppv_four_glms) )
colnames(all_accs_2) <- c("CC-Log", "1-Log", "n-Log")

accs <- -1*as.matrix(all_accs_2)
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.9)
order_methods <- names(out1$means)
accs <- accs[ ,order_methods]
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.95, main=paste("Nemenyi test ranks using PPV "))
out1$fpval


all_accs_2 <- cbind.data.frame(as.vector(npv_td_stages), as.vector(npv_glm_stages), as.vector(npv_four_glms) )
colnames(all_accs_2) <- c("SC-Log", "1-Log", "n-Log")

accs <- -1*as.matrix(all_accs_2)
out1 <- tsutils::nemenyi(accs, plottype = "none", conf.level=0.9)
out1$fpval



