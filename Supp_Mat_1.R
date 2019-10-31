#  THIS FILE CONTAINS EXPERIMENTS USING SYNTHETIC DATA
# -----------------------------------------------------------------------------------
# Experiment 1 - Synthetic data - CC-Log, 1-Log and n-Log
# -----------------------------------------------------------------------------------

library("abind")
library("AtmRay")
library("eventstream")
library("pROC")
library("ggplot2")

# -----------------------------------------------------------------------------------
# Experiment 1 - Synthetic data - CC-Log, 1-Log and n-Log
# -----------------------------------------------------------------------------------

nnn <- 5
accuracy_td_stages <- matrix(0, nrow=nnn, ncol=4)
accuracy_glm_stages <- matrix(0, nrow=nnn, ncol=4)
four_glm_accuracy_stages <- matrix(0, nrow=nnn, ncol=4)



for(i in 1:nnn){
  
  # Generate data
  # sd corresponds to the seed
  # Each iteration
  out <- gen_stream(8, sd=i, vis=FALSE)
  feat_all <- extract_event_ftrs(out$data, supervised=TRUE, details=out$details, thres=0.95, step_size = 8, folder="None", vis=FALSE, tt=8)
  class_col <- dim(feat_all)[2]
  
  # Extract time, other features and class labels from feat_all
  tt <- feat_all[,"width",]
  Y <- feat_all[, class_col, ]
  
  # Remove unnecessary variables, cluster_id, centroid_x, centroid_y, width=age_of_event and class
  X0 <- feat_all[,-c(1, 7, 8, 14:class_col) , ] 
  X0[is.na(X0)] <- 0
  rm_inds <- which(apply(tt,1,function(x) sum(is.na(x))>0 ))
  if(length(rm_inds)>0){
    X0 <- X0[-rm_inds, , ]
    tt <- tt[-rm_inds, ]
    Y <- Y[-rm_inds, ]
  }
  
  # change format to suit model
  train_t <- c(tt[,1],tt[,2], tt[,3],tt[,4])
  Y_tr <- c(Y[,1], Y[,2], Y[,3], Y[,4])
  all_data <- rbind.data.frame( X0[,,1], X0[,,2], X0[,,3], X0[,,4])
  
  
  
  # train model
  model_td <- td_logistic(train_t, all_data, Y_tr, quad=TRUE, interact=FALSE, lambda=0.05 ) # 
  
  # generate test data
  test_out <- gen_stream(2, sd=(i+10))
  
  feat_test <- extract_event_ftrs(test_out$data, supervised=TRUE, details=test_out$details, thres=0.95, step_size = 8,save=FALSE, folder="None", vis=FALSE, tt=8)
  
  # Extract time, other features and class labels from feat_all
  tt2 <- feat_test[,"width",]
  Y2 <- feat_test[, class_col, ]
  Xt <- feat_test[, -c(1, 7, 8,  14:class_col), ]    
  
  # change format to suit model
  test_data <- rbind.data.frame(Xt[,,1], Xt[,,2], Xt[,,3], Xt[,,4])
  test_t <- c(tt2[,1],tt2[,2],tt2[,3],tt2[,4])
  test_Y <- c(Y2[,1],Y2[,2],Y2[,3],Y2[,4])
  
  # predict using model
  pr <- predict_tdl(model_td ,test_t, test_data, probs = TRUE)
  prr <- pr
  pr[prr>=0.5] <- 1
  pr[prr< 0.5] <- 0
  
  # accuracy at each stage
  q_1 <- length(pr)/4
  accuracy_td_stages[i,1] <- sum(pr[1:q_1]==test_Y[1:q_1])/q_1
  accuracy_td_stages[i,2] <- sum(pr[(q_1+1):(2*q_1)]==test_Y[(q_1+1):(2*q_1)])/q_1
  accuracy_td_stages[i,3] <- sum(pr[(2*q_1+1):(3*q_1)]==test_Y[(2*q_1+1):(3*q_1)])/q_1
  accuracy_td_stages[i,4] <- sum(pr[(3*q_1+1):(4*q_1)]==test_Y[(3*q_1+1):(4*q_1)])/q_1
  

  # For glm
  quadpara = TRUE
  dat <- rbind.data.frame(all_data,test_data)
  if(quadpara){
    dat1 <- cbind.data.frame(dat, dat^2)
    colnames(dat1) <- c( colnames(dat), paste(colnames(dat), "_SQ", sep="") )
    dat <- dat1
  }
  Y <- c(Y_tr, test_Y)
  ll <- length(Y_tr)
  en <- length(Y)
  df <- cbind.data.frame(dat,Y)
  xyglm_tr <- df[1:ll,]
  
  mod_glm <- glm(Y~., data=df[1:ll,], family="binomial")
  preds <- predict(mod_glm, newdata = df[(ll+1):en, ])
  prds <- rep(0, length(preds))
  prds[preds>=0.5] <- 1
  
  
  
  accuracy_glm_stages[i,1] <- sum(prds[1:q_1]==test_Y[1:q_1])/q_1
  accuracy_glm_stages[i,2] <- sum(prds[(q_1+1):(2*q_1)]==test_Y[(q_1+1):(2*q_1)])/q_1
  accuracy_glm_stages[i,3] <- sum(prds[(2*q_1+1):(3*q_1)]==test_Y[(2*q_1+1):(3*q_1)])/q_1
  accuracy_glm_stages[i,4] <- sum(prds[(3*q_1+1):(4*q_1)]==test_Y[(3*q_1+1):(4*q_1)])/q_1
  
  

  # TRAIN AND TEST FOUR SEPARATE GLMS
  # USE THE SAME TRAIN AND TEST DATA AS BEFORE
  q1_tr <- dim(all_data)[1]/4
  age_1_tr <- 1:q1_tr
  age_2_tr <- (q1_tr+1):(2*q1_tr)
  age_3_tr <- (2*q1_tr+1):(3*q1_tr)
  age_4_tr <- (3*q1_tr +1):4*q1_tr
  
  model_glm_t1 <- glm(Y ~ ., data=xyglm_tr, subset = age_1_tr, family="binomial") 
  model_glm_t2 <- glm(Y ~ ., data=xyglm_tr, subset = age_2_tr, family="binomial") 
  model_glm_t3 <- glm(Y ~ ., data=xyglm_tr, subset = age_3_tr, family="binomial") 
  model_glm_t4 <- glm(Y ~ ., data=xyglm_tr, subset = age_4_tr,  family="binomial")
  
  q1_ts <- dim(test_data)[1]/4
  age_1_ts <- 1:q1_ts
  age_2_ts <- (q1_ts+1):(2*q1_ts)
  age_3_ts <- (2*q1_ts+1):(3*q1_ts)
  age_4_ts <- (3*q1_ts+1):(4*q1_ts)
  
  xyglm <- df
  num <- ll
  ccol <- dim(df)[2]
  pr_glm_t1 <- predict(model_glm_t1, newdata = xyglm[(num+age_1_ts), -c(ccol)], type="response")
  pr_glm_t2 <- predict(model_glm_t2, newdata = xyglm[(num+age_2_ts), -c(ccol)], type="response" )
  pr_glm_t3 <- predict(model_glm_t3, newdata = xyglm[(num+age_3_ts), -c(ccol)], type="response" )
  pr_glm_t4 <- predict(model_glm_t4, newdata = xyglm[(num+age_4_ts), -c(ccol)], type="response" )
  pr_glm_t1r <- pr_glm_t1
  pr_glm_t2r <- pr_glm_t2
  pr_glm_t3r <- pr_glm_t3
  pr_glm_t4r <- pr_glm_t4
  
  pr_glm_t1 <- ifelse(pr_glm_t1>0.5,1,0)
  pr_glm_t2 <- ifelse(pr_glm_t2>0.5,1,0)
  pr_glm_t3 <- ifelse(pr_glm_t3>0.5,1,0)
  pr_glm_t4 <- ifelse(pr_glm_t4>0.5,1,0)
  
  four_glm_accuracy_stages[i, 1] <- sum(pr_glm_t1==test_Y[age_1_ts])/length(age_1_ts)
  four_glm_accuracy_stages[i, 2] <- sum(pr_glm_t2==test_Y[age_2_ts])/length(age_2_ts)
  four_glm_accuracy_stages[i, 3] <- sum(pr_glm_t3==test_Y[age_3_ts])/length(age_3_ts)
  four_glm_accuracy_stages[i, 4] <- sum(pr_glm_t4==test_Y[age_4_ts])/length(age_4_ts)
  
}

# accuracy_td_stages
apply(accuracy_td_stages,2,mean)
apply(accuracy_td_stages,2,sd)

#accuracy_glm_stages
apply(accuracy_glm_stages,2,mean)
apply(accuracy_glm_stages,2,sd)

#four_glm_accuracy_stages
apply(four_glm_accuracy_stages,2, mean)
apply(four_glm_accuracy_stages,2, sd)



all_accs_2 <- cbind.data.frame(as.vector(accuracy_td_stages), as.vector(accuracy_glm_stages), as.vector(four_glm_accuracy_stages) )
colnames(all_accs_2) <- c("CC-Log", "1-Log", "n-Log")

accs <- -1*as.matrix(all_accs_2)
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.95)
order_methods <- names(out1$means)
accs <- accs[ ,order_methods]
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.9, main="Nemenyi test ranks for synthetic data")
out1$fpval

