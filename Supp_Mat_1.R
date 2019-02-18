#  THIS FILE CONTAINS EXPERIMENTS USING SYNTHETIC DATA
# -----------------------------------------------------------------------------------
# Experiment 1.1 - Synthetic data - SAVEC & glm
# Experiment 1.2 - Synthetic data - DAVEC
# Output 1.1 - Synthetic data associated graphs
# -----------------------------------------------------------------------------------


library("dma")
library("abind")
library("AtmRay")
library("eventstream")
library("pROC")
library("reshape2")
library("ggplot2")

# -----------------------------------------------------------------------------------
# Experiment 1.1 - Synthetic data - SAVEC & glm
# -----------------------------------------------------------------------------------

accuracy_td_stages <- matrix(0, nrow=5, ncol=4)
accuracy_glm_stages <- matrix(0, nrow=5, ncol=4)

for(i in 1:5){

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
  model_td <- td_logistic(train_t, all_data, Y_tr, interact=TRUE, lambda=0.05)

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
  pr <- predict_tdl(model_td ,test_t, test_data)

  # accuracy at each stage
  q_1 <- length(pr)/4
  accuracy_td_stages[i,1] <- sum(pr[1:q_1]==test_Y[1:q_1])/q_1
  accuracy_td_stages[i,2] <- sum(pr[(q_1+1):(2*q_1)]==test_Y[(q_1+1):(2*q_1)])/q_1
  accuracy_td_stages[i,3] <- sum(pr[(2*q_1+1):(3*q_1)]==test_Y[(2*q_1+1):(3*q_1)])/q_1
  accuracy_td_stages[i,4] <- sum(pr[(3*q_1+1):(4*q_1)]==test_Y[(3*q_1+1):(4*q_1)])/q_1

  # For glm
  dat <- rbind.data.frame(all_data,test_data)
  Y <- c(Y_tr, test_Y)
  ll <- length(Y_tr)
  en <- length(Y)
  df <- cbind.data.frame(dat,Y)

  mod_glm <- glm(Y~., data=df[1:ll,], family="binomial")
  preds <- predict(mod_glm, newdata = df[(ll+1):en, ])
  prds <- rep(0, length(preds))
  prds[preds>=0.5] <- 1



  accuracy_glm_stages[i,1] <- sum(prds[1:q_1]==test_Y[1:q_1])/q_1
  accuracy_glm_stages[i,2] <- sum(prds[(q_1+1):(2*q_1)]==test_Y[(q_1+1):(2*q_1)])/q_1
  accuracy_glm_stages[i,3] <- sum(prds[(2*q_1+1):(3*q_1)]==test_Y[(2*q_1+1):(3*q_1)])/q_1
  accuracy_glm_stages[i,4] <- sum(prds[(3*q_1+1):(4*q_1)]==test_Y[(3*q_1+1):(4*q_1)])/q_1
}


accuracy_td_stages
apply(accuracy_td_stages,2,mean)
apply(accuracy_td_stages,2,sd)


accuracy_glm_stages
apply(accuracy_glm_stages,2,mean)
apply(accuracy_glm_stages,2,sd)



# -----------------------------------------------------------------------------------
# Experiment 1.2 - Synthetic data - DAVEC
# -----------------------------------------------------------------------------------
accuracy_dma_stages <- matrix(0, ncol=4, nrow=5)
colnames(accuracy_dma_stages) <- c("dma_T1", "dma_T2", "dma_T3", "dma_T4")

for(j in 1:5){
  out <- gen_stream(10, sd=i, vis=FALSE)
  feat_all <- extract_event_ftrs(out$data, supervised=TRUE, details=out$details, thres=0.95, step_size = 8, folder="None", vis=FALSE, tt=8)

  class_col <- dim(feat_all)[2]
  Y <- feat_all[, class_col, ]
  X0 <- feat_all[, -c(1, 7, 8, 14:class_col), ]  

  X1.1 <- scale(X0[,,1], center=TRUE, scale=TRUE)
  X1.2 <- scale(X0[,,2], center=TRUE, scale=TRUE)
  X1.3 <- scale(X0[,,3], center=TRUE, scale=TRUE)
  X1.4 <- scale(X0[,,4], center=TRUE, scale=TRUE)

  X1  <- abind(X1.1,X1.2,X1.3,X1.4, along=3)

  mmat1 <- matrix(  rep(1, (dim(X1)[2])),  nrow=1, byrow=TRUE)
  tt <- floor(4*dim(X1)[1]/5)
  en <- dim(X1)[1]

  mod1 <- logdma.init(X1[1:tt,,1],Y[1:tt,1], mmat1)
  pred.mod.1 <- logdma.predict(mod1, X1[1:tt,,1])

  B <- cbind(X1[1:tt, , 2], pred.mod.1[,1])
  mmat <- matrix(  rep(1, (dim(X1)[2] +1 )),  nrow=1, byrow=TRUE)
  mod2 <- logdma.init(B,Y[1:tt,2], mmat1)
  pred.mod.2 <- logdma.predict(mod2, B[1:tt,])

  B <- cbind(X1[1:tt, , 3], pred.mod.2[,1])
  mod3 <- logdma.init(B,Y[1:tt,3], mmat1)
  pred.mod.3 <- logdma.predict(mod3, B[1:tt,])

  B <- cbind(X1[1:tt, , 4], pred.mod.3[,1])
  mod4 <- logdma.init(B,Y[1:tt,4], mmat1)

  # Prediction and model updates
  preds <- array(0, dim =c(en,4))
  preds[1:tt,1] <- t(mod1$yhatmodel)
  preds[1:tt,2] <- t(mod2$yhatmodel)
  preds[1:tt,3] <- t(mod3$yhatmodel)
  preds[1:tt,4] <- t(mod4$yhatmodel)
  for(i in (tt+1):en){
    ## Prediction
    preds[i,1] <- logdma.predict(mod1,X1[i, , 1])
    X2 <- c(X1[i, , 2], preds[i,1] )
    preds[i,2] <- logdma.predict(mod2,X2)
    X3 <- c(X1[i, , 3], preds[i,2] )
    preds[i,3] <- logdma.predict(mod3,X3)
    X4 <- c(X1[i, , 4], preds[i,3] )
    preds[i,4] <- logdma.predict(mod4,X4)
    # Update Models
    mod1 <- logdma.update(mod1, X1[i, , 1], Y[i,1], lambda = 0.99, autotune=FALSE)
    mod2 <- logdma.update(mod2, X2, Y[i,2], lambda = 0.99, autotune=FALSE)
    mod3 <- logdma.update(mod3, X3, Y[i,3], lambda = 0.99, autotune=FALSE)
    mod4 <- logdma.update(mod4, X4, Y[i,4], lambda = 0.99, autotune=FALSE)
  }

  preds.onezero.1 <- array(0, dim=c(en-tt, 4))
  preds.onezero.1 <-   ifelse(preds[(tt+1):en,]>0.5,1,0)

  accuracy_dma_stages[j,1] <- sum(preds.onezero.1[,1]==Y[(tt+1):en,1])/length(Y[(tt+1):en,1])

  accuracy_dma_stages[j,2] <- sum(preds.onezero.1[,2]==Y[(tt+1):en,2])/length(Y[(tt+1):en,2])

  accuracy_dma_stages[j,3] <- sum(preds.onezero.1[,3]==Y[(tt+1):en,3])/length(Y[(tt+1):en,3])

  accuracy_dma_stages[j,4] <- sum(preds.onezero.1[,4]==Y[(tt+1):en,4])/length(Y[(tt+1):en,4])
}
accuracy_dma_stages
apply(accuracy_dma_stages,2,mean)
apply(accuracy_dma_stages,2,sd)


# -----------------------------------------------------------------------------------
# Output 1.1 - Synthetic data associated graphs
# -----------------------------------------------------------------------------------
accuracy_td_stages1 <- accuracy_td_stages

accuracy_td_stages <- cbind.data.frame(1:5, accuracy_td_stages)
colnames(accuracy_td_stages) <- c("Fold","T1", "T2", "T3", "T4")
accuracy_td_stages <- as.data.frame(accuracy_td_stages)
df1 <- melt(id="Fold",accuracy_td_stages)
df1 <- cbind.data.frame(df1, rep("SAVEC", dim(df1)[1]))
df1
colnames(df1) <- c("Fold", "Age", "Accuracy", "Classifier")
df1$Accuracy <- df1$Accuracy*100

accuracy_glm_stages1 <- accuracy_glm_stages

accuracy_glm_stages <- cbind.data.frame(1:5, accuracy_glm_stages)
colnames(accuracy_glm_stages) <- c("Fold","T1", "T2", "T3", "T4")
accuracy_glm_stages <- as.data.frame(accuracy_glm_stages)
df2 <- melt(id="Fold",accuracy_glm_stages)
df2 <- cbind.data.frame(df2, rep("Logistic Reg.", dim(df2)[1]))
df2
colnames(df2) <- c("Fold", "Age", "Accuracy", "Classifier")
df2$Accuracy <- df2$Accuracy*100

accuracy_dma_stages1 <- accuracy_dma_stages

accuracy_dma_stages <- cbind.data.frame(1:5, accuracy_dma_stages)
colnames(accuracy_dma_stages) <- c("Fold","T1", "T2", "T3", "T4")
accuracy_dma_stages <- as.data.frame(accuracy_dma_stages)
df3 <- reshape2::melt(id="Fold",accuracy_dma_stages)
df3 <- cbind.data.frame(df3, rep("DAVEC", dim(df3)[1]))
df3
colnames(df3) <- c("Fold", "Age", "Accuracy", "Classifier")
df3$Accuracy <- df3$Accuracy*100

df <- rbind(df1, df3, df2)

pp <- ggplot(df, aes(Age,Accuracy)) + geom_point(aes(color=Classifier))
pp +  geom_line(data=df1, aes(color=Classifier, group=Fold)) + geom_line(data=df2, aes(color=Classifier, group=Fold) ) + geom_line(data=df3, aes(color=Classifier, group=Fold)) + facet_grid(.~Classifier)  +  xlab("Event Age")  + theme_bw()
## Save as 3_Classifiers_1.pdf


pp <- ggplot(df, aes(Age,Accuracy)) + geom_point(aes(color=Classifier))
pp +  geom_line(data=df1, aes(color=Classifier, group=Fold)) + geom_line(data=df2, aes(color=Classifier, group=Fold) ) + geom_line(data=df3, aes(color=Classifier, group=Fold))  + facet_grid(Fold~.)  +  xlab("Event Age")  + theme_bw()
## Save as 3_Classifiers_2.pdf

