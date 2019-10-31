library("raster")
library("maps")
library("eventstream")
library("ggplot2")
library("pROC")
library("latex2exp")
library("tsutils")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# CLASSIFICATION OF NO2 CLUSTERS
# ------------------------------------------------------------------------------
data(NO2_2010)
data(NO2_2011)
data(NO2_2012)
data(NO2_2013)
data(NO2_2014)
data(NO2_2015)
data(NO2_2016)
data(NO2_2017)
data(NO2_2018)
data(NO2_2019)
ftrs_2010 <- extract_event_ftrs(NO2_2010, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=TRUE, tt=1, vis=FALSE) 
ftrs_2011 <-  extract_event_ftrs(NO2_2011, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE) 
ftrs_2012 <-  extract_event_ftrs(NO2_2012, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2013 <-  extract_event_ftrs(NO2_2013, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2014 <-  extract_event_ftrs(NO2_2014, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2015 <-  extract_event_ftrs(NO2_2015, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2016 <-  extract_event_ftrs(NO2_2016, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2017 <-  extract_event_ftrs(NO2_2017, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2018 <-  extract_event_ftrs(NO2_2018, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)
ftrs_2019 <-  extract_event_ftrs(NO2_2019, thres=0.97, epsilon = 2, miniPts = 20, win_size=4, step_size=1, rolling=FALSE, tt=1, vis=FALSE)


scc_accuracy_stages <- matrix(0, nrow=10, ncol=4)
four_glm_accuracy_stages <- one_glm_accuracy_stages <- matrix(0, nrow=10, ncol=4)
colnames(four_glm_accuracy_stages) <- colnames(scc_accuracy_stages) <- colnames(one_glm_accuracy_stages) <- c("T1_Accuracy", "T2_Accuracy","T3_Accuracy","T4_Accuracy" )

for(jj in 1:10){
  if(jj==1){
    test_ftrs <- ftrs_2010
    train_ftrs <-  abind::abind( ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==2){
    test_ftrs <- ftrs_2011
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==3){
    test_ftrs <- ftrs_2012
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==4){
    test_ftrs <- ftrs_2013
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==5){
    test_ftrs <- ftrs_2014
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==6){
    test_ftrs <- ftrs_2015
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2016, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==7){
    test_ftrs <- ftrs_2016
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2017, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==8){
    test_ftrs <- ftrs_2017
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2018, ftrs_2019, along=1)
  }else if(jj==9){
    test_ftrs <- ftrs_2018
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2019, along=1)
  }else if(jj==10){
    test_ftrs <- ftrs_2019
    train_ftrs <-  abind::abind( ftrs_2010, ftrs_2011, ftrs_2012, ftrs_2013, ftrs_2014, ftrs_2015, ftrs_2016, ftrs_2017, ftrs_2018, along=1)
  }
  

  train_labs1 <- ifelse( (train_ftrs[ ,11, 4] - train_ftrs[ ,11, 2] )>0, 1, 0)
  train_labs <- train_labs1 
  train_ftrs[ ,17, ] <- train_labs
  

  test_labs1 <- ifelse( ( test_ftrs[ ,11, 4] - test_ftrs[ ,11, 2] )>0, 1, 0)
  test_labs2 <- ifelse( ( (test_ftrs[ ,2, 4]-test_ftrs[ ,2, 3])/test_ftrs[ ,2, 2] )>=1, 1, 0)
  test_labs <-  ifelse( (test_labs1==1)|(test_labs2==1), 1, 0)
  test_ftrs[ ,17, ] <- test_labs
  
  X0 <-  train_ftrs[ ,-c(1,8,9,10,17),] 
  train_t <- c(X0[,2,1], X0[,2,2], X0[,2,3], X0[,2,4])
  Y_tr <- rep(train_labs,4)
  all_data <- rbind.data.frame( X0[,,1], X0[,,2], X0[,,3], X0[,,4])
  
  # Remove NAs
  na_sums <- apply(X0, 1, function(x) sum(is.na(x)))
  rm_rows <- which(na_sums >0 )
  if(length(rm_rows)>0){
    train_t <- train_t[-rm_rows]
    all_data <- all_data[-rm_rows, ]
    Y_tr <- Y_tr[-rm_rows]
  }
  
  # TRAIN CC-Log MODEL
  model_td <- td_logistic(train_t, all_data, Y_tr, quad=TRUE, lambda=0.05) 
  
  
  # TEST CC-Log MODEL
  X0 <-  test_ftrs[ ,-c(1,8,9,10,17),] 
  test_t <- c(X0[,2,1], X0[,2,2], X0[,2,3], X0[,2,4])
  Y_t <- rep(test_labs,4)
  test_data <- rbind.data.frame( X0[,,1], X0[,,2], X0[,,3], X0[,,4])
  
  # Remove NAs
  na_sums <- apply(test_data, 1, function(x) sum(is.na(x)))
  rm_rows <- which(na_sums >0 )
  if(length(rm_rows)>0){
    test_t <- test_t[-rm_rows]
    test_data <- test_data[-rm_rows, ]
    Y_t <- Y_t[-rm_rows]
  }
  
  
  prpr <- predict_tdl(model_td ,test_t, test_data,  probs=TRUE)
  pr2 <- prpr
  pr2[prpr >=0.5] <- 1
  pr2[prpr < 0.5] <- 0
  age_1 <- which(test_data$length==1)
  scc_accuracy_stages[jj,1] <- sum(pr2[age_1]==Y_t[age_1])/length(age_1)
  
  age_2 <- which(test_data$length==2)
  scc_accuracy_stages[jj,2] <- sum(pr2[age_2]==Y_t[age_2])/length(age_2)
  
  age_3 <- which(test_data$length==3)
  scc_accuracy_stages[jj,3] <- sum(pr2[age_3]==Y_t[age_3])/length(age_3)
  
  age_4 <- which(test_data$length==4)
  scc_accuracy_stages[jj,4] <- sum(pr2[age_4]==Y_t[age_4])/length(age_4)
  

  # TRAIN AND TEST ONE GLM
  quadpara <- TRUE
  train_data <- all_data
  test_data <- test_data
  train_labels <- Y_tr
  test_labels <- Y_t
  xglm <- rbind(train_data, test_data)
  yglm <- c(train_labels, test_labels)
  
  if(quadpara){
    dat1 <- cbind.data.frame(xglm, xglm^2)
    colnames(dat1) <- c( colnames(xglm), paste(colnames(xglm), "_SQ", sep="") )
    xglm <- dat1
  }
  
  xyglm <- cbind.data.frame(xglm, yglm)
  num <- dim(train_data)[1]
  xyglm_tr <- xyglm[1:num, ]
  en <- dim(xyglm)[1]
  ccol <- dim(xyglm)[2]
  
  
  model_glm <- glm(yglm ~ ., data=xyglm_tr, family="binomial")
  pr_glm <- predict(model_glm , newdata = xyglm[(num+1):en,-ccol], type="response" )
  pred_glm <- rep(0,length(pr_glm))
  pred_glm[pr_glm > 0.5] <- 1
  pred_glm[pr_glm <= 0.5] <- 0
  pred_glm
  test_labels
  
  age_1 <- which(test_data$length==1)
  one_glm_accuracy_stages[jj, 1] <- sum(pred_glm[age_1]==Y_t[age_1])/length(age_1)
  
  age_2 <- which(test_data$length==2)
  one_glm_accuracy_stages[jj, 2] <-sum(pred_glm[age_2]==Y_t[age_2])/length(age_2)
  
  age_3 <- which(test_data$length==3)
  one_glm_accuracy_stages[jj, 3] <-sum(pred_glm[age_3]==Y_t[age_3])/length(age_3)
  
  age_4 <- which(test_data$length==4)
  one_glm_accuracy_stages[jj, 4] <-sum(pred_glm[age_4]==Y_t[age_4])/length(age_4)
  

  # TRAIN AND TEST FOUR SEPARATE GLMS
  # USE THE SAME TRAIN AND TEST DATA AS BEFORE
  age_1_tr <- which(train_data$length==1)
  age_2_tr <- which(train_data$length==2)
  age_3_tr <- which(train_data$length==3)
  age_4_tr <- which(train_data$length==4)
  
  model_glm_t1 <- glm(yglm ~ ., data=xyglm_tr[ ,-c(2,9,10,11,14,21:23)], subset = age_1_tr, family="binomial") 
  model_glm_t2 <- glm(yglm ~ ., data=xyglm_tr[ ,-c(2,10,11,14,22,23)], subset = age_2_tr, family="binomial", na.action="na.exclude") 
  model_glm_t3 <- glm(yglm ~ ., data=xyglm_tr[ ,-c(2,10,14,22)], subset = age_3_tr, family="binomial", na.action="na.exclude") 
  model_glm_t4 <- glm(yglm ~ ., data=xyglm_tr[ ,-c(2,10,14,22)], subset = age_4_tr,  family="binomial", na.action="na.exclude")
  
  age_1_ts <- which(test_data$length==1)
  age_2_ts <- which(test_data$length==2)
  age_3_ts <- which(test_data$length==3)
  age_4_ts <- which(test_data$length==4)
  
  pr_glm_t1r <- predict(model_glm_t1, newdata = xyglm[(num+age_1_ts), -c(2,9,10,11,14,21:23,ccol)], type="response")
  pr_glm_t2r <- predict(model_glm_t2, newdata = xyglm[(num+age_2_ts), -c(2,10,11,14,22,23, ccol)], type="response" )
  pr_glm_t3r <- predict(model_glm_t3, newdata = xyglm[(num+age_3_ts), -c(2,10,14,22, ccol)], type="response" )
  pr_glm_t4r <- predict(model_glm_t4, newdata = xyglm[(num+age_4_ts), -c(2,10,14,22, ccol)], type="response" )
  
  pr_glm_t1 <- ifelse(pr_glm_t1r>0.5,1,0)
  pr_glm_t2 <- ifelse(pr_glm_t2r>0.5,1,0)
  pr_glm_t3 <- ifelse(pr_glm_t3r>0.5,1,0)
  pr_glm_t4 <- ifelse(pr_glm_t4r>0.5,1,0)
  
  four_glm_accuracy_stages[jj, 1] <- sum(pr_glm_t1==Y_t[age_1])/length(age_1)
  four_glm_accuracy_stages[jj, 2] <- sum(pr_glm_t2==Y_t[age_2])/length(age_2)
  four_glm_accuracy_stages[jj, 3] <- sum(pr_glm_t3==Y_t[age_3])/length(age_3)
  four_glm_accuracy_stages[jj, 4] <- sum(pr_glm_t4==Y_t[age_4])/length(age_4)
  

}
apply(scc_accuracy_stages, 2, mean)
apply(one_glm_accuracy_stages, 2, mean)
apply(four_glm_accuracy_stages, 2, mean)


apply(scc_accuracy_stages, 2, sd)
apply(one_glm_accuracy_stages, 2, sd)
apply(four_glm_accuracy_stages, 2, sd)


all_accs_1 <- cbind.data.frame(as.vector(scc_accuracy_stages), as.vector(one_glm_accuracy_stages), as.vector(four_glm_accuracy_stages) )
colnames(all_accs_1) <- c("CC-Log", "1-Log", "n-Log")

accs <- -1*as.matrix(all_accs_1)
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.95)
order_methods <- names(out1$means)
accs <- accs[ ,order_methods]
out1 <- tsutils::nemenyi(accs, plottype = "matrix", conf.level=0.95, main=TeX("Nemenyi test ranks for $NO_2$ data"))
out1$fpval
friedman.test(as.matrix(all_accs_1))



