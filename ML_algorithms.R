#REQUIREMENTS
library("igraph")
library("dplyr")
library("ggplot2")
library("tidyr")
library("enrichR")
library("plyr")
#set pathway
setwd("~/R-4.2.1/predict_complexes_ML")
options(scipen=999)

#import_files
stats_pdbs_net<-read.csv("pdbs_net.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
positive_cases<-read.csv("final_positive_cases.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
#all_possible_triangles<-readRDS("onen_triangles_net_new_ids.RData")
edges_net_hippie<-read.csv("edges_new_ids.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
nodes_net_hippie<-read.csv("features_nodes_new_ids.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
ML_data<-read.csv("randomized_v1v2.csv", fileEncoding = "UTF-8-BOM", header=TRUE)

#statistics/plots
hist_pdb_num_prot_net<-ggplot(stats_pdbs_net, aes(x = number_proteins)) +
  geom_histogram(binwidth = 1, fill = "#F8766D", color = "black", alpha = 1) +
  labs(x = "Number of proteins/PDB",
       y = "Frequency")+theme_bw()+
  scale_x_continuous(breaks = seq(3,52, by = 2))
hist_pdb_num_prot_net

hist_pdb_num_ot_pdb<-ggplot(stats_pdbs_net, aes(x = number_ot)) +
  geom_histogram(binwidth = 14, fill = "#F8766D", color = "black", alpha = 7) +
  labs(x = "Number of open triangles",
       y = "Frequency")+theme_bw()
hist_pdb_num_ot_pdb

pdb_stats<-ggarrange(hist_pdb_num_prot_net,hist_pdb_num_ot_pdb,
                        ncol = 2,labels = c("A", "B"), nrow = 1)
pdb_stats

#ML analysis
library("randomForest")
require(caTools)
library('caret')
library("ROSE")
library("pROC")

df_ml<-ML_data%>%
  dplyr::select(1,5:46) %>%
  tidyr::drop_na()

summary(df_ml)
names(df_ml) <- make.names(names(df_ml))
df_ml$positive<-as.factor(df_ml$positive)
table(df_ml$positive)

#split in train and test set
set.seed(452016)
sample_df_ml = sample.split(df_ml$positive, SplitRatio = .70)
train_df_ml = subset(df_ml, sample_df_ml == TRUE)
test_df_ml  = subset(df_ml, sample_df_ml == FALSE)
dim(train_df_ml)
dim(test_df_ml)
table(train_df_ml$positive)

#undersampling
data_balanced_under <- ovun.sample(positive~.,data = train_df_ml,method = "under",N=810)$data
data_balanced_under_replace <- data_balanced_under %>%
  dplyr::mutate(positive = ifelse(positive == 0,"No","Yes"))
test_df_ml_replace<-test_df_ml %>%
  dplyr::mutate(positive = ifelse(positive == 0,"No","Yes"))

#random forest
rf_under<- train(positive ~ ., data = data_balanced_under,
                 method = 'rf',metric="Accuracy",
                 trControl=trainControl(method="repeatedcv",
                                        repeats = 10,
                                        number=5,
                                        returnResamp='final',savePredictions =  TRUE))
rf_under

pred_class_rf= predict(rf_under,newdata=test_df_ml[-1],type='raw')
#pred_class_rf_prob= predict(rf_under,newdata=test_df_ml[-1],type='prob')
cm_rf = as.data.frame(table(test_df_ml[,1], as.numeric(pred_class_rf)))
cm_rf<-caret::confusionMatrix(pred_class_rf,as.factor(test_df_ml$positive),
                              mode="everything")
cm_rf

#SVM
svm_under<- train(positive ~ ., data = data_balanced_under_replace,
                    method = "svmLinear",metric="Accuracy",
                    trControl=trainControl(method="repeatedcv",
                                           repeats = 10,
                                           number=5,returnResamp='final',
                                           savePredictions =  TRUE,classProbs =  TRUE))
svm_under
pred_test_svm_class= predict(svm_under,newdata=test_df_ml_replace[-1],type='raw')
#pred_test_svm_prob= predict(svm_under,newdata=test_df_ml_replace[-1],type='prob')
cm_svm = as.data.frame(table(test_df_ml_replace[,1], as.numeric(pred_test_svm_class)))
cm_svm<-caret::confusionMatrix(pred_test_svm_class,as.factor(test_df_ml_replace$positive),
                               mode="everything")
cm_svm

#logistic regression
lg_under<- train(positive ~ ., data = data_balanced_under,
                 method = 'glmnet',metric="Accuracy",
                 trControl=trainControl(method="repeatedcv",
                                        repeats = 10,
                                        number=5,
                                        returnResamp='final',savePredictions =  TRUE))
lg_under

pred_class_lg= predict(lg_under,newdata=test_df_ml[-1],type='raw')
#pred_prob_lg= predict(lg_under,newdata=test_df_ml[-1],type='prob')
cm_lg = as.data.frame(table(test_df_ml[,1], as.numeric(pred_class_lg)))
cm_lg<-caret::confusionMatrix(pred_class_lg,as.factor(test_df_ml$positive),
                              mode="everything")
cm_lg

#decission_tree
dt_under<- train(positive ~ ., data = data_balanced_under,
                 method = 'rpart',metric="Accuracy",
                 trControl=trainControl(method="repeatedcv",
                                        repeats = 10,
                                        number=5,
                                        returnResamp='final',savePredictions =  TRUE))
dt_under

pred_class_dt= predict(dt_under,newdata=test_df_ml[-1],type='raw')
#pred_prob_dt= predict(dt_under,newdata=test_df_ml[-1],type='prob')
cm_dt = as.data.frame(table(test_df_ml[,1], as.numeric(pred_class_dt)))
cm_dt<-caret::confusionMatrix(pred_class_dt,as.factor(test_df_ml$positive),
                              mode="everything")
cm_dt
#knn
knn_under<- train(positive ~ ., data = data_balanced_under,
                 method = 'knn',metric="Accuracy",
                 trControl=trainControl(method="repeatedcv",
                                        repeats = 10,
                                        number=5,
                                        returnResamp='final',savePredictions =  TRUE))
knn_under

pred_class_knn= predict(knn_under,newdata=test_df_ml[-1],type='raw')
#pred_prob_knn= predict(knn_under,newdata=test_df_ml[-1],type='prob')
cm_knn = as.data.frame(table(test_df_ml[,1], as.numeric(pred_class_knn)))
cm_knn<-caret::confusionMatrix(pred_class_knn,as.factor(test_df_ml$positive),
                              mode="everything")
cm_knn

#apply the model to the dataset
predictions_triplets= predict(rf_under,newdata=df_ml,type='prob')
final_predictions_triplets<-cbind(ML_data,predictions_triplets)

final_predictions<-final_predictions_triplets%>%
  dplyr::select(1:4,48)%>%
  dplyr::rename("score"=5)

cooperative<-final_predictions[final_predictions$score >= 0.9,]
competitive<-final_predictions[final_predictions$score <= 0.1,]

hist(final_predictions$score)
final_predictions_top_score<-final_predictions[final_predictions$score >= 0.5,]

#feature importance
i_scores_caret <- varImp(rf_under,scale = FALSE)
print(i_scores_caret)
feat_imp<-plot(i_scores_caret)
feat_imp
importance_df<-as.data.frame(i_scores_caret$importance)
importance_df$feature <- rownames(importance_df)
rownames(importance_df) <- c()
importance_df_or<-importance_df[order(-importance_df$Overall),]
library("forcats")
feat_impo<-ggplot(importance_df_or, aes(x=fct_reorder(feature, Overall),y=Overall)) +
  geom_bar(stat="identity",color="black",fill="#00BFC4",alpha=0.7,width = 0.6)+
  theme_bw()+coord_flip()+xlab("")+
  ylab("Feature importance")
feat_impo

#roc curve
library("pROC")
roc_score=roc(test_df_ml[,5], as.numeric(pred_down_inside[[2]]))
auc <- roc_score$auc

roc_plot<-ggroc(roc_score, colour = '#00BFC4', size = 0.8,legacy.axes = T) +
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),color="darkgrey")
roc_plot

plot_model<-ggarrange(accu_mod,roc_plot,
                      ncol = 2,labels = c("A", "B"), nrow = 1)
plot_model

