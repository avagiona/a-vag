#REQUIREMENTS
library("igraph")
library("dplyr")
library("ggplot2")
library("tidyr")
library("writexl")
library("plyr")
library("tibble")
#set pathway
setwd("~/R-4.2.1/ML_PTMs")
options(scipen=999)

#ML analysis
library("randomForest")
require("caTools")
library('caret')
library("ROSE")
library("pROC")

final_df_ML_feat<-read.csv("final_df_ML_feat.csv",fileEncoding = "UTF-8-BOM", header=TRUE)
df_ml<-final_df_ML_feat%>%
  dplyr::select(c(4,5,8:20)) %>%
  tidyr::drop_na()

summary(df_ml)
names(df_ml) <- make.names(names(df_ml))
df_ml$enzymatic<-as.factor(df_ml$enzymatic)
table(df_ml$enzymatic)

#split in train and test set
set.seed(32658)
sample_df_ml = sample.split(df_ml$enzymatic, SplitRatio = .70)
train_df_ml = subset(df_ml, sample_df_ml == TRUE)
test_df_ml  = subset(df_ml, sample_df_ml == FALSE)
dim(train_df_ml)
dim(test_df_ml)

#hyperparameters
#train the model without spliting train and test
set.seed(32658)
down_inside<- train(enzymatic ~ ., data = train_df_ml,
                    method = 'rf',metric='Accuracy',
                    trControl=trainControl(method="repeatedcv",
                                           repeats = 10,
                                           number=5,sampling='down',
                                           returnResamp='final',savePredictions =  TRUE))
down_inside

#plot accuracies of the models from cross validation
df_accu_models_cross<-(down_inside$resample)
split<-data.frame(do.call("rbind", strsplit(as.character(df_accu_models_cross$Resample), ".", fixed = TRUE)))
df_accu_models_cross<-cbind(df_accu_models_cross,split)

accu_mod<-ggplot(df_accu_models_cross, aes(x=X1, y=Accuracy))+  
  geom_line(aes(colour = X2,group = X2),linewidth=0.8)+ 
  xlab("Number of folds")+theme_bw()+ labs(color='Repeats')+
  theme(legend.position = "none") 
accu_mod

#apply the model
set.seed(32658)
pred_down_inside= predict(down_inside,newdata=test_df_ml[-5],type='prob')#remove the last column you want to predict
pred_down_inside_class= predict(down_inside,newdata=test_df_ml[-5],type='raw')
cm_down_inside = table(test_df_ml[,5], as.numeric(pred_down_inside[[2]]))
cm_down_inside_df<-as.data.frame(cm_down_inside)
cm<-caret::confusionMatrix(pred_down_inside_class,as.factor(test_df_ml$enzymatic))

#feature importance
i_scores_caret <- varImp(down_inside,scale = FALSE)
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
  #ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))+
  theme_bw()+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),color="darkgrey")
roc_plot

plot_model<-ggarrange(accu_mod,roc_plot,
                      ncol = 2,labels = c("A", "B"), nrow = 1)
plot_model
