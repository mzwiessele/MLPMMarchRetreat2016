do.SVM<-function(work_dir,data,design_file,verbose=F,output_folder,microarray){
  save.folder<-paste0(work_dir,"/",output_folder)
  dir.create(save.folder,recursive=T,showWarnings=F)
  
  library(e1071)  
  rnaseq_design<-read.table(paste0(work_dir,"/",design_file),header=T,sep=",",row.names=1,stringsAsFactors=F)
  #head(rnaseq_design)
  #str(rnaseq_design)
  
  rnaseq<-read.table(paste0(work_dir,"/",data),header=T,sep=",",row.names=1,stringsAsFactors=F)
  #head(rnaseq)
  #class(rnaseq[1,1])
  
  thresh<-quantile(apply(rnaseq,1,var))[2]
  no_variance<-which(apply(rnaseq,1,var)<=thresh)
  if(length(no_variance)>0){
    rnaseq<-rnaseq[-no_variance,]
  }
  
  
  
  if(microarray){  
    idx_train<-which(rnaseq_design$Train.Test=="Training Set")
  }else{
    idx_train<-which(rnaseq_design$Train.Test=="TRAINING SET")
  }
  
  
  ## Train Set
  trainset<-rnaseq[,idx_train]
  trainlabel<-rnaseq_design$Effect[idx_train]
  train_df<-data.frame(t(trainset),"Labels"=trainlabel,stringsAsFactors=T)
  
  # Tes Set
  testset<-rnaseq[,-idx_train]
  testlabel<-rnaseq_design$Effect[-idx_train]
  test_df<-data.frame("Labels"=testlabel,t(testset),stringsAsFactors=T)
  
  
  # tune.out<-tune.svm(Labels~.,data=train_df, kernel='radial', cost=2^seq(-5,5,len=5), gamma=10^seq(-3,6,len=5),tunecontrol=tune.control(cross=5))
  # summary(tune.out)
  # svmmodel<-svm(Labels ~., data=train_df, kernel='radial', type="C-classification",cross=10,gamma=tune.out$best.parameters$gamma,cost=tune.out$best.parameters$cost) # for factor C-class...
  
  svmmodel<-svm(Labels ~., data=train_df, kernel='radial', type="C-classification",cross=10)
  svm_prediction <- predict(svmmodel,test_df[,-1],type="decision")
  
  pred_tab_subpath<-table(pred=svm_prediction, true=test_df$Labels)
  # print(confusionMatrix(pred_tab_subpath))
  # Compute accuracy
  accuracy<-sum(as.character(svm_prediction)==as.character(test_df$Labels))/length(test_df$Labels)
  data_name<-strsplit(data,"[.]")[[1]][1]
  save.image(paste0(save.folder,"/",data_name,".RData"))
  if(verbose){
    cat(accuracy," ",data_name,"\n")
  }
  return(accuracy)
}
