rm(list=ls())

# Download data from https://drive.google.com/drive/folders/0B7hBcCIRGNHXTk8teURpSUthUk0

# Example for the arguments
# work_dir= "/Data/microarray/" or "/Data/rnaseq/"
# design_file= "microarray_design.csv" or "rnaseq_design.csv"
# output_folder= a folder name as you wish
# microarray= if the data orgin is microarray set TRUE otherwise FALSE

########### utils ###########

source("../retreat_utils.R")

########### for microarray ########

work_dir<-"/home/cankut/Desktop/other_projects/mlpm_retreat/Analysis_Retreat_challenge_cankut/microarray/"
files<-list.files(work_dir,pattern=".csv")
idx_design_file<-grep("design",files)
files<-files[-idx_design_file]

acc_list<-list()


for(data in files){
  acc<-do.SVM(work_dir, data, design_file="microarray_design.csv", verbose=T, output_folder="result_array", microarray=T)
  acc_list[[data]]<-acc
}
array_res<-do.call("rbind",acc_list)


########### for RNAseq ###########

work_dir<-"/home/cankut/Desktop/other_projects/mlpm_retreat/Analysis_Retreat_challenge_cankut/rnaseq/"
files<-list.files(work_dir,pattern=".csv")
idx_design_file<-grep("design",files)
files<-files[-idx_design_file]

acc_list<-list()

for(data in files){
  acc<-do.SVM(work_dir, data, design_file="rnaseq_design.csv", verbose=T, output_folder="result_rnaseq", microarray=F)
  acc_list[[data]]<-acc
}
rnaseq_res<-do.call("rbind",acc_list)


################################
all_res<-cbind(array_res,rnaseq_res)
colnames(all_res)<-c("microarray","rnaseq")
rownames(all_res)<-gsub(".csv","",gsub("microarray_","",rownames(all_res)))

write.table(all_res,"../all_res.tsv",,quote = F,sep = "\t",col.names = T,row.names = T)

################ END ##########

