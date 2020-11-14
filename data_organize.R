library(stringr)
library(tidyverse)


setwd('D:/Mirror/Documents/my study/PITT/biostatistics/class/2069 Statistical Methods for omics/final project')

#sample sheet, file id is the folder name, from file name can know type (htseq.counts/FPKM/FPKM-UQ)
sample_sheet<-read.delim('gdc_sample_sheet.2020-11-04.tsv')
htseq<-sample_sheet[str_detect(sample_sheet[,2],'htseq.counts'),]
FPKM<-sample_sheet[str_detect(sample_sheet[,2],'FPKM.txt'),]
FPKM_UQ<-sample_sheet[str_detect(sample_sheet[,2],'FPKM-UQ.txt'),]

#check if genes in all samples are same
genes<-list()
for(i in 1:nrow(htseq)){
htseq_file<-gzfile(paste0('gdc_download_20201104_160434.677213/',htseq[i,1],'/',htseq[i,2]),open='r')
genes[[i]]<-read.table(htseq_file)[,1]
close(htseq_file)
}
length(unique(genes)) #[1] 1 all equal
rm(genes)
#generate gene exp df
genes<-read.table(gzfile(paste0('gdc_download_20201104_160434.677213/',htseq[i,1],'/',htseq[i,2]),open='r'))[,1]
htseq_df<-data.frame(Gene=genes)
for(i in 1:nrow(htseq)){
  htseq_file<-gzfile(paste0('gdc_download_20201104_160434.677213/',htseq[i,1],'/',htseq[i,2]),open='r')
  htseq_i<-read.table(htseq_file)
  close(htseq_file)
  colnames(htseq_i)<-c('Gene',as.character(htseq[i,7]))
  htseq_df<-merge(htseq_df,htseq_i,by='Gene')
}
write.table(htseq_df,file='htseq_df.txt',quote=F,sep='\t',row.names = F)

genes<-list()
for(i in 1:nrow(FPKM)){
  FPKM_file<-gzfile(paste0('gdc_download_20201104_160434.677213/',FPKM[i,1],'/',FPKM[i,2]),open='r')
  genes[[i]]<-read.table(FPKM_file)[,1]
  close(FPKM_file)
}
length(unique(genes)) #[1] 1 all equal
rm(genes)
#generate gene exp df
genes<-read.table(gzfile(paste0('gdc_download_20201104_160434.677213/',FPKM[i,1],'/',FPKM[i,2]),open='r'))[,1]
FPKM_df<-data.frame(Gene=genes)
for(i in 1:nrow(FPKM)){
  FPKM_file<-gzfile(paste0('gdc_download_20201104_160434.677213/',FPKM[i,1],'/',FPKM[i,2]),open='r')
  FPKM_i<-read.table(FPKM_file)
  close(FPKM_file)
  colnames(FPKM_i)<-c('Gene',as.character(FPKM[i,7]))
  FPKM_df<-merge(FPKM_df,FPKM_i,by='Gene')
}
write.table(FPKM_df,file='FPKM_df.txt',quote=F,sep='\t',row.names = F)


rm(list = ls())
library('XML')
library('methods')
getwd()
sample_sheet<-read.delim('gdc_sample_sheet.2020-11-13.tsv',stringsAsFactors = F)
setwd('D:/Mirror/Documents/my study/PITT/biostatistics/class/2069 Statistical Methods for omics/final project/gdc_download_20201114_002432.893894')
file_path=paste0('./',sample_sheet[,1],'/',sample_sheet[,2])
cl <- lapply( file_path
            ,function(x){
              result = xmlParse(x)  #简单说就是把文件丢给‘XML’包去处理了
              rootnode = xmlRoot(result)
              xmldataframe = xmlToDataFrame(rootnode[2])  #构建矩阵，但是用的是rootnood【2】，因为【1】是目录文件，【2】才是有实质内容的
              return(t(xmldataframe))
            })

cl_df <- t(do.call(cbind.data.frame,cl)) 
write.table(cl_df,file = '../GDC_TCGA_GBM_clinical_df.txt',quote = F,row.names=T,sep='\t')  
