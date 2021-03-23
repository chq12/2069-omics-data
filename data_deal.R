# setwd('D:/Mirror/Documents/my study/PITT/biostatistics/class/2069 Statistical Methods for omics/final project/')
# 
# group<-read.table('group.txt',sep='\t',header = T)
library(stringr)
library(tidyverse)
# FPKM_df<-read.table('FPKM_df.txt',header=T,sep='\t')
# rownames(FPKM_df)<-FPKM_df[,1]
# FPKM_df<-FPKM_df[,-1]
# colnames(FPKM_df)<- str_sub(colnames(FPKM_df),start = 1,end=12)
# colnames(FPKM_df)<- str_replace_all(colnames(FPKM_df),'\\.','-')
# colnames(FPKM_df)  %>% unique()
# group$bcr_patient_barcode %>% unique()
# table(group$bcr_patient_barcode %in% colnames(FPKM_df))
# 
# FPKM_df_mut<-FPKM_df[,colnames(FPKM_df) %in% group$bcr_patient_barcode]
# group<-group[group$bcr_patient_barcode %in% colnames(FPKM_df_mut),]
# 
# write.table(FPKM_df_mut,file='samples160_FPKM.txt',quote=T,sep='\t',row.names = T)
# write.table(group,file='samples160_group.txt',quote=T,sep='\t',row.names = F)

FPKM_df_mut<-read.table(file='samples160_FPKM.txt',sep='\t',header = T)
group<-read.table(file='samples160_group.txt',header = T,sep='\t')

colnames(FPKM_df_mut)<-str_sub(colnames(FPKM_df_mut),1,12)
colnames(FPKM_df_mut)<- str_replace_all(colnames(FPKM_df_mut),'\\.','-')
dup_sample<-colnames(FPKM_df_mut)[duplicated(colnames(FPKM_df_mut))]

for(i in 1:length(dup_sample)){
FPKM_df_mut[,which(colnames(FPKM_df_mut)==dup_sample[i])] <- rowMeans(FPKM_df_mut[,which(colnames(FPKM_df_mut)==dup_sample[i])])
}
colnames(FPKM_df_mut)
FPKM_df_mut.1<-FPKM_df_mut[,!duplicated(colnames(FPKM_df_mut))]
FPKM_df_mut.1<-FPKM_df_mut.1[,match(group$bcr_patient_barcode,colnames(FPKM_df_mut.1))]

rowMeans(FPKM_df_mut.1) %>% median()
FPKM_df_mut.2<-FPKM_df_mut.1[rowMeans(FPKM_df_mut.1)>0.0652,]
FPKM_df_mut.3<-FPKM_df_mut.2[apply(FPKM_df_mut.2,1,function(x){sd(x)>3}),]

library(preprocessCore)
a<-log(FPKM_df_mut.3+0.01)
rnames<-rownames(a)
cnames<-colnames(a)
a<-normalize.quantiles(as.matrix(a))
rownames(a)<-rnames
colnames(a)<-cnames
plot(density(a),main='density plot')
boxplot(a)

#samr
library(samr)
sam.data<-list(x=as.matrix(a),y=as.numeric(factor(group$oneyear)),geneid= 1:nrow(a), genenames= row.names(a), logged2=T)

#Perform SAM
samr.obj<-samr(sam.data, resp.type="Two class unpaired", nperms=100, random.seed=12345)
par(mfrow=c(1,1))
samr.plot(samr.obj,del=0.4)
samr.plot(samr.obj,del=2)

delta.table <- samr.compute.delta.table(samr.obj)
siggenes.table<-samr.compute.siggenes.table(samr.obj, del=0, sam.data, delta.table)

up<-siggenes.table$genes.up
up<-up[as.numeric(up[,8])<5, ]
nrow(up) #Total number of up-regulated genes: 383; FDR=1%

lo<-siggenes.table$genes.lo
lo<-lo[as.numeric(lo[,8])<5, ]
nrow(lo) #Total number of down-regulated genes: 239; FDR=1

#Now let us draw a heatmap of the upregulated and down regulated genes. First, generate the matrix.
heatmap.data<-sam.data$x[as.numeric(c(up[,3], lo[,3])),]
#annotation_col = data.frame(survival=recode_factor(group$survival,`1`='less than 1 Y', `2`="1 ~ 5 Y",`3`="moer than 5 Y"))
annotation_col = data.frame(oneyear = recode_factor(group$oneyear,`0`='LESS THAN 1 YEAR',`1`='MORE THAN 1 YEAR'))
rownames(annotation_col)<-colnames(heatmap.data)
pheatmap(heatmap.data,scale='row',cluster_rows = F,cellheight = 5,fontsize_row = 6,fontsize_col=3 ,annotation_col = annotation_col)


group$survival[group$days_to_death < 365] <- 1
group$survival[group$days_to_death >= 365 & group$days_to_death < 1825] <- 2
group$survival[group$days_to_death >= 1825 | is.na(group$days_to_death)] <- 3
group$oneyear[group$days_to_death < 365] <- 0
group$oneyear[group$days_to_death >= 365 | is.na(group$days_to_death)] <- 1
group$fiveyear[group$days_to_death < 1825] <- 0
group$fiveyear[group$days_to_death >= 1825 | is.na(group$days_to_death)] <- 1
#t.test(as.numeric(FPKM_df_mut.2[14850,])~group$oneyear)
#t.test(as.numeric(FPKM_df_mut.2[14850,])~group$fiveyear)

a<-heatmap.data

#Dimension reduction
library(umap)
set.seed(123)
data.umap = umap(t(a))
plot(data.umap$layout[,1], data.umap$layout[,2],type="n", xlab="UMAP_1", ylab="UMAP_2")
text(data.umap$layout[group$survival==1,1], data.umap$layout[group$survival==1,2], '<1Y',col=1)
text(data.umap$layout[group$survival==2,1], data.umap$layout[group$survival==2,2], '1~5Y',col=2)
text(data.umap$layout[group$survival==3,1], data.umap$layout[group$survival==3,2], '>5Y',col=3)

text(data.umap$layout[group$oneyear==0,1], data.umap$layout[group$oneyear==0,2], '<1Y',col=1)
text(data.umap$layout[group$oneyear==1,1], data.umap$layout[group$oneyear==1,2], '>1Y',col=2)

library(Rtsne)
set.seed(123)
data.tsne<-Rtsne(t(a), dim=2, perplexity=10)
plot(data.tsne$Y[,1],data.tsne$Y[,2], type="n", xlab="tSNE_1", ylab="tSNE_2")
text(data.tsne$Y[group$survival==1,1], data.tsne$Y[group$survival==1,2], '<1Y',col=1)
text(data.tsne$Y[group$survival==2,1], data.tsne$Y[group$survival==2,2], '1~5Y',col=2)
text(data.tsne$Y[group$survival==3,1], data.tsne$Y[group$survival==3,2], '>5Y',col=3)

text(data.tsne$Y[group$oneyear==0,1], data.tsne$Y[group$oneyear==0,2], '<1Y',col=1)
text(data.tsne$Y[group$oneyear==1,1], data.tsne$Y[group$oneyear==1,2], '>1Y',col=2)


#cluster
library(NbClust)
best.K = NbClust(data=as.matrix(t(a)),min.nc = 2, max.nc = 10,method="kmeans",index="hartigan")$Best.nc[1]
# perform K means using selected K
set.seed(123)
data.kmeans<-kmeans(as.matrix(t(a)),4 )
table(as.factor(group$oneyear), data.kmeans$cluster)
# 1  2  3
# 0 13 19 28
# 1 47 14 39

adjustedRandIndex(group$oneyear, data.kmeans$cluster) #0.02759485
adjustedRandIndex(group$survival, data.kmeans$cluster) #0.01666278

library(mclust)
BIC<-mclustBIC(t(a), 1:10, "EII") # Model name EII : spherical, equal volume
plot(BIC)
summary(BIC, t(a))
data.emclust2<-mclustBIC(t(a), X=BIC) 
cluster<-summary(data.emclust2, t(a))$classification
table(as.factor(group$oneyear), as.factor(cluster))
# 1  2  3  4  5  6  7  8
# 0  9  8 13  7  5 12  5  1
# 1 25  8  9  2 14 19  9 14
# 1  2  3
# 0 16 31 13
# 1 60 34  6
adjustedRandIndex(group$oneyear, cluster) #0.09768906
adjustedRandIndex(group$survival, cluster) #0.05358887


library(sparcl)
km.perm <- KMeansSparseCluster.permute(t(a),K=3,wbounds=seq(2,20,len=20),nperms=5)
print(km.perm)
# run sparse k-means
km.out <- KMeansSparseCluster(t(a),K=3,wbounds=km.perm$bestw)
print(km.out)
table(as.factor(group$oneyear), as.factor(km.out[[1]]$Cs))
adjustedRandIndex(group$oneyear, km.out[[1]]$Cs) #0.01725721

# 1  2  3
# 0 16 24 20
# 1 28 25 47




#-------------------htseq   -------------------
  
group<-read.table('group.txt',sep='\t',header = T)
library(stringr)
library(tidyverse)
htseq_df<-read.table('htseq_df.txt',header=T,sep='\t')
rownames(htseq_df)<-htseq_df[,1]
htseq_df<-htseq_df[,-1]
colnames(htseq_df)<- str_sub(colnames(htseq_df),start = 1,end=12)
colnames(htseq_df)<- str_replace_all(colnames(htseq_df),'\\.','-')
colnames(htseq_df)  %>% unique()
group$bcr_patient_barcode %>% unique()
table(group$bcr_patient_barcode %in% colnames(htseq_df))

htseq_df_mut<-htseq_df[,colnames(htseq_df) %in% group$bcr_patient_barcode]
group<-group[group$bcr_patient_barcode %in% colnames(htseq_df_mut),]

write.table(htseq_df_mut,file='samples160_htseq.txt',quote=T,sep='\t',row.names = T)
write.table(group,file='samples160_group.txt',quote=T,sep='\t',row.names = F)

htseq_df_mut<-read.table(file='samples160_htseq.txt',sep='\t',header = T)
group<-read.table(file='samples160_group.txt',header = T,sep='\t')

colnames(htseq_df_mut)<-str_sub(colnames(htseq_df_mut),1,12)
colnames(htseq_df_mut)<- str_replace_all(colnames(htseq_df_mut),'\\.','-')
dup_sample<-colnames(htseq_df_mut)[duplicated(colnames(htseq_df_mut))]

for(i in 1:length(dup_sample)){
  htseq_df_mut[,which(colnames(htseq_df_mut)==dup_sample[i])] <- rowMeans(htseq_df_mut[,which(colnames(htseq_df_mut)==dup_sample[i])])
}
colnames(htseq_df_mut)
htseq_df_mut.1<-htseq_df_mut[,!duplicated(colnames(htseq_df_mut))]
htseq_df_mut.1<-htseq_df_mut.1[,match(group$bcr_patient_barcode,colnames(htseq_df_mut.1))]

write.table(htseq_df_mut.1,file='samples160_htseq2.txt',quote=T,sep='\t',row.names = T)


