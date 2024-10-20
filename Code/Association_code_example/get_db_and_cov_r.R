library('stringr')
library(BEDMatrix)
library('org.Hs.eg.db')
library('data.table')
library(RSQLite)

main_dir <- '~/MABM/generate_db_and_cov'
annotation_file_name='~/MABM/generate_db_and_cov/gencode.v27.GRCh37.txt'
plink_file_name='~/MABM/generate_db_and_cov/eas_QC_ins0_all_pruned_0.8_snp'
expression_file_name='~/MABM/generate_db_and_cov/eas_protein.txt'
output_file_name='MABM_ukbeas_db_and_cov'

eas_genotype<-BEDMatrix('~/MABM/generate_db_and_cov/eas_QC_ins0_all_pruned_0.8_snp.bed')
colnames(eas_genotype)<-str_split_fixed(colnames(eas_genotype),"_",n=2)[,1]
eas_genotype<-as.matrix(eas_genotype)
#fill in genotype missing values with mean
for (i in 1:dim(eas_genotype)[2]){
  eas_genotype[is.na(eas_genotype[, i]), i] <- mean(eas_genotype[, i], na.rm=T)
}

best_cutoff = read.table("~/MABM/5-fold-cv/best_cv.txt",header=T)

genename<-str_split_fixed(best_cutoff$name,"_",n=2)[,2]
geneid<-mapIds(org.Hs.eg.db, keys = genename, keytype = "SYMBOL", column="ENSEMBL")  #get genes ENSEMBL ID
best_cutoff$geneid<-geneid



#--------------get weights and covariance-----------------------------

#mkdir
if(!dir.exists(paste0(main_dir,'/gene/'))){
  dir.create(paste0(main_dir,'/gene/'))
}
if(!dir.exists(paste0(main_dir,'/weights/'))){
  dir.create(paste0(main_dir,'/weights/'))
}
if(!dir.exists(paste0(main_dir,'/cov/'))){
  dir.create(paste0(main_dir,'/cov/'))
}

for (i in 1:dim(best_cutoff)[1]) {
  if (best_cutoff[i,3]=='TransLasso' & best_cutoff[i,2]>0.01) {
    protein = best_cutoff[i,1]
    weights <- read.table(file=paste("~/prediction_model/trans_lasso_result/", protein, ".txt", sep=""), colClasses = c("character","character","numeric"))
  } else if (best_cutoff[i,3]=='Lasso' & best_cutoff[i,2]>0.01) {
    protein = best_cutoff[i,1]
    weights <- read.table(file=paste("~/prediction_model/lasso_result/", protein, ".txt", sep=""), colClasses = c("character","character","numeric"))
  } else if (best_cutoff[i,3]=='HEAT' & best_cutoff[i,2]>0.01) {
    protein = best_cutoff[i,1]
    weights <- read.table(file=paste("~/prediction_model/HEAT_result/", protein, ".txt", sep=""), colClasses = c("character","character","numeric"))
} else {next}
  #generate weight df
  weights_df <- data.frame('rsid'=weights$V1, 'weights'=weights$V3, 'effect_allele'=weights$V2, stringsAsFactors=F)
  #performance
  r2 <- best_cutoff[i, 2]
  weights_df$r2<-r2
  if (nrow(weights_df)>0 & r2>0.01){
    #generate covariance matrix
    cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=1)) 
    colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
    #calculate covariance
    if(nrow(weights_df)>1){
      snp_list<-weights_df$rsid
      o_i=1
      for (k in 1:length(snp_list)){
        for (j in k:length(snp_list)){
          cov_df[o_i,2]<-snp_list[k]
          cov_df[o_i,3]<-snp_list[j]
          cov_df[o_i,4]<-round(cov(eas_genotype[,snp_list[k]],eas_genotype[,snp_list[j]]),3)
          o_i=o_i+1
        }
      }
    }else{
      cov_df[,'RSID1']=cov_df[,'RSID2']=weights_df[1,'rsid']
      cov_df[,'VALUE']=round(cov(eas_genotype[,weights_df$rsid],eas_genotype[,weights_df$rsid]),3)
    }
    cov_df[,1]<-geneid[i]

    #output cov
    write.table(cov_df,paste0(main_dir,'/cov/',geneid[i],'.txt'),sep='\t',quote = F,row.names = F)
    #output weights
    write.table(weights_df,paste0(main_dir,'/weights/',geneid[i],'.txt'),sep='\t',quote = F,row.names = F)
  }else{
    print('INFO failed to build the model')
  } 
} 


#-----------------get .db and .cov----------------------------------
#mkdir
if(!dir.exists(paste0(main_dir,'/output/'))){
  dir.create(paste0(main_dir,'/output/'))
}

#load annotation
anno<-read.table(paste0(annotation_file_name),header = T,stringsAsFactors = F)

#gene list
gene_list<-sub('....$','',dir(paste0(main_dir,'/cov/')))


#load weights
print('INFO loading weights...')
print(paste0('INFO totally ',length(gene_list),' imputable genes'))
db_weights<-list()
n_rows=0
for (i in 1:length(gene_list)){
  if(i %in% c(seq(1,1000)*100,length(gene_list))){
    print(paste0('     ',i,' loaded'))
  }
  geneid<-gene_list[i]
  db_weights[[i]]<-read.table(paste0(main_dir,'/weights/',geneid,'.txt'),header = T,stringsAsFactors = F)
  n_rows=n_rows+nrow(db_weights[[i]])
}

#db.weight
print('INFO generating db file...')
weights<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows))
colnames(weights)<-c('rsid','gene','weight','eff_allele')

#db.sample_info
n_sample<-read.table(paste0(expression_file_name),header = T,nrows=1,stringsAsFactors = F,check.names = F)
sample_info<-data.frame('n.samples'=ncol(n_sample)-1)

#db.extra
extra<-as.data.frame(matrix(data=NA,ncol=4,nrow=length(gene_list)))
colnames(extra)<-c('gene','genename','pred.perf.R2','n.snps.in.model') 

n_rows=0
for (i in 1:length(gene_list)){
  #extra.gene
  extra[i,1]<-gene_list[i]
  #extra.genename
  extra[i,2]<-genename[which(best_cutoff$geneid==gene_list[i])]
  #extra.r2
  extra[i,3]<-db_weights[[i]][1,4]
  #extra.n.snps
  extra[i,4]<-nrow(db_weights[[i]]) 
  
  weights[(n_rows+1):(n_rows+nrow(db_weights[[i]])),]<-as.matrix(data.frame('rsid'=db_weights[[i]]$rsid,'gene'=gene_list[i],'weight'=db_weights[[i]]$weights,'effect_allele'=db_weights[[i]]$effect_allele))
  n_rows=n_rows+nrow(db_weights[[i]])  
}

#load bim to get ref allele
bim<-read.table(paste0(plink_file_name,'.bim'),header = F,stringsAsFactors = F)
bim<-bim[bim[,2]!='.',];bim<-bim[!duplicated(bim[,2]),]
bim<-bim[which(bim[,2] %in% weights[,1]),-1]

weights<-merge(weights,bim,by=1)
weights$ref_allele<-ifelse(weights$eff_allele==weights$V5,weights$V6,weights$V5)
weights<-weights[,c(1,2,3,9,4)]

if(file.exists(paste0(main_dir,"/output/",output_file_name,".db"))){
  file.remove(paste0(main_dir,"/output/",output_file_name,".db"))
}

db<-dbConnect(RSQLite::SQLite(), paste0(main_dir,"/output/",output_file_name,".db"))
dbWriteTable(db, "weights", weights)
dbWriteTable(db, "extra", extra)
dbWriteTable(db, "sample_info", sample_info)
dbDisconnect(db)

#load cov
print('INFO loading cov files...')
cov<-list()
n_rows=0
for (i in 1:length(gene_list)){
  if(i %in% c(seq(1,1000)*100,length(gene_list))){
    print(paste0('     ',i,' loaded'))
  }
  geneid<-gene_list[i]
  cov[[i]]<-read.table(paste0(main_dir,'/cov/',geneid,'.txt'),header = T,stringsAsFactors = F)
  n_rows=n_rows+nrow(cov[[i]])
}
cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows)) 
colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
#combine cov
print('INFO combining cov file...')
n_rows=0
for (i in 1:length(gene_list)){
  cov_df[(n_rows+1):(n_rows+nrow(cov[[i]])),]<-cov[[i]]
  n_rows=n_rows+nrow(cov[[i]])
}

#output cov
print('INFO writing cov file...')
write.table(cov_df,paste0(main_dir,'/output/',output_file_name,'.cov'),sep='\t',row.names = F,quote = F)

print('INFO finished')
