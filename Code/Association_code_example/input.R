#------------------------------get expression file---------------------------------------

library('data.table')
library(BEDMatrix)
library('org.Hs.eg.db')
library('stringr')

path <- "~/olink_all_sample/prune/QC_ins0_all_pruned_0.8_snp.bed"
genotype <- BEDMatrix(path)
eas_protein = fread("~/ukb_protein/olink_instance_0_eas.tsv",sep="\t")
rnt_eas <- eas_protein[,1:2817]
rnt_eas[, sample:=paste(eas_protein[,FID],eas_protein[,IID],sep="_")]

#matching genotype and protein samples
sample <- intersect(rnt_eas$sample, rownames(genotype))


##get expression file
gene_list<-str_split_fixed(colnames(rnt_eas)[-ncol(rnt_eas)],'_',n=2)[,2]
gene_list_ENSEMBL<-mapIds(org.Hs.eg.db, keys = gene_list, keytype = "SYMBOL", column="ENSEMBL")
gene_list<-gene_list[-which(is.na(gene_list_ENSEMBL))]
gene_list_ENSEMBL<-gene_list_ENSEMBL[-which(is.na(gene_list_ENSEMBL))]
rnt_gene_list<-paste0('rnt','_',gene_list)
eas_protein<-t(rnt_eas[,..rnt_gene_list])
colnames(eas_protein)<-rnt_eas$sample
eas_protein<-cbind.data.frame(gene_list_ENSEMBL,eas_protein)
rownames(eas_protein)<-NULL
colnames(eas_protein)<-c('geneid',str_split_fixed(rnt_eas$sample,'_',n=2)[,2])
write.table(eas_protein,file='~/MABM/generate_db_and_cov/eas_protein.txt',row.names=F,col.names=T,quote=F)


#----------------------------------get ukb_eas plink files----------------------------------
sample_list<-data.frame(str_split_fixed(sample,'_',n=2)[,1],str_split_fixed(sample,'_',n=2)[,2])
write.table(sample_list,"~/MABM/generate_db_and_cov/sample_list.txt",row.names=F,col.names=F,quote=F)

##------------------------------------plink-------------------------------------------------
plink --bfile ~/olink_all_sample/prune/QC_ins0_all_pruned_0.8_snp --keep ~/MABM/generate_db_and_cov/sample_list.txt --make-bed --out ~/MABM/generate_db_and_cov/eas_QC_ins0_all_pruned_0.8_snp
