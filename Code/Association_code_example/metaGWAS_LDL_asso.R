library('data.table')
library('RSQLite')
library('stringr')


#get gene list
a = list.files('~/MABM/generate_db_and_cov/weights')
gene_list = str_split_fixed(a,".txt",n=2)[,1]
write.table(gene_list, '~/MABM/gene_list.txt',row.names=F,col.names=F,quote=F)

db_path = '~/MABM/generate_db_and_cov/output/MABM_ukbeas_db_and_cov.db'
cov_path = '~/MABM/generate_db_and_cov/output/MABM_ukbeas_db_and_cov.cov'

bim_path = '~/MABM/generate_db_and_cov/eas_QC_ins0_all_pruned_0.8_snp.bim'
gene_list_path = '~/MABM/gene_list.txt'
gwas_path = '~/LDLC_gwas/GCST90278665.h.tsv.gz'
gwas_variant_col = 'rsid'
gwas_beta_col = 'beta'
gwas_se_col = 'standard_error'
gwas_p_col = 'p_value'
gwas_eff_allele_col = 'effect_allele'
gwas_ref_allele_col = 'other_allele'
asso_out_path = '~/MABM/metaGWAS_LDL_asso_result.txt'


cat('INFO loading weight file ...\n')
con <- dbConnect(RSQLite::SQLite(), dbname=db_path) #establish connections
weights = dbReadTable(con,"weights")
colnames(weights) = c('rsid', 'gene', 'weight', 'model_ref_allele', 'model_eff_allele')
extra_info = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

cat('INFO loading covariance ...\n')
cov_all = data.frame(fread(cov_path))
cat('INFO loading GWAS results ...\n')
gwas <- data.frame(fread(gwas_path))

#convert beta, or, se, pval
gwas <- gwas[!duplicated(gwas[,gwas_variant_col]),]
if(is.na(gwas_beta_col)){
  gwas$beta = log(gwas[,gwas_or_col],2.718)
  gwas_beta_col = 'beta'
}
if(is.na(gwas_se_col)){
  if(is.na(gwas_zscore_col)){
    gwas$z_score = qnorm(1-gwas[,gwas_p_col]/2)
    gwas$se = abs(gwas[,gwas_beta_col] / gwas[,'z_score'])
    gwas_se_col = 'se'
    
  }else{
    gwas$se = gwas[,gwas_beta_col] / gwas[,gwas_zscore_col]
    gwas_se_col = 'se'
  }
}

#upper allele
gwas[,gwas_eff_allele_col] = toupper(gwas[,gwas_eff_allele_col])
gwas[,gwas_ref_allele_col] = toupper(gwas[,gwas_ref_allele_col])

#get the list of snps both available from gwas and prediction model
weights <- merge(weights, gwas, by.x = 'rsid', by.y = gwas_variant_col)

#strand check
cat('INFO strand check and allele flipping ...\n')
weights$allele_idx1 = weights[,'model_eff_allele'] == weights[,gwas_eff_allele_col]
weights$allele_idx2 = weights[,'model_eff_allele'] == weights[,gwas_ref_allele_col]

same_strand_pos = which(weights$allele_idx1 | weights$allele_idx2)
num_of_variants_have_same_strand = length(same_strand_pos)

if (num_of_variants_have_same_strand / nrow(weights) <0.9){
  stop('More than 10% of the variants in GWAS summary statistics have different strand, please check')
}

weights = weights[same_strand_pos,]

#flip alleles
weights[,gwas_beta_col] = ifelse(weights$allele_idx1, weights[,gwas_beta_col], weights[,gwas_beta_col]*-1)


gene_list <- unique(weights$gene)

if(!is.na(gene_list_path)){
  specified_gene_list = read.table(gene_list_path, header = T,stringsAsFactors = F)
  gene_list = intersect(gene_list, specified_gene_list[,1])
}


cat(paste0('INFO Association test for ',length(gene_list),' imputable genes ...\n'))

asso<-function(gene_id){
  
  df = weights[weights$gene == gene_id, ]
  df$weight = as.numeric(df$weight)
  df$id = seq(1,nrow(df))
  
  cov_tmp = cov_all[cov_all$GENE == gene_id, ]
  cov_tmp = merge(cov_tmp, df[,c('rsid','id')], by.x = 'RSID1', by.y = 'rsid')
  cov_tmp = merge(cov_tmp, df[,c('rsid','id')], by.x = 'RSID2', by.y = 'rsid')
  
  cov_matrix = matrix(data = NA, nrow = nrow(df), ncol = nrow(df))
  
  for(i in 1:nrow(cov_tmp)){
    cov_matrix[cov_tmp[i,'id.x'], cov_tmp[i,'id.y']] = cov_tmp[i, 'VALUE']
    cov_matrix[cov_tmp[i,'id.y'], cov_tmp[i,'id.x']] = cov_tmp[i, 'VALUE']
  }
  
  #var(GReX)
  var_GReX = as.numeric(t(df$weight) %*% cov_matrix %*% df$weight)
  
  #beta
  effect_size = sum(df$weight * df[,gwas_beta_col] * diag(cov_matrix) / var_GReX)
  
  #z-score
  zscore = sum(df$weight * sqrt(diag(cov_matrix)) * df[,gwas_beta_col] / sqrt(var_GReX) / df[,gwas_se_col], na.rm = T)
  ans = list()
  ans[['effect_size']] = effect_size
  ans[['zscore']] = zscore
  ans[['n_snps_used']] = nrow(df)
  return(ans)
  
}

result <- as.data.frame(t(sapply(gene_list, function(x) asso(x))))
result$effect_size = as.numeric(result$effect_size)
result$zscore = as.numeric(result$zscore)
result$n_snps_used = as.numeric(result$n_snps_used)

#output
result = as.data.frame(result, stringsAsFactors = F)
result$gene = gene_list
result$pvalue = (1-pnorm(abs(result$zscore)))*2  #two sided
result = result[which(result$effect_size * result$zscore>0),]
result = merge(result, extra_info, by='gene')

result = result[,c('gene','genename','zscore','effect_size','pvalue','pred.perf.R2','n_snps_used','n.snps.in.model')]
colnames(result) = c('gene','genename','zscore','effect_size','pvalue','pred_perf_R2','n_snps_used','n_snps_in_model')

result = result[order(abs(result$zscore),decreasing=T),]
write.table(result, asso_out_path, quote = F, row.names = F, sep = '\t')

cat('INFO done\n')



