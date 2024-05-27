# Ensembl's GO 
# Setting- Up 
betaVMPsig <- read.csv("/home/gaytri/Downloads/Cov_Files/hetero_VMPs_bumblebee_data_FDR_SIG.csv", na.strings=c("","NA"))
ensemblass <- read.csv("/home/gaytri/Downloads/Cov_Files/EnsemblBetaVMPsSIG_final.csv", na.strings=c("","NA"))

# Combine both dataframes 
mother <- cbind(betaVMPsig,ensemblass)

# Genes & GO Terms only as vector
Gene <- na.omit(mother$ensembl_gene_id)
GO <- na.omit(mother$go_id)

# Count for gene
genelist <- unique(Gene)
genelist_count <- numeric(length(genelist))
for(i in 1:length(genelist)){
  for(j in 1:length(Gene)){
    cur_gene <- genelist[i]
    cur_dum <- Gene[j]
    if(cur_gene == cur_dum){
      genelist_count[i] <- genelist_count[i] + 1
    }
  }
}
Gene_res <- cbind(genelist,genelist_count)
write.csv(Gene_res, "/home/gaytri/Downloads/Cov_Files/Gene_res_final.csv", row.names=FALSE)


# Count for GO Terms 
GOterms <- unique(GO)
GOterms_count <- numeric(length(GOterms))
for(i in 1:length(GOterms)){
  for(j in 1:length(GO)){
    cur_GO <- GOterms[i]
    cur_dum <- GO[j]
    if(cur_GO == cur_dum){
      GOterms_count[i] <- GOterms_count[i] + 1
    }
  }
}
GO_res <- cbind(GOterms, GOterms_count)
write.csv(GO_res, "/home/gaytri/Downloads/Cov_Files/GO_res.csv", row.names=FALSE)
