### Making sense of the VMPs
# Setting- Up
betaVMPs <- read.csv("/home/gaytri/Downloads/Cov_Files/hetero_VMPs_bumblebee_data_.csv")

# P-Value Adjustment
P.value_adj_beta <- p.adjust(betaVMPs$P.value,method="holm")
betaVMPs <- cbind(betaVMPs,P.value_adj_beta)

# Extracting Significant 0.05
betaVMPsSIG <- betaVMPs[betaVMPs$P.value_adj_beta < 0.05, ]
#49,012 survived the correction

# Exporting Adjusted
write.csv(betaVMPsSIG, "/home/gaytri/Downloads/Cov_Files/hetero_VMPs_bumblebee_data_FDR_SIG.csv", row.names=FALSE)
write.csv(betaVMPs, "/home/gaytri/Downloads/Cov_Files/hetero_VMPs_bumblebee_data_FDR_ALL.csv", row.names=FALSE)

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("biomaRt")

# Setting up Ensembl database
library(biomaRt)
ensembl_metazoa <- useEnsemblGenomes(biomart = "metazoa_mart")
searchDatasets(ensembl_metazoa, pattern = "Nasonia")
searchDatasets(ensembl_metazoa, pattern = "Bombus")
ensembl = useEnsemblGenomes(biomart="metazoa_mart", dataset="btgca910591885v2_eg_gene") 
attributes <- listAttributes(ensembl) #all the columns you can get
emsemblgenes <- getBM(attributes=c('ensembl_gene_id',
                                   "description",
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "strand",
                                   "go_id"), 
                      mart=ensembl)
ensemblgenes <- as.data.frame(emsemblgenes)
## N.B.: GO Terms are NOT specific to BP.

# Reformatting/ Splitting ID
library(stringr)
boundID <- betaVMPsSIG$ID
sepID_matrix <- str_split_fixed(boundID, " ", 3)
colnames(sepID_matrix) <- c("chr","start","end")
betaVMPsSIG <- cbind(betaVMPsSIG, sepID_matrix)

# Changing chr name so it's compatible with Ensembl data
unique(betaVMPsSIG$chr) #Just trying to find out the variations (see notes)

Ensembl_chr <- betaVMPsSIG$chr
for (i in 1:nrow(betaVMPsSIG)) {
  if (Ensembl_chr[i] == "NC_015762") {
    Ensembl_chr[i] <- "1"
    } else if (Ensembl_chr[i] == "NC_015763") {
    Ensembl_chr[i] <- "2"
    } else if (Ensembl_chr[i] == "NC_015764") {
    Ensembl_chr[i] <- "3"
    } else if (Ensembl_chr[i] == "NC_015765") {
    Ensembl_chr[i] <- "4"
    } else if (Ensembl_chr[i] == "NC_015766") {
    Ensembl_chr[i] <- "5"
    } else if (Ensembl_chr[i] == "NC_015767") {
    Ensembl_chr[i] <- "6"
    } else if (Ensembl_chr[i] == "NC_015768") {
    Ensembl_chr[i] <- "7"
    } else if (Ensembl_chr[i] == "NC_015769") {
    Ensembl_chr[i] <- "8"
    } else if (Ensembl_chr[i] == "NC_015770") {
    Ensembl_chr[i] <- "9"
    } else if (Ensembl_chr[i] == "NC_015771") {
    Ensembl_chr[i] <- "10"
    } else if (Ensembl_chr[i] == "NC_015772") {
    Ensembl_chr[i] <- "11"
    } else if (Ensembl_chr[i] == "NC_015773") {
    Ensembl_chr[i] <- "12"
    } else if (Ensembl_chr[i] == "NC_015774") {
    Ensembl_chr[i] <- "13"
    } else if (Ensembl_chr[i] == "NC_015775") {
    Ensembl_chr[i] <- "14"
    } else if (Ensembl_chr[i] == "NC_015776") {
    Ensembl_chr[i] <- "15"
    } else if (Ensembl_chr[i] == "NC_015777") {
    Ensembl_chr[i] <- "16"
    } else if (Ensembl_chr[i] == "NC_015778") {
    Ensembl_chr[i] <- "17"
    } else if (Ensembl_chr[i] == "NC_015779") {
    Ensembl_chr[i] <- "18"
    } else {
    Ensembl_chr[i] <- ""
    }
  }
betaVMPsSIG <- cbind(betaVMPsSIG, Ensembl_chr)

# Annotating significant VMPs with Ensembl database
ensemblcompare <- betaVMPsSIG[, c(8,7,6)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(betaVMPsSIG)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

#Export
write.csv(ensemblres, "/home/gaytri/Downloads/Cov_Files/EnsemblBetaVMPsSIG_final.csv", row.names=FALSE)
