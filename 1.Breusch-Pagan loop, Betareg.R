#Required packages
library(betareg)
library(lmtest)
library(vegan)
library(dplyr)

# Bruesch Pagan Loop
Meth_data <- read.delim("/home/gaytri/Downloads/Cov_Files/percentage_meth_bumblebee_data_decitabine.txt")
Meth_data %>% mutate(across(control1:treatment3, ~ .x/100))->Meth_data # turning into betas
Meth_data$id<-paste(Meth_data$chr,Meth_data$start,Meth_data$end) #Added the ID at end
vmp_data <- dplyr::select(Meth_data, 1:6,10) # Selecting the first 6 columns
treatment<-as.factor(c("cont","cont","cont","deci","deci","deci"))#0=control, 1= decitabine

# Extremes 0s and 1s
n <- 6
new_zeros <- (0*(n-1)+0.5)/n
new_ones <- (1*(n-1)+0.5)/n
vmp_data[vmp_data == 0] <- new_zeros
vmp_data[vmp_data == 1] <- new_ones
vmp_data %>%
  mutate_if(is.numeric, round, digits=4) %>%
  filter(if_any(control1:treatment3, ~ .x != control1))->vmp_data # removing rows that have all the same values (crashes beta regression)

# Looping - Beta Regression Model ONLY
out_p_store_beta <- numeric(length=nrow(vmp_data))
out_bp_store_beta <- numeric(length=nrow(vmp_data))

for (i in c(1:nrow(vmp_data))){
  # Just extracting one column & binding metadata
  oneSample <- rbind(vmp_data[i,c(1:6)],treatment)
  rownames(oneSample) <- c("beta","treatment")
  oneSample <- t(oneSample)
  oneSample <- as.data.frame(oneSample)
  
  # Fit Beta Regression
  betaModel <- betareg(oneSample$beta ~ oneSample$treatment)

  # Breusch Pagan Test
  one_bp <- bptest(betaModel, 
                   varformula = ~ fitted.values(betaModel),
                   studentize = FALSE)
  
  # Outputing for storage
  out_p_store_beta[i] <- one_bp$p.value
  out_bp_store_beta[i] <- one_bp$statistic 
  
   #Print progress 
  cat("\r",round(i/nrow(vmp_data)*100),'% done')
}


# Making vector meaningful into data frame
out_store_bound_beta <- cbind(out_p_store_beta, out_bp_store_beta, vmp_data$id)
colnames(out_store_bound_beta) <- c("P-value","BP Stats","ID")
out_store_bound_beta <- as.data.frame(out_store_bound_beta)

#Export dataframe out
write.csv(out_store_bound_beta, "/home/gaytri/Downloads/Cov_Files/hetero_VMPs_bumblebee_data_.csv", row.names=FALSE)

