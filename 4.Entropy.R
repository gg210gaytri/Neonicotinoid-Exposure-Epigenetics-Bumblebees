library(tidyverse)
library(tidyr)
library(viridis)
library(ggpubr)
library(ggfortify)
library(emmeans)
set.seed(123)

#Getting in the cpgs
#kristi_data <- read.delim("/home/emb3/Dropbox/Projects/Gaytri_files/percentage_meth_bumblebee_data_decitabine.txt")
kristi_data <- read.delim("/Users/emb3/Dropbox/Projects/Gaytri_files/percentage_meth_bumblebee_data_decitabine.txt")#mac
kristi_data %>% mutate(across(control1:treatment3, ~ .x/100))->kristi_data # turning into betas
kristi_data$cpg<- 1:56948
kristi_data$id<-paste(kristi_data$chr,kristi_data$start,kristi_data$end)


# Wrangling the data #
entropy_data_loci<-select(kristi_data,c(1:6,11)) # just diff cpgs

#entropy_data_loci<-select(kristi_data,c(23,5:21,22)) # all meth cpgs
entropy_data<- gather(entropy_data_loci, sample, betavalue, control1:treatment3, factor_key=TRUE)
entropy_data$treatment<- substr(entropy_data$sample,1,1)

#entropy_data$sex<- as.factor(entropy_data$sex)
#entropy_data$chrono<- substr(entropy_data$sample,2,2)
#entropy_data$chrono<-as.numeric(recode(entropy_data$chrono, "1" = "16"))
beta_normalize <- function(x) {
  x_ <- ((x - min(x)) / (max(x) - min(x)))
  (x_ * (length(x_) - 1) + 0.5) / length(x_)
} #https://static1.squarespace.com/static/58a7d1e52994ca398697a621/t/5a2ebc43e4966b0fab6b02de/1513012293857/betareg_politics.pdf

entropy_data$betavalue<-beta_normalize(entropy_data$betavalue)


######Calculating entropy
entropy_factor<-(1/(length(unique(entropy_data$id))))*log2(0.5)
#entropy_factor<-0.1
resid_entropy_values = entropy_data %>% group_by(sample)  %>%
  summarise(entropy_pre = 
            sum(
              (betavalue*log2(betavalue))
              +((1-betavalue)*log2((1-betavalue))
              )),
            .groups = 'drop') # equation from https://www-sciencedirect-com.ezproxy3.lib.le.ac.uk/science/article/pii/S1097276512008933#sec3
# details on summarise https://sparkbyexamples.com/r-programming/group-by-summarise-in-r/
resid_entropy_values$entropy_measure<-entropy_factor*resid_entropy_values$entropy_pre
resid_entropy_values$treatment<-factor(c("Control","Control","Control","Decitabine","Decitabine","Decitabine"))

#Beta regression
library(betareg)
entropy_model<-betareg(resid_entropy_values$entropy_measure~resid_entropy_values$treatment)
summary(entropy_model)

# Producing figure for paper
ggboxplot(resid_entropy_values, x = "treatment", y = "entropy_measure",
          color = "treatment",  palette = "viridis",
         add = "jitter", ylab = "Entropy", xlab = "Treatment") + theme(legend.title = element_blank(),legend.position = c(0.1, 0.8))  

ggsave("/Users/emb3/Dropbox/Projects/Gaytri_files/gaytri_entropy_scatter.pdf")
