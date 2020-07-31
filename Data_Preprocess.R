## Pre-process Data for Predictive Model

# Load in Packages
library(dplyr)
library(tidyr)
library(ROSE) # for oversampling

#import data
gene_expression <-read.table('/home/ruth/METABRIC/data_mRNA_median_Zscores.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL)
clinical<-read.csv('/home/ruth/METABRIC/data_clinical_patient.txt', header=T, skip = 4, sep='\t', na.strings=c(""," ","NA"))


#check dimension
dim(clinical) # 2509   21
dim(gene_expression) # 18534  1906

clinical <- na.omit(clinical)

dim(clinical) # 1519 21

ten_year_surv <- clinical

over_ten <- filter(ten_year_surv, ten_year_surv$OS_MONTHS >120) # 735
under_ten <- filter(ten_year_surv, ten_year_surv$OS_MONTHS <120) # 784

# only include those who died from disease before 10 years
under_ten <- filter(under_ten, under_ten$VITAL_STATUS == "Died of Disease") # 390

# add vector for over 10
over_ten_vec <- rep("overten", nrow(over_ten))
over_ten$status <- over_ten_vec

# add vector for under 10
under_ten_vec <- rep("underten", nrow(under_ten))
under_ten$status <- under_ten_vec

# combine under and over ten
ten_year_surv_clinical <- dplyr::full_join(under_ten, over_ten, by = c("PATIENT_ID", "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "CELLULARITY", "CHEMOTHERAPY", "ER_IHC", "HER2_SNP6", "HORMONE_THERAPY", "INFERRED_MENOPAUSAL_STATE", "INTCLUST", "AGE_AT_DIAGNOSIS", "OS_MONTHS", "OS_STATUS", "CLAUDIN_SUBTYPE", "THREEGENE", "VITAL_STATUS", "LATERALITY", "RADIO_THERAPY", "HISTOLOGICAL_SUBTYPE", "BREAST_SURGERY", "status", "COHORT"))
ten_year_surv_clinical$status <- as.factor(ten_year_surv_clinical$status)
ten_year_data <- ten_year_surv_clinical

## For Gene expression Data
driver=c("MAP2K4","ARID1A", "PIK3CA", "TBX3", "MAP3K1", "CBFB", "TP53", "KMT2C",
         "AKT1", "GATA3", "RUNX1", "PTEN", "CDH1", "NF1", "PIK3R1", "RB1", "CDKN1B",
         "NCOR1", "CDKN2A", "ERBB2", "ERBB3", "FOXO3", "SMAD4", "KRAS", "BRCA2",
         "BAP1", "GPS2", "AGTR2", "ZFP36L1", "MEN1", "CHEK2", "SF3B1", "AHNAK2",
         "SYNE1", "MUC16")
length(driver) #35

gene_expression=gene_expression[,-2]#remove Entrez_Gene_Id column
gene_expression = gene_expression %>% filter (gene_expression$Hugo_Symbol %in% driver) # filter for driver genes
rownames(gene_expression) = gene_expression$Hugo_Symbol # set Hugo Symbol as row names
gene_expression=gene_expression[,-1] #remove Hugo_Symbol column
gene_expression = as.matrix(gene_expression)
#intersect of samples between RNA_seq and clinical_exp
clin_exprs_intersect=intersect(colnames(gene_expression), clinical$PATIENT_ID)
length(clin_exprs_intersect) #1519
clin_exprs_match<-match(clin_exprs_intersect,colnames(gene_expression))
gene_expression = gene_expression[,clin_exprs_match]
clin_exprs_match2<-match(clin_exprs_intersect,clinical$PATIENT_ID)
clinical_exp <- clinical[clin_exprs_match2,]

gene_expression[1:5,1:5]
clinical_exp[1:5,1:4]

### Corellation with survival

library('survival')
library('Hmisc')
library('survivalAnalysis')
# create event vector for RNASeq data
#>median is up-regulated genes and <median is down regulated genes
event_rna <- apply(gene_expression,1, function(x) ifelse(x > median(x),"up","down"))
event_rna <- as.data.frame(event_rna) #should be as data frame

#create new column ‘status’ with survival event as binary
clinical_exp$status <- clinical_exp$OS_STATUS
clinical_exp <- mutate(clinical_exp, status = recode(status, '0:LIVING'= 0, "1:DECEASED" = 1))
library(survivalAnalysis)

#add time and event columns of clinical_exp to event_rna
event_rna <- cbind(event_rna,clinical_exp$OS_MONTHS) #time
event_rna <- cbind(event_rna,clinical_exp$status) #event
colnames(event_rna)[36:37] <- c("time", "event")
event_rna[1:5,1:5] 

set.seed(420)
### this works for extraction of results
univ_formulas <- sapply(driver,
                        function(x) as.formula(paste('Surv(time, event)~', as.factor(x))))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = event_rna)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
res$p.value <- as.numeric(as.character(res$p.value))
surv_sig_genes <- filter(res, res$p.value <= 0.05) # filter for significant genes (pvalue<0.05)

# vector of 14 significant genes
sig_genes <- rownames(surv_sig_genes)

# this is the filtered event rna dataframe with 14 sig genes
filtered_event_rna <- select(event_rna, one_of(sig_genes))

# 10 year filtered clinical data 
rownames(ten_year_data) <- ten_year_data$PATIENT_ID

# combined genomic & clinical by row names
ten_year_combined <- merge(ten_year_data, filtered_event_rna, by = 0)

ten_year_combined <- droplevels(ten_year_combined)

## Random Oversampling
n_new <- 735 / (1 - 0.50)
oversampling_result <- ovun.sample(formula = status ~ ., data = ten_year_combined,
                                   method = "over", N = n_new, seed = 2018)

oversampled_ten_combined <- oversampling_result$data

write.csv(oversampled_ten_combined, "combined_data.csv") # export as .csv file
