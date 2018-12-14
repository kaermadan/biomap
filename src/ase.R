#!/user/bin/Rscript

require(tidyverse)
require(reshape2)
require(magrittr)
library(wesanderson)
require(RColorBrewer)
require(ggsci)
library(scales)
require(qvalue)
library(scales)#use percent in numbers
library(viridis)#put plots together
require(cowplot)
require(readxl)
library(purrr)
library(broom)
library(tibble)

options(dplyr.width = Inf)



print(paste("Job -- ase_newData_allMale -- begin running!!! on", date()))

ase_new = readRDS('~/projects/biomap/data/rearrange_merg_ase_nRC181101.rds')
test1= ase_new %>% filter(ncft/(n0+n1) < 0.1, (n0+n1) >15, (nRC_female+nRC_male) >15, Female !='Q381')
test= test1

print(paste("Executed dataset row number after filtering is", nrow(test1)))
################# Two binom test and one Fisher test#######################################
##keep male/female order consistent across binom and Fisher!! especially for p!=0.5
##181214 when do binom in parents, for endosperm tissue, p should set as 0.5

get_test <- function(test){

	#1 parent binom
	test_parent_binom <- apply(test, 1, function( x ) {

		if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] !='Endosperm') {
			model_binom <- binom.test( x = as.numeric( x[8] ),
                             n = as.numeric( x[8])+as.numeric( x[9] ), 
                             p = 0.5, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(Parent.binom.p = as.numeric(model_binom$p.value))) 
            #Parent.binom.p_CI95_low = as.numeric(model_binom$conf.int[1]),
            #Parent.binom.p_CI95_high = as.numeric(model_binom$conf.int[2])))

		} else if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] =='Endosperm') {
			model_binom <- binom.test( x = as.numeric( x[8] ),
                             n = as.numeric( x[8])+as.numeric( x[9] ), 
                             p = 0.5, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(Parent.binom.p = as.numeric(model_binom$p.value))) 

		}else {
			return( c(Parent.binom.p = as.numeric(NA)))
            #Parent.binom.p_CI95_low = as.numeric(NA),
            #Parent.binom.p_CI95_high =as.numeric(NA)))
		}			      
	})

	#F1 binom
	test_hyb_binom <- apply(test, 1, function( x ) {

		if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] !='Endosperm') {
			model_binom <- binom.test( x = as.numeric( x[5] ),
                             n = as.numeric( x[4])+as.numeric( x[5] ), 
                             p = 0.5, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(Hybrid.binom.p = as.numeric(model_binom$p.value))) 
            #Hybrid.binom.p_CI95_low = as.numeric(model_binom$conf.int[1]),
            #Hybrid.binom.p_CI95_high = as.numeric(model_binom$conf.int[2])))

		} else if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] =='Endosperm') {

			model_binom <- binom.test( x = as.numeric( x[5] ),
                             n = as.numeric( x[4])+as.numeric( x[5] ), 
                             p = 0.6666666, 
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(Hybrid.binom.p = as.numeric(model_binom$p.value))) 

		} else {
			return( c(Hybrid.binom.p = as.numeric(NA))) 
            #Hybrid.binom.p_CI95_low = as.numeric(NA),
            #Hybrid.binom.p_CI95_high =as.numeric(NA)))
		}			      
	})

	##fisher， note the order of th efour value in fisher test!!!!! and endosperm allele in F1 for female was divided by 2!
	test_fisher<- apply(test, 1, function( x ) {

		if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] !='Endosperm') {
			model_fisher <- fisher.test( x = matrix(c(as.numeric( x[5] ),
                            as.numeric( x[4]), as.numeric( x[8]), as.numeric( x[9])),ncol=2), 
                           
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(trans.p = as.numeric(model_fisher$p.value))) 
            #Hybrid.binom.p_CI95_low = as.numeric(model_binom$conf.int[1]),
            #Hybrid.binom.p_CI95_high = as.numeric(model_binom$conf.int[2])))

		} else if ((as.numeric(x[6])/(as.numeric(x[5])+as.numeric(x[4])) < 0.1) & ((as.numeric(x[5])+as.numeric(x[4])) >15) & ((as.numeric(x[8])+as.numeric(x[9])) >15) & x[7] =='Endosperm') {

			model_fisher <- fisher.test( x = matrix(c(round(as.numeric( x[5])*0.5, 0),
                            as.numeric( x[4]), as.numeric( x[8]), as.numeric( x[9])),ncol=2), 
                           
                             alternative = "two.sided", 
                             conf.level = 0.95)
			return( c(trans.p = as.numeric(model_fisher$p.value))) 

		} else {
			return( c(trans.p = as.numeric(NA))) 
            #Hybrid.binom.p_CI95_low = as.numeric(NA),
            #Hybrid.binom.p_CI95_high =as.numeric(NA)))
		}			      
	})

	#cbind results from three test
	cbind( Parent.binom.p = test_parent_binom,Hybrid.binom.p = test_hyb_binom, trans.p = test_fisher )	
}


#do parallel computing, not working!
#y = bplapply(1:nrow(test), my, BPPARAM = bpparam)
p_value = get_test(test)

##combine original dataset and the three p-values	
test_result <- do.call(cbind, list(test, p_value))


##################################### qvalue ################################################
#fisher test may output some values as 1, which is actually greater than 1 in R system, so do some trim here
test_result$Parent.binom.p <- ifelse(test_result$Parent.binom.p > 1, 1, test_result$Parent.binom.p)
test_result$Hybrid.binom.p <- ifelse(test_result$Hybrid.binom.p > 1, 1, test_result$Hybrid.binom.p)
test_result$trans.p <- ifelse(test_result$trans.p > 1, 1, test_result$trans.p)

#convert factor col to character to aovid NA generating in subsequent rbind  step
factor_col= names(test_result)[c(2,3,7)]
test_result[factor_col]<- sapply(test_result[factor_col], as.character)

##qvalue must be calculated based on each sample, not all sample together!

sample_all=NULL
for (i in unique(test_result$SampleID)){
	
	sample= test_result[which(test_result$SampleID==i),]
	sample_q= sample %>% mutate(
	                             Parent.binom.qvalue = ifelse(is.na(Parent.binom.p) == FALSE, qvalue(Parent.binom.p)$qvalues, 'NA'),
	                             Hybrid.binom.qvalue = ifelse(is.na(Hybrid.binom.p) == FALSE, qvalue(Hybrid.binom.p)$qvalues, 'NA'),
	                             trans.qvalue = ifelse(is.na(trans.p) == FALSE, qvalue(trans.p)$qvalues, 'NA'))

	sample_all= rbind(sample_all, sample_q)
}

################# assign genes to regulatory categories ######################################
ELSE=TRUE
test_result_f= sample_all %>% mutate(Reg.Cat= case_when(
          (Parent.binom.qvalue < 0.05 & Hybrid.binom.qvalue < 0.05 & trans.qvalue > 0.05) ~ 'Cis',
          (Parent.binom.qvalue < 0.05 & Hybrid.binom.qvalue > 0.05 & trans.qvalue < 0.05) ~ 'Trans',
          (Parent.binom.qvalue < 0.05 & Hybrid.binom.qvalue < 0.05 & trans.qvalue < 0.05 & (n1-n0)/(n1+n0)*(nRC_male-nRC_female)/(nRC_male+nRC_female) > 0 ) ~ 'Cis+Trans',
          (Parent.binom.qvalue < 0.05 & Hybrid.binom.qvalue < 0.05 & trans.qvalue < 0.05 & (n1-n0)/(n1+n0)*(nRC_male-nRC_female)/(nRC_male+nRC_female) < 0 ) ~ 'CisxTrans',
          (Parent.binom.qvalue > 0.05 & Hybrid.binom.qvalue < 0.05 & trans.qvalue < 0.05) ~ 'Compensatory',
          (Parent.binom.qvalue > 0.05 & Hybrid.binom.qvalue > 0.05 ) ~ 'Conserved',
              
       	   ELSE ~ 'Ambiguous'))


#write.csv(test_result_f,'/home/hirschc1/lixx5447/PAV/data/count/newData/ase_newData_allMale.csv')

saveRDS(test_result_f, "~/projects/biomap/data/ase_newData181214.rds")
######
head(test_result_f)
print('ase_newData_allMale completed!!!,
	Results were saved in "ase_newData_allMale.rds"
	 and "ase_newData_allMale.csv"')

print(paste("Job -- ase_newData_allMale -- Ending!!! on", date()))


