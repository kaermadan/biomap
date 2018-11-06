#!/usr/bin/Rscript

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

#updataed 181101
load('~/PAV/data/count/newData/update181101/10.rc.ase.rda')
#samples that retain T2 (LH156xMo17_I_T2, Mo17xPHG29_L_T2,Mo44xMo17_L_T2,NC230xB73_E_T2)
#the corresponding T1 are ('bm237','bm278','bm287','bm310')
#bm252 LH93xB73 root is clustered wrong, removed
T2 =c('bm040','bm067','bm163','bm166','bm228','bm234','bm237','bm240','bm265','bm278','bm287','bm310','bm390','bm443','bm454','bm465','bm252')

th %<>% select(-sizeFactor, -libSize, -normFactor) %>% filter(! SampleID %in% T2)

merg = left_join(tm, th, by=c('SampleID')) %>% select(-ReadCount,-rCPM,-rFPKM,-FPKM,-Replicate)
ase = left_join(ta, th, by=c('SampleID')) %>% select(-Replicate)  %>% filter(inbred =='FALSE',Genotype != 'W64AxMo17xB73') %>% separate(Genotype, c('Female','Male'), sep='x') %>% select(-SampleID, -inbred)

head(merg)
head(ase)

inb = merg %>% filter(inbred =='TRUE') %>% select(-inbred)
hyb = merg %>% filter(inbred =='FALSE', Genotype != 'W64AxMo17xB73') %>% separate(Genotype, c('Female','Male'), sep='x') %>% select(-inbred)

hyb %<>% dplyr::rename(nRC_hybrid = nRC, CPM_hybrid = CPM) 
inb %<>% dplyr::rename(nRC_female = nRC, CPM_female = CPM, Female = Genotype) %>% select(-SampleID)
inb$Female <- as.character(inb$Female)

male = inb %>% filter(Female  %in% c("Mo17","PH207","B73","Oh43","PHG29")) %>% dplyr::rename(Male= Female) %>% dplyr::rename(nRC_male = nRC_female, CPM_male = CPM_female)

head(hyb)
head(inb)
head(male)

######
###### test######
#j1= left_join(hyb[1:100000,], male, by=c('gid','Tissue','Male'))
#j2= left_join(j1, inb[1:100000,], by=c('gid','Tissue','Female'))
#head(j2) %>% as.data.frame()

#j3= left_join(j2, ase[1:100000,],by=c('gid','Tissue','Female','Male'))
#head(j3) %>% as.data.frame()

###real run, will be killed in console
j1= left_join(hyb, male, by=c('gid','Tissue','Male'))
j2= left_join(j1, inb, by=c('gid','Tissue','Female')) 
j3= left_join(j2, ase, by=c('gid','Tissue','Female','Male'))

head(j3) %>% as.data.frame()
nrow(j3)

#rearrange
df1=j3[,c(7,6,1,13,12,14,5,10,8,3,4,9,11,2)]

df1$nRC_female =  round(df1$nRC_female,0)
df1$nRC_male =  round(df1$nRC_male,0)
df1$Tissue = as.character(df1$Tissue)
head(df1) %>% as.data.frame()

saveRDS(df1, "~/PAV/data/count/newData/update181101/rearrange_merg_ase_nRC181101.rds")


