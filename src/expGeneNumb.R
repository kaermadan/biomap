#!/user/bin/Rscript

#just use to  test expGNumb use different CPM
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

options(dplyr.width = Inf)#print all col in tibble
###@@@@@@@@@@@@@@@@@@@@@@@@@@@ Make some Plot @@@@@@@@@@@@@@@@@@@@@@@@@@@put here for convinent
pavFilter = readRDS("~/PAV/data/count/newData/rearrange_merg_ase_nRC_addGeneType_addPAV_addPatterns_addSync_addHybExp_spOff_rmDup_F.rds")

#!!!based on previously PCA results, LH93xB73 -Root should be Leaf, or removed
pavFilter %<>% mutate(tmp = paste(Female,'x', Male,Tissue, sep='')) %>% filter(tmp !='LH93xB73Root') %>% select(-tmp) %>% filter(!Female %in% c('Q381','PHG83'))

#pavFilter %<>% mutate(tmp = paste(Female,'x', Male,Tissue, sep='')) %>% filter(tmp !='LH93xB73Root') %>% select(-tmp) %>% filter(!Female %in% c('Q381','PHG83'))



##################################################################

load('~/PAV/data/count/newData/10.norm.RData')
load('~/PAV/biomap/analysis/11.ase.RData')
objects()

head(tl)
print(paste("Job -- -- begin running!!! on", date()))
###Mo44 only have R,S tissue in inbred, but have S,L,R,I tissue for both xB73, and xMo17 hybrids, so the other four hybrids will not have SPE data

#samples that retain T2 (LH156xMo17_I_T2, Mo17xPHG29_L_T2,Mo44xMo17_L_T2,NC230xB73_E_T2)
#the corresponding T1 are ('bm237','bm278','bm287','bm310')
T2 =c('bm040','bm067','bm163','bm166','bm228','bm234','bm237','bm240','bm265','bm278','bm287','bm310','bm390','bm443','bm454','bm465')

#remove these sample in further analysis
rmv = grep(paste('Q381*','PHG83*',sep='|') ,tl$Genotype,value=T)
rmv1 = tl %>% filter(Genotype %in% rmv)

#based on previously PCA results, LH93xB73 -Root should be Leaf, or removed, which is bm252, and also remove triple hybrid, bm466
rmv2 = c(unique(rmv1$SampleID),'bm252','bm466')
rmv3 = unique(c(rmv2, T2))

##look at libsize summmary 
tl1 = tl %>% filter(! SampleID %in% rmv3) %>% summarise_at(vars(libSize), funs(median, mean, max, min, sd))

#in summary 180830
##430 samples, 270 hybrids, 160 inbreds; 135 inbreds--- 270 hybrids;

tl %<>% select(-sizeFactor, -libSize, -normFactor) %>% filter(! SampleID %in% T2)

tl$Genotype <- as.character(tl$Genotype)
tl$Tissue <- as.character(tl$Tissue)

merg = left_join(tm, tl, by=c('SampleID')) %>% select(-ReadCount,-rCPM,-rFPKM,-FPKM,-Treatment)
#ase = left_join(t_ase, tl, by=c('sid' = 'SampleID')) %>% select(-Treatment) %>% dplyr::rename(SampleID=sid) %>% filter(inbred =='FALSE',Genotype != 'W64AxMo17xB73') %>% separate(Genotype, c('Female','Male'), sep='x') %>% select(-SampleID, -inbred)


inb = merg %>% filter(inbred =='TRUE') %>% select(-inbred)
hyb = merg %>% filter(inbred =='FALSE', Genotype != 'W64AxMo17xB73') %>% separate(Genotype, c('Female','Male'), sep='x') %>% select(-inbred)

hyb %<>% dplyr::rename(nRC_hybrid = nRC, CPM_hybrid = CPM) 
inb %<>% dplyr::rename(nRC_female = nRC, CPM_female = CPM, Female = Genotype) %>% select(-SampleID)
inb$Female <- as.character(inb$Female)
inb$Tissue <- as.character(inb$Tissue)
hyb$Tissue <- as.character(hyb$Tissue)

male = inb %>% filter(Female  %in% c("Mo17","PH207","B73","Oh43","PHG29")) %>% dplyr::rename(Male= Female) %>% dplyr::rename(nRC_male = nRC_female, CPM_male = CPM_female)

print('male done')
head(male)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   p002    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
##############>>>p002 summary expressed gene numbers
#updated 181112, corr of gene numb diff and genetic distance, etc
inb = pavFilter %>% filter(gene_type=='protein_coding') %>% filter(CPM_female >  0.2) %>% group_by(Male, Female, Tissue) %>% dplyr::count() %>% dplyr::rename(inb = n)

hyb = pavFilter %>% filter(gene_type=='protein_coding') %>% filter(CPM_hybrid >  0.2) %>% group_by(Male, Female, Tissue) %>% dplyr::count() %>% dplyr::rename(hyb = n)

inb %<>% mutate(state='Inbred') %>% dplyr::rename(ExpGeneNumb = inb)

hyb %<>% mutate(state='Hybrid') %>% dplyr::rename(ExpGeneNumb = hyb)

comb = rbind(hyb, inb)
comb

## expressed gene numbers, keep PHG83 here(remove PHG83, and make average between male and Female for inbred expGene, 180825)

#add male number
#male is from line 460
ml = male %>% filter(CPM_male > 0.2) %>% group_by(Tissue, Male) %>% dplyr::count() 
ml1 = ml %>% filter(! Male %in% c('Oh43','PHG29'))
ml2 = ml %>% filter( Male %in% c('Oh43','PHG29'))

# mean(male, female)
comb1 = comb %>% filter(state =='Inbred') 
combH = comb %>% filter(state =='Hybrid') 
combI = left_join(comb1, ml, by= c('Tissue','Male')) %>% dplyr::rename(fmlExpGNumb = ExpGeneNumb, mlExpGNumb = n) %>% mutate(ExpGeneNumb = (fmlExpGNumb + mlExpGNumb)*0.5) %>%select(-fmlExpGNumb, -mlExpGNumb)

combF = rbind(combH, combI)



## summary expression genes/hybrid-inbred comparison

e1= left_join(combH,combI, by=c('Male','Female','Tissue')) %>% filter(!Female %in% c('Q381','PHG83')) %>% mutate(diff = ExpGeneNumb.x -ExpGeneNumb.y)

write.table(e1,'~/projects/biomap/data/paperRevise/expGeneNumb.txt', sep='\t', quote=FALSE)


