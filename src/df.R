#!/user/bin/Rscript

source('/home/hirschc1/lixx5447/projects/src/pkgs.R')

dir= '/home/hirschc1/lixx5447/projects/biomap/graph'
ase= readRDS('/home/hirschc1/lixx5447/projects/biomap/data/ase_newData181214.rds')
#ase= readRDS('/home/hirschc1/lixx5447/projects/biomap/data/ase_newData181218_filter20.rds')
ase %<>% as.tibble() 
ase1 = ase %>% group_by(Male, Female, Tissue) %>% dplyr::count() %>% filter(n > 5000) %>% mutate(tmp =paste0(Female, Male, Tissue))
ase2 = ase %>% mutate(tmp = paste0(Female, Male, Tissue)) %>% filter(tmp %in% ase1$tmp) %>% select(-tmp)
#
a1= ase %>% group_by(Reg.Cat,Female, Male,Tissue) %>% dplyr::count()
a2= a1 %>% ungroup() %>% group_by(Male, Female, Tissue) %>% mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop=n/sum*100) %>% filter(sum> 5000)

a2 %>% group_by(Reg.Cat) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))
deg= c('Cis','Trans','CisxTrans','Cis+Trans')

z2= a2 %>% filter(Reg.Cat %in% deg)
z2 %>% group_by(Reg.Cat) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))
z2 %>% group_by(Reg.Cat, Tissue) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))
