#!/user/bin/Rscripr

source('/home/hirschc1/lixx5447/projects/src/pkgs.R')

ht = read_excel('/panfs/roc/groups/14/hirschc1/lixx5447/PAV/graph/paper/heterotic_group.xlsx')
ht = ht[,c(1:2)]
colnames(ht)<- c('geno','heterotic')


rmhtp = grep('*UK*', unique(q12$htp), value=T)

q12 %<>% filter(!htp %in% rmhtp)

q12 %>% group_by(htp, Tissue) %>% dplyr::count()

q12$htp<- factor(q12$htp, levels=c("SSSxNSS", "NSSxSSS","IODENTxSSS","SSSxIODENT","IODENTxNSS","NSSxIODENT","SSSxSSS", "NSSxNSS","IODENTxIODENT"), labels=c("SSSxNSS", "NSSxSSS","IODENTxSSS","SSSxIODENT","IODENTxNSS","NSSxIODENT","SSSxSSS", "NSSxNSS","IODENTxIODENT"))

q12$Tissue<- factor(q12$Tissue, levels=c('Endosperm','Leaf','Root','Seedling','Internode'), labels=c('Endosperm','Leaf','Root','Seedling','Internode'))
