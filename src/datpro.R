#!/user/bin/Rscript

source('/home/hirschc1/lixx5447/projects/src/pkgs.R')

dir= '/home/hirschc1/lixx5447/projects/biomap/graph'
ase= readRDS('/home/hirschc1/lixx5447/projects/biomap/data/ase_newData181214.rds')
ase %<>% as.tibble()

a1= ase %>% group_by(Reg.Cat,Female, Male,Tissue) %>% dplyr::count() 
a2= a1 %>% ungroup() %>% group_by(Male, Female, Tissue) %>% mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop=n/sum*100) %>% filter(sum> 5000) 

a2 %>% group_by(Reg.Cat) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))
deg= c('Cis','Trans','CisxTrans','Cis+Trans')

z2= a2 %>% filter(Reg.Cat %in% deg)
z2 %>% group_by(Reg.Cat) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))
z2 %>% group_by(Reg.Cat, Tissue) %>% summarise_at(vars(prop), funs(mean, max,min,median,sd))

#heterotic
ht = read_excel('/panfs/roc/groups/14/hirschc1/lixx5447/PAV/graph/paper/heterotic_group.xlsx')
ht = ht[,c(1:2)]
colnames(ht)<- c('geno','heterotic')
# 
q12= left_join(a2, ht, by=c('Female'= 'geno')) %>% left_join(.,ht, by=c('Male'= 'geno')) %>% mutate(htp = paste(heterotic.x, 'x', heterotic.y, sep=''))
rmhtp = grep('*UK*', unique(q12$htp), value=T)
q12 %<>% filter(!htp %in% rmhtp)
q12 %>% group_by(htp, Tissue) %>% dplyr::count()
q12 %>% group_by(Reg.Cat,htp) %>% summarise_at(vars(prop),funs(mean,max,min,median,sd))

q12$htp<- factor(q12$htp, levels=c("SSSxNSS", "NSSxSSS","IODENTxSSS","SSSxIODENT","IODENTxNSS","NSSxIODENT","SSSxSSS", "NSSxNSS","IODENTxIODENT"), labels=c("SSSxNSS", "NSSxSSS","IODENTxSSS","SSSxIODENT","IODENTxNSS","NSSxIODENT","SSSxSSS", "NSSxNSS","IODENTxIODENT"))
q12$Tissue<- factor(q12$Tissue, levels=c('Endosperm','Leaf','Root','Seedling','Internode'), labels=c('Endosperm','Leaf','Root','Seedling','Internode'))

fp01= file.path(dir, paste('reg.cat.htp.prop.bar','.png', sep=''))
dk2 = brewer.pal(n = 6, name = "Paired")
add3= c('gray40','gray60','gray80')
dk3 = c(dk2, add3)
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x)))  
}
p1= q12 %>% ggplot(aes(x=Reg.Cat, y=prop, fill=htp))+ 
            #geom_bar(stat='identity',position= position_dodge(width=0.85), width=0.7) +
            geom_boxplot(notch = F, outlier.size = 0.13, width = .65,lwd=0.12, fatten=1.0, position=position_dodge(width=0.75))+ 
            scale_x_discrete(name = '', expand=c(0,0.5))+
            scale_y_continuous(name = 'Proportion (%)',limits = c(0,45), expand = c(0,1.05))+ 
            scale_fill_manual(values = dk3)+
            theme_bw() +
            theme(plot.margin = unit(c(.0,.0,.0,0), "lines")) +
            theme(legend.key.size = unit(0.8,'lines'),
                  legend.position = c(0.7,0.94),
                  legend.direction = 'horizontal',
                  legend.text = element_text(size=8),
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,-8,0))+
            theme(panel.grid.major = element_blank())+
            guides(fill = guide_legend(title = ''))+
            stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2.2, vjust= -2.1, hjust=0.5)+ 
            ggsave(filename = fp01, width=9, height=6)


