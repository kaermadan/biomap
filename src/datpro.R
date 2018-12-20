#!/user/bin/Rscript

source('/home/hirschc1/lixx5447/projects/src/pkgs.R')

dir= '/home/hirschc1/lixx5447/projects/biomap/graph'
ase= readRDS('/home/hirschc1/lixx5447/projects/biomap/data/ase_newData181214.rds')
#ase= readRDS('/home/hirschc1/lixx5447/projects/biomap/data/ase_newData181218_filter20.rds')
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
            theme(axis.text.x= element_text( size =10),
                  axis.text.y= element_text(size=10),
                  axis.title.y= element_text(size=10))+
            theme(legend.key.size = unit(0.8,'lines'),
                  legend.position = c(0.7,0.94),
                  legend.direction = 'horizontal',
                  legend.text = element_text(size=8),
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,-8,0))+
            theme(panel.grid.major = element_blank())+
            guides(fill = guide_legend(title = ''))+
            stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2.0, vjust= -2.1, hjust=0.5)+ 
            ggsave(filename = fp01, width=7, height=6)

#facet tissues
fp01 = file.path(dir, paste0('reg.cat.het.tissues','.png'))
p1= q12 %>% ggplot(aes(x=Reg.Cat, y=prop, fill=htp))+ 
            #geom_bar(stat='identity',position= position_dodge(width=0.85), width=0.7) +
            geom_boxplot(notch = F, outlier.size = 0.13, width = .65,lwd=0.12, fatten=1.0, position=position_dodge(width=0.75))+
            facet_wrap(~Tissue, ncol=1)+ 
            scale_x_discrete(name = '', expand=c(0,0.5))+
            scale_y_continuous(name = 'Proportion (%)',limits = c(0,45), expand = c(0,1.05))+ 
            scale_fill_manual(values = dk3)+
            theme_bw() +
            theme(plot.margin = unit(c(.0,.0,.0,0), "lines")) +
            theme(axis.text.x= element_text( size =10),
                  axis.text.y= element_text(size=10),
                  axis.title.y= element_text(size=10))+
            theme(legend.key.size = unit(0.8,'lines'),
                  legend.position = c(0.7,0.98),
                  legend.direction = 'horizontal',
                  legend.text = element_text(size=6),
                  legend.margin = margin(0,0,0,0),
                  legend.box.margin = margin(0,0,-8,0))+
            theme(panel.grid.major = element_blank())+
            guides(fill = guide_legend(title = ''))+
            stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2.0, vjust= -2.1, hjust=0.5)+ 
            ggsave(filename = fp01, width=7, height=10)

# ase in genotypes across tissue
#convert tibble to matrix
matrix.please <- function(x){
   x<- as.data.frame(x)
   m<- as.matrix(x[,-1])
   rownames(m)<- x[,1]
   m
}

tis = as.character(unique(q12$Tissue))
cln = unique(q12$Reg.Cat)[c(2,7,3,4,5,6)]
for ( i in tis ) { 
	endo = q12 %>% ungroup() %>% filter(Tissue == i, Reg.Cat != 'Ambiguous' ) %>% mutate(Genotype= paste0(Female, 'x', Male)) %>%  select(Genotype, prop, Reg.Cat) %>% spread(., Reg.Cat, prop)

	e1 = matrix.please(endo)
	#mycol = colorRampPalette(brewer.pal(8, "PiYG"))(25)
	mycol = terrain.colors(12)
	range_1<- extendrange(e1, r = range(e1, na.rm = TRUE), f = 0.01)
        #bre=c(seq(range_1[1], range_1[2],length=15))# change length=? to achieve the best effect
        bre = c(seq(0,35, length = 13))
	png(filename = paste0(dir, '/heatmap.endosperm.reg.cat_', i, '.png'), height = 5, width = 4, res = 500,units = 'in')
	p2 = heatmap.2(e1, Rowv = FALSE, Colv = cln, col = mycol, 
	  dendrogram = 'none',
	  key = T, trace = 'none',
	  margins = c(4.6,4.6), 
	  cexRow = 0.5, cexCol = 0.8,
	  na.rm = TRUE, na.color = 'gray23',
	  key.title = i,
	  key.xlab = '',
	  tracecol = 'gray31',
	  symm = F,symkey = F,symbreaks = T,
	  offsetRow = -0.3, offsetCol = -0.3, 
	  scale = "none", breaks = bre,
          keysize=1.6) 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          #key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          #lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))
	 dev.off()
}


