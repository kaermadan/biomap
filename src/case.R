#!/user/bin/Rscript
source('/home/hirschc1/lixx5447/projects/biomap/src/df.R')
#Tissue specific
#demo
b7= ase %>% filter(Male=='B73')
mo= b7 %>% filter(Female=='Mo17')
mt = mo[,c(3,7,ncol(mo))]
mt1 = spread(mt, Tissue, Reg.Cat)
#mt2 = mt1[complete.cases(mt1),]
mt3 = mt1[complete.cases(mt1[,2:6]),]
#
female = unique(ase2$Female)
male = unique(ase2$Male)
otp <- NULL
tisOut<- NULL
tisg <- NULL
#try faster method for append result
#otp1 = list()
for (i in male ){
	for (j in female){
		gnp = ase2 %>% filter(Male ==i, Female == j) %>% select(gid, Tissue, Reg.Cat) %>% spread(., Tissue, Reg.Cat)
		if (nrow(gnp) > 0){
                        if (ncol(gnp) >5){
				#shg = gnp[complete.cases(gnp[,2:ncol(gnp)]),]
				#shg1 = gnp[complete.cases(gnp[,3:ncol(gnp)]),]
				#shg2 = gnp[complete.cases(gnp[,4:ncol(gnp)]),]
				#shg3 = gnp[complete.cases(gnp[,5:ncol(gnp)]),]
                                #tissue constitutive
                                shg50 = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na <1) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling) %>% filter(is.na(Reg.Cat) == FALSE) 
                                shg5 = shg50 %>% group_by(Tissue, Reg.Cat) %>%dplyr::count() %>% ungroup() %>% group_by(Tissue) %>%  mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop = round(n/sum*100,1)) %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop) %>% mutate(type = 'Tissue Constitutive')
				#intermedia express tissue
				shg40 = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na <2 & na >0) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling) %>% filter(is.na(Reg.Cat) == FALSE) 
                                shg4 = shg40 %>% group_by(Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 1)) %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop) %>% mutate(type = 'shared4')
				shg30 = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na <3 & na >1) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling) %>% filter(is.na(Reg.Cat) == FALSE) 
                                shg3 = shg30 %>% group_by(Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 1)) %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop) %>% mutate(type = 'shared3')
				shg20 = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na <4 & na >2) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling) %>% filter(is.na(Reg.Cat) == FALSE) 
                                shg2 = shg20 %>% group_by(Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 1)) %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop) %>% mutate(type = 'shared2')
                                #tissue specific
                                shg10 = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na >3) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling)  %>% filter(is.na(Reg.Cat) == FALSE) 
                                shg1 = shg10 %>% group_by(Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 1)) %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop) %>% mutate(type = 'Tissue Specific')
                                #
                                shgx0 = rbind(shg50, shg40, shg30, shg20, shg10) %>% mutate(id = paste0(j, 'x', i)) %>% select(id, Tissue, na, gid, Reg.Cat)
                                tisg = rbind(tisg, shgx0)
                                tsp = rbind(shg5, shg4, shg3, shg2, shg1)
                                tisOut = rbind(tisOut, tsp)
                                #calculate prop
				total = nrow(gnp)
				shared5 = nrow(shg5)
				shared4 = nrow(shg4)
				shared3 = nrow(shg3)
				shared2 = nrow(shg2) 
                                shared1 = nrow(shg1)
				prop5 = round(shared5/total*100,1)
				prop4 = round(shared4/total*100,1)
				prop3 = round(shared3/total*100,1)
				prop2 = round(shared2/total*100,1)
				prop1 = round(shared1/total*100,1)
				df = data.frame(Male = i, Female = j, nGene = total, shared5, prop5, shared4, prop4, shared3, prop3, shared2, prop2, shared1, prop1)
				otp = rbind(otp,df)
			}
		}
	}
}
#
give.n <- function(x){
  return(c(y = median(x)*1.00, label = length(x)))
}
#p01 = tisOut %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Tissue, y = prop, fill = Reg.Cat))+
#        geom_boxplot(notch=F, outlier.size = 0.080, width = 0.7, position = position_dodge(0.8), lwd = 0.1, fatten = 3) +
#    scale_x_discrete( name ='', expand = c(0,.5)) +
#    scale_y_continuous(name = 'Proportion (%)', limits = c(0,50), expand = c(0,1.05)) +
#    scale_fill_npg() +
#    #coord_flip() +
#    theme_bw() +
#    #labs(x= '', y='No. of PAV Genes')+
#    theme(plot.margin = unit(c(.0,.0,.0,0), "lines")) +
#    theme(axis.text.x= element_text(angle=0, size =10),
#          axis.text.y= element_text(size=10),
#          axis.title.y = element_text(size =10),
#           axis.title.x = element_text(size =10))+
#    theme(#legend.justification=c(0.1,0.05),
#           legend.position='top',
#           legend.text=element_text(size=8),
#           legend.direction = 'horizontal',
#           legend.key.size=unit(0.7,'lines'),
#           legend.margin = margin(0,0,0,0),
#           legend.box.margin= margin(0,0,-8,0))+
#    theme(panel.grid.major.x =  element_blank())+
#
#    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
#                  position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
#    guides(fill=guide_legend(title="",nrow=1))+
#    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
#        ggsave(filename = fp01,width =5,height=3.0)

#plot constittutive and tissue specific together
#test sig
cbd0 = tisOut %>% mutate(pro = prop*0.01)
cbd = cbd0 %>% filter(!type %in% c('shared4', 'shared3', 'shared2')) 
library(psych)
library(betareg)
library(emmeans)
library(multcompView)
tis = unique(cbd$Tissue)
reg = unique(cbd$Reg.Cat)
btp <- NULL
for (i in tis){
	for (j in reg){
		cb1 = cbd %>% filter(Tissue == i, Reg.Cat == j)
		model = betareg(pro ~ type, data = cb1)
		joint_tests(model)
		marginal = emmeans(model,
				   ~ type)
		cmp = pairs(marginal,
		      adjust="tukey")
		Sum = CLD(marginal,
			  alpha   = 0.05,
			  Letters = letters,         ###  Use lowercase letters for .group
			  adjust  = "tukey")         ###  Tukey-adjusted comparisons
		cmp %<>% data.frame(.) %>% mutate(Tissue = i, Reg.Cat = j)
		btp = rbind(btp, cmp)
	}
}
btp %<>% select(-contrast, -df) %>% arrange(p.value) %>% select(Reg.Cat, Tissue, p.value, estimate:z.ratio)
btp

#plot
fp01 = file.path(dir, paste0('tissue.constitutive-specific.reg.cat.box','.png'))
p01 = cbd %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.080, width = 0.5, position = position_dodge(0.8), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1)+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', limits = c(0,40), expand = c(0,1.05)) +
    scale_color_brewer(palette = 'Paired') +
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 10, hjust = 1),
          axis.text.y= element_text(size=10),
          axis.title.y = element_text(size =10),
           axis.title.x = element_text(size =10))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
    	  strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+

    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
        ggsave(filename = fp01,width = 4, height=6.0)

#try plot all together
cbd.tem = cbd0 %>% mutate(type = ifelse(type == 'Tissue Constitutive', 'shared5', type)) %>% mutate(type = ifelse(type == 'Tissue Specific', 'shared1', type))
cbd.tem %>% group_by(type, Reg.Cat) %>% summarise_at(vars(n), funs(max, min, mean, median, sd))
dk2 = brewer.pal(n =8, name ='Greys')
dk2 = dk2[4:8]
fp01 = file.path(dir, paste0('tissue.shared1-5.reg.cat.box','.png'))
p01 = cbd.tem %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.040, width = 0.5, position = position_dodge(0.68), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1)+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', limits = c(0,45), expand = c(0,1.05)) +
    scale_fill_manual(values = dk2) +
    #scale_color_brewer(palette = 'Paired')+
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0.1), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 9, hjust = 1),
          axis.text.y= element_text(size=9),
          axis.title.y = element_text(size =9),
           axis.title.x = element_text(size =9))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
    	  strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+
    #stat_summary(fun.data = give.n, geom = "text", fun.y = median,
     #             position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
        ggsave(filename = fp01,width = 4.5, height = 6)

#filter out n< 20 categories
cbd.tem %>% filter(n >=20 ) %>% group_by(type, Reg.Cat) %>% summarise_at(vars(n), funs(max, min, mean, median, sd)) %>% data.frame()
fp01 = file.path(dir, paste0('tissue20.shared1-5.reg.cat.box','.png'))
p01 = cbd.tem %>% filter( n >= 20) %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.040, width = 0.5, position = position_dodge(0.68), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1)+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', limits = c(0,45), expand = c(0,1.05)) +
    scale_fill_manual(values = dk2) +
    #scale_color_brewer(palette = 'Paired')+
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0.1), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 9, hjust = 1),
          axis.text.y= element_text(size=9),
          axis.title.y = element_text(size =9),
           axis.title.x = element_text(size =9))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
          strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+
    #stat_summary(fun.data = give.n, geom = "text", fun.y = median,
     #             position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
        ggsave(filename = fp01,width = 4.5, height = 6)

###Genotype specific
#demo
leaf = ase2 %>% filter(Tissue == 'Leaf')
leaf1 = leaf[,c(1:3,7,ncol(leaf))] %>% mutate(id = paste0(Female,'x', Male))  
leaf2 = leaf1 %>% select(-Male, -Female, -Tissue) %>% spread(., id, Reg.Cat)
leaf3 = leaf2 %>% mutate(na = rowSums(is.na(.))) 
##
tis = unique(ase2$Tissue)
geo <- NULL
geog <- NULL
for (i in tis){
	gsc = ase2 %>% filter(Tissue == i) 
	gsc1 = gsc[, c(1:3, 7, ncol(gsc))] %>% mutate(id = paste0(Female, 'x', Male)) 
	gsc2 = gsc1 %>% select(-Male, -Female, -Tissue) %>% spread(., id, Reg.Cat) %>% mutate(na = rowSums(is.na(.)))  
	#how many genes are expressed in how many genotypes
	nog = gsc1 %>% group_by(gid) %>% dplyr::count() %>% ungroup() %>% group_by(n) %>% dplyr::count()
	#80%
	gsc80 = gsc2 %>% filter(na <= (ncol(.)-2)*0.2) %>% gather(., id, Reg.Cat, contains('x')) %>% filter(is.na(Reg.Cat) == FALSE) 
	gsc801 = gsc80 %>% group_by(id, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(id) %>%  mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop = round(n/sum*100,1)) %>% mutate(Tissue = i) %>% mutate(type = 'Genotype Constitutive') 
	#20-80%
	gsc50 = gsc2 %>% filter(na > (ncol(.)-2)*0.2, na < (ncol(.)-2)*0.8) %>% gather(., id, Reg.Cat, contains('x')) %>% filter(is.na(Reg.Cat) == FALSE) 
	gsc501 = gsc50 %>% group_by(id, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(id) %>%  mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop = round(n/sum*100,1)) %>% mutate(Tissue = i) %>% mutate(type = 'Intermedia Frequency') 
	#20%
	gsc20 = gsc2 %>% filter(na >= (ncol(.)-2)*0.8) %>% gather(., id, Reg.Cat, contains('x')) %>% filter(is.na(Reg.Cat) == FALSE) 
	gsc201 = gsc20 %>% group_by(id, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(id) %>%  mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop = round(n/sum*100,1)) %>% mutate(Tissue = i) %>% mutate(type = 'Genotype Specific') 
        #
        gtemg = rbind(gsc80, gsc50, gsc20) %>% mutate(Tissue = i) %>% select(id, Tissue, na, gid, Reg.Cat)
        geog = rbind(geog, gtemg)
        gtem = rbind(gsc801, gsc501, gsc201)
        geo = rbind(geo, gtem)
}
geo %<>% mutate(pro = prop*0.01)

#sig test
tis = unique(geo$Tissue)
reg = unique(geo$Reg.Cat)
gsig<- NULL
for (i in tis){
        for (j in reg){
                cb1 = geo %>% filter(type != 'Intermedia Frequency', pro >0) %>% filter(Tissue == i, Reg.Cat == j)
                model = betareg(pro ~ type, data = cb1)
                joint_tests(model)
                marginal = emmeans(model,
                                   ~ type)
                cmp = pairs(marginal,
                      adjust="tukey")
                Sum = CLD(marginal,
                          alpha   = 0.05,
                          Letters = letters,         ###  Use lowercase letters for .group
                          adjust  = "tukey")         ###  Tukey-adjusted comparisons
                cmp %<>% data.frame(.) %>% mutate(Tissue = i, Reg.Cat = j)
                gsig = rbind(gsig, cmp)
        }
}
gsig %<>% select(-contrast, -df) %>% arrange(p.value) %>% select(Reg.Cat, Tissue, p.value, estimate:z.ratio)
gsig

#plot genotype constitutive and genotype specific
geo28 = geo %>% filter(type != 'Intermedia Frequency') 
fp01 = file.path(dir, paste0('genotype.constitutive-specific.reg.cat.box','.png'))
p01 = geo28 %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.080, width = 0.5, position = position_dodge(0.8), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1)+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', expand = c(0,1.05)) +
    scale_color_brewer(palette = 'Paired') +
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 10, hjust = 1),
          axis.text.y= element_text(size=10),
          axis.title.y = element_text(size =10),
           axis.title.x = element_text(size =10))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
          strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+
    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
    ggsave(filename = fp01,width = 4, height=6.0)

#plot genotype specific all together
geo$type <- factor(geo$type, levels = c('Genotype Specific', 'Intermedia Frequency','Genotype Constitutive'), labels = c('shared20', 'shared20-80', 'shared80'))
dk2 = brewer.pal(n =8, name ='Greys')
dk2 = dk2[4:8]
fp01 = file.path(dir, paste0('genotype.specific.reg.cat.box','.png'))
p01 = geo %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.040, width = 0.5, position = position_dodge(0.68), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1, scale = 'free_y')+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', expand = c(0,1.05)) +
    scale_fill_manual(values = dk2) +
    #scale_color_brewer(palette = 'Paired')+
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0.1), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 9, hjust = 1),
          axis.text.y= element_text(size=9),
          axis.title.y = element_text(size =9),
           axis.title.x = element_text(size =9))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
    	  strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+
    #stat_summary(fun.data = give.n, geom = "text", fun.y = median,
     #             position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
    ggsave(filename = fp01,width = 4.5, height=6.50)

#filter out n< 20 samples
geo20 = geo %>% filter(n >= 20)
geo20  %>% group_by(type, Reg.Cat) %>% summarise_at(vars(n), funs(max, min, mean, median, sd)) %>% data.frame()
fp01 = file.path(dir, paste0('genotype.specific_20.reg.cat.box','.png'))
p01 = geo20 %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = type))+
        geom_boxplot(notch=F, outlier.size = 0.040, width = 0.5, position = position_dodge(0.68), lwd = 0.1, fatten = 3) +
    facet_wrap(~Tissue, ncol=1, scale = 'free_y')+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', expand = c(0,1.05)) +
    scale_fill_manual(values = dk2) +
    #scale_color_brewer(palette = 'Paired')+
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0.1), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 9, hjust = 1),
          axis.text.y= element_text(size=9),
          axis.title.y = element_text(size =9),
           axis.title.x = element_text(size =9))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
          strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank())+
    #stat_summary(fun.data = give.n, geom = "text", fun.y = median,
     #             position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
    ggsave(filename = fp01,width = 4.5, height=6.50)

#explore overlap of genotype-specific and genotype specific
#demo
#constitutive
dgc = geog %>% filter(na <= 8, id =='Mo17xB73', Tissue=='Leaf')
dtc= tisg %>% filter(na == 0, id =='Mo17xB73', Tissue =='Leaf')
dolpc = inner_join(dgc, dtc, by=c('gid')) 
dchk = length(intersect(dgc$gid, dtc$gid))
if (dchk == nrow(dolpc)){
	cat('ok, inner_join and intersect return the same result')
} else {
	cat('warning!!, inner_join and intersect results are inconsistent!')
}
#specific
dgs = geog %>% filter(na >= 35, id =='Mo17xB73', Tissue=='Leaf')
dts= tisg %>% filter(na == 4, id =='Mo17xB73', Tissue =='Leaf')
dolps = inner_join(dgs, dts, by=c('gid'))
dchk = length(intersect(dgs$gid, dts$gid))
if (dchk == nrow(dolps)){
        cat('ok, inner_join and intersect return the same result')
} else {
        cat('warning!!, inner_join and intersect results are inconsistent!')
}

#run
tis = unique(geog$Tissue)
genotype = unique(tisg$id)
geneList <-NULL
smyResult <- NULL
warnings <- NULL
for (i in tis){
	gtis = geog %>% filter(Tissue == i)
	for (j in genotype){
		ntis = tisg %>% filter(id == j) 

		#constitutive
		rgc = gtis %>% filter(id == j, na <= (max(gtis$na)+1)*0.2)
		rtc = ntis %>% filter(Tissue == i, na < (max(ntis$na)+1)*0.2) 
		rolpc = inner_join(rgc, rtc, by = c('gid', 'Tissue', 'Reg.Cat','id')) %>% mutate(state = 'Constitutive_overlap')
		rchkc = length(intersect(rgc$gid, rtc$gid))
		if (rchkc == nrow(rolpc)){
			cmsg = 'ok, inner_join and intersect return the same result'
		} else {
			cmsg = 'warning!!, inner_join and intersect results are inconsistent!'
		}
		rgc_ratio = round(nrow(rolpc)/nrow(rgc)*100, 2) 
		rtc_ratio = round(nrow(rolpc)/nrow(rtc)*100, 2)

		#specific
		rgs = gtis %>% filter(id == j, na >= (max(gtis$na)+1)*0.8)
		rts = ntis %>% filter(Tissue == i, na >= (max(ntis$na)+1)*0.8)
		rolps = inner_join(rgs, rts, by = c('gid', 'Tissue', 'Reg.Cat','id')) %>% mutate(state = 'Specific_overlap')
		rchks = length(intersect(rgs$gid, rts$gid))
		if (rchks == nrow(rolps)){
			smsg = 'ok, inner_join and intersect return the same result'
		} else {
			smsg = 'warning!!, inner_join and intersect results are inconsistent!'
		}
		rgs_ratio = round(nrow(rolps)/nrow(rgs)*100, 2)
		rts_ratio = round(nrow(rolps)/nrow(rts)*100,2)
		#
		dfOut = data.frame(id = j, Tissue = i, geno_constitutive = nrow(rgc), tis_constitutive = nrow(rtc), overlap_c = nrow(rolpc), prop_gc = rgc_ratio, prop_tc = rtc_ratio, geno_specific = nrow(rgs), tis_specific = nrow(rts), overlap_s = nrow(rolps), prop_gs = rgs_ratio, prop_ts = rts_ratio) 

		smyResult = rbind(smyResult, dfOut)
		#
		msg = data.frame(id = j, Tissue = i, state_c = cmsg, state_s = smsg) 
		warnings = rbind(warnings, msg)
		#
		geneOut = rbind(rolpc, rolps)
		geneList = rbind(geneOut, geneList)
	}
}

#
csShared = geneList %>% group_by(id, state,Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue, id, state) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 2))
cs1 = csShared %>% group_by(Tissue, state, Reg.Cat) %>% summarise_at(vars(n), funs(N=n(), max,  min, mean, median, sd)) %>% data.frame() %>% mutate_at(vars(mean, sd), funs(round(.,1)))
s1= cs1 %>% filter(state=='Specific_overlap', Reg.Cat != 'Ambiguous') %>% arrange(Reg.Cat)
c1= cs1 %>% filter(state=='Constitutive_overlap', Reg.Cat != 'Ambiguous') %>% arrange(Reg.Cat)

#plot
fp01 = file.path(dir, paste0('shared_constitutive.specific.reg.cat.box','.png'))
p01 = csShared %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Reg.Cat, y = prop, fill = state))+
        geom_boxplot(notch=F, outlier.size = 0.040, width = 0.6, position = position_dodge(0.68), lwd = 0.1, fatten = 3) +
    facet_grid(state ~ Tissue)+
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', expand = c(0,1.05)) +
    #scale_fill_manual(values = dk2) +
    #scale_color_brewer(palette = 'Paired')+
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.1,.0,0.1), "lines")) +
    theme(axis.text.x= element_text(angle=45, size = 9, hjust = 1),
          axis.text.y= element_text(size=9),
          axis.title.y = element_text(size =9),
           axis.title.x = element_text(size =9))+
    theme(#legend.justification=c(0.1,0.05),
           legend.position='top',
           legend.text=element_text(size=8),
           legend.direction = 'horizontal',
           legend.key.size=unit(0.7,'lines'),
           legend.margin = margin(0,0,0,0),
           legend.box.margin= margin(0,0,-8,0))+
    theme(strip.background = element_rect(colour = 'gray88', fill = 'gray88', size=0.45),
          strip.text.x = element_text(size=9, color='black', margin = margin(0.13,0,0.13,0, "cm")))+
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank())+
    #stat_summary(fun.data = give.n, geom = "text", fun.y = median,
     #             position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    #geom_vline(xintercept = c(1.5,2.5,3.5,4.5, 5.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
    ggsave(filename = fp01,width = 5.2, height=4.500)






