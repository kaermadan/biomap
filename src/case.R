#!/user/bin/Rscript
source('/home/hirschc1/lixx5447/projects/biomap/src/df.R')

#demo
b7= ase %>% filter(Male=='B73')
mo= b7 %>% filter(Female=='Mo17')
mt = mo[,c(3,7,ncol(mo))]
mt1 = spread(mt, Tissue, Reg.Cat)
mt2 = mt1[complete.cases(mt1),]
mt3 = mt1[complete.cases(mt1[,2:6]),]

#
female = unique(ase2$Female)
male = unique(ase2$Male)
otp <- NULL
smy <- NULL
tisOut<- NULL
#try faster method for append result
#otp1 = list()
for (i in male ){
	for (j in female){
		gnp = ase2 %>% filter(Male ==i, Female == j) %>% select(gid, Tissue, Reg.Cat) %>% spread(., Tissue, Reg.Cat)
		if (nrow(gnp) > 0){
                        if (ncol(gnp) >5){
				shg = gnp[complete.cases(gnp[,2:ncol(gnp)]),]
				shg1 = gnp[complete.cases(gnp[,3:ncol(gnp)]),]
				shg2 = gnp[complete.cases(gnp[,4:ncol(gnp)]),]
				shg3 = gnp[complete.cases(gnp[,5:ncol(gnp)]),]
				total = nrow(gnp)
				shared5 = nrow(shg)
				shared4 = nrow(shg1)
				shared3 = nrow(shg2)
				shared2 = nrow(shg3) 
				ntis = ncol(gnp)-1
				prop5 = round(shared5/total*100,1)
				prop4 = round(shared4/total*100,1)
				prop3 = round(shared3/total*100,1)
				prop2 = round(shared2/total*100,1)
				df = data.frame(Male = i, Female = j, nGene = total, nTissue = ntis, shared5, prop5, shared4, prop4, shared3, prop3, shared2,prop2)
				otp = rbind(otp,df)
                                
                                smy1 = shg %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling) %>% group_by(Tissue, Reg.Cat) %>%dplyr::count() %>% ungroup() %>% group_by(Tissue) %>%  mutate_at(vars(n), funs(sum=sum)) %>% mutate(prop = round(n/sum*100,1))
                                smy2 = smy1 %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop)
                                smy = rbind(smy, smy2)

                                #tissue specific
                                tsg = gnp %>% mutate(na = rowSums(is.na(.))) %>% filter(na >3) %>% gather(., Tissue, Reg.Cat, Endosperm:Seedling)
				tsp = tsg %>% filter(is.na(Reg.Cat) == FALSE) %>% group_by(Tissue, Reg.Cat) %>% dplyr::count() %>% ungroup() %>% group_by(Tissue) %>% mutate_at(vars(n), funs(sum = sum)) %>% mutate(prop = round(n/sum*100, 1))
                                tsp1 = tsp %>% mutate(Male = i, Female = j) %>% select(Male, Female, Tissue:prop)
                                tisOut = rbind(tisOut, tsp1)
			}
		}
	}
}
otp %<>% arrange(nTissue)
#try ## = do.call(rbind, otp1) 
#resul#t %<>% arrange(Proportion)
#or re#sult = dplyr::bind_rows(otp1)
#or re#sult = data.table::rbindlist(otp1)

fp01 = file.path(dir, paste0('shared5.tis.reg.cat.box','.png'))
give.n <- function(x){
  return(c(y = median(x)*1.00, label = length(x)))  
}
p01 = smy %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Tissue, y = prop, fill = Reg.Cat))+
        geom_boxplot(notch=F, outlier.size = 0.080, width = 0.7, position = position_dodge(0.8), lwd = 0.1, fatten = 3) +
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', limits = c(0,50), expand = c(0,1.05)) +
    scale_fill_npg() +
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+ 
    theme(plot.margin = unit(c(.0,.0,.0,0), "lines")) +
    theme(axis.text.x= element_text(angle=0, size =10),
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
    theme(panel.grid.major.x =  element_blank())+
  
    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
        ggsave(filename = fp01,width =5,height=3.0)

#tissue specific
fp01 = file.path(dir, paste0('tissue.specific.reg.cat.box','.png'))
give.n <- function(x){
  return(c(y = median(x)*1.00, label = length(x)))
}
p01 = tisOut %>% filter( Reg.Cat != 'Ambiguous') %>% ggplot(aes(x = Tissue, y = prop, fill = Reg.Cat))+
        geom_boxplot(notch=F, outlier.size = 0.080, width = 0.7, position = position_dodge(0.8), lwd = 0.1, fatten = 3) +
    scale_x_discrete( name ='', expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (%)', limits = c(0,50), expand = c(0,1.05)) +
    scale_fill_npg() +
    #coord_flip() +
    theme_bw() +
    #labs(x= '', y='No. of PAV Genes')+
    theme(plot.margin = unit(c(.0,.0,.0,0), "lines")) +
    theme(axis.text.x= element_text(angle=0, size =10),
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
    theme(panel.grid.major.x =  element_blank())+

    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                  position = position_dodge(width = 0.75), size = 2,vjust= -2.47, hjust=0.5)+
    guides(fill=guide_legend(title="",nrow=1))+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype= 'dashed',color = 'gray38', size = 0.1)+
        ggsave(filename = fp01,width =5,height=3.0)

#plot constittutive and tissue specific together
tisOut %<>% mutate(type = 'Tissue Constitutive')
smy %<>% mutate(type = 'Tissue Specific')
cbd = rbind(tisOut, smy)
#test sig
cbd %<>% mutate(pro = prop*0.01)
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
		Sum = cld(marginal,
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
        ggsave(filename = fp01,width = 4, height=7.0)





