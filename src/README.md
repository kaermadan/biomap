#
look at how many genes with ase data are shared across tissue in single triplet

#next
how about across different triplet?
#tissue specific
For 32 triplets, 9143-18973 genes are covered by SNPs and can be xxx five tissues. Among them, xxx-xxx genes are spcificlly expressed in only one tissue, and xxx-xxx genes are constitutively expressed in all five tissue
cbd.tem %>% group_by(type) %>% summarise_at(vars(sum),funs(max,min,mean,sd))
# A tibble: 5 x 5
  type      max   min  mean    sd
  <chr>   <dbl> <dbl> <dbl> <dbl>
1 shared1  2574   169  524.  397.
2 shared2  1690    94  713.  334.
3 shared3  2531   130 1181.  534.
4 shared4  5746   469 2092.  817.
5 shared5  9634  3748 8087. 1459. 
#just to make clear, for each tissue of each triplet, the total number of genes in it can be first divided into five types: shared5-shared1, then for each type, it can be further divided into 7 Reg.Cat: Cis, Trans....  

#Genotype specific
for leaf tissue, in total, there are 22886 genes with SNP across 43 triplets. Among them, 2954 genes are expressed in over 80% of the genotypes (35): 1077-2855 genes in different genotypes. `gsc80 %>% group_by(id) %>% count()`
Such genes in each genotype can be further classified into seven regulatory categories: Cis, Trans..., plot can be made from the distribution of numbers and proportions of them genes in each triplet in each tissue.
