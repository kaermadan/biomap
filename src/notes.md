#
look at how many genes with ase data are shared across tissue in single triplet

#next
how about across different triplet?
#tissue specific
#just to make clear:
For 32 triplets that have data for all five tissue, regulatory divergence analysis can be conducted for 9143-18973 genes which are covered by SNP markers. For each tissue of each triplet, the total number of genes in it can be first divided into five types: shared5-shared1, Among them, 169-2574 genes are spcificlly expressed in only one tissue, and 3748-9634 genes are constitutively expressed in all five tissues. Then for each type, it can be further classified into 7 regulatory categories: Cis, Trans....
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

#ideas
- Look at known genes in our dataset, like tb1, ZmCCT, ZCN8, gt1, etc., to check if they 
  have the same regulation pattern as here.
- Look at regulation pattern of homeolog gene pairs, see [this paper](/https://onlinelibrary.wiley.com/doi/epdf/10.1111/tpj.14228)

