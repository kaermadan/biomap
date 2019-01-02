R data file (/home/hirschc1/lixx5447/projects/biomap/data/10.norm.RData.) containing raw read count tables, normalized expression values and allele-specific read counts:

tm - tibble for biomap expression data
gid: Gene ID (AGP_v4, Ensembl Plants v37, 46,117 in total)
SampleID: bm001 - bm467
ReadCount: raw read count
nRC: normalized read count (nRC = ReadCount / sizeFactor)
rCPM: raw CPM (adds up to 1,000,000 for each sample/library)
rFPKM: raw FPKM calculated using rCPM and gene exon length
CPM: CPM calculated by edgeR (CPM = rCPM / normFactor)
FPKM: FPKM calculated using CPM and gene exon length

tl - tibble for library (sample), with columns:
SampleID: bm001 - bm467
Tissue: Leaf, Internode, Root, etc.
Genotype: B73, Mo17xPH207, etc
Treatment: replicate 1 or 2
inbred: whether this is the inbred parent (TRUE) or hybrid (FALSE)
sizeFactor, libSize: library size and normalization factor calculated using the median log ratio approach by DESeq2, accounts for library size
normFactor: library normalzation factor computed by edgeR using the TMM approach, does NOT account for library size

#Important notes: Q381 is equal to PH207, PHG83 is equal to B73; based on previously PCA results, LH93xB73 -Root(bm252) should be removed; bm466 is a triple-crossing hybrid.
