## ----style, echo=FALSE, results="asis", message=FALSE, warnings = FALSE----
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE,fig.width=6, fig.height=6)

## ----Load package, eval=TRUE, echo=TRUE, include=TRUE--------------------
library(NoRCE)

## ---- Enrichment analysis based on gene neighbourhood when input set is consist of genes, eval=TRUE, echo = FALSE----
#GO enrichment results based on mRNAs that exon regions fall into the upstream and downstream of the input non-coding genes
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene')

#GO enrichment results based on mRNAs that any region falls into the upstream and downstream of the input miRNA genes
mirGO<-mirnaGOEnricher(gene = brain_mirna, hg='hg19', near=TRUE,downstream = 20000, upstream = 500, searchRegion = "all")

## ----Convert txt file to bed formatted data, eval=TRUE-------------------
#txt formatted data
data("ncRegion")

#Write data to a txt file
write.table(ncRegion,paste("ncRegion.txt"),sep = '\t',row.names = FALSE,col.names = FALSE)

#Convert txt file that just created to bed formatted data  
regionNC<-readbed("ncRegion.txt")


## ----Enrichment analysis based on gene neighbourhood when input set is consist of regions, eval=TRUE----
data("ncRegion")

#Convert dataframe to bed formatted data
regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)


#Perform enrichment based on neighourhood of non-coding RNA region
regionGO<-geneRegionGOEnricher(region = regionNC, hg= 'hg38', near = TRUE)



## ---- Enrichment of targeted nearest genes for the input, eval=TRUE------
#Intersection of the nearest genes of the input gene set and the potential target set is carries out for enrichment analysis
mirGO<-mirnaGOEnricher(gene = brain_mirna, hg='hg19', near=TRUE, target=TRUE)


## ---- Enrichment based on TAD, eval=TRUE---------------------------------
#Intersection of the nearest genes of the input gene set, the potential target set and overlap of the genes that are in TAD regions and the input list is carries out for enrichment analysis
mirGO<-mirnaGOEnricher(gene = brain_mirna, hg='hg19', near=TRUE, target=TRUE, isTADSearch = TRUE)

#Only genes that are targeted by the input list is put through the enrichment analysis
mirGO<-mirnaGOEnricher(gene = brain_mirna, hg='hg19', near=FALSE, target=FALSE, isTADSearch = TRUE)


## ---- Enrichment based on TAD cellline, eval=TRUE------------------------
#Retrieve list of cell-line 
a<-listTAD(TADName = tad_hg19)

#Intersection of the nearest genes of the input gene set, the potential target set and overlap of the genes that are in TAD regions for A549, AdrenalGland cell-line and the input list is carries out for enrichment analysis
mirGO<-mirnaGOEnricher(gene = brain_mirna, hg='hg19', near=TRUE, target=TRUE, isTADSearch = TRUE, TAD = tad_hg19, cellline = c('A549', 'AdrenalGland'))


## ---- GO enrichment based on custom TAD, eval=TRUE-----------------------

#Bed formatted txt file that contains TAD regions
cus_TAD<-system.file("extdata", "DER-18_TAD_adultbrain.txt", package = "NoRCE")

#Convert TAD regions to bed format
tad_custom <- readbed(cus_TAD, isText = TRUE)

ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, isTADSearch = TRUE, TAD = tad_custom)


## ---- Enrichment based on correlation analysis for the predefined correlations, eval=FALSE----
#  ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', express = TRUE,cancer = 'BRCA',minAbsCor = 0.3,databaseFile = 'C://Users/gulden olgun/Downloads/miRCancer.db/miRCancer.db')
#  

## ---- Enrichment based on correlation analysis, eval=FALSE---------------
#  nc4 <- mirnaGOEnricher(gene = brain_mirna, hg = 'hg19',target=TRUE, express = TRUE, exp1 = brain_mirnathres, exp2 = brainmrna,label1 = brain_mirnagene, label2 = brainmrna_genes, isCustomExp = TRUE, minAbsCor = 0.3)
#  

## ---- Enrichment based on correlation analysis without headers, eval=TRUE----
nc1 <- mirnaGOEnricher(gene = brain_mirna, hg = 'hg19',target=TRUE, express = TRUE, exp1 = mirna, exp2 = mrna, isCustomExp = TRUE, minAbsCor = 0.1)

## ---- Pathway gene enrichment, eval=TRUE---------------------------------
#Pathway enrichment based on the gen sets that falls into the TAD regions
ncRNAPathway<-genePathwayEnricher(gene = brain_disorder_ncRNA[1:500,], hg='hg19', pathwayType = 'kegg', isTADSearch = TRUE,TAD = tad_hg19, genetype = 'Ensembl_gene', min=2)

## ---- Pathway enrichment using custom database, eval=FALSE---------------
#  #Pathway enrichment using custom database
#  nc2 <- genePathwayEnricher(gene = brain_disorder_ncRNA, hg = 'hg19', pathwayType = 'other', near = TRUE, genetype = 'Ensembl_gene',  min = 2, gmtName = 'Human_AllPathways_March_01_2019_entrezgene.gmt', isSymbol = FALSE)
#  

## ---- Gene Enrichment Analysis, eval=FALSE-------------------------------
#  a <-mirnaPathwayEnricher(gene = brain_mirna, hg = "hg19", near = TRUE, gmtName = brain_gene, isSymbol = TRUE, isGeneEnrich = TRUE, pathwayType = "other")

## ----Biotype filtering, eval=FALSE, echo=TRUE----------------------------
#  #Load differentially expressed non-coding genes
#  data("brain_disorder_ncRNA")
#  
#  biotypes <- c('lincRNA','antisense')
#  
#  #It is suggested to use Zip format of the GENCODE gtf file
#  lncGenes <-filterBiotype(gtfFile = "pathTotheGTFFile//gtf", biotypes = biotypes)
#  
#  #All of the lincRNA and antisense genes are carried out to the GO enrichment analysis
#  lncEnrich<-geneGOEnricher(gene = unique(lncGenes), hg = "hg19", genetype = "Ensembl_gene")
#  
#  #Find lincRNA and antisense lncRNAs in the given input list and perform enrichment
#  lncEnrich<-geneGOEnricher(gene = intersect(lncGenes$gene,brain_disorder_ncRNA$V1), hg = "hg19", genetype = "Ensembl_gene")
#  

## ----Correlation Analysis for the given miRNA, eval=FALSE, echo=TRUE-----
#  #Load miRNA gene set
#  data("brain_mirna")
#  
#  database<- "filePathOfTheDatabase//miRCancer.db"
#  
#  brcaCorr<- corrbased(mirnagene = brain_mirna, cancer = 'BRCA', minAbsCor = 0.3, database)
#  

## ----Correlation analysis for the given mRNA, eval=FALSE, echo=TRUE------
#  data("breastmRNA")
#  
#  database<- "filePathOfTheDatabase//miRCancer.db"
#  
#  breast_corr<-corrbasedMrna(mRNAgene = as.data.frame(breastmRNA),cancer = 'BRCA', minAbsCor = 0.3, database)
#  
#  corr_miRNAEnrichment<-mirnaGOEnricher(gene=breast_corr$mirna_base,hg='hg38',near = TRUE,target = FALSE)

## ----Custom correlation analysis, eval=TRUE, echo=TRUE-------------------
#Assume that mirnanorce and mrnanorce are custom patient by gene data

dataCor <- calculateCorr(exp1 =  mirna, exp2 = mrna)


## ----Write tabular format of the results, eval=TRUE----------------------
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, pAdjust = "none")
writeEnrichment(mrnaObject = ncGO, fileName = 'enrichment.txt', sept = '\t',type = 'pvalue')

#Top 10 enrichment results are written
writeEnrichment(mrnaObject = ncGO, fileName = 'enrichment.txt', sept = '\t',type = 'pvalue', n=10)


## ---- Draw dot plot for the GO enrichment, eval=TRUE, echo=TRUE----------
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, pAdjust = "none")

drawDotPlot(mrnaObject = ncGO, type = 'pvalue', n = 20)

## ----GO:mRNA network, eval=TRUE, echo=TRUE-------------------------------
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, pAdjust = "none")

createNetwork(mrnaObject = ncGO, type = 'pvalue', n =3)

createNetwork(mrnaObject = ncGO, type = 'pvalue', n = 3, isNonCode = FALSE,takeID = TRUE)

## ----GO:ncRNA network, eval=TRUE, echo=TRUE------------------------------
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, pAdjust = "none")

createNetwork(mrnaObject = ncGO, type = 'pvalue', n = 3, isNonCode = TRUE,takeID = TRUE)

## ----GO-term DAG, eval=TRUE----------------------------------------------
ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE, genetype = 'Ensembl_gene', searchRegion = "exon", slim = FALSE, pAdjust = "none")

getGoDag(mrnaObject = ncGO, type = 'pvalue', n =3, filename = 'gonetwork.png',imageFormat = 'png')

## ----Create KEGG pathway map, eval=FALSE---------------------------------
#  ncRNAPathway<-mirnaPathwayEnricher(gene = brain_mirna[1:100,], hg='hg19',target = TRUE)
#  
#  #For the enriched pathway hsa04742, pathway map is created
#  getKeggDiagram(mrnaObject = ncRNAPathway, hg = 'hg19',pathway = ncRNAPathway@ID[1])
#  

## ----Create Reactome pathway map, eval=FALSE-----------------------------
#  miPath<-mirnaPathwayEnricher(gene=brain_mirna,hg='hg38',target=TRUE, pathwayType = 'reactome')
#  
#  getReactomeDiagram(mrnaObject = miPath, pathway = miPath@ID[1], imageFormat = 'png')
#  

