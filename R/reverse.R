#Type of the input gene list. Possible values are "Entrez", "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol "NCBI" value, for Entrez gene id "Entrez", for mirbase id "mirna" is used.
commonGene <- function(mrnaobject, hg, downstream, upstream,inputGene, inGeneType){
  a<-unique(unlist(mrnaobject@geneList))
  aa<-convertGeneID(genelist = a,genetype = 'NCBI',hg = hg)

  a1 <- convertGeneID(genelist = inputGene,genetype = inGeneType,hg = hg)

  big_islands <-
    resize(a1, width = upstream + width(a1), fix = "start")
  hits <- findOverlaps(big_islands,aa, ignore.strand=TRUE)
  tmp <- big_islands[queryHits(hits),]
  tmp1 <-aa[subjectHits(hits),]
  pairss <-data.frame(tmp$gene,tmp1$gene)

  big_islands <-
    resize(a1, width = downstream + width(a1), fix = "end")
  hits <- findOverlaps(big_islands,aa, ignore.strand=TRUE)
  tmp <- big_islands[queryHits(hits),]
  tmp1 <-aa[subjectHits(hits),]

  pairss<-rbind(pairss, data.frame(tmp$gene,tmp1$gene))

  getNoncode<-function(x){a <- pairss[which(pairss[,2] %in% unlist(mrnaobject@geneList[[x]])),1]
  a[!duplicated(a)] }

  ab <- lapply(1:length(mrnaobject@Term), getNoncode)

  ab[isEmpty(ab)]<-'NA'
  if(length(mrnaobject@geneList)==1)
    ab<-list(ab)
  return(ab)
}

