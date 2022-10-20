commonGene <- function(mrnaobject,
                       org_assembly,
                       downstream,
                       upstream,
                       inputGene,
                       inGeneType) {
  a <- unique(unlist(mrnaobject@geneList))
  aa <-
    convertGeneID(genelist = a,
                  genetype = 'NCBI',
                  org_assembly = org_assembly)
  
  a1 <- convertGeneID(genelist = inputGene,
                      genetype = inGeneType,
                      org_assembly = org_assembly)
  
  big_islands <-
    resize(a1, width = upstream + width(a1), fix = "start")
  hits <- findOverlaps(big_islands, aa, ignore.strand = TRUE)
  tmp <- big_islands[S4Vectors::queryHits(hits), ]
  tmp1 <- aa[S4Vectors::subjectHits(hits), ]
  pairss <- data.frame(tmp$gene, tmp1$gene)
  
  big_islands <-
    resize(a1, width = downstream + width(a1), fix = "end")
  hits <- findOverlaps(big_islands, aa, ignore.strand = TRUE)
  tmp <- big_islands[S4Vectors::queryHits(hits), ]
  tmp1 <- aa[S4Vectors::subjectHits(hits), ]
  
  pairss <- unique(rbind(pairss, data.frame(tmp$gene, tmp1$gene)))
  
  getNoncode <- function(x) {
    a <-
      pairss[which(pairss[, 2] %in% unlist(mrnaobject@geneList[[x]])), 1]
    a[!duplicated(a)]
  }
  
  ab <- lapply(seq_along(mrnaobject@Term), getNoncode)
  
  ab[IRanges::isEmpty(ab)] <- 'NA'
  if (length(mrnaobject@geneList) == 1)
    ab <- list(ab)
  return(ab)
}

#' @importFrom GenomeInfoDb seqlevels
commonGeneRegion <- function(mrnaobject,
                             org_assembly,
                             downstream,
                             upstream,
                             inRegion) {
  a <- unique(unlist(mrnaobject@geneList))
  aa <-
    convertGeneID(genelist = a,
                  genetype = 'NCBI',
                  org_assembly = org_assembly)
  
  regions <- paste(sub('chr','',seqlevels(inRegion)), start(inRegion), end(inRegion), sep=":")
  if(org_assembly == 'hg19' | org_assembly == 'hg38'){
    results <- getBM(attributes = c( "chromosome_name", "start_position","end_position", 'strand', "hgnc_symbol"),
                     filters = c("chromosomal_region"),
                     values=regions,
                     mart=pkg.env$mart)
  }else{
    results <- getBM(attributes = c( "chromosome_name", "start_position","end_position", 'strand', "external_gene_name"),
                     filters = c("chromosomal_region"),
                     values=regions,
                     mart=pkg.env$mart)
  }
  
  colnames(results)[5] = 'hgnc_symbol'
  
  for(i in 1 : dim(results)){
    if(results$hgnc_symbol[i] == ''){
      results$hgnc_symbol[i] = paste0('chr',results$chromosome_name, ':', results$start_position,'-', results$end_position)
    }}
  
  file1 <-
    with(results, GRanges(
      paste0("chr", chromosome_name),
      IRanges::IRanges(start_position, end_position),
      strand,
      hgnc_symbol
    ))
  
  big_islands <-
    resize(file1, width = upstream + width(inRegion), fix = "start")
  hits <- findOverlaps(big_islands, aa, ignore.strand = TRUE)
  tmp <- big_islands[S4Vectors::queryHits(hits), ]
  tmp1 <- aa[S4Vectors::subjectHits(hits), ]
  pairss <- data.frame(tmp$hgnc_symbol, tmp1$gene)
  
  big_islands <-
    resize(file1, width = downstream + width(inRegion), fix = "end")
  hits <- findOverlaps(big_islands, aa, ignore.strand = TRUE)
  tmp <- big_islands[S4Vectors::queryHits(hits), ]
  tmp1 <- aa[S4Vectors::subjectHits(hits), ]
  
  pairss <- unique(rbind(pairss, data.frame(tmp$hgnc_symbol, tmp1$gene)))
  
  getNoncode <- function(x) {
    a <-
      pairss[which(pairss[, 2] %in% unlist(mrnaobject@geneList[[x]])), 1]
    a[!duplicated(a)]
  }
  
  ab <- lapply(seq_along(mrnaobject@Term), getNoncode)
  
  ab[IRanges::isEmpty(ab)] <- 'NA'
  if (length(mrnaobject@geneList) == 1)
    ab <- list(ab)
  return(ab)
}
