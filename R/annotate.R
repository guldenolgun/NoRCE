#' Annotate the set of genes with the GO terms for a given species and assembly
#'
#' @param genes List of mRNA genes. Supported format for genes is Hugo.
#' @param GOtype Hierarchical category of the GO ontology. Possible
#'     values are 'BP', 'CC', 'MF'.
#' @param org_assembly Genome assembly of interest. Possible assemblies are
#'     'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit
#'     fly, 'ce11' for worm, 'hg19' and 'hg38' for human
#'
#' @return data frame of the GO term annotation of the genes
#'
#' @import GO.db
#' @import org.Hs.eg.db
#'
annotate <- function(genes,GOtype = c("BP", "CC", "MF"),
                     org_assembly = c("hg19", "hg38", "mm10", "dre10", "rn6",
                                      "dm6", "ce11", "sc3")) {
  if (missing(genes)) {
    message("Genes are missing. Expected input: FOXP2 SOX4 HOXC6")
  }
  if (missing(GOtype)) {
    message("GO type is missing. Please select among 'BP', 'CC', 'MF'.")
  }
  if (missing(org_assembly)) {
    message("Genome assembly is missing.")
    assembly(org_assembly)
  }
  
  genes <- data.frame(genes)
  goData <- AnnotationDbi::select(
    GO.db,
    keys = GOtype ,
    columns = c("GOID",
                "TERM",
                "ONTOLOGY"),
    keytype = "ONTOLOGY"
  )
  gogene <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = goData$GOID,
    columns = c("GO", "SYMBOL"),
    keytype = "GO"
  )
  annot <- merge(goData, gogene, by.x = "GOID", by.y = "GO")
  
  
  annot <- unique(annot[, c(1, 6, 3)])
  colnames(annot) <- c('GOID', 'Gene', 'GOTerm')
  a <- annot[which(annot$Gene %in% genes[, 1]), ]
  
  return(list(annot, a))
}
