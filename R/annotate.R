#' Annotate the set of genes with the GO terms for a given species and 
#' assembly
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
#' 
#'
annGO <- function(genes,GOtype = c("BP", "CC", "MF"),
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
  
  types <-
    rbind(
      c("hg19","Hs.eg.db"), c("hg38","Hs.eg.db"), c("mm10","Mm.eg.db"),
      c("dre10","Dr.eg.db"), c("rn6","Rn.eg.db"), c("sc3","Sc.sgd.db"),
      c("dm6","Dm.eg.db"), c("ce11","Ce.eg.db") )
  index = which(org_assembly == types[, 1])
  
  genes <- data.frame(genes)
  goData <- AnnotationDbi::select(
    GO.db,
    keys = GOtype ,
    columns = c("GOID",
                "TERM",
                "ONTOLOGY"),
    keytype = "ONTOLOGY"
  )
  
  z <- paste0("org.",types[index,2])
  if ( !requireNamespace(z, quietly = TRUE))
    stop("Install package ",td," in order to use this function.")
  else
    lapply(z, require, character.only = TRUE)
  
  gogene <- AnnotationDbi::select(
         eval(as.name(z)),
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
