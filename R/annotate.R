#' Annotate the set of genes with the GO terms for a given species and assembly
#'
#' @param genes List of mRNA genes. Supported format for genes is Hugo.
#' @param GOtype Hierarchical category of the GO ontology. Possible values are 'BP', 'CC', 'MF'.
#' @param hg Genome assembly of interest. Possible assemblies are 'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit fly, 'ce11' for worm, 'hg19' and 'hg38' for human
#' @param slim Boolean value stating whether set of annotation should be performed for high level GO terms (GO slim)
#'
#' @return data frame of the GO term annotation of the genes
#'
#' @examples
#' \dontrun{
#' FOXP2annot<-annotate(genes = 'FOXP2', GOtype = 'BP', hg = 'hg19')
#' write.table(FOXP2annot,'foxp2.txt',row.names=FALSE)
#' }
#'
annotate <- function(genes,
                     GOtype,
                     hg,
                     slim = FALSE) {
  if (missing(genes)) {
    message("Genes are missing. Expected input: FOXP2 SOX4 HOXC6")
  }
  if (missing(GOtype)) {
    message("GO type is missing. Please select among 'BP', 'CC', 'MF'.")
  }
  if (missing(hg)) {
    message(
      "Genome assembly is missing. Possible assemblies are 'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit fly, 'ce11' for worm, 'hg19' and 'hg38' for human."
    )
    assembly(hg)
  }
  if (hg == 'mm10') {
    if (GOtype == "BP") {
      go <- NoRCE::mouse_p
    }
    if (GOtype == "CC") {
      go <- NoRCE::mouse_c
    }
    if (GOtype == "MF") {
      go <- NoRCE::mouse_f
    }
  }
  else if (hg == 'dre10') {
    if (GOtype == "BP") {
      go <- NoRCE::zebra_p
    }
    if (GOtype == "CC") {
      go <- NoRCE::zebra_c
    }
    if (GOtype == "MF") {
      go <- NoRCE::zebra_f
    }
  }
  else if (hg == 'rn6') {
    if (GOtype == "BP") {
      go <- NoRCE::rat_p
    }
    if (GOtype == "CC") {
      go <- NoRCE::rat_c
    }
    if (GOtype == "MF") {
      go <- NoRCE::rat_f
    }

  }
  else if (hg == 'dm6') {
    if (GOtype == "BP") {
      go <- NoRCE::fly_p
    }
    if (GOtype == "CC") {
      go <- NoRCE::fly_c
    }
    if (GOtype == "MF") {
      go <- NoRCE::fly_f
    }
  }
  else if (hg == 'ce11') {
    if (GOtype == "BP") {
      go <- NoRCE::worm_p
    }
    if (GOtype == "CC") {
      go <- NoRCE::worm_c
    }
    if (GOtype == "MF") {
      go <- NoRCE::worm_f
    }
  }
  else if (hg == 'sc3') {
    if (slim) {
      if (GOtype == "BP") {
        go <- NoRCE::goslim_yeast_p
      }
      if (GOtype == "CC") {
        go <- NoRCE::goslim_yeast_cc
      }
      if (GOtype == "MF") {
        go <- NoRCE::goslim_yeast_fp
      }
    }
    else{
      if (GOtype == "BP") {
        go <- NoRCE::yeast_p
      }
      if (GOtype == "CC") {
        go <- NoRCE::yeast_c
      }
      if (GOtype == "MF") {
        go <- NoRCE::yeast_f
      }
    }
  }
  else{
    if (slim) {
      if (GOtype == "BP") {
        go <- NoRCE::slimb
      }
      if (GOtype == "CC") {
        go <- NoRCE::slimc
      }
      if (GOtype == "MF") {
        go <- NoRCE::slimf
      }
    }
    else{
      if (GOtype == "BP") {
        go <- NoRCE::h_p
      }
      if (GOtype == "CC") {
        go <- NoRCE::h_c
      }
      if (GOtype == "MF") {
        go <- NoRCE::h_f
      }
    }
  }
  assign('go', go, envir = .GlobalEnv)
  geneMatrix <- as.matrix(genes)
  annot <-
    data.frame(GOID = character(),
               Gene = character(),
               GOTerm = character())
  for (i in 1:dim(geneMatrix)[1])
  {
    if (length(which(geneMatrix[i] == go$X2)) > 0)
    {
      t <- go[which(geneMatrix[i] == go$X2), ]
      tmp <- data.frame(t$X1, t$X2, t$X5)
      colnames(tmp) <- c('GOID', 'Gene', 'GOTerm')
      annot <- rbind(annot, tmp)
    }
  }
  return(annot)
}
