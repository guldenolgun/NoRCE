#' Differentially expressed human brain data
#' @format Not Available
#' @source \url{http://resource.psychencode.org/}
#' @examples
#' data(brain_mirna)
"brain_mirna"

#' Protein coding genes that are differentially expressed in TCGA
#' breast cancer RNAseq data.
#' @format Not Available
#' @source \url{https://portal.gdc.cancer.gov/}
#' @examples
#' data(breastmRNA)
"breastmRNA"

#' Differentially expressed non-coding gene regions
#' @format Not Available
#' @source \url{http://resource.psychencode.org/}
#' @examples
#' data(ncRegion)
"ncRegion"

#' Differentially expressed non-coding gene
#' @format Not Available
#' @source \url{http://resource.psychencode.org/}
#' @examples
#' data(brain_disorder_ncRNA)
"brain_disorder_ncRNA"

#' TAD regions for the fly
#' @format Not Available
#' @source \url{http://chorogenome.ie-freiburg.
#' mpg.de/data_sources.html#hi-c_datasets }
#' @examples
#' data(tad_dmel)
"tad_dmel"

#' TAD regions for human hg19 assembly
#' @format Not Available
#' @source \url{http://promoter.bx.psu.edu/hi-c/publications.html}
#' @examples
#' data(tad_hg19)
"tad_hg19"

#' TAD regions for human hg38 assembly
#' @format Not Available
#' @source \url{http://promoter.bx.psu.edu/hi-c/publications.html}
#' @examples
#' data(tad_hg38)
"tad_hg38"

#' TAD regions for mouse
#' @format Not Available
#' @source \url{http://promoter.bx.psu.edu/hi-c/publications.html}
#' @examples
#' data(tad_mm10)
"tad_mm10"

#' Brain miRNA expression retrieved from the TCGA
#' @format Not Available
#' @source \url{https://www.gencodegenes.org/}
#' @examples
#' data(mirna)
"mirna"

#' Brain mRNA expression retrieved from the TCGA
#' @format Not Available
#' @source \url{https://www.gencodegenes.org/}
#' @examples
#' data(mrna)
"mrna"

conn <- mirna_base <- feature<- count <- . <- NULL
Pvalue <- EnrichGeneNumber <- PAdjust <- geneLoc <- study <- NULL
globalVariables(c("tad_hg19", "tad_hg38", "tad_mm10", "tad_dmel"))
