#' Get the biotype of the non-coding genes. It is suitable for the
#' GENCODE gtf files
#'
#' @param gtfFile Path of the input gtf file which contains biotype
#'     information. The gtf file must be provided from the Ensembl
#'     or Gencode site. For space efficiency, gft files should be
#'     in a zip format.
#'
#' @return Tabular form of the gtf file with the required features
#'   such as gene id and biotypes
#'
#' @import readr
#'
#' @examples
#' fileImport<-system.file("extdata", "temp.gtf", package = "NoRCE")
#' gtf <- extractBiotype(gtfFile = fileImport)
#'
#' @export
#'
extractBiotype <- function(gtfFile) {
  mydata <- read_table(gtfFile, comment = '#', col_names = FALSE)
  temp <- strsplit(as.character(mydata$X1), "[;\t]+")
  
  #r<-lapply(temp,function(x) grepl("type|transcript_id|gene_id",x))
  r <- lapply(temp, function(x)
    grepl("gene_type|gene_id", x))
  r <- lapply(r, as.logical)
  
  # a <- list()
  # for (i in seq_along(r)) {
  #   tr <- which(r[[i]] == 'TRUE')
  #   a[[i]] <- temp[[i]][tr]
  # }
  #
  # mat <- t(lapply(a,
  #                 function(x, m)
  #                   c(x, rep(NA, m - length(
  #                     x
  #                   ))),
  #                 max(rapply(a, length))))
  # mat<-matrix(unlist(lapply(mat, `[[`, 1)))
  
  a <- data.frame()
  for (i in seq_along(r)) {
    tr <- which(r[[i]] == 'TRUE')
    a[i, seq_len(2)] <- temp[[i]][tr]
  }
  
  #gtf <- gsub("^.* ", "", mat, perl = TRUE)
  gtf <- apply(a, 2, function(x)
    (gsub("^.* ", "", x)))
  
  gtf <- as.data.frame(gsub("\"", "", gtf))
  return(gtf)
}

#' Extract the genes that have user provided biotypes. This method is useful
#' when input gene list is mixed or when research of the interest is only
#' focused on specific group of genes.
#'
#' @param gtfFile Input gtf file for the genes provided by the extractBiotype
#'     function
#' @param biotypes Selected biotypes for the genes
#'
#' @return Table format of genes with a given biotypes
#'
#' @examples
#' biotypes <- c('unprocessed_pseudogene','transcribed_unprocessed_pseudogene')
#' fileImport<-system.file("extdata", "temp.gtf", package = "NoRCE")
#' extrResult <- filterBiotype(fileImport, biotypes)
#'
#' @export
filterBiotype <- function(gtfFile, biotypes) {
  gtf <- extractBiotype(gtfFile = gtfFile)
  all <- data.frame(gene = character(), stringsAsFactors = FALSE)
  for (i in seq_along(biotypes)) {
    index <- which(gtf == biotypes[i], arr.ind = TRUE)
    all <- rbind(all, as.data.frame(gtf[index[, 1], 1]))
  }
  all <-
    as.data.frame(
      unlist(apply(unique(all), 2, strsplit, '[.]'))[c(TRUE, FALSE)])
  colnames(all) <- 'gene'
  return(all)
}
