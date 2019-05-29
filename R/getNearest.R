op <- options(warn = (-1))
options(readr.num_columns = 0)

#' Convert txt file or data frame to bed format. First three columns should be chr number, start and end, respectively.
#'
#' @param dm_file Bed formated txt file or data frame
#' @param isText Boolean value that holds whether input data is txt file. If it is TRUE, input file has to be txt.
#'
#' @return Bed file
#'
#' @importFrom readr read_table
#'
#' @examples
#'
#' #txt formatted data
#' data("ncRegion")
#' \dontrun{
#'
#' #Write data to a txt file
#' write.table(ncRegion,paste("ncRegion.txt"),sep = '\t',row.names = FALSE,col.names = FALSE)
#'
#' #Convert txt file that just created to bed formatted data
#' regionNC <-readbed("ncRegion.txt")
#'
#' #Directly convert data frame to bed format
#' regionNC <- readbed(dm_file = ncRegion,isText = FALSE)
#'
#' }
#' @export
readbed <- function(dm_file, isText=TRUE) {
  if (missing(dm_file)) {
    message("Please provide path of the bed file")
    dm_file <- readline("New path: ")
    readbed(dm_file)
  }
  if(isText){
  data <- read.table(dm_file, header = FALSE)
  colnames(data) <- c('chr', 'start', 'end')
  bedfile <- with(data, GRanges(chr, IRanges(start, end)))
  }
  else
  {
    colnames(dm_file)<- c('chr', 'start', 'end')
    bedfile <- with(dm_file, GRanges(chr, IRanges(start, end)))
  }
  return(bedfile)
}

#' @importFrom TxDb.Celegans.UCSC.ce11.refGene TxDb.Celegans.UCSC.ce11.refGene
#' @importFrom TxDb.Dmelanogaster.UCSC.dm6.ensGene TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' @importFrom TxDb.Drerio.UCSC.danRer10.refGene TxDb.Drerio.UCSC.danRer10.refGene
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
#' @importFrom TxDb.Rnorvegicus.UCSC.rn6.refGene TxDb.Rnorvegicus.UCSC.rn6.refGene
#' @importFrom biomaRt getBM useEnsembl useMart
assembly <- function(hg) {
  if (hg == "hg19") {
    data <- NoRCE::ucsc_hg19
    colnames(data) <-
      c('gene', 'chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Hsapiens.UCSC.hg19.knownGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart = useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      host = "grch37.ensembl.org",
      path = "/biomart/martservice",
      dataset = "hsapiens_gene_ensembl"
    )
    assign('mart', mart, envir = .GlobalEnv)
  }
  if (hg == "hg38")
  {
    data <- NoRCE::ucsc_hg38
    colnames(data) <-
      c('gene', 'chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Hsapiens.UCSC.hg38.knownGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)
  }
  if (hg == 'mm10') {
    #mm10
    data <- NoRCE::ucsc_mm10
    colnames(data) <-
      c('gene', 'chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Mmusculus.UCSC.mm10.knownGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)
  }
  if (hg == 'dre10') {
    #zebrafish

    data <- NoRCE::ucsc_dre10
    colnames(data) <-
      c('gene', 'chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Drerio.UCSC.danRer10.refGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useMart(host = "useast.ensembl.org",
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "drerio_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)

  }
  if (hg == 'rn6') {
    #rat

    data <- NoRCE::ucsc_rn6
    colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Rnorvegicus.UCSC.rn6.refGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)
  }
  if (hg == 'sc3') {
    data <- NoRCE::sc3
    colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useMart("ensembl", dataset="scerevisiae_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)

  }
  if (hg == 'dm6') {
    #Fruit fly
    data <- NoRCE::ucsc_dm6
    colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useMart("ensembl",
              dataset = "dmelanogaster_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)
  }
  if (hg == 'ce11') {
    #Worm
    data <- NoRCE::ucsc_ce11
    colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
    ucsc <-
      with(data, GRanges(chr, IRanges(start, end), strand, symbol))
    assign('ucsc', ucsc, envir = .GlobalEnv)

    genomee <- TxDb.Celegans.UCSC.ce11.refGene
    assign('genomee', genomee, envir = .GlobalEnv)

    mart <-
      useMart(host = "useast.ensembl.org",
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "celegans_gene_ensembl")
    assign('mart', mart, envir = .GlobalEnv)
  }
}

#' Get nearest genes for the window of the upstream/downstream region.
#'
#' When downstream = 0 / upstream = 0, function converts bed formated regions to HUGO genes
#'
#'
#' @param bedfile Bed formated input gene regions
#' @param upstream Maximum upstream distance from the transcription start region of the input gene
#' @param downstream Maximum downstream distance from the transcription end region of the input gene
#' @param hg genomee assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#'
#'
#' @return genes
#'
#' @import GenomicFeatures
#' @importFrom GenomicRanges as.data.frame duplicated end findOverlaps intersect match merge order pintersect resize split start strand width
#'
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#'  \dontrun{
#'
#' neighbour <- getUCSC(bedfile = regionNC, upstream = 1000, downstream = 1000,hg = 'hg19')
#'
#'}
#'@export
getUCSC <-
  function(bedfile, upstream, downstream, hg) {
    if (!exists("mart") | !exists("ucsc")) {
      assembly(hg)
    }
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }
    if (missing(upstream)) {
      message("Upstream information is missing?")
    }
    if (missing(downstream)) {
      message("Downstream information is missing?")
    }
    if (missing(hg)) {
      message(
        "genomee assembly version is missing. Possible assemblies are 'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit fly, 'ce11' for worm, 'hg19' and 'hg38' for human."
      )
    }

    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    rt1 <- subsetByOverlaps(ucsc, unstrand(big_islands))

    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    rt2 <- subsetByOverlaps(ucsc, unstrand(big_islands))
    results2 <- rt2$symbol
    results1 <- rt1$symbol
    result <- as.data.frame(unique(c(results1, results2)))
    return(result)
  }

#' Get only those neighbouring genes that fall within exon region
#'
#' @param bedfile Input bed formated file
#' @param upstream Maximum upstream distance from the transcription start position
#' @param downstream Maximum downstream distance from the transcription end position
#' @param hg genomee assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#'
#' @return genes
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' \dontrun{
#' r<-getNearToExon(bedfile = regionNC,upstream = 1000, downstream = 2000,hg = 'hg19')
#' }
#'
#' @export
getNearToExon <-
  function(bedfile, upstream, downstream, hg) {
    if (!exists("mart") | !exists("ucsc") | !exists("genomee")) {
      assembly(hg)
    }
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }
    if (missing(upstream)) {
      message("Upstream information is missing?")
    }
    if (missing(downstream)) {
      message("Downstream information is missing?")
    }

    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    exons = subsetByOverlaps(exons(genomee), big_islands)
    region1 <- subsetByOverlaps(ucsc, exons)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    exons = subsetByOverlaps(exons(genomee), big_islands)
    region2 <- subsetByOverlaps(ucsc, exons)
    results2 <- region2$symbol
    results1 <- region1$symbol
    result <- as.data.frame(unique(c(results1, results2)))
    return(result)
  }

#' Get only those neighbouring genes that fall within intron region
#'
#' @param bedfile Bed file
#' @param upstream upstream distance
#' @param downstream downstream distance
#' @param hg genomee assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#'
#' @return genes
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#'
#' \dontrun{
#'
#' r<-getNearToExon(bedfile = regionNC,upstream = 1000, downstream = 2000,hg = 'hg19')
#'}
#'
#' @importFrom GenomicFeatures as.list intronsByTranscript
#'
#' @export
getNearToIntron <-
  function(bedfile, upstream, downstream, hg) {
    if (missing(hg)) {
      message("Assembly is missing?")
    }
    if (!exists("ucsc") ) {
      assembly(hg)
    }
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }
    if (missing(upstream)) {
      message("Upstream information is missing?")
    }
    if (missing(downstream)) {
      message("Downstream information is missing?")
    }
    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    intron = subsetByOverlaps(intronsByTranscript(genomee), big_islands)

    region1 <- subsetByOverlaps(ucsc, intron)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    intron = subsetByOverlaps(intronsByTranscript(genomee), bedfile)
    region2 <- subsetByOverlaps(ucsc, intron)
    results2 <- region2$symbol
    results1 <- region1$symbol
    result <- as.data.frame(unique(c(results1, results2)))
    return(result)
  }

#' For given region of interest, overlapped genes in the TAD regions are found. Results can be filtered according to the available cell lines.
#'
#' @param bedfile Region of interest
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or any new TAD regions can be used for the analysis. TAD regions must be formated as GRanges object. Predefined TAD regions are \'tad_hg19\', \'tad_hg38\', \'tad_mm10\', \'tad_dmel\' for hg19, hg38, mm9 and dm6 assembly, respectively.
#' @param cellline Cell lines for TAD regions.
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be considered in the analysis
#' @param upstream Holds upstream distance from the transcription start position
#' @param downstream Holds downstream distance from the transcription end position
#'
#' @return List of protein coding genes that falls into the TAD regions
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' \dontrun{
#' r<-getTADOverlap(bedfile = regionNC,TAD = tad_hg38, hg = 'hg19',cellline = 'HUVEC')
#' }
#'
#' @export
getTADOverlap <-
  function(bedfile,
           hg,
           TAD = NoRCE::tad_hg19, near = FALSE,  upstream = 10000,
           downstream = 10000,
           cellline = 'all') {
    if (missing(hg)) {
      message("Assembly is missing?")
    }
    if (!exists("ucsc") ) {
      assembly(hg)
    }
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }

    if (cellline == 'all') {
      tad <- TAD
    }
    else{
      temp <- which(TAD$celline == cellline)
      tad <-  TAD[temp, ]
    }

    if(near){
      treg <- Reduce(subsetByOverlaps, list(tad, bedfile, resize(bedfile, width = downstream + width(bedfile), fix = "end")))
      region2 <- subsetByOverlaps(ucsc, treg)
      treg1 <- Reduce(subsetByOverlaps, list(tad, bedfile, resize(bedfile, width = upstream + width(bedfile), fix = "start")))
      region1 <- subsetByOverlaps(ucsc, treg1)
      region <- data.frame(symbol = unique(c(region1$symbol,region2$symbol)))
    }
     else
     {
       tadRegion = subsetByOverlaps(bedfile, tad)
       region <- subsetByOverlaps(ucsc, tadRegion)
     }

    return(unique(as.data.frame(region$symbol)))
  }

#' Convert gene ids according to the gene type
#'
#' @param genetype Type of the input gene list. Possible values are "Entrez", "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol "NCBI" value, for Entrez gene id "Entrez", for mirbase id "mirna" is used.
#' @param genelist Input gene list
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#'
#' @return GRange object of the given input
#'
#' @importFrom biomaRt getBM
#'
#' @examples
#'
#' convGene <-convertGeneID(genetype = "mirna", genelist = brain_mirna[1:100,],hg = 'hg19')
#'
#'
#' @export
convertGeneID <- function(genetype, genelist, hg) {
  if (missing(hg)) {
    message(
      "genomee assembly version is missing. Possible assemblies are 'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit fly, 'ce11' for worm, 'hg19' and 'hg38' for human."
    )
  }
  if (missing(genetype)) {
    message(
      "Format of the gene is missing. Possible values are 'Entrez', 'mirna', 'Ensembl_gene', 'Ensembl_trans', 'NCBI'."
    )
  }
  if (missing(genelist)) {
    message("List of gene is missing.")
  }
  if (!exists("ucsc") ) {
    assembly(hg)
  }
  attributes <-
    c("chromosome_name",
      "start_position",
      "end_position",
      "strand")
  if (genetype == "Entrez") {
    output <-
      getBM(
        attributes = c( "entrezgene", attributes),
        filters = "entrezgene",
        values = genelist,
        mart = mart
      )
  }

  else if (genetype == "mirna") {
    output <-
      getBM(
        attributes = c("mirbase_id",attributes),
        filters = "mirbase_id",
        values = apply(as.data.frame(genelist),2,tolower),
        mart = mart
      )
  }
  else if (genetype == "Ensembl_gene") {
    output <-
      getBM(
        attributes = c("ensembl_gene_id", attributes),
        filters = "ensembl_gene_id",
        values = genelist,
        mart = mart
      )
  }
  else if (genetype == "Ensembl_trans") {
    output <-
      getBM(
        attributes = c("ensembl_transcript_id",attributes),
        filters = "ensembl_transcript_id",
        values = genelist,
        mart = mart
      )
  }
  else if (genetype == "NCBI") {
    output <-
      getBM(
        attributes = c("hgnc_symbol",attributes),
        filters = "hgnc_symbol",
        values = genelist,
        mart = mart
      )
  }
  else if (genetype == "mgi_symbol")
  {
    output <-
      getBM(
        attributes = c("mgi_symbol",attributes),
        filters = "mgi_symbol",
        values = genelist,
        mart = mart
      )
  }
  else if (genetype == "external_gene_name")
  {
    output <-
      getBM(
        attributes = c("external_gene_name",attributes),
        filters = "external_gene_name",
        values = genelist,
        mart = mart
      )
  }

  colnames(output)<- c("gene","chromosome_name",
                       "start_position",
                       "end_position",
                       "strand")

  if(hg == 'hg19' | hg == 'hg38'){
    geneList<-as.data.frame(genelist)
    a <- gencode %>% dplyr::filter(gencode$gene %in% geneList[,1])
    #a <- cbind(gencode %>% dplyr::filter(gencode$gene %in% geneList[,1]),'*')
    if(dim(a)[1] > 0){
      a<-cbind(a,'*')
      colnames(a) <-colnames(output)
      output <- unique(rbind(output,a))
      colnames(a) <-colnames(output)
      output <- unique(rbind(output,a))
    }
  }

  file1 <-
    with(output, GRanges(
      paste0("chr", chromosome_name),
      IRanges(start_position, end_position),'*',gene
    ))


  return(file1)
}

#' List cell line of the given topological domain regions
#'
#' @param TADName input TAD regions
#'
#' @return cell line of the input tad data
#'
#' @examples
#'
#' listTAD(TADName = tad_hg19)
#'
#'
#' @export
#'
listTAD<-function(TADName){
  return(unique(TADName$celline))
}
