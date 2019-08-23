op <- options(warn = (-1))
options(readr.num_columns = 0)

pkg.env <- new.env()
pkg.env$ucsc <- GRanges(c(
  seqnames = NULL,
  ranges = NULL,
  strand = NULL
))
pkg.env$genomee <- data.frame()
pkg.env$mart <- data.frame()

#' Convert txt file or data frame to bed format. First three columns should
#' be chr number, start and end, respectively.
#'
#' @param dm_file Bed formated txt file or data frame
#' @param isText Boolean value that holds whether input data is txt file. If
#'    it is TRUE, input file has to be txt.
#'
#' @return Bed file
#'
#' @importFrom readr read_table
#'
#' @examples
#'
#' #txt formatted data
#' data("ncRegion")
#'
#' #Directly convert data frame to bed format
#' regionNC <- readbed(dm_file = ncRegion,
#'                     isText = FALSE)
#'
#' @export
readbed <- function(dm_file, isText = TRUE) {
  if (missing(dm_file)) {
    message("Please provide path of the bed file")
    dm_file <- readline("New path: ")
    readbed(dm_file)
  }
  if (isText) {
    data <- read.table(dm_file, header = FALSE)
    colnames(data) <- c('chr', 'start', 'end')
    bedfile <- with(data, GRanges(chr, IRanges(start, end)))
  }
  else
  {
    colnames(dm_file) <- c('chr', 'start', 'end')
    bedfile <- with(dm_file, GRanges(chr, IRanges(start, end)))
  }
  return(bedfile)
}

#' @importFrom biomaRt getBM useEnsembl useMart
assembly <- function(org_assembly = c("hg19",
                                      "hg38",
                                      "mm10",
                                      "dre10",
                                      "rn6",
                                      "dm6",
                                      "ce11",
                                      "sc3")) {
  myses <- browserSession()
  
  if (org_assembly == "hg19") {
    if (length(packageCheck(c(
      "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"
    ))) == 0) {
      genome(myses) <- "hg19"
      data <-
        getTable(ucscTableQuery(myses, track = "wgEncodeGencodeV31lift37"))
      data <- data[, c(3, 4, 5, 6, 13)]
      colnames(data) <-
        c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Hsapiens.UCSC.hg19.knownGene
      pkg.env$genomee <- genomee
      
      mart = useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = "grch37.ensembl.org",
        path = "/biomart/martservice",
        dataset = "hsapiens_gene_ensembl"
      )
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
                    )))
    }
  }
  
  if (org_assembly == "hg38")
  {
    if (length(packageCheck(c(
      "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"
    ))) == 0) {
      genome(myses) <- "hg38"
      data <-
        getTable(ucscTableQuery(myses, track = "wgEncodeGencodeV31"))
      data <- data[, c(3, 4, 5, 6, 13)]
      colnames(data) <-
        c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Hsapiens.UCSC.hg38.knownGene
      pkg.env$genomee <- genomee
      
      mart <-
        useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
                    )))
    }
  }
  if (org_assembly == 'mm10') {
    if (length(packageCheck(c(
      "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"
    ))) == 0) {
      genome(myses) <- "mm10"
      data <-
        getTable(ucscTableQuery(myses, track = "wgEncodeGencodeVM22"))
      data <- data[, c(3, 4, 5, 6, 13)]
      colnames(data) <-
        c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Mmusculus.UCSC.mm10.knownGene
      pkg.env$genomee <- genomee
      
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db")
                    )))
    }
  }
  if (org_assembly == 'dre10') {
    #zebrafish
    if (length(packageCheck(c(
      "TxDb.Drerio.UCSC.danRer10.refGene", "org.Dr.eg.db"
    ))) == 0) {
      genome(myses) <- "danRer11"
      a <- getTable(ucscTableQuery(myses, track = "ensGene"))
      a1 <-
        getTable(ucscTableQuery(myses, track = "ensGene",
                                table = "ensemblToGeneName"))
      data <- merge(a1, a)
      data <- data[, c(4, 5, 6, 7, 2)]
      
      colnames(data) <-
        c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Drerio.UCSC.danRer10.refGene
      pkg.env$genomee <- genomee
      
      mart <-
        useMart(host = "useast.ensembl.org",
                biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "drerio_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Drerio.UCSC.danRer10.refGene", "org.Dr.eg.db")
                    )))
    }
  }
  if (org_assembly == 'rn6') {
    #rat
    if (length(packageCheck(c(
      "TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db"
    ))) == 0) {
      genome(myses) <- "rn6"
      a <- getTable(ucscTableQuery(myses, track = "ensGene"))
      a1 <-
        getTable(ucscTableQuery(myses, track = "ensGene",
                                table = "ensemblToGeneName"))
      data <- merge(a1, a)
      data <- data[, c(4, 5, 6, 7, 2)]
      
      colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Rnorvegicus.UCSC.rn6.refGene
      pkg.env$genomee <- genomee
      
      mart <-
        useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db")
                    )))
    }
  }
  if (org_assembly == 'sc3') {
    if (length(packageCheck(
      c(
        "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
        "org.Sc.sgd.db"
      )
    )) == 0) {
      genome(myses) <- "sacCer3"
      data <-
        getTable(ucscTableQuery(myses, track = "refSeqComposite"))
      data <- data[, c(3, 4, 5, 6, 13)]
      colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
      pkg.env$genomee <- genomee
      
      mart <-
        useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c(
                        "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
                        "org.Sc.sgd.db"
                      )
                    )))
    }
  }
  if (org_assembly == 'dm6') {
    #Fruit fly
    if (length(packageCheck(
      c("TxDb.Dmelanogaster.UCSC.dm6.ensGene", "org.Dm.eg.db")
    )) == 0) {
      genome(myses) <- "dm6"
      a <- getTable(ucscTableQuery(myses, track = "ensGene"))
      a1 <- getTable(ucscTableQuery(myses, track = "ensGene",
                                    table = "ensemblToGeneName"))
      data <- merge(a1, a)
      data <- data[, c(4, 5, 6, 7, 2)]
      colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
      pkg.env$genomee <- genomee
      
      mart <-
        useMart("ensembl",
                dataset = "dmelanogaster_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Dmelanogaster.UCSC.dm6.ensGene", "org.Dm.eg.db")
                    )))
    }
  }
  if (org_assembly == 'ce11') {
    #Worm
    if (length(packageCheck(c(
      "TxDb.Celegans.UCSC.ce11.refGene", "org.Ce.eg.db"
    ))) == 0) {
      genome(myses) <- "ce11"
      a <- getTable(ucscTableQuery(myses, track = "ensGene"))
      a1 <-
        getTable(ucscTableQuery(myses, track = "ensGene",
                                table = "ensemblToGeneName"))
      data <- merge(a1, a)
      data <- data[, c(4, 5, 6, 7, 2)]
      colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
      ucsc <-
        with(data, GRanges(chr, IRanges(start, end), strand, symbol))
      pkg.env$ucsc <- ucsc
      
      genomee <- TxDb.Celegans.UCSC.ce11.refGene
      pkg.env$genomee <- genomee
      
      mart <-
        useMart(host = "useast.ensembl.org",
                biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "celegans_gene_ensembl")
      pkg.env$mart <- mart
    }
    else{
      message(paste("Please install",
                    packageCheck(
                      c("TxDb.Celegans.UCSC.ce11.refGene", "org.Ce.eg.db")
                    )))
    }
  }
}

#' Get nearest genes for the window of the upstream/downstream region.
#'
#' When downstream = 0 / upstream = 0, function converts bed formated regions
#' to HUGO genes
#'
#'
#' @param bedfile Bed formated input gene regions
#' @param upstream Maximum upstream distance from the transcription start
#'      region of the input gene
#' @param downstream Maximum downstream distance from the transcription end
#'      region of the input gene
#' @param org_assembly genomee assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#'
#'
#' @return genes
#'
#' @import GenomicFeatures
#' @importFrom GenomicRanges as.data.frame duplicated end findOverlaps
#' @importFrom GenomicRanges intersect match merge order pintersect resize
#' @importFrom GenomicRanges split start strand width
#'
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#'
#' neighbour <- getUCSC(bedfile = regionNC,
#'                      upstream = 1000,
#'                      downstream = 1000,
#'                      org_assembly = 'hg19')
#'
#'@export
getUCSC <-
  function(bedfile,
           upstream,
           downstream,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3")) {
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
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
    if (missing(org_assembly)) {
      message("genomee assembly version is missing.")
    }
    
    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    rt1 <- subsetByOverlaps(pkg.env$ucsc, unstrand(big_islands))
    
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    rt2 <- subsetByOverlaps(pkg.env$ucsc, unstrand(big_islands))
    results2 <- rt2$symbol
    results1 <- rt1$symbol
    result <- data.frame(unique(results1, results2))
    return(result)
  }

#' Get only those neighbouring genes that fall within exon region
#'
#' @param bedfile Input bed formated file
#' @param upstream Maximum upstream distance from the TSS position
#' @param downstream Maximum downstream distance from the TES position
#' @param org_assembly genomee assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'     "hg38" for human
#'
#' @return genes
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' r<-getNearToExon(bedfile = regionNC,
#'                  upstream = 1000,
#'                  downstream = 2000,
#'                  org_assembly = 'hg19')
#'
#' @export
getNearToExon <-
  function(bedfile,
           upstream,
           downstream,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3")) {
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
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
    exons = subsetByOverlaps(exons(pkg.env$genomee), big_islands)
    region1 <- subsetByOverlaps(pkg.env$ucsc, exons)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    exons = subsetByOverlaps(exons(pkg.env$genomee), big_islands)
    region2 <- subsetByOverlaps(pkg.env$ucsc, exons)
    results2 <- region2$symbol
    results1 <- region1$symbol
    result <- data.frame(unique(results1, results2))
    return(result)
  }

#' Get only those neighbouring genes that fall within intron region
#'
#' @param bedfile Bed file
#' @param upstream upstream distance
#' @param downstream downstream distance
#' @param org_assembly genomee assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#'
#' @return genes
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#'
#' r<-getNearToExon(bedfile = regionNC,
#'                  upstream = 1000,
#'                  downstream = 2000,
#'                  org_assembly = 'hg19')
#'
#' @importFrom GenomicFeatures as.list intronsByTranscript
#'
#' @export
getNearToIntron <-
  function(bedfile,
           upstream,
           downstream,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3")) {
    if (missing(org_assembly)) {
      message("Assembly is missing?")
    }
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
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
    intron = subsetByOverlaps(intronsByTranscript(pkg.env$genomee),
                              big_islands)
    
    region1 <- subsetByOverlaps(pkg.env$ucsc, intron)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    intron = subsetByOverlaps(intronsByTranscript(pkg.env$genomee), bedfile)
    region2 <- subsetByOverlaps(pkg.env$ucsc, intron)
    results2 <- region2$symbol
    results1 <- region1$symbol
    result <- data.frame(unique(results1, results2))
    return(result)
  }

#' For given region of interest, overlapped genes in the TAD regions are found.
#' Results can be filtered according to the available cell lines.
#'
#' @param bedfile Region of interest
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'      assembly, respectively.
#' @param cellline Cell lines for TAD regions.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'       considered in the analysis
#' @param upstream Holds upstream distance from the transcription start 
#'     position
#' @param downstream Holds downstream distance from the transcription end
#'      position
#'
#' @return List of protein coding genes that falls into the TAD regions
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' r<-getTADOverlap(bedfile = regionNC,
#'                  TAD = tad_hg38,
#'                  org_assembly = 'hg19',
#'                  cellline = 'HUVEC')
#'
#' @export
getTADOverlap <-
  function(bedfile,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           TAD = c(NoRCE::tad_hg19,
                   NoRCE::tad_dmel,
                   NoRCE::tad_hg38,
                   NoRCE::tad_mm10),
           near = FALSE,
           upstream = 10000,
           downstream = 10000,
           cellline = 'all') {
    if (missing(org_assembly)) {
      message("Assembly is missing?")
    }
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
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
    
    if (near) {
      treg <-
        Reduce(subsetByOverlaps, list(
          tad,
          bedfile,
          resize(bedfile, width = downstream + width(bedfile), fix = "end")
        ))
      region2 <- subsetByOverlaps(pkg.env$ucsc, treg)
      treg1 <-
        Reduce(subsetByOverlaps, list(
          tad,
          bedfile,
          resize(bedfile, width = upstream + width(bedfile), fix = "start")
        ))
      region1 <- subsetByOverlaps(pkg.env$ucsc, treg1)
      region <-
        data.frame(symbol = unique(region1$symbol, region2$symbol))
    }
    else
    {
      tadRegion = subsetByOverlaps(bedfile, tad)
      region <- subsetByOverlaps(pkg.env$ucsc, tadRegion)
    }
    
    return(unique(as.data.frame(region$symbol)))
  }

#' Convert gene ids according to the gene type
#'
#' @param genetype Type of the input gene list. Possible values are "Entrez",
#'      "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol
#'      "NCBI" value, for Entrez gene id "Entrez", for mirbase id "mirna" is
#'      used.
#' @param genelist Input gene list
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#'
#' @return GRange object of the given input
#'
#' @importFrom biomaRt getBM
#'
#' @examples
#'
#' convGene <-convertGeneID(genetype = "mirna",
#'                          genelist = brain_mirna[1:30,],
#'                          org_assembly = 'hg19')
#'
#'
#' @export
convertGeneID <-
  function(genetype = c("Entrez",
                        "mirna",
                        "Ensembl_gene",
                        "Ensembl_trans",
                        "NCBI"),
           genelist,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3")) {
    if (missing(org_assembly)) {
      message("genomee assembly version is missing.")
    }
    if (missing(genetype)) {
      message("Format of the gene is missing.")
    }
    if (missing(genelist)) {
      message("List of gene is missing.")
    }
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
    }
    attributes <-
      c("chromosome_name",
        "start_position",
        "end_position",
        "strand")
    if (genetype == "Entrez") {
      output <-
        getBM(
          attributes = c("entrezgene_id", attributes),
          filters = "entrezgene_id",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    
    else if (genetype == "mirna") {
      output <-
        getBM(
          attributes = c("mirbase_id", attributes),
          filters = "mirbase_id",
          values = apply(as.data.frame(genelist), 2, tolower),
          mart = pkg.env$mart
        )
    }
    else if (genetype == "Ensembl_gene") {
      output <-
        getBM(
          attributes = c("ensembl_gene_id", attributes),
          filters = "ensembl_gene_id",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    else if (genetype == "Ensembl_trans") {
      output <-
        getBM(
          attributes = c("ensembl_transcript_id", attributes),
          filters = "ensembl_transcript_id",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    else if (genetype == "NCBI") {
      output <-
        getBM(
          attributes = c("hgnc_symbol", attributes),
          filters = "hgnc_symbol",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    else if (genetype == "mgi_symbol")
    {
      output <-
        getBM(
          attributes = c("mgi_symbol", attributes),
          filters = "mgi_symbol",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    else if (genetype == "external_gene_name")
    {
      output <-
        getBM(
          attributes = c("external_gene_name", attributes),
          filters = "external_gene_name",
          values = genelist,
          mart = pkg.env$mart
        )
    }
    
    colnames(output) <- c("gene",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand")
    
    
    file1 <-
      with(output, GRanges(
        paste0("chr", chromosome_name),
        IRanges(start_position, end_position),
        '*',
        gene
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
listTAD <- function(TADName) {
  return(unique(TADName$celline))
}


#' Check the package availability for the given assembly
packageCheck <- function(pkg)
{
  notInstalled <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  return(notInstalled)
}
