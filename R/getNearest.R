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

pkg.env$upstream <- 10000
pkg.env$downstream <- 10000
pkg.env$searchRegion <-
  "all" #searchRegion = c('all',"exon","intron")
pkg.env$GOtype = "BP"  #c("BP", "CC", "MF")
pkg.env$pCut = 0.05
pkg.env$pAdjCut = 0.05
#c("holm","hochberg","hommel","bonferroni", "BH", "BY","fdr", "none")
pkg.env$pAdjust = "none"
#c("hyper", "binom", "fisher", "chi")
pkg.env$enrichTest = "hyper"
pkg.env$varCutoff = 0.0025
pkg.env$minAbsCor = 0.3
pkg.env$pcut = 0.05
pkg.env$conf = 0.95
pkg.env$min = 5
pkg.env$cellline = 'all'
pkg.env$corrMethod = "pearson" #c("pearson", "kendall", "spearman")
pkg.env$alternate = "two.sided" #c('greater', "two.sided", "less")
pkg.env$pathwayType = 'kegg' #c('kegg', 'reactome','wiki','other'),
pkg.env$isSymbol = FALSE

#' Get the required information for the given assembly
#'
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38"
#'     for human
#'
#' @return setting required information
#'
#' @importFrom biomaRt getBM useEnsembl useMart
#' @importFrom rtracklayer browserSession genome getTable ucscTableQuery
#' @importFrom rtracklayer genome<-
#' @importFrom IRanges IRanges
#'
#' @examples
#'
#' assembly('hg19')
#'
#' @export
assembly <- function(org_assembly = c("hg19",
                                      "hg38",
                                      "mm10",
                                      "dre10",
                                      "rn6",
                                      "dm6",
                                      "ce11",
                                      "sc3")) {
  myses <- browserSession()
  types <-
    rbind(
      c(
        "Hsapiens",
        "hg19",
        "hg19",
        "wgEncodeGencodeV31lift37",
        3,
        4,
        5,
        6,
        13,
        "knownGene",
        "hsapiens_gene_ensembl",
        "Hs.eg.db"
      ),
      c(
        "Hsapiens",
        "hg38",
        "hg38",
        "wgEncodeGencodeV31",
        3,
        4,
        5,
        6,
        13,
        "knownGene",
        "hsapiens_gene_ensembl",
        "Hs.eg.db"
      ),
      c(
        "Mmusculus",
        "mm10",
        "mm10",
        "wgEncodeGencodeVM22",
        3,
        4,
        5,
        6,
        13,
        "knownGene",
        "mmusculus_gene_ensembl",
        "Mm.eg.db"
      ),
      c(
        "Drerio",
        "danRer10",
        "dre10",
        "ensGene",
        4,
        5,
        6,
        7,
        2,
        "refGene",
        "drerio_gene_ensembl",
        "Dr.eg.db"
      ),
      c(
        "Rnorvegicus",
        "rn6",
        "rn6",
        "ensGene",
        4,
        5,
        6,
        7,
        2,
        "refGene",
        "rnorvegicus_gene_ensembl",
        "Rn.eg.db"
      ),
      c(
        "Scerevisiae",
        "sacCer3",
        "sc3",
        "refSeqComposite",
        3,
        4,
        5,
        6,
        13,
        "sgdGene",
        "scerevisiae_gene_ensembl",
        "Sc.sgd.db"
      ),
      c(
        "Dmelanogaster",
        "dm6",
        "dm6",
        "ensGene",
        4,
        5,
        6,
        7,
        2,
        "ensGene",
        "dmelanogaster_gene_ensembl",
        "Dm.eg.db"
      ),
      c(
        "Celegans",
        "ce11",
        "ce11",
        "ensGene",
        4,
        5,
        6,
        7,
        2,
        "refGene",
        "celegans_gene_ensembl",
        "Ce.eg.db"
      )
    )
  
  index = which(org_assembly == types[, 3])
  
  genome(myses) <- types[index, 2]
  
  if (index == 4 & index == 5 & index == 7 & index == 8) {
    a <- getTable(ucscTableQuery(myses, track = "ensGene"))
    a1 <-
      getTable(ucscTableQuery(myses, track = "ensGene",
                              table = "ensemblToGeneName"))
    data <- merge(a1, a)
    data <- data[, c(4, 5, 6, 7, 2)]
  }
  else{
    data <- getTable(ucscTableQuery(myses, track = types[index, 4]))
    data <- data[, as.double(types[index, 5:9])]
  }
  colnames(data) <- c('chr', 'strand', 'start', 'end', 'symbol')
  ucsc <-
    with(data, GRanges(chr, IRanges::IRanges(start, end), strand, symbol))
  pkg.env$ucsc <- ucsc
  td <-
    paste0("TxDb.", types[index, 1], ".UCSC.", types[index, 2], ".",
           types[index, 10])
  yy <- paste0("org.", types[index, 12])
  x <- list(yy, td)
  
  if (!requireNamespace(td, quietly = TRUE) |
      !requireNamespace(yy, quietly = TRUE))
    stop("Install package ", td, " or ", yy, " in order to use this function.")
  else
    lapply(x, require, character.only = TRUE)
  
  genomee <- eval(as.name(td))
  pkg.env$genomee <- genomee
  
  if (index == 1) {
    mart = useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      host = "grch37.ensembl.org",
      path = "/biomart/martservice",
      dataset = "hsapiens_gene_ensembl"
    )
  }
  else{
    mart <-
      useEnsembl(biomart = "ensembl", dataset = types[index, 11])
  }
  pkg.env$mart <- mart
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
#' @importFrom IRanges subsetByOverlaps
#' @import zlibbioc
#'
#' @examples
#'
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
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
    
    assembly(org_assembly)
    
    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    rt1 <-
      IRanges::subsetByOverlaps(pkg.env$ucsc, 
                                BiocGenerics::unstrand(big_islands))
    
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    rt2 <-
      IRanges::subsetByOverlaps(pkg.env$ucsc, 
                                BiocGenerics::unstrand(big_islands))
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
#' regions <- system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#'
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
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }
    if (missing(org_assembly)) {
      message("Assembly is missing?")
    }
    if (missing(upstream)) {
      message("Upstream information is missing?")
    }
    if (missing(downstream)) {
      message("Downstream information is missing?")
    }
    assembly(org_assembly)
    
    big_islands <-
      resize(bedfile, width = downstream + width(bedfile), fix = "end")
    exons = IRanges::subsetByOverlaps(exons(pkg.env$genomee), big_islands)
    region1 <- IRanges::subsetByOverlaps(pkg.env$ucsc, exons)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    exons = IRanges::subsetByOverlaps(exons(pkg.env$genomee), big_islands)
    region2 <- IRanges::subsetByOverlaps(pkg.env$ucsc, exons)
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
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
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
    assembly(org_assembly)
    
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
    intron = IRanges::subsetByOverlaps(intronsByTranscript(pkg.env$genomee),
                                       big_islands)
    
    region1 <- IRanges::subsetByOverlaps(pkg.env$ucsc, intron)
    big_islands <-
      resize(bedfile, width = upstream + width(bedfile), fix = "start")
    intron = IRanges::subsetByOverlaps(intronsByTranscript(pkg.env$genomee),
                                       bedfile)
    region2 <- IRanges::subsetByOverlaps(pkg.env$ucsc, intron)
    results2 <- region2$symbol
    results1 <- region1$symbol
    result <- data.frame(unique(results1, results2))
    return(result)
  }

#' For given region of interest, overlapped genes in the TAD regions are
#' found. Results can be filtered according to the available cell lines.
#'
#' @param bedfile Region of interest
#' @param tad TAD genomic regions for the species. Predefined TAD regions or
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
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#'
#' r<-getTADOverlap(bedfile = regionNC,
#'                  tad = tad_hg19,
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
           tad = c(tad_hg19,
                   tad_dmel,
                   tad_hg38,
                   tad_mm10),
           near = FALSE,
           upstream = 10000,
           downstream = 10000,
           cellline = 'all') {
    if (missing(org_assembly)) {
      message("Assembly is missing?")
    }
    
    assembly(org_assembly)
    
    if (missing(bedfile)) {
      message("Bed file is missing?")
    }
    
    if (cellline != 'all') {
      temp <- which(tad$celline == cellline)
      tad <-  tad[temp,]
    }
    
    if (near) {
      treg <-
        Reduce(IRanges::subsetByOverlaps, list(
          tad,
          bedfile,
          resize(bedfile, width = downstream + width(bedfile), fix = "end")
        ))
      region2 <- IRanges::subsetByOverlaps(pkg.env$ucsc, treg)
      treg1 <-
        Reduce(IRanges::subsetByOverlaps, list(
          tad,
          bedfile,
          resize(bedfile, width = upstream + width(bedfile), fix = "start")
        ))
      region1 <- IRanges::subsetByOverlaps(pkg.env$ucsc, treg1)
      region <-
        data.frame(symbol = unique(region1$symbol, region2$symbol))
    }
    else
    {
      tadRegion = IRanges::subsetByOverlaps(bedfile, tad)
      region <- IRanges::subsetByOverlaps(pkg.env$ucsc, tadRegion)
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
#' @importFrom IRanges IRanges
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
      message("Genome assembly version is missing.")
    }
    if (missing(genetype)) {
      message("Format of the gene is missing.")
    }
    if (missing(genelist)) {
      message("List of gene is missing.")
    }
    
    assembly(org_assembly)
    
    attributes <-
      c("chromosome_name",
        "start_position",
        "end_position",
        "strand")
    
    ifelse(
      genetype == "Entrez",
      output <-
        getBM(
          attributes = c("entrezgene_id", attributes),
          filters = "entrezgene_id",
          values = genelist,
          mart = pkg.env$mart
        ),
      ifelse(
        genetype == "mirna" ,
        output <-
          getBM(
            attributes = c("mirbase_id", attributes),
            filters = "mirbase_id",
            values = apply(as.data.frame(genelist), 2, tolower),
            mart = pkg.env$mart
          ),
        ifelse(
          genetype == "Ensembl_gene",
          output <-
            getBM(
              attributes = c("ensembl_gene_id", attributes),
              filters = "ensembl_gene_id",
              values = genelist,
              mart = pkg.env$mart
            ),
          ifelse(
            genetype == "Ensembl_trans",
            output <-
              getBM(
                attributes = c("ensembl_transcript_id", attributes),
                filters = "ensembl_transcript_id",
                values = genelist,
                mart = pkg.env$mart
              ),
            ifelse(
              (genetype == "NCBI" & (org_assembly == 'hg19' | org_assembly == 'hg38')),
              output <-
                getBM(
                  attributes = c("hgnc_symbol", attributes),
                  filters = "hgnc_symbol",
                  values = genelist,
                  mart = pkg.env$mart),
              ifelse(
                (genetype == "NCBI" & org_assembly == 'mm10'),
                output <-
                  getBM(
                    attributes = c("mgi_symbol", attributes),
                    filters = "mgi_symbol",
                    values = genelist,
                    mart = pkg.env$mart
                  ),
                output <-
                  getBM(
                    attributes = c("external_gene_name", attributes),
                    filters = "external_gene_name",
                    values = genelist,
                    mart = pkg.env$mart
                  )
              )
            )
          )
        )
      )
    )
    
    
    colnames(output) <- c("gene",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand")
    
    
    file1 <-
      with(output, GRanges(
        paste0("chr", chromosome_name),
        IRanges::IRanges(start_position, end_position),
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
#' @param pkg Required packages
#'
#' @return return install packages
#'
#' @importFrom utils installed.packages
#'
packageCheck <- function(pkg)
{
  notInstalled <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  return(notInstalled)
}

#' Set the parameters
#'
#' Parameters:
#' upstream: Upstream distance from the transcription start position
#' downstream: Downstream distance from the transcription end position
#' searchRegion: Search space of the cis-region. Possible values are
#'     "all", "exon", "intron"
#' GOtype: Hierarchical category of the GO ontology. Possible values
#'     are "BP", "CC", "MF"
#' pCut: Threshold value for the pvalue. Default value is 0.05
#' pAdjCut: Cutoff value for the adjusted p-values using one of given
#'     method. Default value is 0.05.
#' pAdjust: Methods of the adjusted p-values. Possible methods are
#'     "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' min: Minimum number of genes that are required for enrichment. By
#'     default, this value is set to 5.
#' cellline: Cell lines for TAD regions.
#' corrMethod Correlation coeffient method that will be used for
#'     evaluation. Possible values are "pearson", "kendall", "spearman"
#' varCutoff: Variance cutt off that genes have less variance than this
#'     value will be trimmed
#' pcut: P-value cut off for the correlation values
#' alternate: Holds the alternative hypothesis and "two.sided", "greater"
#'     or "less" are the possible values.
#' conf: Confidence level for the returned confidence interval. It is
#'     only used for the Pearson correlation coefficient if there are at
#'     least 4 complete pairs of observations.
#' minAbsCor: Cut-off value for the Pearson correlation coefficient of
#'     the miRNA-mRNA
#' pathwayType: Pathway database for enrichment. Possible values are
#'     'reactome' for Reactome, 'kegg' for KEGG, 'wiki' for WikiPathways,
#'     'other' for custom database
#' enrichTest: Types of enrichment methods to perform enrichment
#'      analysis. Possible values are "hyper"(default), "binom", "fisher",
#'      "chi".
#' isSymbol: Boolean variable that hold the gene format of the gmt file.
#'      If it is set as TRUE, gene format of the gmt file should be symbol.
#'      Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
#'
#' @param type List of parameter names
#' @param value New values for the parameters. Value and the parameter names
#'  must be in the same order.
#'
#' @return changed parameters
#'
#' @examples
#'
#' type <- c('downstream','upstream')
#'
#' value <- c(2000,30000)
#'
#' setParameters(type,value)
#'
#' @export
setParameters <- function(type, value) {
  for (i in seq_along(type)) {
    eval(parse(text = paste0("pkg.env$", type[i], " <- value[i]")))
  }
}