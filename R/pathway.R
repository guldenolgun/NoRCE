#' KEGG pathway enrichment
#'
#' @param genes Input genes
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#'
#'
#' @return KEGG pathway enrichment results
#'
#' @importFrom KEGGREST keggGet
#'
#' @examples
#' subsetGene <- breastmRNA[1:100,]
#'
#' \dontrun{
#' br_enr<-KeggEnrichment(genes = subsetGene,hg='hg19')
#' getKeggDiagram(mrnaObject = br_enr,pathway ="hsa04060" ,hg='hg19')
#'
#'}
#' @export
#'
KeggEnrichment <-
  function(genes,
           hg,
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           min = 5) {
    if (missing(genes)) {
      message("Gene is missing.")
    }
    if (missing(hg)) {
      message(
        "Assembly version is missing. Possible assemblies are \"mm10\" for mouse, \"dre10\" for zebrafish, \"rn6\" for rat, \"dm6\" for fruit fly, \"ce11\" for worm, \"hg19\" and \"hg38\" for human."
      )
    }
    assembly(hg)
    pathTable <- unique(keggPathwayDB(hg))
    genes <- as.data.frame(genes)
    colnames(genes) <- 'g'
    annot <- pathTable[which(pathTable$symbol %in% genes$g),]


    pathfreq <- as.data.frame(table(annot$pathway))
    pathfreq <- pathfreq[which(pathfreq$Freq > 0),]


    geneSize = length(unique(pathTable$symbol))

    bckfreq <- as.data.frame(table(pathTable$pathway))
    notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
    freq<-merge(pathfreq,notGene, by = "Var1")
    found <- freq$Freq.x
    M <- freq$Freq.y
    n <- rep(length(unique(annot$symbol)), length(M))

    pvalues <-
      phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)

    pAdjust1 <- p.adjust(pvalues, method = pAdjust)

    GeneRatio <-
      apply(data.frame(found, n), 1, function(x)
        paste(x[1], "/", x[2], sep = "", collapse = ""))

    BgRatio <-
      apply(data.frame(M, geneSize), 1, function(x)
        paste(x[1], "/", x[2], sep = "", collapse = ""))

    enrich <-
      which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
    pathT <- as.character(freq$Var1[enrich])
    gratio <- GeneRatio[enrich]
    bgratio <- BgRatio[enrich]
    padj <- pAdjust1[enrich]
    pval <- pvalues[enrich]
    r <- annot[annot$pathway %in% pathT,]
    enric <- list()
    for (i in 1:length(pathT))
    {
      if (length(which(pathT[i] == r$pathway)) > 0)
      {
        enric <-
          c(enric, setNames(list(as.character(r[which(pathT[i] == r$pathway),]$symbol)), paste(pathT[i])))
      }
    }

    pathways <- data.frame(unique(pathT))
    tmp <- character(length(pathT))
    if (dim(pathways)[1] > 0) {
      tmp <-
        sapply(1:dim(pathways)[1], function(x)
          tmp[x] <- try(keggGet(pathT[x])[[1]]$NAME)
        )
    }
    return(
      new(
        "NoRCE",
        ID = pathT,
        Term = tmp,
        geneList = enric,
        pvalue = pval,
        pAdj = padj,
        GeneRatio = gratio,
        BckRatio = bgratio
      )
    )
  }


#' Reactome pathway enrichment
#'
#' @param genes Input genes
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#'
#'
#' @return Reactome pathway enrichment results
#'
#'
#' @importFrom org.Ce.eg.db org.Ce.egPATH2EG org.Ce.egSYMBOL
#' @importFrom org.Dr.eg.db org.Dr.egPATH2EG org.Dr.egSYMBOL
#' @importFrom org.Hs.eg.db org.Hs.egPATH2EG org.Hs.egSYMBOL
#' @importFrom org.Mm.eg.db org.Mm.egPATH2EG org.Mm.egSYMBOL
#' @importFrom org.Rn.eg.db org.Rn.egPATH2EG org.Rn.egSYMBOL
#' @importFrom org.Sc.sgd.db org.Sc.sgdGENENAME org.Sc.sgdPATH2ORF
#' @importFrom reactome.db reactomePATHID2EXTID reactomePATHID2NAME
#'
#' @examples
#' data(breastmRNA)
#'
#' \dontrun{
#'br_enr<-reactomeEnrichment(genes = breastmRNA,hg='hg19')
#' getReactomeDiagram(mrnaObject = br_enr, pathway = br_enr@ID[1], imageFormat = 'svg')
#' }
#'
#' @export
reactomeEnrichment <-
  function(genes,
           hg,
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           min = 5) {
    if (missing(genes)) {
      message("Gene is missing.")
    }
    if (missing(hg)) {
      message(
        "Assembly version is missing. Possible assemblies are \"mm10\" for mouse, \"dre10\" for zebrafish, \"rn6\" for rat, \"dm6\" for fruit fly, \"ce11\" for worm, \"hg19\" and \"hg38\" for human."
      )
    }
    assembly(hg)
    pathTable <- unique(reactomePathwayDB(hg))
    genes <- as.data.frame(genes)
    colnames(genes) <- 'g'
    annot <- pathTable[which(pathTable$symbol %in% genes$g),]


    pathfreq <- as.data.frame(table(annot$pathway))
    pathfreq <- pathfreq[which(pathfreq$Freq > 0),]

    geneSize = length(unique(pathTable$symbol))

    bckfreq <- as.data.frame(table(pathTable$pathway))
    notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
    freq<-merge(pathfreq,notGene, by = "Var1")
    found <- freq$Freq.x
    M <- freq$Freq.y

    n <- rep(length(unique(annot$symbol)), length(M))

    pvalues <-
      phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)

    pAdjust1 <- p.adjust(pvalues, method = pAdjust)

    GeneRatio <- apply(data.frame(found, n), 1, function(x)
      paste(x[1], "/", x[2], sep = "", collapse = ""))

    BgRatio <- apply(data.frame( M, geneSize), 1, function(x)
      paste(x[1], "/", x[2], sep = "", collapse = ""))

    enrich <-
      which (pvalues < pCut & pAdjust1 < pAdjCut & found >= min)
    pathT <- as.character(freq$Var1[enrich])
    gratio <- GeneRatio[enrich]
    bgratio <- BgRatio[enrich]
    padj <- pAdjust1[enrich]
    pval <- pvalues[enrich]
    r <- annot[annot$pathway %in% pathT,]
    rt <- unique(r[, (2:3)])

    enric <- list()
    for (i in 1:length(pathT))
    {
      if (length(which(pathT[i] == r$pathway)) > 0)
      {
        enric <-
          c(enric, setNames(list(as.character(r[which(pathT[i] == r$pathway),]$symbol)), paste(pathT[i])))
      }
    }
    return(
      new(
        "NoRCE",
        ID = pathT,
        Term = as.character(rt[order(match(rt$pathway, pathT)), ]$name),
        geneList = enric,
        pvalue = pval,
        pAdj = padj,
        GeneRatio = gratio,
        BckRatio = bgratio
      )
    )
  }

reactomePathwayDB <- function(hg) {
  xx <- as.list(reactomePATHID2EXTID)
  table1 <- data.frame(pathway = rep(names(xx), lapply(xx, length)),
                       gene = unlist(xx))
  pn <- as.list(reactomePATHID2NAME)
  pn <- data.frame(pathway = rep(names(pn), lapply(pn, length)),
                   name = unlist(pn))
  if (hg == 'hg19' | hg == 'hg38') {
    ty <- table1[grepl("^R-HSA", table1$pathway),]
    pn1 <- pn[grepl("^R-HSA", pn$pathway),]
    symb <- as.data.frame(org.Hs.egSYMBOL)
  }
  if (hg == 'mm10') {
    ty <- table1[grepl("^R-MMU", table1$pathway),]
    pn1 <- pn[grepl("^R-MMU", pn$pathway),]
    symb <- as.data.frame(org.Mm.egSYMBOL)
  }
  if (hg == 'dre10') {
    ty <- table1[grepl("^R-DRE", table1$pathway),]
    pn1 <- pn[grepl("^R-DRE", pn$pathway),]
    symb <- as.data.frame(org.Dr.egSYMBOL)
  }
  if (hg == 'rn6') {
    ty <- table1[grepl("^R-RNO", table1$pathway),]
    pn1 <- pn[grepl("^R-RNO", pn$pathway),]
    symb <- as.data.frame(org.Rn.egSYMBOL)
  }
  if (hg == 'ce11') {
    ty <- table1[grepl("^R-CEL", table1$pathway),]
    pn1 <- pn[grepl("^R-CEL", pn$pathway),]
    symb <- as.data.frame(org.Ce.egSYMBOL)
  }
  colnames(symb) <- c("gene", "symbol")
  merge1 <- merge(x = pn1,
                  y = ty,
                  by = "pathway",
                  all.y = TRUE)
  path <- merge(merge1, symb, by = "gene")

  return(path)
}

keggPathwayDB <- function(hg) {
  if (hg == 'hg19' | hg == 'hg38') {
    kegg <- org.Hs.egPATH2EG
    x <- org.Hs.egSYMBOL
    prefix <- 'hsa'
  }
  if (hg == 'mm10') {
    kegg <- org.Mm.egPATH2EG
    x <- org.Mm.egSYMBOL
    prefix <- 'mmu'
  }
  if (hg == 'dre10') {
    kegg <- org.Dr.egPATH2EG
    x <- org.Dr.egSYMBOL
    prefix <- 'dre'
  }
  if (hg == 'rn6') {
    kegg <- org.Rn.egPATH2EG
    x <- org.Rn.egSYMBOL
    prefix <- 'rno'
  }
  if (hg == 'ce11') {
    kegg <- org.Ce.egPATH2EG
    x <- org.Ce.egSYMBOL
    prefix <- 'cel'
  }
  if (hg == 'sc3') {
    kegg <- org.Sc.sgdPATH2ORF
    x <- org.Sc.sgdGENENAME
    prefix <- 'sce'
  }
  mapped <- mappedkeys(kegg)
  kegg2 <- as.list(kegg[mapped])
  pathTable <-
    data.frame(pathway = paste0(prefix, rep(names(kegg2), lapply(kegg2, length))),
               gene = unlist(kegg2))
  x <- as.data.frame(x)
  colnames(x) <- c("gene", "symbol")
  pathTable <- merge(pathTable, x, by = "gene")
  return(pathTable)
}

WikiPathwayDB <- function(hg) {
  if (hg == 'hg19' | hg == 'hg38')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens", format = "gmt")
  if (hg == 'mm10')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt")
  if (hg == 'dre10')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Danio rerio", format = "gmt")
  if (hg == 'rn6')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Rattus norvegicus", format = "gmt")
  if (hg == 'dre10')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Drosophila melanogaster", format = "gmt")
  if (hg == 'ce11')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Caenorhabditis elegans", format = "gmt")
  if (hg == 'sc3')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "	Saccharomyces cerevisiae", format = "gmt")

  gmtFile <- readGMT(gmtName = wp.gmt)
  tmp <-
    do.call(rbind, strsplit(as.character(gmtFile$pathTerm), '%'))
  gmtFile <-
    unique(
      data.frame(
        Entrez = gmtFile$Entrez,
        gene = gmtFile$symbol,
        pathID = tmp[, 3],
        pathTerm = tmp[, 1]
      )
    )
  return(gmtFile)
}

#' WikiPathways Enrichment
#'
#' @param genes Input genes
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#'
#' @return Wiki Pathway Enrichment
#'
#' @importFrom rWikiPathways downloadPathwayArchive
#'
#' @examples
#'
#' geneInterest <- breastmRNA[1:100,]
#' \dontrun{
#' br_enr<-WikiEnrichment(genes = geneInterest,hg='hg19')
#'}
#' @export
WikiEnrichment <- function(genes,
                           hg,
                           pCut = 0.05,
                           pAdjCut = 0.05,
                           pAdjust = "BH",
                           min = 5) {
  if (missing(genes)) {
    message("Gene is missing.")
  }
  if (missing(hg)) {
    message(
      "Assembly version is missing. Possible assemblies are \"mm10\" for mouse, \"dre10\" for zebrafish, \"rn6\" for rat, \"dm6\" for fruit fly, \"ce11\" for worm, \"hg19\" and \"hg38\" for human."
    )
  }
  assembly(hg)
  pathTable <- unique(WikiPathwayDB(hg))
  genes <- as.data.frame(genes)
  colnames(genes) <- 'g'
  annot <- pathTable[which(pathTable$gene %in% genes$g),]

  pathfreq <- as.data.frame(table(annot$pathID))
  pathfreq <- pathfreq[which(pathfreq$Freq > 0),]

  geneSize = length(unique(pathTable$gene))
  bckfreq <- as.data.frame(table(pathTable$pathID))
  notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
  freq<-merge(pathfreq,notGene, by = "Var1")
  found <- freq$Freq.x
  M <- freq$Freq.y

  n <- rep(length(unique(annot$gene)), length(M))

  pvalues <-
    phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)

  pAdjust1 <- p.adjust(pvalues, method = pAdjust)
  GeneRatio <-
    apply(data.frame(found, n), 1, function(x)
      paste(x[1], "/", x[2], sep = "", collapse = ""))

  BgRatio <-
    apply(data.frame(M, geneSize), 1, function(x)
      paste(x[1], "/", x[2], sep = "", collapse = ""))

  enrich <-
    which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
  pathT <- as.character(freq$Var1[enrich])
  gratio <- GeneRatio[enrich]
  bgratio <- BgRatio[enrich]
  padj <- pAdjust1[enrich]
  pval <- pvalues[enrich]
  r <- annot[annot$pathID %in% pathT,]
  enric <- list()
  pathTerms <- as.character(r$pathTerm[match(pathT, r$pathID)])
  for (i in 1:length(pathT))
  {
    if (length(which(pathT[i] == r$pathID)) > 0)
    {
      enric <-
        c(enric, setNames(list(as.character(r[which(pathT[i] == r$pathID),]$gene)), paste(pathT[i])))
    }
  }

  rm(ucsc,go,mart)

  return(
    new(
      "NoRCE",
      ID = pathT,
      Term = pathTerms,
      geneList = enric,
      pvalue = pval,
      pAdj = padj,
      GeneRatio = gratio,
      BckRatio = bgratio
    )
  )
}

#' For a given gmt file of a specific pathway database, pathway enrichment can be performed. Function supports Entrez ID and symbol based gmt file.
#'
#' @param genes Input genes
#' @param gmtFile File path of the gmt file
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param isSymbol Boolean value that controls the gene formats. If it is TRUE, gene format of the gmt file should be symbol. Otherwise, gene format must be ENTREZ ID.
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#' @param isGeneEnrich Boolean value whether gene enrichment should be performed
#'
#' @return Pathway Enrichment
#'
#' @examples
#' tmpGene <- breastmRNA[1:100,]
#'\dontrun{
#' br_enr<-pathwayEnrichment(genes = tmpGene,hg='hg19',
#'     gmtFile = 'Human_AllPathways_February_01_2019_symbol.gmt',isSymbol = TRUE)
#'}
#' @export
pathwayEnrichment <- function(genes,
                              gmtFile,
                              hg,
                              pCut = 0.05,
                              pAdjCut = 0.05,
                              pAdjust = "BH",
                              isSymbol,
                              min = 5, isGeneEnrich = FALSE) {
  if (missing(genes))
    message("Gene is missing.")

  if (missing(hg)){
    message(
      "Assembly version is missing. Possible assemblies are \"mm10\" for mouse, \"dre10\" for zebrafish, \"rn6\" for rat, \"dm6\" for fruit fly, \"ce11\" for worm, \"hg19\" and \"hg38\" for human."
    )
    assembly(hg)
  }
  if (missing(gmtFile))
    message("GMT file is missing.")

  if (missing(isSymbol))
    message("GMT gene format is missing.")

  assembly(hg)

  genes <- as.data.frame(genes)
  colnames(genes) <- 'g'

  if(!isGeneEnrich){
    pathTable <- unique(readGMT(gmtName = gmtFile, isSymbol = isSymbol))
  }
  else{
    pathTable <- geneListEnrich(f = gmtFile, isSymbol = isSymbol)
  }
  annot <- pathTable[which(pathTable$symbol %in% genes$g),]
  pathfreq <- as.data.frame(table(annot$pathTerm))
  pathfreq <- pathfreq[which(pathfreq$Freq > 0),]


  if(!isGeneEnrich)
    geneSize = length(unique(pathTable$symbol))
  else
    geneSize = length(ucsc)


  bckfreq <- as.data.frame(table(pathTable$pathTerm))

  notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
  freq<-merge(pathfreq,notGene, by = "Var1")
  found <- freq$Freq.x
  M <- freq$Freq.y

  n <- rep(length(unique(annot$symbol)), length(M))

  pvalues <-
    phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)

  pAdjust1 <- p.adjust(pvalues, method = pAdjust)

  GeneRatio <- apply(data.frame(found, n), 1, function(x)
    paste(x[1], "/", x[2], sep = "", collapse = ""))

  BgRatio <- apply(data.frame(M, geneSize), 1, function(x)
    paste(x[1], "/", x[2], sep = "", collapse = ""))

  enrich <-
    which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
  pathT <- as.character(freq$Var1[enrich])
  gratio <- GeneRatio[enrich]
  bgratio <- BgRatio[enrich]
  padj <- pAdjust1[enrich]
  pval <- pvalues[enrich]

  r <- annot[annot$pathTerm %in% pathT,]
  enric <- list()
  pathTerms <- as.character(r$pathTerm[match(pathT, r$pathID)])

  for (i in 1:length(pathT))
  {
    if (length(which(pathT[i] == r$pathTerm)) > 0)
      enric <-
        c(enric, setNames(list(as.character(r[which(pathT[i] == r$pathTerm),]$symbol)), paste(pathT[i])))
  }
  rm(ucsc,go,mart)
  return(
    new(
      "NoRCE",
      ID = pathT,
      Term = pathTerms,
      geneList = enric,
      pvalue = pval,
      pAdj = padj,
      GeneRatio = gratio,
      BckRatio = bgratio
    )
  )
}

#' Import gmt formatted pathway file to the NoRCE
#'
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file. If it is set as TRUE, gene format of the gmt file should be symbol. Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
#'
#' @return return gmt file
#'
#'
#' @importFrom biomaRt getBM
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gmtFile <- readGMT(gmtName = wp.gmt)
#'}
#'
readGMT <- function(gmtName,
                    isSymbol = FALSE) {
  x <- scan(gmtName, what = "", sep = "\n")
  x <- strsplit(x, '\t')
  max.col <- max(sapply(x, length))
  cn <- paste("V", 1:max.col, sep = "")

  x <-
    read.table(
      gmtName,
      fill = TRUE,
      col.names = cn,
      sep = "\t",
      quote = NULL
    )
  x <- x[,-c(2)]
  f <- melt(x, id = c('V1'))
  f <- f[,-c(2)]
  if (length(which(f$value %in% '')) > 0)
    f <- f[-c(which(f$value %in% '')), ]
  if (!isSymbol) {
    output <-
      getBM(
        attributes = c('entrezgene', 'hgnc_symbol'),
        filters = "entrezgene",
        values = f$value,
        mart = mart
      )
    f[, 3] <- output[match(f$value, output$entrezgene), 2]
    colnames(f) <- c('pathTerm', 'Entrez', 'symbol')
  }
  else{
    output <-
      getBM(
        attributes = c('entrezgene', 'hgnc_symbol'),
        filters = "hgnc_symbol",
        values = f$value,
        mart = mart
      )
    f[, 3] <- output[match(f$value, output$hgnc_symbol), 1]
    colnames(f) <- c('pathTerm', 'symbol', 'Entrez')
  }
  return(f)
}
geneListEnrich <- function(f, isSymbol=FALSE)
{
  if(dim(f)[2]==1)
    f <- cbind("NotAva",f)

  colnames(f) <- c("Annot","value")

  if (!isSymbol) {
    output <-
      getBM(
        attributes = c('entrezgene', 'hgnc_symbol'),
        filters = "entrezgene",
        values = f$value,
        mart = mart
      )
    f[, 3] <- output[match(f$value, output$entrezgene), 2]
    colnames(f) <- c('pathTerm', 'Entrez', 'symbol')
  }
  else{
    output <-
      getBM(
        attributes = c('entrezgene', 'hgnc_symbol'),
        filters = "hgnc_symbol",
        values = f$value,
        mart = mart
      )
    f[, 3] <- output[match(f$value, output$hgnc_symbol), 1]
    colnames(f) <- c('pathTerm', 'symbol', 'Entrez')
  }
  return(f)
}

