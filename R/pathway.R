#' KEGG pathway enrichment
#'
#' @param genes Input genes
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment.
#'      By default, it is set to 5.
#' @param gmtFile File path of the gmt file
#' @param isSymbol Boolean value that controls the gene formats. If it is TRUE,
#'      gene format of the gmt file should be symbol. Otherwise, gene format
#'      must be ENTREZ ID.
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#' @return KEGG pathway enrichment results
#'
#' @examples
#' subsetGene <- breastmRNA[1:30,]
#'
#' br_enr<-KeggEnrichment(genes = subsetGene,
#'                        org_assembly='hg19')
#'
#' @export
#'
KeggEnrichment <-
  function(genes,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = c("holm",
                       "hochberg",
                       "hommel",
                       "bonferroni",
                       "BH",
                       "BY",
                       "fdr",
                       "none"),
           min = 5,
           gmtFile = '',
           isSymbol = '',
           isGeneEnrich = '') {
    if (missing(genes)) {
      message("Gene is missing.")
    }
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }
      assembly(org_assembly)
    
    pathTable <- unique(keggPathwayDB(org_assembly))
    genes <- as.data.frame(genes)
    colnames(genes) <- 'g'
    annot <- pathTable[which(pathTable$symbol %in% genes$g),]
    
    
    pathfreq <- as.data.frame(table(annot$pathway))
    pathfreq <- pathfreq[which(pathfreq$Freq > 0),]
    
    
    geneSize = length(unique(pathTable$symbol))
    
    bckfreq <- as.data.frame(table(pathTable$pathway))
    notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
    freq <- merge(pathfreq, notGene, by = "Var1")
    found <- freq$Freq.x
    M <- freq$Freq.y
    n <- rep(length(unique(annot$symbol)), length(M))
    
    pvalues <-
      phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
    
    pAdjust1 <- p.adjust(pvalues, method = pAdjust)
    
    GeneRatio <-
      apply(data.frame(found, n), 1, function(x)
        file.path(x[1], x[2]))
    
    BgRatio <-
      apply(data.frame(M, geneSize), 1, function(x)
        file.path(x[1], x[2]))
    
    enrich <-
      which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
    pathT <- as.character(freq$Var1[enrich])
    gratio <- GeneRatio[enrich]
    bgratio <- BgRatio[enrich]
    padj <- pAdjust1[enrich]
    pval <- pvalues[enrich]
    r <- annot[annot$pathway %in% pathT,]
    enric <- list()
    for (i in seq_along(pathT))
    {
      if (length(which(pathT[i] == r$pathway)) > 0)
      {
        enric <-
          c(enric,
            setNames(list(
              as.character(r[which(pathT[i] == r$pathway),]$symbol)),
                     paste(pathT[i])))
      }
    }
    
    pathways <- data.frame(unique(pathT))
    tmp <- character(length(pathT))
    if (nrow(pathways) > 0) {
      tmp <-
        unlist(lapply(seq_len(nrow(pathways)), function(x)
          tmp[x] <- try(KEGGREST::keggGet(pathT[x])[[1]]$NAME)
        ))
    }
    return(
      methods::new(
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
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'     "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By
#'      default, it is set to 5.
#' @param gmtFile File path of the gmt file
#' @param isSymbol Boolean value that controls the gene formats. If it is TRUE,
#'      gene format of the gmt file should be symbol. Otherwise, gene format
#'      must be ENTREZ ID.
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#'
#' @return Reactome pathway enrichment results
#'
#'
#' @examples
#' br_enr<-reactomeEnrichment(genes = breastmRNA,org_assembly='hg19')
#'
#' @export
reactomeEnrichment <-
  function(genes,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = c("holm",
                       "hochberg",
                       "hommel",
                       "bonferroni",
                       "BH",
                       "BY",
                       "fdr",
                       "none"),
           min = 5,
           gmtFile = '',
           isSymbol = '',
           isGeneEnrich = '') {
    if (missing(genes)) {
      message("Gene is missing.")
    }
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }

      assembly(org_assembly)
    
    pathTable <- unique(reactomePathwayDB(org_assembly))
    genes <- as.data.frame(genes)
    colnames(genes) <- 'g'
    annot <- pathTable[which(pathTable$symbol %in% genes$g),]
    
    
    pathfreq <- as.data.frame(table(annot$pathway))
    pathfreq <- pathfreq[which(pathfreq$Freq > 0),]
    
    geneSize = length(unique(pathTable$symbol))
    
    bckfreq <- as.data.frame(table(pathTable$pathway))
    notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
    freq <- merge(pathfreq, notGene, by = "Var1")
    found <- freq$Freq.x
    M <- freq$Freq.y
    
    n <- rep(length(unique(annot$symbol)), length(M))
    
    pvalues <-
      phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
    
    pAdjust1 <- p.adjust(pvalues, method = pAdjust)
    
    GeneRatio <- apply(data.frame(found, n), 1, function(x)
      file.path(x[1], x[2]))
    
    BgRatio <- apply(data.frame(M, geneSize), 1, function(x)
      file.path(x[1], x[2]))
    
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
    for (i in seq_along(pathT))
    {
      if (length(which(pathT[i] == r$pathway)) > 0)
      {
        enric <-
          c(enric,
            setNames(
              list(as.character(r[which(pathT[i] == r$pathway),]$symbol)),
                     paste(pathT[i])))
      }
    }
    return(
      methods::new(
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


reactomePathwayDB <- function(org_assembly = c("hg19",
                                               "hg38",
                                               "mm10",
                                               "dre10",
                                               "rn6",
                                               "dm6",
                                               "ce11",
                                               "sc3")) {
  xx <- as.list(reactome.db::reactomePATHID2EXTID)
  table1 <- data.frame(pathway = rep(names(xx), lapply(xx, length)),
                       gene = unlist(xx))
  pn <- as.list(reactome.db::reactomePATHID2NAME)
  pn <- data.frame(pathway = rep(names(pn), lapply(pn, length)),
                   name = unlist(pn))
  if (org_assembly == 'hg19' | org_assembly == 'hg38') {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("Install package org.Hs.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-HSA", table1$pathway),]
    pn1 <- pn[grepl("^R-HSA", pn$pathway),]
    symb <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
  }
  if (org_assembly == 'mm10') {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop("Install package org.Mm.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-MMU", table1$pathway),]
    pn1 <- pn[grepl("^R-MMU", pn$pathway),]
    symb <- as.data.frame(org.Mm.eg.db::org.Mm.egSYMBOL)
  }
  if (org_assembly == 'dre10') {
    if (!requireNamespace("org.Dr.eg.db", quietly = TRUE))
      stop("Install package org.Dr.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-DRE", table1$pathway),]
    pn1 <- pn[grepl("^R-DRE", pn$pathway),]
    symb <- as.data.frame(org.Dr.eg.db::org.Dr.egSYMBOL)
  }
  if (org_assembly == 'rn6') {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE))
      stop("Install package org.Rn.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-RNO", table1$pathway),]
    pn1 <- pn[grepl("^R-RNO", pn$pathway),]
    symb <- as.data.frame(org.Rn.eg.db::org.Rn.egSYMBOL)
  }
  if (org_assembly == 'ce11') {
    if (!requireNamespace("org.Ce.eg.db", quietly = TRUE))
      stop("Install package org.Ce.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-CEL", table1$pathway),]
    pn1 <- pn[grepl("^R-CEL", pn$pathway),]
    symb <- as.data.frame(org.Ce.eg.db::org.Ce.egSYMBOL)
  }
  if (org_assembly == 'sc3') {
    if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE))
      stop("Install package org.Sc.sgd.db in order to use this function.")
    
    ty <- table1[grepl("^R-SCE", table1$pathway),]
    pn1 <- pn[grepl("^R-SCE", pn$pathway),]
    symb <-
      AnnotationDbi::select(
        org.Sc.sgd.db::org.Sc.sgd.db,
        keys = as.character(ty$gene),
        columns = c("ENTREZID", "GENENAME"),
        keytype = 'ENTREZID'
      )
  }
  if (org_assembly == 'dm6') {
    if (!requireNamespace("org.Dm.eg.db", quietly = TRUE))
      stop("Install package org.Dm.eg.db in order to use this function.")
    
    ty <- table1[grepl("^R-DME", table1$pathway),]
    pn1 <- pn[grepl("^R-DME", pn$pathway),]
    symb <- as.data.frame(org.Dm.eg.db::org.Dm.egSYMBOL)
  }
  colnames(symb) <- c("gene", "symbol")
  merge1 <- merge(x = pn1,
                  y = ty,
                  by = "pathway",
                  all.y = TRUE)
  path <- merge(merge1, symb, by = "gene")
  
  return(path)
}

#' @importFrom AnnotationDbi mappedkeys
keggPathwayDB <- function(org_assembly = c("hg19",
                                           "hg38",
                                           "mm10",
                                           "dre10",
                                           "rn6",
                                           "dm6",
                                           "ce11",
                                           "sc3")) {
  if (org_assembly == 'hg19' | org_assembly == 'hg38') {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("Install package org.Hs.eg.db in order to use this function.")
    
    kegg <- org.Hs.eg.db::org.Hs.egPATH2EG
    x <- org.Hs.eg.db::org.Hs.egSYMBOL
    prefix <- 'hsa'
  }
  if (org_assembly == 'mm10') {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop("Install package org.Mm.eg.db in order to use this function.")
    
    kegg <- org.Mm.eg.db::org.Mm.egPATH2EG
    x <- org.Mm.eg.db::org.Mm.egSYMBOL
    prefix <- 'mmu'
  }
  if (org_assembly == 'dre10') {
    if (!requireNamespace("org.Dr.eg.db", quietly = TRUE))
      stop("Install package org.Dr.eg.db in order to use this function.")
    
    kegg <- org.Dr.eg.db::org.Dr.egPATH2EG
    x <- org.Dr.eg.db::org.Dr.egSYMBOL
    prefix <- 'dre'
  }
  if (org_assembly == 'rn6') {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE))
      stop("Install package org.Rn.eg.db in order to use this function.")
    
    kegg <- org.Rn.eg.db::org.Rn.egPATH2EG
    x <- org.Rn.eg.db::org.Rn.egSYMBOL
    prefix <- 'rno'
  }
  if (org_assembly == 'ce11') {
    if (!requireNamespace("org.Ce.eg.db", quietly = TRUE))
      stop("Install package org.Ce.eg.db in order to use this function.")
    
    kegg <- org.Ce.eg.db::org.Ce.egPATH2EG
    x <- org.Ce.eg.db::org.Ce.egSYMBOL
    prefix <- 'cel'
  }
  if (org_assembly == 'sc3') {
    if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE))
      stop("Install package org.Sc.sgd.db in order to use this function.")
    
    kegg <- org.Sc.sgd.db::org.Sc.sgdPATH2ORF
    x <- org.Sc.sgd.db::org.Sc.sgdGENENAME
    prefix <- 'sce'
  }
  if (org_assembly == 'dm6') {
    if (!requireNamespace("org.Dm.eg.db", quietly = TRUE))
      stop("Install package org.Dm.eg.db in order to use this function.")
    
    kegg <- org.Dm.eg.db::org.Dm.egPATH2EG
    x <- org.Dm.eg.db::org.Dm.egGENENAME
    prefix <- 'dme'
  }
  mapped <- mappedkeys(kegg)
  kegg2 <- as.list(kegg[mapped])
  pathTable <-
    data.frame(pathway = paste0(prefix, rep(names(kegg2),
                                            lapply(kegg2, length))),
               gene = unlist(kegg2))
  x <- as.data.frame(x)
  colnames(x) <- c("gene", "symbol")
  pathTable <- merge(pathTable, x, by = "gene")
  return(pathTable)
}

WikiPathwayDB <- function(org_assembly = c("hg19",
                                           "hg38",
                                           "mm10",
                                           "dre10",
                                           "rn6",
                                           "dm6",
                                           "ce11",
                                           "sc3")) {
  if (org_assembly == 'hg19' | org_assembly == 'hg38')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens",
                                            format = "gmt")
  if (org_assembly == 'mm10')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Mus musculus",
                                            format = "gmt")
  if (org_assembly == 'dre10')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Danio rerio",
                                            format = "gmt")
  if (org_assembly == 'rn6')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(organism = "Rattus norvegicus",
                                            format = "gmt")
  if (org_assembly == 'dm6')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(
        organism = "Drosophila melanogaster",format = "gmt")
  if (org_assembly == 'ce11')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(
        organism = "Caenorhabditis elegans",format = "gmt")
  if (org_assembly == 'sc3')
    wp.gmt <-
      rWikiPathways::downloadPathwayArchive(
        organism = "Saccharomyces cerevisiae",format = "gmt")
  
  gmtFile <- convertGMT(gmtName = wp.gmt,org_assembly = org_assembly,isSymbol = FALSE)
  tmp <-
    do.call(rbind, strsplit(as.character(gmtFile$pathTerm), '%'))
  gmtFile <-
    unique(
      data.frame(
        Entrez = gmtFile$Entrez,
        gene = gmtFile[,3],
        pathID = tmp[, 3],
        pathTerm = tmp[, 1]
      )
    )
  return(gmtFile)
}

#' WikiPathways Enrichment
#'
#' @param genes Input genes
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By
#'      default, it is set to 5.
#' @param gmtFile File path of the gmt file
#' @param isSymbol Boolean value that controls the gene formats. If it is TRUE,
#'      gene format of the gmt file should be symbol. Otherwise, gene format
#'      must be ENTREZ ID.
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#' @return Wiki Pathway Enrichment
#' 
#' @export
#' 
#' 
WikiEnrichment <- function(genes,
                           org_assembly = c("hg19",
                                            "hg38",
                                            "mm10",
                                            "dre10",
                                            "rn6",
                                            "dm6",
                                            "ce11",
                                            "sc3"),
                           pCut = 0.05,
                           pAdjCut = 0.05,
                           pAdjust = c("holm",
                                       "hochberg",
                                       "hommel",
                                       "bonferroni",
                                       "BH",
                                       "BY",
                                       "fdr",
                                       "none"),
                           min = 5,
                           gmtFile = '',
                           isSymbol = '',
                           isGeneEnrich = '') {
  if (missing(genes)) {
    message("Gene is missing.")
  }
  if (missing(org_assembly)) {
    message("Assembly version is missing. ")
  }
  pathTable <- unique(WikiPathwayDB(org_assembly))
  genes <- as.data.frame(genes)
  colnames(genes) <- 'g'
  annot <- pathTable[which(pathTable$gene %in% genes$g),]
  
  pathfreq <- as.data.frame(table(annot$pathID))
  pathfreq <- pathfreq[which(pathfreq$Freq > 0),]
  
  geneSize = length(unique(pathTable$gene))
  bckfreq <- as.data.frame(table(pathTable$pathID))
  notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
  freq <- merge(pathfreq, notGene, by = "Var1")
  found <- freq$Freq.x
  M <- freq$Freq.y
  
  n <- rep(length(unique(annot$gene)), length(M))
  
  pvalues <-
    phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
  
  pAdjust1 <- p.adjust(pvalues, method = pAdjust)
  GeneRatio <-
    apply(data.frame(found, n), 1, function(x)
      file.path(x[1], x[2]))
  
  BgRatio <-
    apply(data.frame(M, geneSize), 1, function(x)
      file.path(x[1], x[2]))
  
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
  for (i in seq_along(pathT))
  {
    if (length(which(pathT[i] == r$pathID)) > 0)
    {
      enric <-
        c(enric, setNames(
          list(as.character(r[which(pathT[i] == r$pathID),]$gene)),
                          paste(pathT[i])))
    }
  }
  
  return(
    methods::new(
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

#' For a given gmt file of a specific pathway database, pathway enrichment
#' can be performed. Function supports Entrez ID and symbol based gmt file.
#'
#' @param genes Input genes
#' @param gmtFile File path of the gmt file
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param isSymbol Boolean value that controls the gene formats. If it is TRUE,
#'      gene format of the gmt file should be symbol. Otherwise, gene format
#'      must be ENTREZ ID.
#' @param min Minimum number of genes that are required for enrichment. By
#'      default, it is set to 5.
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#' @return Pathway Enrichment
#' 
pathwayEnrichment <- function(genes,
                              gmtFile,
                              org_assembly = c("hg19",
                                               "hg38",
                                               "mm10",
                                               "dre10",
                                               "rn6",
                                               "dm6",
                                               "ce11",
                                               "sc3"),
                              pCut = 0.05,
                              pAdjCut = 0.05,
                              pAdjust = c("holm",
                                          "hochberg",
                                          "hommel",
                                          "bonferroni",
                                          "BH",
                                          "BY",
                                          "fdr",
                                          "none"),
                              isSymbol,
                              min = 5,
                              isGeneEnrich = FALSE) {
  if (missing(genes))
    message("Gene is missing.")
  
  if (missing(org_assembly)) {
    message("Assembly version is missing.")
  }
  if (missing(gmtFile))
    message("GMT file is missing.")
  
  if (missing(isSymbol))
    message("GMT gene format is missing.")

    assembly(org_assembly)
  
  genes <- as.data.frame(genes)
  colnames(genes) <- 'g'
  
  if (!isGeneEnrich) {
    pathTable <-
      unique(convertGMT(gmtName = gmtFile, isSymbol = isSymbol,
                        org_assembly =org_assembly))
  }
  else{
    pathTable <- geneListEnrich(f = gmtFile, isSymbol = isSymbol)
  }
  annot <- pathTable[which(pathTable$symbol %in% genes$g),]
  pathfreq <- as.data.frame(table(annot$pathTerm))
  pathfreq <- pathfreq[which(pathfreq$Freq > 0),]
  
  
  if (!isGeneEnrich)
    geneSize = length(unique(pathTable$symbol))
  else
    geneSize = length(pkg.env$ucsc)
  
  
  bckfreq <- as.data.frame(table(pathTable$pathTerm))
  
  notGene <- bckfreq[bckfreq$Var1 %in% pathfreq$Var1,]
  freq <- merge(pathfreq, notGene, by = "Var1")
  found <- freq$Freq.x
  M <- freq$Freq.y
  
  n <- rep(length(unique(annot$symbol)), length(M))
  
  pvalues <-
    phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
  
  pAdjust1 <- p.adjust(pvalues, method = pAdjust)
  
  GeneRatio <- apply(data.frame(found, n), 1, function(x)
    file.path(x[1], x[2]))
  
  BgRatio <- apply(data.frame(M, geneSize), 1, function(x)
    file.path(x[1], x[2]))
  
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
  
  for (i in seq_along(pathT))
  {
    if (length(which(pathT[i] == r$pathTerm)) > 0)
      enric <-
        c(enric, setNames(
          list(as.character(r[which(pathT[i] == r$pathTerm),]$symbol)),
                          paste(pathT[i])))
  }
  return(
    methods::new(
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

#' Convert gmt formatted pathway file to the Pathway ID, Entrez, symbol
#'  formatted data frame
#'
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file.
#'     If it is set as TRUE, gene format of the gmt file should be symbol.
#'     Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#'
#' @return return data frame
#'
#'
#' @importFrom biomaRt getBM
#' @importFrom reshape2 melt
#' 
#' @export
#' 
convertGMT <- function(gmtName, org_assembly,
                       isSymbol = FALSE) {
  types <-
    rbind(
      c("hg19","Hs.eg.db"),
      c("hg38", "Hs.eg.db"),
      c( "mm10", "Mm.eg.db"),
      c("dre10","Dr.eg.db"),
      c("rn6", "Rn.eg.db" ),
      c( "sc3", "Sc.sgd.db"),
      c("dm6","Dm.eg.db"),
      c("ce11", "Ce.eg.db")
    )
  index = which(org_assembly == types[, 1])
  yy <- as.name(paste0("org.", types[index, 2]))
  
  if (!requireNamespace(yy, quietly = TRUE))
    stop("Install package ", yy, " in order to use this function.")
  else
    lapply(deparse(substitute(yy)), require, character.only = TRUE)
  
  x <- scan(gmtName, what = "", sep = "\n")
  x <- strsplit(x, '\t')
  max.col <- max(vapply(x, length, FUN.VALUE = integer(1)))
  cn <- paste("V", seq_len(max.col), sep = "")
  
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
    # output <-
    #   getBM(
    #     attributes = c('entrezgene_id', 'hgnc_symbol'),
    #     filters = "entrezgene_id",
    #     values = f$value,
    #     mart = pkg.env$mart
    #   )
   output<- AnnotationDbi::select(eval(yy), keys =as.character(f$value),
                                  columns=c("ENTREZID", "SYMBOL"),
                                  keytype="ENTREZID")
   # f[, 3] <- output[match((f$value), output$entrezgene_id), 2]
    f[, 3] <- output[match((f$value), output$ENTREZID), 2]
    colnames(f) <- c('pathTerm', 'Entrez', 'symbol')
  }
  else{
     # output <-
     #   getBM(
     #     attributes = c('entrezgene_id', 'hgnc_symbol'),
     #     filters = "hgnc_symbol",
     #     values = f$value,
     #     mart = pkg.env$mart
     #   )
    output<- AnnotationDbi::select(eval(yy), keys =as.character(f$value),
                                  columns=c("ENTREZID", "SYMBOL"),
                                 keytype="SYMBOL")
    #f[, 3] <- output[match((f$value), output$hgnc_symbol), 1]
    f[, 3] <- output[match((f$value), output$SYMBOL), 2]
    colnames(f) <- c('pathTerm', 'symbol', 'Entrez')
  }
  f <- na.omit(f)
  return(f)
}

geneListEnrich <- function(f, isSymbol = FALSE)
{
  if (ncol(f) == 1)
    f <- cbind("NotAva", f)
  
  colnames(f) <- c("Annot", "value")
  
  if (!isSymbol) {
    output <-
      getBM(
        attributes = c('entrezgene_id', 'hgnc_symbol'),
        filters = "entrezgene_id",
        values = f$value,
        mart = pkg.env$mart
      )
    f[, 3] <- output[match(f$value, output$entrezgene_id), 2]
    colnames(f) <- c('pathTerm', 'Entrez', 'symbol')
  }
  else{
    output <-
      getBM(
        attributes = c('entrezgene_id', 'hgnc_symbol'),
        filters = "hgnc_symbol",
        values = f$value,
        mart = pkg.env$mart
      )
    f[, 3] <- output[match(f$value, output$hgnc_symbol), 1]
    colnames(f) <- c('pathTerm', 'symbol', 'Entrez')
  }
  return(f)
}
