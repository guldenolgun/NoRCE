#' An S4 class to represent enrichment
#'
#' @slot ID factor
#' @slot Term factor
#' @slot geneList factor
#' @slot ncGeneList factor
#' @slot pvalue factor
#' @slot pAdj factor
#' @slot GeneRatio factor
#' @slot BckRatio factor
#'
#'
#'
#' @export
setClass(
  Class = "NoRCE",
  slots = c(
    ID = "character",
    Term = "character",
    geneList = "list",
    pvalue = "numeric",
    pAdj = "numeric",
    GeneRatio = "character",
    BckRatio = "character",
    ncGeneList = "list"
  )
)
#' Perform enrichment analysis of the given genes
#'
#' @param genes Set of input genes. Supported format HUGO.
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param GOtype Hierarchical category of the GO ontology. Possible values are "BP"(default), "CC", "MF".
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of gene that are required for enrichment. By default, it is set to 5
#' @param enrichTest Types of enrichment methods to perform enrichment analysis. Possible values are "hyper"(default), "binom", "fisher", "chi".
#' @param slim Boolean value stating whether set of annotation should be performed for high level GO terms (GO slim)
#' @param backG The set of genes that tested against to the input(background gene)
#' @param backGType Type of the background gene. If miRNA gene set is used for background gene, backGType should be set to the 'mirna'
#'
#' @return GO enrichment results
#'
#' @examples
#' subsetGene <- breastmRNA[1:100,]
#' \dontrun{
#' breastEnr <- goEnrichment(genes = subsetGene,hg = 'hg19',GOtype = 'MF')
#'
#' #Enriched genes with Generic GO-term
#' breastEnr <- goEnrichment(genes = subsetGene,hg = 'hg19',GOtype = 'MF', slim = TRUE,min = 3)
#' }
#' @importFrom stats chisq.test cor cor.test fisher.test na.omit p.adjust pbinom phyper reorder setNames var
#'
#' @export
goEnrichment <-
  function(genes,
           hg,
           GOtype = "BP",
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           min = 5,
           backG = '',
           backGType = 'pc_gene',
           enrichTest = "hyper",
           slim = FALSE) {
    if (missing(hg)) {
      message(
        "Genome assembly version is missing. Possible assemblies are 'mm10' for mouse, 'dre10' for zebrafish, 'rn6' for rat, 'dm6' for fruit fly, 'ce11' for worm, 'hg19' and 'hg38' for human."
      )
      assembly(hg)
    }
    if (missing(genes)) {
      message("Genes are missing. Expected input: FOXP2 SOX4 HOXC6")
    }
    annot <- unique(annotate(genes, GOtype, hg, slim = slim))
    uniqueGO <- annot$GOID[!duplicated(annot$GOID)]


    gofreq <- as.data.frame(table(annot$GOID))
    notGene <-
      getBackGenes(
        backgroundGene = backG,
        GOtype = GOtype,
        gofreq = gofreq,
        hg = hg,
        type = backGType
      )

    if (dim(gofreq)[1] > dim(notGene)[1]) {
      gofreq <- gofreq[gofreq$Var1 %in% notGene$Var1, ]
    }

    freq<-merge(gofreq,notGene, by = "Var1")
    found <- freq$Freq.x
    if (backG == '')
      geneSize = length(unique(go$X2))
    else
      geneSize = length(unique(backG))


    M <- freq$Freq.y
    n <- rep(length(unique(annot$Gene)), length(M))


    if (enrichTest == "binom") {
      pvalues <- 2 * (1 - pbinom(found, n, pCut))
    }
    if (enrichTest == "fisher") {
      pvalues <-
        fisher.test(matrix(c(
          found, (M - found), (n - found), (geneSize - M - n + found)
        ), 2, 2), alternative = 'greater')$p.value
    }
    if (enrichTest == "chi") {
      pvalues <-
        chisq.test(matrix(c(
          found, (M - found), (n - found), (geneSize - M - n + found)
        )))$p.value
    }
    else {
      pvalues <- phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
    }


    pAdjust1 <- p.adjust(pvalues, method = pAdjust)



    GeneRatio <-
      apply(data.frame(found, n), 1, function(x)
        paste(x[1], "/", x[2], sep = "", collapse = ""))

    BgRatio <-
      apply(data.frame(M, geneSize), 1, function(x)
        paste(x[1], "/", x[2], sep = "", collapse = ""))

    enrich <-
      which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
    goT <- as.character(freq$Var1[enrich])
    gratio <- GeneRatio[enrich]
    bgratio <- BgRatio[enrich]
    padj <- pAdjust1[enrich]
    pval <- pvalues[enrich]
    r <- annot[annot$GOID %in% goT, 1:2]

    enric <- list()
    for (i in 1:length(goT))
    {
      if (length(which(goT[i] == r$GOID)) > 0)
      {
        enric <-
          c(enric, setNames(list(as.character(r[which(goT[i] == r$GOID), ]$Gene)), paste(goT[i])))
      }
    }
    goTe <-
      as.character(annot[match(goT, annot[, 'GOID']), ]$GOTerm)

    return(
      new(
        "NoRCE",
        ID = goT,
        Term = goTe,
        geneList = enric,
        pvalue = pval,
        pAdj = padj,
        GeneRatio = gratio,
        BckRatio = bgratio
      )
    )
  }

getBackGenes <-
  function(backgroundGene,
           GOtype,
           gofreq,
           hg,
           type = 'pc_gene') {
    if (backgroundGene == '') {
      bckfreq <- as.data.frame(table(go$X1))
    }
    else{
      backgroundGene = as.data.frame(backgroundGene)
      colnames(backgroundGene) = 'bg'

      if (type == 'mirna') {
        a <-
          as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "", backgroundGene$bg))
        colnames(a) <- 'gene'
        geneTargetLoc <-
          convertGeneID(genetype = "mirna",
                        genelist = a$gene,
                        hg = hg)
        backgroundGene <- getUCSC(geneTargetLoc, 10000, 10000, hg)
        colnames(backgroundGene) = 'bg'
      }
      annot <- unique(annotate(backgroundGene$bg, GOtype, hg))
      uniqueGO <- annot$GOID[!duplicated(annot$GOID)]

      bckfreq <- as.data.frame(table(annot$GOID))
    }
    bb <- bckfreq[bckfreq$Var1 %in% gofreq$Var1, ]

    return(bb)
  }
