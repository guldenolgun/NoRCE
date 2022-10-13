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
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param GOtype Hierarchical category of the GO ontology. Possible values are
#'      "BP"(default), "CC", "MF".
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are
#'      "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of gene that are required for enrichment. By
#'       default, it is set to 5
#' @param enrichTest Types of enrichment methods to perform enrichment
#'       analysis. Possible values are "hyper"(default), "binom", "fisher",
#'       "chi".
#' @param backG The set of genes that tested against to the input
#'       (background gene)
#' @param backGType Type of the background gene. If miRNA gene set is used for
#'        background gene, backGType should be set to the 'mirna'
#'
#' @return GO enrichment results
#'
#' @examples
#' subsetGene <- breastmRNA[1:30,]
#' breastEnr <- goEnrichment(genes = subsetGene,
#'                           org_assembly = 'hg19',
#'                           GOtype = 'MF',
#'                           min = 2)
#'
#' @importFrom stats chisq.test cor cor.test fisher.test na.omit p.adjust
#' @importFrom stats pbinom phyper reorder setNames var
#'
#' @export
goEnrichment <-
  function(genes,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           GOtype = c("BP", "CC", "MF"),
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
           backG = '',
           backGType = 'pc_gene',
           enrichTest = c("hyper", "binom", "fisher", "chi")) {
    if (missing(org_assembly)) {
      message("Genome assembly version is missing.")
    }
    if (missing(genes)) {
      message("Genes are missing. Expected input: FOXP2 SOX4 HOXC6")
    }
    assembly(org_assembly)
    #  annot <- unique(annGO(genes, GOtype, org_assembly))
    goData <- annGO(genes, GOtype, org_assembly)
    annot <- goData[[2]]
    uniqueGO <- annot$GOID[!duplicated(annot$GOID)]
    
    if(is.null(dim(backG)[1])){
      gofreq <- as.data.frame(table(annot$GOID))
    }else{
      
    }
   
    notGene <-
      getBackGenes(
        backgroundGene = backG,
        all =  goData[[1]],
        GOtype = GOtype,
        gofreq = gofreq,
        org_assembly = org_assembly,
        type = backGType
      )
    
    if (nrow(gofreq) > nrow(notGene)) {
      gofreq <- gofreq[gofreq$Var1 %in% notGene$Var1,]
    }
    
    freq <- merge(gofreq, notGene, by = "Var1")
    found <- freq$Freq.x
    
    ifelse(is.null(dim(backG)[1]),
           geneSize <- length(unique(goData[[1]]$Gene)),
           geneSize <- dim(unique(backG))[1])
    
    M <- freq$Freq.y
    n <- rep(length(unique(annot$Gene)), length(M))
    
    ifelse(
      pkg.env$enrichTest == "binom",
      pvalues <-
        2 * (1 - pbinom(found, n, pCut)),
      ifelse(
        pkg.env$enrichTest == "fisher",
        pvalues <-
          fisher.test(matrix(c(
            found, (M - found), (n - found), (geneSize - M - n + found)
          ), 2, 2), alternative = 'greater')$p.value,
        ifelse(
          pkg.env$enrichTest == "chi",
          pvalues <-
            chisq.test(matrix(c(
              found, (M - found), (n - found), (geneSize - M - n + found)
            )))$p.value,
          pvalues <-
            phyper(found - 1, M, geneSize - M, n, lower.tail = FALSE)
        )
      )
    )
    
    pAdjust1 <- p.adjust(pvalues, method = pAdjust)
    
    
    
    GeneRatio <-
      apply(data.frame(found, n), 1, function(x)
        file.path(x[1], x[2]))
    
    BgRatio <-
      apply(data.frame(M, geneSize), 1, function(x)
        file.path(x[1], x[2]))
    
    enrich <-
      which (pvalues <= pCut & pAdjust1 <= pAdjCut & found >= min)
    goT <- as.character(freq$Var1[enrich])
    gratio <- GeneRatio[enrich]
    bgratio <- BgRatio[enrich]
    padj <- pAdjust1[enrich]
    pval <- pvalues[enrich]
    r <- annot[annot$GOID %in% goT, seq_len(2)]
    
    enric <- list()
    for (i in seq_along(goT))
    {
      if (length(which(goT[i] == r$GOID)) > 0)
      {
        enric <-
          c(enric,
            setNames(list(as.character(r[which(goT[i] == r$GOID), ]$Gene)),
                     paste(goT[i])))
      }
    }
    goTe <-
      as.character(annot[match(goT, annot[, 'GOID']),]$GOTerm)
    
    return(
      methods::new(
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
           all,
           GOtype = c("BP", "CC", "MF"),
           gofreq,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           type = c('pc_gene', 'mirna')) {
    backgroundGene = as.data.frame(backgroundGene)
    if (length(backgroundGene[,1]) == 1) {
      bckfreq <- as.data.frame(table(all$GOID))
    }
    else{
      colnames(backgroundGene) = 'bg'
      
      if (type == 'mirna') {
        a <-
          as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "",
                             backgroundGene$bg))
        colnames(a) <- 'gene'
        geneTargetLoc <-
          convertGeneID(
            genetype = "mirna",
            genelist = a$gene,
            org_assembly = org_assembly
          )
        backgroundGene <-
          getUCSC(geneTargetLoc, 10000, 10000, org_assembly)
        colnames(backgroundGene) = 'bg'
      }
      goData <-annGO(backgroundGene$bg, GOtype, org_assembly)
      annot <- goData[[2]]
      
      bckfreq <- as.data.frame(table(annot$GOID))
    }
    bb <- bckfreq[bckfreq$Var1 %in% gofreq$Var1,]
    
    return(bb)
  }