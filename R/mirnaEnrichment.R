op <- options(warn = (-1))
options(readr.num_columns = 0)

#' GO term enrichments of the microRNA genes with mRNAs that fall in the
#' given upstream/downstream regions of the microRNA genes
#'
#' @param gene Input microRNA gene. It supports both pre-miRNA and mature
#'     miRNA, however, when target prediction is performed (target= TRUE),
#'     miRNA genes should be mature.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38"
#'     for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'     considered in the analysis
#' @param target Boolean value shows whether miRNA target prediction should
#'     be performed
#' @param backGenes The set of genes that tested against to the input
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'     performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'     any new TAD regions can be used for the analysis. TAD regions must be
#'     formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'     'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly,
#'     respectively.
#' @param backGType Type of the background gene. If miRNA gene set is used for
#'     background gene, backGType should be set to the 'mirna'
#' @param express Boolean variable whether co-expression analysis is performed.
#'     If this option is set to TRUE, co-expression analysis will be performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with
#'     custom data will be performed. When this option is set, exp1 and exp2
#'     parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for
#'     correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC, CHOL,
#'     COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LGG,
#'     LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES,
#'     TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows
#'     must be patients. If gene names are provided as header, no need to
#'     redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows
#'     must be patients. If gene names are provided as header, no need to
#'     redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not
#'     provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not
#'     provided, column name of the exp2 data will be taken.
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'     of the co-expression analysis and the other analysis should be
#'     considered
#'
#' @return MiRNA GO term enrichment object for the given input
#' 
#' @examples
#' \dontrun{
#' subsetGene <- brain_mirna[1:30,]
#'
#' miGO <-mirnaGOEnricher(gene=subsetGene,
#'                        org_assembly='hg19',
#'                        near = TRUE,
#'                        target = FALSE) }
#' @export mirnaGOEnricher
mirnaGOEnricher <-
  function(gene,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           near = FALSE,
           target = FALSE,
           backGenes = '',
           backGType = 'pc_gene',
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           isUnionCorGene = FALSE,
           databaseFile = '') {
    if (missing(gene)) {
      message("Gene is missing.")
    }
    
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }
    if (!is.data.frame(gene) &
        !is.character(gene) & !is.factor(gene))
      message("Type of the gene should be data.frame or character")
    
    
    assembly(org_assembly)
    
    
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    
    if (target) {
      targetResult <- predictmiTargets(gene = gene$genes,
                                       type = "mirna",
                                       org_assembly = org_assembly)
      if (is.null(targetResult))
      {
        message("There is no target!")
        return(NULL)
      }
      targetResult <- unique(targetResult)
      
      geneTargetLoc <-
        convertGeneID(genetype = "Ensembl_gene",
                      genelist = targetResult[, 3],
                      org_assembly = org_assembly)
    }
    
    a <-
      as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "",
                         gene$genes))
    colnames(a) <- 'genes'
    a <- unique(rbind(a, gene$genes))
    geneLoc <-
      convertGeneID(
        genetype = "mirna",
        genelist = a$genes,
        org_assembly = org_assembly
      )
    
    
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        miNearGene_temp <-
          getUCSC(
            geneLoc,
            pkg.env$upstream,
            pkg.env$downstream,
            org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          miNearGene_temp <-
            getNearToExon(
              geneLoc,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            ),
          miNearGene_temp <-
            getNearToIntron(
              geneLoc,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            )
        )
      )
      
      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      if (target) {
        geneL <- IRanges::findOverlapPairs(geneLoc_temp, geneTargetLoc)
        geneLo <- pintersect(geneL, ignore.strand = TRUE)
        miNearGene <-
          getUCSC(
            bedfile = geneLo,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
      else{
        miNearGene <- miNearGene_temp
      }
    }
    
    else{
      if (target) {
        miNearGene <- targetResult[, 2]
      }
      else{
        miNearGene <- gene
      }
    }
    
    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          tad = TAD,
          cellline = pkg.env$cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = pkg.env$upstream,
          downstream = pkg.env$downstream
        )
      if (near | target)
        miNearGene <-
          as.data.frame(intersect(unlist(miNearGene), unlist(tadGene)))
      else
        miNearGene <- tadGene
    }
    
    if (express) {
      if (!isCustomExp) {
        nearG <- corrbased(
          mirnagene = a$genes,
          cancer = cancer,
          minAbsCor = pkg.env$minAbsCor,
          databaseFile = databaseFile
        )
        d <- nearG[which(a$genes %in% nearG$mirna_base), ]
        
        if (!isUnionCorGene)
          miNearGene <- intersect(unlist(miNearGene), d$feature)
        else
          miNearGene <- union(unlist(miNearGene), d$feature)
      }
      else{
        nearG <-  calculateCorr(
          exp1 = exp1,
          exp2 = exp2,
          label1 = label1 ,
          label2 = label2,
          corrMethod = pkg.env$corrMethod,
          varCutoff = pkg.env$varCutoff,
          corCutoff = pkg.env$minAbsCor,
          pcut = pkg.env$pcut,
          alternate = pkg.env$alternate,
          conf = pkg.env$conf
        )
        tt <-
          lapply(seq_len(nrow(a)), function(x)
            unlist(which(
              nearG$firstExp %in% tolower(a$genes[x])
            )))
        nearG <- nearG[unlist(tt), ]
        if (!isUnionCorGene)
          miNearGene <-
          intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      methods::new(
        "NoRCE",
        ID = '',
        Term = '',
        geneList = list(),
        pvalue = 0,
        GeneRatio = '',
        BckRatio = ''
      )
    }
    else{
      miEnrich <-
        goEnrichment(
          genes = miNearGene,
          GOtype = pkg.env$GOtype,
          org_assembly = org_assembly,
          pCut = pkg.env$pCut,
          pAdjCut = pkg.env$pAdjCut,
          pAdjust = pkg.env$pAdjust,
          backG = backGenes,
          backGType = backGType,
          min = pkg.env$min,
          enrichTest = pkg.env$enrichTest
        )
      if (length(miEnrich@Term)) {
        miEnrich@ncGeneList <-
          commonGene(
            mrnaobject = miEnrich,
            org_assembly = org_assembly,
            downstream = pkg.env$downstream,
            upstream = pkg.env$upstream,
            inputGene = rbind(a, gene),
            inGeneType = 'mirna'
          )
      }
      return(miEnrich)
    }
  }

#' Pathway enrichments of the microRNA genes with mRNAs that fall in the
#' given upstream/downstream regions of the microRNA genes
#'
#' @param gene Input microRNA gene. It supports both pre-miRNA and mature miRNA,
#'     however, when target prediction is performed(target= TRUE), miRNA genes
#'     should be mature.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'     "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'     considered in the analysis
#' @param target Boolean value shows whether miRNA target prediction should
#'     be performed
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'     performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'     any new TAD regions can be used for the analysis. TAD regions must be
#'     formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'     'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'     assembly, respectively.
#' @param gmtName Custom pathway gmt file
#' @param express Boolean variable whether co-expression analysis is performed.
#'      If this option is set to TRUE, co-expression analysis will be
#'      performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with
#'      custom data will be performed. When this option is set, exp1 and exp2
#'      parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for
#'      correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC,
#'      CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC,
#'      KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM,
#'      STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'      redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'      redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not
#'      provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not
#'      provided, column name of the exp2 data will be taken.
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'      of the co-expression analysis and the other analysis should be
#'      considered
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#' @return MiRNA pathway enrichment object for the given input
#'
#' @examples
#' miPath <- mirnaPathwayEnricher(gene = brain_mirna,
#'                                org_assembly = 'hg19',
#'                                near = TRUE)
#'
#' @export
mirnaPathwayEnricher <-
  function(gene,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           near = FALSE,
           target = FALSE,
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           gmtName = '',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           isUnionCorGene = FALSE,
           databaseFile,
           isGeneEnrich = FALSE) {
    if (missing(gene)) {
      message("Gene is missing.")
    }
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }
    
    assembly(org_assembly)
    
    
    if (!is.data.frame(gene) &
        !is.character(gene) & !is.factor(gene))
      message("Type of the gene should be data.frame or character")
    
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    
    if (target) {
      targetResult <- predictmiTargets(gene = gene$genes,
                                       type = "mirna",
                                       org_assembly = org_assembly)
      if (is.null(targetResult))
      {
        message("There is no target!")
        return(NULL)
      }
      
      targetResult <- unique(targetResult)
      
      geneTargetLoc <-
        convertGeneID(genetype = "Ensembl_trans",
                      genelist = targetResult,
                      org_assembly = org_assembly)
    }
    
    a <-
      as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "",
                         gene$genes))
    colnames(a) <- 'genes'
    a <- unique(rbind(a, gene$genes))
    geneLoc <-
      convertGeneID(
        genetype = "mirna",
        genelist = a$genes,
        org_assembly = org_assembly
      )
    
    
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        miNearGene_temp <-
          getUCSC(
            geneLoc,
            pkg.env$upstream,
            pkg.env$downstream,
            org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          miNearGene_temp <-
            getNearToExon(
              geneLoc,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            ),
          miNearGene_temp <-
            getNearToIntron(
              geneLoc,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            )
        )
      )
      
      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      if (target) {
        geneL <- IRanges::findOverlapPairs(geneLoc_temp, geneTargetLoc)
        geneLo <- pintersect(geneL, ignore.strand = TRUE)
        miNearGene <-
          getUCSC(
            bedfile = geneLo,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
      else{
        miNearGene <- miNearGene_temp
      }
    }
    
    else{
      if (target) {
        miNearGene <- targetResult[, 2]
      }
      else{
        miNearGene <- gene
      }
    }
    
    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          tad = TAD,
          cellline = pkg.env$cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = pkg.env$upstream,
          downstream = pkg.env$downstream
        )
      if (near | target)
        miNearGene <-
          as.data.frame(intersect(unlist(miNearGene), unlist(tadGene)))
      else
        miNearGene <- tadGene
    }
    
    if (express) {
      if (!isCustomExp) {
        nearG <- corrbased(
          mirnagene = a$genes,
          cancer = cancer,
          minAbsCor = pkg.env$minAbsCor,
          databaseFile = databaseFile
        )
        d <- nearG[which(a$genes %in% nearG$mirna_base), ]
        
        if (!isUnionCorGene)
          miNearGene <- intersect(unlist(miNearGene), d$feature)
        else
          miNearGene <- union(unlist(miNearGene), d$feature)
      }
      else{
        nearG <-  calculateCorr(
          exp1 = exp1,
          exp2 = exp2,
          label1 = label1 ,
          label2 = label2,
          corrMethod = pkg.env$corrMethod,
          varCutoff = pkg.env$varCutoff,
          corCutoff = pkg.env$minAbsCor,
          pcut = pkg.env$pcut,
          alternate = pkg.env$alternate,
          conf = pkg.env$conf
        )
        tt <-
          lapply(seq_len(nrow(a)), function(x)
            unlist(which(
              nearG$firstExp %in% tolower(a$genes[x])
            )))
        nearG <- nearG[unlist(tt), ]
        if (!isUnionCorGene)
          miNearGene <-
          intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      methods::new(
        "NoRCE",
        ID = '',
        Term = '',
        geneList = list(),
        pvalue = 0,
        pAdj = 0,
        GeneRatio = '',
        BckRatio = ''
      )
    }
    else{
      ifelse(
        pkg.env$pathwayType == 'kegg',
        pth <- 1,
        ifelse(
          pkg.env$pathwayType == 'reactome',
          pth <- 2,
          ifelse(pkg.env$pathwayType == 'wiki',
                 pth <- 3,
                 pth <- 4)
        )
      )
      
      funclist <- list(KeggEnrichment,
                       reactomeEnrichment,
                       WikiEnrichment,
                       pathwayEnrichment)
      
      miEnrich <- funclist[[pth]](
        genes = miNearGene,
        gmtFile = gmtName,
        org_assembly = org_assembly,
        pCut = pkg.env$pCut,
        pAdjCut = pkg.env$pAdjCut,
        pAdjust = pkg.env$pAdjust,
        isSymbol = pkg.env$isSymbol,
        min = pkg.env$min,
        isGeneEnrich = isGeneEnrich
      )
      
      if (length(miEnrich@Term) > 0)
      {
        miEnrich@ncGeneList <-
          commonGene(
            mrnaobject = miEnrich,
            org_assembly = org_assembly,
            downstream = pkg.env$downstream,
            upstream = pkg.env$upstream,
            inputGene = rbind(a, gene),
            inGeneType = 'mirna'
          )
      }
      
      return(miEnrich)
    }
  }


#' GO enrichments of the microRNA regions with mRNAs that fall in the given
#' upstream/downstream regions of the microRNA genes
#'
#' @param region MiRNA region in a bed format
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param backG The set of genes that tested against to the input
#' @param near Boolean value presents whether cis-neighbourhood should be
#'      considered in the analysis
#' @param target Boolean value shows whether miRNA target prediction should
#'      be performed
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'      performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'      assembly, respectively.
#' @param backGType Type of the background gene. If miRNA gene set is used for
#'      background gene, backGType should be set to the 'mirna'
#' @param express Boolean variable whether co-expression analysis is performed.
#'      If this option is set to TRUE, co-expression analysis will be
#'      performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with
#'      custom data will be performed. When this option is set, exp1 and exp2
#'      parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for
#'      correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC,
#'      CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC,
#'      KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM,
#'      STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'      redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'      redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not
#'      provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not
#'      provided, column name of the exp2 data will be taken.
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'      of the co-expression analysis and the other analysis should be
#'      considered
#'
#' @return MiRNA GO enrichment object for the given input
#' 
#'
#'@examples
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#'
#' a<- mirnaRegionGOEnricher(region = regionNC,
#'                           org_assembly = 'hg19',
#'                           near = TRUE)
#'
#' @export
mirnaRegionGOEnricher <-
  function(region,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           near = FALSE,
           target = FALSE,
           backG = '',
           backGType = 'pc-genes',
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           isUnionCorGene = FALSE,
           databaseFile) {
    if (missing(region)) {
      message("Region of interest is missing.")
    }
    
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }
    
    assembly(org_assembly)
    
    
    if (target) {
      genes <-
        getUCSC(
          bedfile = region,
          downstream = 0,
          upstream = 0,
          org_assembly = org_assembly
        )
      targetResult <- predictmiTargets(gene = genes,
                                       type = "NCBI",
                                       org_assembly = org_assembly)
      if (is.null(targetResult))
      {
        message("There is no target!")
        return(NULL)
      }
      targetResult <- unique(targetResult)
      
      geneTargetLoc <-
        convertGeneID(genetype = "Ensembl_trans",
                      genelist = targetResult,
                      org_assembly = org_assembly)
      
    }
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        miNearGene_temp <-
          getUCSC(
            region,
            pkg.env$upstream,
            pkg.env$downstream,
            org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          miNearGene_temp <-
            getNearToExon(
              region,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            ),
          miNearGene_temp <-
            getNearToIntron(
              region,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            )
        )
      )
      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      
      if (target) {
        geneL <- IRanges::findOverlapPairs(geneLoc_temp, geneTargetLoc)
        geneLo <- pintersect(geneL, ignore.strand = TRUE)
        miNearGene <-
          getUCSC(
            bedfile = geneLo,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
      else{
        miNearGene <- miNearGene_temp
      }
    }
    else{
      if (target) {
        miNearGene <- targetResult[, 2]
      }
      else{
        miNearGene <-
          getUCSC(
            bedfile = region,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
    }
    
    
    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = region,
          tad = TAD,
          cellline = pkg.env$cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = pkg.env$upstream,
          downstream = pkg.env$downstream
        )
      if (near | target)
        miNearGene <-
          as.data.frame(intersect(unlist(miNearGene), unlist(tadGene)))
      else
        miNearGene <- tadGene
    }
    
    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = miNearGene,
            cancer = cancer,
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          miNearGene <- intersect(unlist(miNearGene), nearG$feature)
        else
          miNearGene <- union(unlist(miNearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(
          exp1 = exp1,
          exp2 = exp2,
          label1 = label1 ,
          label2 = label2,
          corrMethod = pkg.env$corrMethod,
          varCutoff = pkg.env$varCutoff,
          corCutoff = pkg.env$minAbsCor,
          pcut = pkg.env$pcut,
          alternate = pkg.env$alternate,
          conf = pkg.env$conf
        )
        if (!isUnionCorGene)
          miNearGene <-
            intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      methods::new(
        "NoRCE",
        ID = '',
        Term = '',
        geneList = list(),
        pvalue = 0,
        pAdj = 0,
        GeneRatio = '',
        BckRatio = ''
      )
    }
    else{
      miEnrich <-
        goEnrichment(
          genes = miNearGene,
          GOtype = pkg.env$GOtype,
          org_assembly = org_assembly,
          pCut = pkg.env$pCut,
          pAdjCut = pkg.env$pAdjCut,
          pAdjust = pkg.env$pAdjust,
          backG = backG,
          backGType = backGType,
          min = pkg.env$min,
          enrichTest = pkg.env$enrichTest
        )
      if (length(miEnrich@Term) > 0)
      {
        miEnrich@ncGeneList <- commonGeneRegion(
          mrnaobject = miEnrich,
          org_assembly = org_assembly,
          downstream = pkg.env$downstream,
          upstream = pkg.env$upstream,
          inRegion =  region
        )
      }
      
      return(miEnrich)
    }
  }

#' Pathway enrichments of the microRNA regions with mRNAs that fall in the
#' given upstream/downstream regions of the microRNA genes
#'
#' @param region MiRNA region in a bed format
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'      considered in the analysis
#' @param target Boolean value shows whether miRNA target prediction should be
#'       performed
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'      performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and
#'      dm6 assembly, respectively.
#' @param gmtName Custom pathway gmt file
#' @param express Boolean variable whether co-expression analysis is performed.
#'       If this option is set to TRUE, co-expression analysis will be
#'       performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with
#'      custom data will be performed. When this option is set, exp1 and exp2
#'      parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for
#'      correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC,
#'      CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC,
#'      KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC,
#'      SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'       redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows
#'      must be patients. If gene names are provided as header, no need to
#'      redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not
#'      provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not
#'      provided, column name of the exp2 data will be taken.
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'       of the co-expression analysis and the other analysis should be
#'       considered
#' @param isGeneEnrich Boolean value whether gene enrichment should be
#'      performed
#'
#' @return miRNA pathway enrichment object for the given input
#' 
#'
#' @examples
#'
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#'
#' a<- mirnaRegionPathwayEnricher(region = regionNC,
#'              org_assembly = 'hg19')
#'
#' @export
mirnaRegionPathwayEnricher <-
  function(region,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           near = FALSE,
           target = FALSE,
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           gmtName = '',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           isUnionCorGene = FALSE,
           databaseFile,
           isGeneEnrich = FALSE) {
    if (missing(region)) {
      message("Region of interest is missing.")
    }
    
    if (missing(org_assembly)) {
      message("Assembly version is missing.")
    }
    
    assembly(org_assembly)
    
    
    if (target) {
      genes <-
        getUCSC(
          bedfile = region,
          downstream = 0,
          upstream = 0,
          org_assembly = org_assembly
        )
      targetResult <- predictmiTargets(gene = genes,
                                       type = "NCBI",
                                       org_assembly = org_assembly)
      if (is.null(targetResult))
      {
        message("There is no target!")
        return(NULL)
      }
      targetResult <- unique(targetResult)
      
      geneTargetLoc <-
        convertGeneID(genetype = "Ensembl_trans",
                      genelist = targetResult,
                      org_assembly = org_assembly)
      
    }
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        miNearGene_temp <-
          getUCSC(
            region,
            pkg.env$upstream,
            pkg.env$downstream,
            org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          miNearGene_temp <-
            getNearToExon(
              region,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            ),
          miNearGene_temp <-
            getNearToIntron(
              region,
              pkg.env$upstream,
              pkg.env$downstream,
              org_assembly
            )
        )
      )
      
      geneLoc_temp <-
        convertGeneID(genetype
                      = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      
      if (target) {
        geneL <- IRanges::findOverlapPairs(geneLoc_temp, geneTargetLoc)
        geneLo <- pintersect(geneL, ignore.strand = TRUE)
        miNearGene <-
          getUCSC(
            bedfile = geneLo,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
      else{
        miNearGene <- miNearGene_temp
      }
    }
    else{
      if (target) {
        miNearGene <- targetResult[, 2]
      }
      else{
        miNearGene <-
          getUCSC(
            bedfile = region,
            upstream = 0,
            downstream = 0,
            org_assembly = org_assembly
          )
      }
    }
    
    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = region,
          tad = TAD,
          cellline = pkg.env$cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = pkg.env$upstream,
          downstream = pkg.env$downstream
        )
      if (near | target)
        miNearGene <-
          as.data.frame(intersect(unlist(miNearGene), unlist(tadGene)))
      else
        miNearGene <- tadGene
    }
    
    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = miNearGene,
            cancer = cancer,
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          miNearGene <- intersect(unlist(miNearGene), nearG$feature)
        else
          miNearGene <- union(unlist(miNearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(
          exp1 = exp1,
          exp2 = exp2,
          label1 = label1 ,
          label2 = label2,
          corrMethod = pkg.env$corrMethod,
          varCutoff = pkg.env$varCutoff,
          corCutoff = pkg.env$minAbsCor,
          pcut = pkg.env$pcut,
          alternate = pkg.env$alternate,
          conf = pkg.env$conf
        )
        if (!isUnionCorGene)
          miNearGene <-
            intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      methods::new(
        "NoRCE",
        ID = '',
        Term = '',
        geneList = list(),
        pvalue = 0,
        pAdj = 0,
        GeneRatio = '',
        BckRatio = ''
      )
    }
    else{
      ifelse(
        pkg.env$pathwayType == 'kegg',
        pth <- 1,
        ifelse(
          pkg.env$pathwayType == 'reactome',
          pth <- 2,
          ifelse(pkg.env$pathwayType == 'wiki',
                 pth <- 3,
                 pth <- 4)
        )
      )
      
      funclist <- list(KeggEnrichment,
                       reactomeEnrichment,
                       WikiEnrichment,
                       pathwayEnrichment)
      
      miEnrich <- funclist[[pth]](
        genes = miNearGene,
        gmtFile = gmtName,
        org_assembly = org_assembly,
        pCut = pkg.env$pCut,
        pAdjCut = pkg.env$pAdjCut,
        pAdjust = pkg.env$pAdjust,
        isSymbol = pkg.env$isSymbol,
        min = pkg.env$min,
        isGeneEnrich = isGeneEnrich
      )
      
      return(miEnrich)
    }
  }

#' Predict the miRNA targets for the miRNA or mRNA genes, which is specified
#' with type parameter
#'
#' @param gene Data frame of miRNA or mRNA gene. Formats should be NCBI gene
#'     name, ENSEMBL gene or transcript id, and mirna
#' @param type Format of the gene, it should be "NCBI" for NCBI gene name,
#'     "Ensembl_gene" for ENSEMBL gene id, "Ensembl_trans" for Ensembl
#'     transcript id and "mirna" for miRNA gene
#' @param org_assembly Analyzed genome assembly. Possible assemblies are
#'     "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for
#'     fruit fly, "ce11" for worm, "hg19" and "hg38" for human
#'
#' @return miRNA:mRNA target sets of the given genes
#'
#' @examples
#'
#' a<- predictmiTargets(gene = brain_mirna[1:100,],
#'                      org_assembly = 'hg19',
#'                      type = "mirna")
#'
#'
#' @export
predictmiTargets <- function(gene, type, org_assembly)
{
  if (missing(gene)) {
    message("Genes are missing.")
  }
  if (missing(type)) {
    message("Format of the gene is missing.")
  }
  if (missing(org_assembly)) {
    message("Genome assembly version is missing.")
  }
  
  if (!exists("targets")) {
    # if(org_assembly == 'mm10') {
    #   targets <- NoRCE::targets_mouse
    # }
    # else if (org_assembly == 'dre10') {
    #   targets <- NoRCE::targets_zebra
    # }
    # else if (org_assembly == 'ce11') {
    #   targets <- NoRCE::targets_worm
    # }
    # else if (org_assembly == 'rn6') {
    #   targets <- NoRCE::targets_rat
    # }
    # else if (org_assembly == 'dm6') {
    #   targets <- NoRCE::targets_fly
    # }
    # else{
    #   targets <- NoRCE::targets_human
    # }
    markk <- 0
    ifelse(
      org_assembly == 'mm10',
      targets <- read.table(paste0("https://raw.githubusercontent.com/",
                                   "guldenolgun/NoRCE-data/master/target/target_mouse.txt")),
      ifelse(
        org_assembly == 'dre10',
        markk <- 2,
        ifelse (
          org_assembly == 'ce11',
          targets <- read.table(paste0("https://raw.githubusercontent.com/",
                                       "guldenolgun/NoRCE-data/master/target/target_worm.txt")),
          ifelse (
            org_assembly == 'rn6',
            markk <- 1,
            ifelse (
              org_assembly == 'dm6',
              targets <- read.table(paste0(
                "https://raw.githubusercontent.com/",
                "guldenolgun/NoRCE-data/master/target/target_fly.txt"),
                skip = 1),
              targets <- read.table(paste0("https://raw.githubusercontent.com/",
                                           "guldenolgun/NoRCE-data/master/target/target_human.txt"))
            )
          )
        )
      )
    )
    colnames(targets) <- c("ens","sym","trans","mir")
  }
  
  if(markk == 1){
    tmp1 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_rat.txt"))
    tmp2 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_rat1.txt"))
    tmp3 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_rat2.txt"))
    tmp4 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_rat3.txt"))
    target <- rbind(tmp1,tmp2,tmp3)
    targets <- merge(target,tmp4, by = 'V1')
    colnames(targets) <- c("ens","mir","sym","trans")
  }
  
  if(markk  == 2){
    tmp1 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_zebra.txt"))
    tmp2 <- read.table(paste0("https://raw.githubusercontent.com/",
                              "guldenolgun/NoRCE-data/master/target/target_zebra1.txt"))
    targets <- cbind(rbind(tmp1,tmp2),"")
    colnames(target) <- c("ens","sym","mir","trans")
  }
  
  gene <- as.data.frame(gene)
  colnames(gene) <- c("genes")
  
  ifelse(type == "NCBI",
         where <-
           targets[which(tolower(targets$sym) %in% tolower(gene$genes)),],
         ifelse (
           type == "mirna",
           where <-
             targets[which(tolower(targets$mir) %in% tolower(gene$genes)),],
           ifelse (type == "Ensembl_gene" ,
                   where <-
                     targets[which(tolower(targets$ens) %in% tolower(gene$genes)),],
                   where <-
                     targets[which(tolower(targets$trans) %in% tolower(gene$genes)),])
         ))
  if (nrow(where) == 0) {
    return(NULL)
  }
  else{
    colnames(where) <- c('genesEns', 'genesHugo', 'geneTrans', 'mirna')
    tmp1 <-
      data.frame(
        trans = unlist(apply(data.frame(where[, 3]), 2, 
                             strsplit, '[.]'))[2 *(seq_len(nrow(where))) - 1])
    tmp2 <-
      data.frame(
        gene =unlist(
          apply(data.frame(
            where[, 1]), 2, strsplit, '[.]'))[2 *(seq_len(nrow(where))) - 1])
    dat <-
      cbind.data.frame(tmp1, where$genesHugo, tmp2, where$mirna)
    return(dat)
  }
}
