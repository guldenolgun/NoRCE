op <- options(warn = (-1))
options(readr.num_columns = 0)

#' GO term enrichments of the microRNA genes with mRNAs that fall in the given
#' upstream/downstream regions of the microRNA genes
#'
#' @param gene Input microRNA gene. It supports both pre-miRNA and mature 
#'     miRNA, however, when target prediction is performed (target= TRUE),
#'     miRNA genes should be mature.
#' @param org_assembly Genome assembly of interest for the analysis. Possible 
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, 
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38"
#'     for human
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are 
#'     "all", "exon", "intron"
#' @param GOtype Hierarchical category of the GO ontology. Possible values 
#'     are "BP", "CC", "MF"
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given 
#'     method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are 
#'     "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By 
#'     default, this value is set to 5.
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
#' @param cellline Cell lines for TAD regions.
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
#' @param corrMethod Correlation coeffient method that will be used for 
#'     evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this 
#'     value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater"
#'     or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is 
#'     only used for the Pearson correlation coefficient if there are at 
#'     least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of
#'     the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output 
#'     of the co-expression analysis and the other analysis should be 
#'     considered
#'
#' @return MiRNA GO term enrichment object for the given input
#'
#' @examples
#' subsetGene <- brain_mirna[1:30,]
#' 
#' miGO <-mirnaGOEnricher(gene=subsetGene,
#'                        org_assembly='hg19',
#'                        near = TRUE,
#'                        target = FALSE, 
#'                        pAdjust = "none")
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
           upstream = 10000,
           downstream = 10000,
           searchRegion = c('all',"exon","intron"),
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
           near = FALSE,
           target = FALSE,
           min = 5,
           backGenes = '',
           backGType = 'pc_gene',
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           cellline = 'all',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = c("pearson","kendall", "spearman"),
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = c('greater',"two.sided", "less"),
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile = '') {
    if (missing(gene)) {
      message("Gene is missing.")
    }

    if (missing(org_assembly)) {
      message(
        "Assembly version is missing."
      )
    }
    if(!is.data.frame(gene) | is.character(gene))
      message("Type of the gene should be data.frame or character")
    
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
    }

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
      convertGeneID(genetype = "mirna",
                    genelist = a$genes,
                    org_assembly = org_assembly)


    if (near) {
      if (searchRegion == 'all')
        miNearGene_temp <-
          getUCSC(geneLoc, upstream, downstream, org_assembly)
      if (searchRegion == 'exon')
        miNearGene_temp <-
          getNearToExon(geneLoc, upstream, downstream, org_assembly)
      if (searchRegion == 'intron')
        miNearGene_temp <-
          getNearToIntron(geneLoc, upstream, downstream, org_assembly)

      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      if (target) {
        geneL <- findOverlapPairs(geneLoc_temp, geneTargetLoc)
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
          TAD = TAD,
          cellline = cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = upstream,
          downstream = downstream
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
          minAbsCor = minAbsCor,
          databaseFile = databaseFile
        )
        d <- nearG[which(a$genes %in% nearG$mirna_base),]

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
          corrMethod = corrMethod,
          varCutoff = varCutoff,
          corCutoff = minAbsCor,
          pcut = pcut,
          alternate = alternate,
          conf = conf
        )
        tt <-
          lapply(seq_len(nrow(a)), function(x) 
            unlist(which(
              nearG$firstExp %in% tolower(a$genes[x])
            )))
        nearG <- nearG[unlist(tt),]
        if (!isUnionCorGene)
          miNearGene <-
          intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      new(
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
          gene = miNearGene,
          GOtype = GOtype,
          org_assembly = org_assembly,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          backG = backGenes,
          backGType = backGType,
          min = min
        )
      if (length(miEnrich@Term)) {
        miEnrich@ncGeneList <-
          commonGene(
            mrnaobject = miEnrich,
            org_assembly = org_assembly,
            downstream = downstream,
            upstream = upstream,
            inputGene = rbind(a, gene),
            inGeneType = 'mirna'
          )
      }
      return(miEnrich)
    }
  }

#' Pathway enrichments of the microRNA genes with mRNAs that fall in the given
#' upstream/downstream regions of the microRNA genes
#'
#' @param gene Input microRNA gene. It supports both pre-miRNA and mature miRNA, 
#'     however, when target prediction is performed(target= TRUE), miRNA genes 
#'     should be mature.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'     assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'     "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and 
#'     "hg38" for human
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are 
#'     "all", "exon", "intron"
#' @param pCut Threshold value for the pvalue. Default value for pCut is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given 
#'     method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "holm", 
#'     "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By 
#'     default, it is set to 5
#' @param pathwayType Pathway database for enrichment. Possible values are 
#'     'reactome' for Reactome, 'kegg' for KEGG, 'wiki' for WikiPathways, 
#'     'other' for custom database
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
#' @param cellline Cell lines for TAD regions.
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file. 
#'      If it is set as TRUE, gene format of the gmt file should be symbol. 
#'      Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
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
#' @param corrMethod Correlation coeffient method that will be used for 
#'      evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this 
#'      value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater"
#'      or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is 
#'      only used for the Pearson correlation coefficient if there are at least
#'      4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of 
#'      the miRNA-mRNA
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
#' subsetGene <- brain_mirna[1:30,]
#'
#' miPath <- mirnaPathwayEnricher(gene = subsetGene,
#'                                org_assembly = 'hg19',
#'                                near = TRUE,
#'                                pAdjust = "none")
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
           upstream = 10000,
           downstream = 10000,
           searchRegion = c('all',"exon","intron"),
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
           pathwayType = c('kegg', 'reactome','wiki','other'),
           near = FALSE,
           target = FALSE,
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           cellline = 'all',
           gmtName = '',
           isSymbol = 'FALSE',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = c("pearson","kendall", "spearman"),
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = c('greater',"two.sided", "less"),
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile,isGeneEnrich = FALSE) {
    if (missing(gene)) {
      message("Gene is missing.")
    }
    if (missing(org_assembly)) {
      message(
        "Assembly version is missing."
      )
    }
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
    }
    
    if(!is.data.frame(gene) | is.character(gene))
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
      convertGeneID(genetype = "mirna",
                    genelist = a$genes,
                    org_assembly = org_assembly)


    if (near) {
      if (searchRegion == 'all')
        miNearGene_temp <-
          getUCSC(geneLoc, upstream, downstream, org_assembly)
      if (searchRegion == 'exon')
        miNearGene_temp <-
          getNearToExon(geneLoc, upstream, downstream, org_assembly)
      if (searchRegion == 'intron')
        miNearGene_temp <-
          getNearToIntron(geneLoc, upstream, downstream, org_assembly)

      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)
      if (target) {
        geneL <- findOverlapPairs(geneLoc_temp, geneTargetLoc)
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
          TAD = TAD,
          cellline = cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = upstream,
          downstream = downstream
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
          minAbsCor = minAbsCor,
          databaseFile = databaseFile
        )
        d <- nearG[which(a$genes %in% nearG$mirna_base),]

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
          corrMethod = corrMethod,
          varCutoff = varCutoff,
          corCutoff = minAbsCor,
          pcut = pcut,
          alternate = alternate,
          conf = conf
        )
        tt <-
          lapply(seq_len(nrow(a)), function(x) 
            unlist(which(
              nearG$firstExp %in% tolower(a$genes[x])
            )))
        nearG <- nearG[unlist(tt),]
        if (!isUnionCorGene)
          miNearGene <-
          intersect(unlist(miNearGene), nearG$SecondExp)
        else
          miNearGene <- union(unlist(miNearGene), nearG$SecondExp)
      }
    }
    if (length(miNearGene) == 0) {
      message("No common gene is found")
      new(
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
      if (pathwayType == 'kegg') {
        miEnrich <-
          KeggEnrichment(
            genes = miNearGene,
            org_assembly = org_assembly,
            pCut = pCut,
            pAdjCut = pAdjCut,
            pAdjust = pAdjust,
            min = min
          )
      }
      else if (pathwayType == 'reactome') {
        miEnrich <-
          reactomeEnrichment(
            genes = miNearGene,
            org_assembly = org_assembly,
            pCut = pCut,
            pAdjCut = pAdjCut,
            pAdjust = pAdjust,
            min = min
          )
      }
      else if (pathwayType == 'wiki') {
        miEnrich <- WikiEnrichment(
          org_assembly = org_assembly,
          genes = miNearGene,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
      }
      else{
        miEnrich <- pathwayEnrichment(
          genes = miNearGene,
          gmtFile = gmtName,
          org_assembly = org_assembly,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          isSymbol = isSymbol,
          min = min, isGeneEnrich = isGeneEnrich
        )
      }
      if (length(miEnrich@Term) > 0)
      {
        miEnrich@ncGeneList <-
          commonGene(
            mrnaobject = miEnrich,
            org_assembly = org_assembly,
            downstream = downstream,
            upstream = upstream,
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
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are 
#'      "all", "exon", "intron"
#' @param GOtype Hierarchical category of the GO ontology. Possible values are
#'      "BP", "CC", "MF"
#' @param pCut Threshold value for the pvalue. Default value for pCut is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given 
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are 
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By 
#'      default, it is set to 5.
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
#' @param cellline Cell lines for TAD regions.
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
#' @param corrMethod Correlation coeffient method that will be used for 
#'      evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this 
#'      value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater" 
#'      or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is only 
#'      used for the Pearson correlation coefficient if there are at least 4 
#'      complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of 
#'      the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'      of the co-expression analysis and the other analysis should be 
#'      considered
#'      
#' @return MiRNA GO enrichment object for the given input
#'
#'@examples
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
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
           upstream = 10000,
           downstream = 10000,
           searchRegion = c('all',"exon","intron"),
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
           near = FALSE,
           target = FALSE,
           min = 5,
           backG = '',
           backGType = 'pc-genes',
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           cellline = 'all',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = c("pearson","kendall", "spearman"),
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = c('greater',"two.sided", "less"),
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile) {
    if (missing(region)) {
      message("Region of interest is missing.")
    }

    if (missing(org_assembly)) {
      message(
        "Assembly version is missing."
      )
    }
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
    }

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
      if (searchRegion == 'all')
        miNearGene_temp <- getUCSC(region, upstream, downstream, org_assembly)
      if (searchRegion == 'exon')
        miNearGene_temp <-
          getNearToExon(region, upstream, downstream, org_assembly)
      if (searchRegion == 'intron')
        miNearGene_temp <-
          getNearToIntron(region, upstream, downstream, org_assembly)

      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)

      if (target) {
        geneL <- findOverlapPairs(geneLoc_temp, geneTargetLoc)
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
          TAD = TAD,
          cellline = cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = upstream,
          downstream = downstream
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
            minAbsCor = minAbsCor,
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
          corrMethod = corrMethod,
          varCutoff = varCutoff,
          corCutoff = minAbsCor,
          pcut = pcut,
          alternate = alternate,
          conf = conf
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
      new(
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
          gene = miNearGene,
          GOtype = GOtype,
          org_assembly = org_assembly,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          backG = backG,
          backGType = backGType,
          min = min
        )

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
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are 
#'      "all", "exon", "intron"
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given 
#'      method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are 
#'      "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#' @param min Minimum number of genes that are required for enrichment. By 
#'      default, it is set to 5.
#' @param pathwayType Pathway database for enrichment. Possible values are 
#'      'reactome' for Reactome, 'kegg' for KEGG, 'wiki' for WikiPathways, 
#'      'other' for custom database
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
#' @param cellline Cell lines for TAD regions
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file. 
#'      If it is set as TRUE, gene format of the gmt file should be symbol. 
#'      Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
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
#' @param corrMethod Correlation coeffient method that will be used for 
#'      evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this
#'      value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater"
#'       or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is 
#'      only used for the Pearson correlation coefficient if there are at 
#'      least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of
#'       the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output
#'       of the co-expression analysis and the other analysis should be 
#'       considered
#' @param isGeneEnrich Boolean value whether gene enrichment should be 
#'      performed
#'
#' @return miRNA pathway enrichment object for the given input
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' 
#' a<- mirnaRegionGOEnricher(region = regionNC, 
#'                           org_assembly = 'hg19', 
#'                           near = TRUE)
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
           upstream = 10000,
           downstream = 10000,
           searchRegion = c('all',"exon","intron"),
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
           pathwayType = c('kegg', 'reactome','wiki','other'),
           near = FALSE,
           target = FALSE,
           isTADSearch = FALSE,
           TAD = c(tad_hg19, tad_dmel, tad_hg38, tad_mm10),
           cellline = 'all',
           gmtName = '',
           isSymbol = FALSE,
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = c("pearson","kendall", "spearman"),
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = c('greater',"two.sided", "less"),
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile, isGeneEnrich = FALSE) {
    if (missing(region)) {
      message("Region of interest is missing.")
    }

    if (missing(org_assembly)) {
      message(
        "Assembly version is missing."
      )
    }
    
    if (class(pkg.env$mart)[1] != "Mart") {
      assembly(org_assembly)
    }

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
      if (searchRegion == 'all')
        miNearGene_temp <- getUCSC(region, upstream, downstream, org_assembly)
      if (searchRegion == 'exon')
        miNearGene_temp <-
          getNearToExon(region, upstream, downstream, org_assembly)
      if (searchRegion == 'intron')
        miNearGene_temp <-
          getNearToIntron(region, upstream, downstream, org_assembly)

      geneLoc_temp <-
        convertGeneID(genetype = "NCBI",
                      genelist = miNearGene_temp,
                      org_assembly = org_assembly)

      if (target) {
        geneL <- findOverlapPairs(geneLoc_temp, geneTargetLoc)
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
          TAD = TAD,
          cellline = cellline,
          org_assembly = org_assembly,
          near = near,
          upstream = upstream,
          downstream = downstream
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
            minAbsCor = minAbsCor,
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
          corrMethod = corrMethod,
          varCutoff = varCutoff,
          corCutoff = minAbsCor,
          pcut = pcut,
          alternate = alternate,
          conf = conf
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
      new(
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
      if (pathwayType == 'kegg') {
        miEnrich <-
          KeggEnrichment(
            genes = miNearGene,
            org_assembly = org_assembly,
            pCut = pCut,
            pAdjCut = pAdjCut,
            pAdjust = pAdjust,
            min = min
          )
      }
      else if (pathwayType == 'reactome') {
        miEnrich <-
          reactomeEnrichment(
            genes = miNearGene,
            org_assembly = org_assembly,
            pCut = pCut,
            pAdjCut = pAdjCut,
            pAdjust = pAdjust,
            min = min
          )
      }
      else if (pathwayType == 'wiki') {
        miEnrich <- WikiEnrichment(
          org_assembly = org_assembly,
          genes = miNearGene,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
      }
      else{
        miEnrich <- pathwayEnrichment(
          genes = miNearGene,
          gmtFile = gmtName,
          org_assembly = org_assembly,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          isSymbol = isSymbol,
          min = min, isGeneEnrich = isGeneEnrich
        )
      }

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
    message(
      "Genes are missing."
    )
  }
  if (missing(type)) {
    message("Format of the gene is missing.")
  }
  if (missing(org_assembly)) {
    message(
      "Genome assembly version is missing."
    )
  }

  if (!exists("targets")) {
    if (org_assembly == 'mm10') {
      targets <- NoRCE::targets_mouse
    }
    else if (org_assembly == 'dre10') {
      targets <- NoRCE::targets_zebra
    }
    else if (org_assembly == 'ce11') {
      targets <- NoRCE::targets_worm
    }
    else if (org_assembly == 'rn6') {
      targets <- NoRCE::targets_rat
    }
    else if (org_assembly == 'dm6') {
      targets <- NoRCE::targets_fly
    }
    else{
      targets <- NoRCE::targets_human
    }
  }

  gene <- as.data.frame(gene)
  colnames(gene) <- c("genes")

  if (type == "NCBI") {
    where <- targets[which(tolower(targets$X2) %in% gene$genes), ]
  }

  else if (type == "mirna") {
    where <-
      targets[which(tolower(targets$X4) %in% tolower(gene$genes)), ]
  }
  else if (type == "Ensembl_gene") {
    where <- targets[which(tolower(targets$X1) %in% gene$genes), ]
  }
  else if (type == "Ensembl_trans") {
    where <- targets[which(tolower(targets$X3) %in% gene$genes), ]
  }
  if (nrow(where) == 0) {
    return(NULL)
  }
  else{
    colnames(where) <- c('genesEns', 'genesHugo', 'geneTrans', 'mirna')
    tmp1 <-
      data.frame(trans = unlist(
        apply((where[, 3]), 2, strsplit, '[.]'))[2 *
                                                   (seq_len(nrow(where))) - 1]) 
    tmp2 <-
      data.frame(gene = 
                  unlist(apply((
                   where[, 1]), 2, strsplit, '[.]'))[2 * 
                                                   (seq_len(nrow(where))) - 1])
    dat <-
      cbind.data.frame(tmp1, where$genesHugo, tmp2, where$mirna)
    return(dat)
  }
}
