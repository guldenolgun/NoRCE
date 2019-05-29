#' Given genes that fall in a given upstream and downstream region of mRNAs of interest, GO term enrichment analysis is carried out
#'
#' @param gene Input genes other than miRNA
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param upstream Holds upstream distance from the transcription start position
#' @param downstream Holds downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are "all", "exon", "intron"
#' @param GOtype Hierarchical category of the GO ontology. Possible values are "BP", "CC", "MF"
#' @param pCut Threshold value for the pvalue. Default value for pCut is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param near Boolean value presents whether cis-neighbourhood should be considered in the analysis
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#' @param slim Boolean value stating whether set of annotation should be performed for high level GO terms (GO slim)
#' @param genetype Type of the input gene list. Possible values are "Entrez", "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol "NCBI" value, for Entrez gene id "Entrez" is used.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or any new TAD regions can be used for the analysis. TAD regions must be formated as GRanges object. Predefined TAD regions are 'tad_hg19', 'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly, respectively.
#' @param isTADSearch Boolean value that shows whether TAD analysis is performed. This value has to be TRUE for TAD analysis.
#' @param cellline Cell lines for TAD regions
#' @param backG The set of genes that tested against to the input(background gene)
#' @param backGType Type of the background gene. If miRNA gene set is used for background gene, backGType should be set to the 'mirna'
#' @param enrichTest Types of enrichment methods to perform enrichment analysis. Possible values are "hyper"(default), "binom", "fisher", "chi".
#' @param express Boolean variable whether co-expression analysis is performed. If this option is set to TRUE, co-expression analysis will be performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with custom data will be performed. When this option is set, exp1 and exp2 parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not provided, column name of the exp2 data will be taken.
#' @param corrMethod Correlation coeffient method that will be used for evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater" or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is only used for the Pearson correlation coefficient if there are at least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output of the co-expression analysis and the other analysis should be considered
#'
#' @return GO term enrichment object for the given input
#'
#'
#' @examples
#' subsetGene <- brain_disorder_ncRNA[1:100,]
#'
#'\dontrun{
#' ncGO<-geneGOEnricher(gene = subsetGene, hg='hg19', near=TRUE, genetype = 'Ensembl_gene',
#'                      searchRegion = "exon", slim = FALSE, pAdjust = "none")
#'
#' writeEnrichment(mrnaObject = ncGO, file='enrichment.txt')
#'
#' getGoDag(mrnaObject = ncGO,filename = 'Dag.png', imageFormat = 'png', type = 'pvalue', n = 3)
#'
#' ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, hg='hg19', near=TRUE,
#'                      genetype = 'Ensembl_gene',
#'                      searchRegion = "exon", slim = FALSE, isTADSearch = TRUE,
#'                      TAD = tad_custom, cancer = 'GBMLGG', express = TRUE,
#'                      minAbsCor = 0.3,databaseFile = database, isUnionCorGene = TRUE)
#'
#' }
#'
#' @export
geneGOEnricher <-
  function(gene,
           hg,
           genetype,
           upstream = 10000,
           downstream = 10000,
           searchRegion = 'all',
           GOtype = "BP",
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           backG = '',
           backGType = 'pc_gene',
           enrichTest = "hyper",
           near = FALSE,
           min = 5,
           slim = FALSE,
           isTADSearch = FALSE,
           TAD = tad_hg19,
           cellline = 'all',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = "pearson",
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = 'greater',
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile) {
    if (missing(gene))
      message("Gene is missing?")

    if (missing(hg))
      message("Assembly version is missing?")

    if (missing(genetype))
      message("Input gene type is missing.")

    assembly(hg)
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    geneLoc <-
      convertGeneID(genetype = genetype,
                    genelist = gene$genes,
                    hg = hg)

    if (near) {
      if (searchRegion == 'all') {
        nearGene <-
          getUCSC(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      }
      else if (searchRegion == 'exon') {
        nearGene <-
          getNearToExon(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      }
      else if (searchRegion == 'intron') {
        nearGene <-
          getNearToIntron(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      }
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = geneLoc,
          upstream = 0,
          downstream = 0,
          hg = hg
        )
    }

    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          TAD = TAD,
          cellline = cellline,
          hg = hg,
          near = near,
          upstream = upstream,
          downstream = downstream
        )
      if(near)
        nearGene <-
          as.data.frame(intersect(unlist(nearGene), unlist(tadGene)))
      else
        nearGene <- tadGene
    }

    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = nearGene,
            cancer = cancer,
            minAbsCor = minAbsCor,
            databaseFile = databaseFile
          )
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(exp1 = exp1,
                           exp2= exp2,
                           label1 = label1 ,
                           label2 = label2,
                           corrMethod = corrMethod,
                           varCutoff = varCutoff,
                           corCutoff = minAbsCor,
                           pcut = pcut,
                           alternate = alternate,
                           conf = conf)
        tt<- sapply(1:dim(gene$genes)[1], function(x) unlist(which(nearG$firstExp %in% (gene$genes[x]))))
        nearG<-nearG[unlist(tt),]
        if(!isUnionCorGene)
        nearGene <- intersect(unlist(nearGene), unlist(nearG$SecondExp))
        else
          nearGene <- union(unlist(nearGene), unlist(nearG$SecondExp))
      }
    }
    if (length(nearGene)==0){
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
    enrichedGene <-
      goEnrichment(
        genes = nearGene,
        GOtype = GOtype,
        hg = hg,
        pCut = pCut,
        pAdjCut = pAdjCut,
        pAdjust = pAdjust,
        backG = backG,
        backGType = backGType,
        enrichTest = enrichTest,
        slim = slim,
        min = min
      )
    if (length(enrichedGene@Term) > 0)
    {
      enrichedGene@ncGeneList <- commonGene(mrnaobject = enrichedGene,hg = hg,downstream = downstream,upstream = upstream,inputGene = gene$genes,inGeneType = genetype)
    }

    objs <- ls(pos = ".GlobalEnv")
    gloVar <- c("mart", "go", "genomee", "ucsc")
    rm(list = objs[which(objs %in% gloVar)], pos = ".GlobalEnv")
    rm(objs,  pos = ".GlobalEnv")
    rm(gloVar,  pos = ".GlobalEnv")
    return(enrichedGene)
    }
  }

#' Given genes that fall in the given upstream and downstream region of mRNAs of interest, pathway enrichment analysis is carried out
#'
#' @param gene Input noncoding genes other than miRNA
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are "all", "exon", "intron"
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#' @param pathwayType Pathway database for enrichment. Possible values are 'reactome' for Reactome, 'kegg' for KEGG, 'wiki' for WikiPathways, 'other' for custom database
#' @param genetype Type of the input gene list. Possible values are "Entrez", "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol "NCBI" value, for Entrez gene id "Entrez", for mirbase id "mirna" is used.
#' @param near Boolean value presents whether cis-neighbourhood should be considered in the analysis
#' @param isTADSearch Boolean value that shows whether TAD analysis is performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or any new TAD regions can be used for the analysis. TAD regions must be formated as GRanges object. Predefined TAD regions are 'tad_hg19', 'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly, respectively.
#' @param cellline Cell lines for TAD regions.
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file. If it is set as TRUE, gene format of the gmt file should be symbol. Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
#' @param express Boolean variable whether co-expression analysis is performed. If this option is set to TRUE, co-expression analysis will be performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with custom data will be performed. When this option is set, exp1 and exp2 parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not provided, column name of the exp2 data will be taken.
#' @param corrMethod Correlation coeffient method that will be used for evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater" or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is only used for the Pearson correlation coefficient if there are at least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output of the co-expression analysis and the other analysis should be considered
#' @param isGeneEnrich Boolean value whether gene enrichment should be performed
#'
#' @return Pathway enrichment object for the given input
#'
#' @examples
#' subsetGene <- brain_disorder_ncRNA[1:100,]
#' \dontrun{
#' #Pathway enrichment based on the gen sets that falls into the TAD regions
#' ncRNAPathway<-genePathwayEnricher(gene = subsetGene , hg='hg19', isTADSearch = TRUE,TAD = tad_hg19,
#'                                   genetype = 'Ensembl_gene', min=2)
#' }
#'
#'
#' @export
genePathwayEnricher <-
  function(gene,
           hg,
           genetype,
           upstream = 10000,
           downstream = 10000,
           searchRegion = 'all',
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           min = 5,
           pathwayType = 'kegg',
           near = TRUE,
           isTADSearch = FALSE,
           TAD = tad_hg19,
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
           corrMethod = "pearson",
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = 'greater',
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile, isGeneEnrich = FALSE) {
    if (missing(gene)) {
      message("Gene is missing?")
    }

    if (missing(hg)) {
      message("Assembly version is missing?")
    }
    assembly(hg)

    if (missing(genetype))
      message("Input gene type is missing.")
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    geneLoc <-
      convertGeneID(genetype = genetype,
                    genelist = gene$genes,
                    hg = hg)

    if (near) {
      if (searchRegion == 'all')
        nearGene <-
          getUCSC(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'exon')
        nearGene <-
          getNearToExon(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'intron')
        nearGene <-
          getNearToIntron(
            bedfile = geneLoc,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = geneLoc,
          upstream = 0,
          downstream = 0,
          hg = hg
        )
    }

    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          TAD = TAD,
          cellline = cellline,
          hg = hg,
          near = near,
          upstream = upstream,
          downstream = downstream
        )
      if(near)
        nearGene <-
          as.data.frame(intersect(unlist(nearGene), unlist(tadGene)))
      else
        nearGene <- tadGene
    }

    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = nearGene,
            cancer = cancer,
            minAbsCor = minAbsCor,
            databaseFile = databaseFile
          )
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(exp1 = exp1,
                                exp2= exp2,
                                label1 = label1 ,
                                label2 = label2,
                                corrMethod = corrMethod,
                                varCutoff = varCutoff,
                                corCutoff = minAbsCor,
                                pcut = pcut,
                                alternate = alternate,
                                conf = conf)
        tt<- sapply(1:dim(gene$genes)[1], function(x) unlist(which(nearG$firstExp %in% (gene$genes[x]))))
        nearG<-nearG[unlist(tt),]
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene)==0){
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
      enrichedGene <-
        KeggEnrichment(
          genes = nearGene,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
    }
    else if (pathwayType == 'reactome') {
      enrichedGene <-
        reactomeEnrichment(
          genes = nearGene,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
    }
    else if (pathwayType == 'wiki') {
      enrichedGene <- WikiEnrichment(
        hg = hg,
        genes = nearGene,
        pCut = pCut,
        pAdjCut = pAdjCut,
        pAdjust = pAdjust,
        min = min
      )
    }
    else{
      enrichedGene <-
        pathwayEnrichment(
          genes = nearGene,
          gmtFile = gmtName,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          isSymbol = isSymbol,
          min = min, isGeneEnrich = isGeneEnrich
        )
    }
    if (length(enrichedGene@Term) > 0)
    {
      enrichedGene@ncGeneList <- commonGene(mrnaobject = enrichedGene,hg = hg,downstream = downstream,upstream = upstream,inputGene = gene$genes,inGeneType = genetype)
    }

      objs <- ls(pos = ".GlobalEnv")
      gloVar <- c("mart", "go", "genomee", "ucsc")
      rm(list = objs[which(objs %in% gloVar)], pos = ".GlobalEnv")
      rm(objs,  pos = ".GlobalEnv")
      rm(gloVar,  pos = ".GlobalEnv")
    return(enrichedGene)
    }
  }

#' Given gene regions that fall in the given upstream and downstream region of mRNAs of interest, GO term enrichment analysis is carried out
#'
#' @param region Bed format of the input gene regions other than miRNA
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are "all", "exon", "intron"
#' @param GOtype Hierarchical category of the GO ontology. Possible values are "BP", "CC", "MF"
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param near Boolean value presents whether cis-neighbourhood should be considered in the analysis
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#' @param slim Boolean value stating whether set of annotation should be performed for high level GO terms (GO slim)
#' @param isTADSearch Boolean value that shows whether TAD analysis is performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or any new TAD regions can be used for the analysis. TAD regions must be formated as GRanges object. Predefined TAD regions are 'tad_hg19', 'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly, respectively.
#' @param cellline Cell lines for TAD regions.
#' @param backG The set of genes that tested against to the input(background gene)
#' @param backGType Type of the background gene. If miRNA gene set is used for background gene, backGType should be set to the 'mirna'
#' @param enrichTest Types of enrichment methods to perform enrichment analysis. Possible values are "hyper"(default), "binom", "fisher", "chi".
#' @param express Boolean variable whether co-expression analysis is performed. If this option is set to TRUE, co-expression analysis will be performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with custom data will be performed. When this option is set, exp1 and exp2 parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not provided, column name of the exp2 data will be taken.
#' @param corrMethod Correlation coeffient method that will be used for evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater" or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is only used for the Pearson correlation coefficient if there are at least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output of the co-expression analysis and the other analysis should be considered

#' @return GO term enrichment object for the given input
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#'
#' \dontrun{
#'  regionGO<-geneRegionGOEnricher(region = regionNC, hg= 'hg19', near = TRUE)
#'
#'  getGoDag(mrnaObject = regionGO, type = 'pvalue', n =3, filename = 'gonetwork.png',
#'          imageFormat = 'png')
#'
#'  createNetwork(mrnaObject = regionGO,type = 'padjust',n = 3)
#' }
#'
#' @export
geneRegionGOEnricher <-
  function(region,
           hg,
           upstream = 10000,
           downstream = 10000,
           searchRegion = 'all',
           GOtype = "BP",
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           near = TRUE,
           backG = '',
           backGType = 'pc_gene',
           enrichTest = "hyper",
           min = 5,
           slim = FALSE,
           isTADSearch = FALSE,
           TAD = tad_hg19,
           cellline = 'all',
           express = FALSE,
           isCustomExp = FALSE,
           cancer,
           exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = "pearson",
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = 'greater',
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile) {
    if (missing(region)) {
      message("Bed file is missing. Please read any bed file is missing?")
    }

    if (missing(hg)) {
      message("Assembly version is missing?")
    }
    assembly(hg)
    if (near) {
      if (searchRegion == 'all')
        nearGene <-
          getUCSC(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'exon')
        nearGene <-
          getNearToExon(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'intron')
        nearGene <-
          getNearToIntron(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = region,
          upstream = 0,
          downstream = 0,
          hg = hg
        )
    }

    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          TAD = TAD,
          cellline = cellline,
          hg = hg,
          near = near,
          upstream = upstream,
          downstream = downstream
        )
      if(near)
        nearGene <-
          as.data.frame(intersect(unlist(nearGene), unlist(tadGene)))
      else
        nearGene <- tadGene
    }
    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = nearGene,
            cancer = cancer,
            minAbsCor = minAbsCor,
            databaseFile = databaseFile
          )
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(exp1 = exp1,
                                exp2= exp2,
                                label1 = label1 ,
                                label2 = label2,
                                corrMethod = corrMethod,
                                varCutoff = varCutoff,
                                corCutoff = minAbsCor,
                                pcut = pcut,
                                alternate = alternate,
                                conf = conf)
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene)==0){
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
    enrichedGene <-
      goEnrichment(
        genes = nearGene,
        GOtype = GOtype,
        hg = hg,
        pCut = pCut,
        pAdjCut = pAdjCut,
        backG = backG,
        backGType = backGType,
        enrichTest = enrichTest,
        pAdjust = pAdjust,
        slim = slim,
        min = min
      )

    objs <- ls(pos = ".GlobalEnv")
    gloVar <- c("mart", "go", "genomee", "ucsc")
    rm(list = objs[which(objs %in% gloVar)], pos = ".GlobalEnv")
    rm(objs,  pos = ".GlobalEnv")
    rm(gloVar,  pos = ".GlobalEnv")
    return(enrichedGene)
    }
  }

#' Given gene regions that fall in the given upstream and downstream region of mRNAs of interest, pathway enrichment analysis is carried out
#'
#' @param region Bed format of input gene regions other than miRNA. Input must be Granges object.
#' @param hg Genome assembly of interest for the analysis. Possible assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and "hg38" for human
#' @param upstream Upstream distance from the transcription start position
#' @param downstream Downstream distance from the transcription end position
#' @param searchRegion Search space of the cis-region. Possible values are "all", "exon", "intron"
#' @param pCut Threshold value for the pvalue. Default value is 0.05
#' @param pAdjCut Cutoff value for the adjusted p-values using one of given method. Default value is 0.05.
#' @param pAdjust Methods of the adjusted p-values. Possible methods are "bonferroni", "holm", "BH"(default)
#' @param min Minimum number of genes that are required for enrichment. By default, it is set to 5.
#' @param pathwayType Pathway database for enrichment. Possible values are 'reactome' for Reactome, 'kegg' for KEGG, 'wiki' for WikiPathways, 'other' for custom database
#' @param near Boolean value presents whether cis-neighbourhood should be considered in the analysis
#' @param isTADSearch Boolean value that shows whether TAD analysis is performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or any new TAD regions can be used for the analysis. TAD regions must be formated as GRanges object. Predefined TAD regions are 'tad_hg19', 'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly, respectively.
#' @param cellline Cell lines for TAD regions.
#' @param gmtName Custom pathway gmt file
#' @param isSymbol Boolean variable that hold the gene format of the gmt file. If it is set as TRUE, gene format of the gmt file should be symbol. Otherwise, gene format should be ENTREZ ID. By default, it is FALSE.
#' @param express Boolean variable whether co-expression analysis is performed. If this option is set to TRUE, co-expression analysis will be performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with custom data will be performed. When this option is set, exp1 and exp2 parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param exp1 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param exp2 Custom expression data matrix. Columns must be genes and rows must be patients. If gene names are provided as header, no need to redefine the headers(labels) of the expression data.
#' @param label1 Gene names of the custom exp1 expression data. If it is not provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is not provided, column name of the exp2 data will be taken.
#' @param corrMethod Correlation coeffient method that will be used for evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cutt off that genes have less variance than this value will be trimmed
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided", "greater" or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is only used for the Pearson correlation coefficient if there are at least 4 complete pairs of observations.
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#' @param isUnionCorGene Boolean value that shows whether union of the output of the co-expression analysis and the other analysis should be considered
#' @param isGeneEnrich Boolean value whether gene enrichment should be performed
#'
#' @return Pathway enrichment object of the given input
#'
#' @examples
#'
#' regionNC <-  readbed(dm_file = ncRegion,isText = FALSE)
#' \dontrun{
#' ncPath<-geneRegionPathwayEnricher(region = regionNC, hg = 'hg38', near = TRUE, upstream =1000,
#'                                  downstream = 0, pathwayType = 'kegg')
#'}
#' @export
geneRegionPathwayEnricher <-
  function(region,
           hg,
           upstream = 10000,
           downstream = 10000,
           searchRegion = 'all',
           pCut = 0.05,
           pAdjCut = 0.05,
           pAdjust = "BH",
           min = 5,
           pathwayType = 'kegg',
           near = FALSE,
           isTADSearch = FALSE,
           TAD = tad_hg19,
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
           corrMethod = "pearson",
           varCutoff = 0.0025,
           minAbsCor = 0.3,
           pcut = 0.05,
           alternate = 'greater',
           isUnionCorGene = FALSE,
           conf = 0.95,
           databaseFile,isGeneEnrich = FALSE) {
    if (missing(region)) {
      message("Bed file is missing. Please read any bed file is missing?")
    }

    if (missing(hg)) {
      message("Assembly version is missing?")
    }
    assembly(hg)

    if (near) {
      if (searchRegion == 'all')
        nearGene <-
          getUCSC(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'exon')
        nearGene <-
          getNearToExon(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
      if (searchRegion == 'intron')
        nearGene <-
          getNearToIntron(
            bedfile = region,
            upstream = upstream,
            downstream = downstream,
            hg = hg
          )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = region,
          upstream = 0,
          downstream = 0,
          hg = hg
        )
    }

    if (isTADSearch) {
      tadGene <-
        getTADOverlap(
          bedfile = geneLoc,
          TAD = TAD,
          cellline = cellline,
          hg = hg,
          near = near,
          upstream = upstream,
          downstream = downstream
        )
      if(near)
        nearGene <-
          as.data.frame(intersect(unlist(nearGene), unlist(tadGene)))
      else
        nearGene <- tadGene
    }

    if (express) {
      if (!isCustomExp) {
        nearG <-
          corrbasedMrna(
            mRNAgene = nearGene,
            cancer = cancer,
            minAbsCor = minAbsCor,
            databaseFile = databaseFile
          )
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
      }
      else{
        nearG <-  calculateCorr(exp1 = exp1,
                                exp2= exp2,
                                label1 = label1 ,
                                label2 = label2,
                                corrMethod = corrMethod,
                                varCutoff = varCutoff,
                                corCutoff = minAbsCor,
                                pcut = pcut,
                                alternate = alternate,
                                conf = conf)
        if(!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene)==0){
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
      enrichedGene <-
        KeggEnrichment(
          genes = nearGene,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
    }
    else if (pathwayType == 'reactome') {
      enrichedGene <-
        reactomeEnrichment(
          genes = nearGene,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          min = min
        )
    }
    else if (pathwayType == 'wiki') {
      enrichedGene <- WikiEnrichment(
        hg = hg,
        genes = nearGene,
        pCut = pCut,
        pAdjCut = pAdjCut,
        pAdjust = pAdjust,
        min = min
      )
    }
    else{
      enrichedGene <-
        pathwayEnrichment(
          genes = nearGene,
          gmtFile = gmtName,
          hg = hg,
          pCut = pCut,
          pAdjCut = pAdjCut,
          pAdjust = pAdjust,
          isSymbol = isSymbol,
          min = min,isGeneEnrich = isGeneEnrich
        )
    }

      objs <- ls(pos = ".GlobalEnv")
      gloVar <- c("mart", "go", "genomee", "ucsc")
      rm(list = objs[which(objs %in% gloVar)], pos = ".GlobalEnv")
      rm(objs,  pos = ".GlobalEnv")
      rm(gloVar,  pos = ".GlobalEnv")
    return(enrichedGene)
    }
  }
