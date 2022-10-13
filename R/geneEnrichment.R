#' Given genes that fall in a given upstream and downstream region of 
#' mRNAs of interest, GO term enrichment analysis is carried out
#'
#' @param gene Input genes other than miRNA
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'      considered in the analysis
#' @param genetype Type of the input gene list. Possible values are "Entrez",
#'      "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol
#'      "NCBI" value, for Entrez gene id "Entrez" is used.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'      assembly, respectively.
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'      performed. This value has to be TRUE for TAD analysis.
#' @param backG The set of genes that tested against to the input
#'      (background gene)
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
#' @return GO term enrichment object for the given input
#'
#'
#' @examples
#' ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, org_assembly='hg19',
#' near=TRUE, genetype = 'Ensembl_gene')
#'
#' @export
geneGOEnricher <-
  function(gene,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           genetype = c("Entrez",
                        "mirna",
                        "Ensembl_gene",
                        "Ensembl_trans",
                        "NCBI"),
           backG = '',
           backGType = 'pc_gene',
           near = FALSE,
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
    if (missing(gene))
      message("Gene is missing?")
    
    if (!is.data.frame(gene) &
        !is.character(gene) & !is.factor(gene))
      message("Type of the gene should be data.frame or character")
    
    if (missing(org_assembly))
      message("Assembly version is missing?")
    
    if (missing(genetype))
      message("Input gene type is missing.")
    
    assembly(org_assembly)
    
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    geneLoc <-
      convertGeneID(
        genetype = genetype,
        genelist = gene$genes,
        org_assembly = org_assembly
      )
    
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        nearGene <-
          getUCSC(
            bedfile = geneLoc,
            upstream = pkg.env$upstream,
            downstream = pkg.env$downstream,
            org_assembly = org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          nearGene <-
            getNearToExon(
              bedfile = geneLoc,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            ),
          nearGene <-
            getNearToIntron(
              bedfile = geneLoc,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            )
        )
      )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = geneLoc,
          upstream = 0,
          downstream = 0,
          org_assembly = org_assembly
        )
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
      if (near)
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
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
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
        tt <- lapply(seq_len(length(gene$genes)), function(x)
          unlist(which(nearG$firstExp %in% (gene$genes[x]))))
        nearG <- nearG[unlist(tt), ]
        if (!isUnionCorGene)
          nearGene <-
          intersect(unlist(nearGene), unlist(nearG$SecondExp))
        else
          nearGene <-
          union(unlist(nearGene), unlist(nearG$SecondExp))
      }
    }
    if (length(nearGene) == 0) {
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
      enrichedGene <-
        goEnrichment(
          genes = nearGene,
          GOtype = pkg.env$GOtype,
          org_assembly = org_assembly,
          pCut = pkg.env$pCut,
          pAdjCut = pkg.env$pAdjCut,
          pAdjust = pkg.env$pAdjust,
          backG = backG,
          backGType = backGType,
          enrichTest = pkg.env$enrichTest,
          min = pkg.env$min
        )
      if (length(enrichedGene@Term) > 0)
      {
        enrichedGene@ncGeneList <- commonGene(
          mrnaobject = enrichedGene,
          org_assembly = org_assembly,
          downstream = pkg.env$downstream,
          upstream = pkg.env$upstream,
          inputGene = gene$genes,
          inGeneType = genetype
        )
      }
      
      return(enrichedGene)
    }
  }

#' Given genes that fall in the given upstream and downstream region of 
#' mRNAs of interest, pathway enrichment analysis is carried out
#'
#' @param gene Input noncoding genes other than miRNA
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param genetype Type of the input gene list. Possible values are "Entrez",
#'      "mirna", "Ensembl_gene", "Ensembl_trans", "NCBI". For HUGO gene symbol
#'      "NCBI" value, for Entrez gene id "Entrez", for mirbase id "mirna" is
#'      used.
#' @param near Boolean value presents whether cis-neighbourhood should be
#'      considered in the analysis
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'      performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'      assembly, respectively.
#' @param gmtName Custom pathway gmt file
#' @param express Boolean variable whether co-expression analysis is performed.
#'      If this option is set to TRUE, co-expression analysis will be
#'      performed.
#' @param isCustomExp Boolean variable whether co-expression analysis with
#'      custom data will be performed. When this option is set, exp1 and exp2
#'      parameters must be defined.
#' @param cancer Defines the name of the TCGA project code such as 'BRCA' for
#'      correlation analysis. Possible cancer types ACC, BLCA, BRCA, CESC,
#'      CHOL,COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN, KIRC, KIRP,
#'      LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, STES,
#'      TGCT, THCA, THYM, UCEC, UCS, UVM, LGG
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
#' @return Pathway enrichment object for the given input
#' 
#'
#' @examples
#' #Pathway enrichment based on the gen sets that falls into the TAD regions
#' ncRNAPathway<-genePathwayEnricher(gene = brain_disorder_ncRNA ,
#'                                   org_assembly='hg19',
#'                                   isTADSearch = TRUE,
#'                                   TAD = tad_hg19,
#'                                   genetype = 'Ensembl_gene')
#'
#'
#' @export
genePathwayEnricher <-
  function(gene,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           genetype = c("Entrez",
                        "mirna",
                        "Ensembl_gene",
                        "Ensembl_trans",
                        "NCBI"),
           near = TRUE,
           isTADSearch = FALSE,
           TAD = tad_hg19,
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
      message("Gene is missing?")
    }
    if (!is.data.frame(gene) &
        !is.character(gene) & !is.factor(gene))
      message("Type of the gene should be data.frame or character")
    
    if (missing(org_assembly)) {
      message("Assembly version is missing?")
    }
    assembly(org_assembly)
    
    
    if (missing(genetype))
      message("Input gene type is missing.")
    gene <- as.data.frame(gene)
    colnames(gene) <- c("genes")
    geneLoc <-
      convertGeneID(
        genetype = genetype,
        genelist = gene$genes,
        org_assembly = org_assembly
      )
    
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        nearGene <-
          getUCSC(
            bedfile = geneLoc,
            upstream = pkg.env$upstream,
            downstream = pkg.env$downstream,
            org_assembly = org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          nearGene <-
            getNearToExon(
              bedfile = geneLoc,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            ),
          nearGene <-
            getNearToIntron(
              bedfile = geneLoc,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            )
        )
      )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = geneLoc,
          upstream = 0,
          downstream = 0,
          org_assembly = org_assembly
        )
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
      if (near)
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
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
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
        tt <- lapply(seq_len(nrow(gene$genes)), function(x)
          unlist(which(nearG$firstExp %in% (gene$genes[x]))))
        nearG <- nearG[unlist(tt), ]
        if (!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene) == 0) {
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
      
      enrichedGene <- funclist[[pth]](
        genes = nearGene,
        gmtFile = gmtName,
        org_assembly = org_assembly,
        pCut = pkg.env$pCut,
        pAdjCut = pkg.env$pAdjCut,
        pAdjust = pkg.env$pAdjust,
        isSymbol = pkg.env$isSymbol,
        min = pkg.env$min,
        isGeneEnrich = isGeneEnrich
      )
      
      if (length(enrichedGene@Term) > 0)
      {
        enrichedGene@ncGeneList <- commonGene(
          mrnaobject = enrichedGene,
          org_assembly = org_assembly,
          downstream = pkg.env$downstream,
          upstream = pkg.env$upstream,
          inputGene = gene$genes,
          inGeneType = genetype
        )
      }
      
      return(enrichedGene)
    }
  }

#' Given gene regions that fall in the given upstream and downstream region
#' of mRNAs of interest, GO term enrichment analysis is carried out
#'
#' @param region Bed format of the input gene regions other than miRNA
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'      considered in the analysis
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'      performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'      any new TAD regions can be used for the analysis. TAD regions must be
#'      formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'      'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6
#'      assembly, respectively.
#' @param backG The set of genes that tested against to the input
#'      (background gene)
#' @param backGType Type of the background gene. If miRNA gene set is used for
#'      background gene, backGType should be set to the 'mirna'
#' @param express Boolean variable whether co-expression analysis is
#'      performed. If this option is set to TRUE, co-expression analysis will
#'      be performed.
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
#'
#' @return GO term enrichment object for the given input
#' 
#' @importFrom rtracklayer import
#' @examples
#'
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#' regionGO<-geneRegionGOEnricher(region = regionNC, org_assembly= 'hg19',
#'                                near = TRUE)
#'
#' @export
geneRegionGOEnricher <-
  function(region,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3"),
           near = TRUE,
           backG = '',
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
           databaseFile) {
    if (missing(region)) {
      message("Bed file is missing. Please read any bed file is missing?")
    }
    
    if (missing(org_assembly)) {
      message("Assembly version is missing?")
    }
    assembly(org_assembly)
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        nearGene <-
          getUCSC(
            bedfile = region,
            upstream = pkg.env$upstream,
            downstream = pkg.env$downstream,
            org_assembly = org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          nearGene <-
            getNearToExon(
              bedfile = region,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            ),
          nearGene <-
            getNearToIntron(
              bedfile = region,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            )
        )
      )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = region,
          upstream = 0,
          downstream = 0,
          org_assembly = org_assembly
        )
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
      if (near)
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
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
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
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene) == 0) {
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
      enrichedGene <-
        goEnrichment(
          genes = nearGene,
          GOtype = pkg.env$GOtype,
          org_assembly = org_assembly,
          pCut = pkg.env$pCut,
          pAdjCut = pkg.env$pAdjCut,
          backG = backG,
          backGType = backGType,
          enrichTest = pkg.env$enrichTest,
          pAdjust = pkg.env$pAdjust,
          min = pkg.env$min
        )
      if (length(enrichedGene@Term) > 0)
      {
        enrichedGene@ncGeneList <- commonGeneRegion(
          mrnaobject = enrichedGene,
          org_assembly = org_assembly,
          downstream = pkg.env$downstream,
          upstream = pkg.env$upstream,
          inRegion =  region
        )
      }
      
      return(enrichedGene)
    }
  }

#' Given gene regions that fall in the given upstream and downstream region
#' of mRNAs of interest, pathway enrichment analysis is carried out
#'
#' @param region Bed format of input gene regions other than miRNA. Input must
#'      be Granges object.
#' @param org_assembly Genome assembly of interest for the analysis. Possible
#'      assemblies are "mm10" for mouse, "dre10" for zebrafish, "rn6" for rat,
#'      "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast, "hg19" and
#'      "hg38" for human
#' @param near Boolean value presents whether cis-neighbourhood should be
#'     considered in the analysis
#' @param isTADSearch Boolean value that shows whether TAD analysis is
#'     performed. This value has to be TRUE for TAD analysis.
#' @param TAD TAD genomic regions for the species. Predefined TAD regions or
#'     any new TAD regions can be used for the analysis. TAD regions must be
#'     formated as GRanges object. Predefined TAD regions are 'tad_hg19',
#'     'tad_hg38', 'tad_mm10', 'tad_dmel' for hg19, hg38, mm9 and dm6 assembly,
#'     respectively.
#' @param gmtName Custom pathway gmt file
#' @param express Boolean variable whether co-expression analysis is
#'      performed. If this option is set to TRUE, co-expression analysis will
#'      be performed.
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
#' @return Pathway enrichment object of the given input
#'
#' @examples
#'
#' regions<-system.file("extdata", "ncRegion.txt", package = "NoRCE")
#' regionNC <- rtracklayer::import(regions, format = "BED")
#' ncPath<-geneRegionPathwayEnricher(region = regionNC,
#'                                   org_assembly = 'hg19',
#'                                   near = TRUE)
#' @export
geneRegionPathwayEnricher <-
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
           isTADSearch = FALSE,
           TAD = tad_hg19,
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
      message("Bed file is missing. Please read any bed file is missing?")
    }
    
    if (missing(org_assembly)) {
      message("Assembly version is missing?")
    }
    
    assembly(org_assembly)
    
    
    if (near) {
      ifelse(
        pkg.env$searchRegion == 'all',
        nearGene <-
          getUCSC(
            bedfile = region,
            upstream = pkg.env$upstream,
            downstream = pkg.env$downstream,
            org_assembly = org_assembly
          ),
        ifelse(
          pkg.env$searchRegion == 'exon',
          nearGene <-
            getNearToExon(
              bedfile = region,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            ),
          nearGene <-
            getNearToIntron(
              bedfile = region,
              upstream = pkg.env$upstream,
              downstream = pkg.env$downstream,
              org_assembly = org_assembly
            )
        )
      )
    }
    else{
      nearGene <-
        getUCSC(
          bedfile = region,
          upstream = 0,
          downstream = 0,
          org_assembly = org_assembly
        )
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
      if (near)
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
            minAbsCor = pkg.env$minAbsCor,
            databaseFile = databaseFile
          )
        if (!isUnionCorGene)
          nearGene <- intersect(unlist(nearGene), nearG$feature)
        else
          nearGene <- union(unlist(nearGene), nearG$feature)
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
          nearGene <- intersect(unlist(nearGene), nearG$SecondExp)
        else
          nearGene <- union(unlist(nearGene), nearG$SecondExp)
      }
    }
    if (length(nearGene) == 0) {
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
      
      enrichedGene <- funclist[[pth]](
        genes = nearGene,
        gmtFile = gmtName,
        org_assembly = org_assembly,
        pCut = pkg.env$pCut,
        pAdjCut = pkg.env$pAdjCut,
        pAdjust = pkg.env$pAdjust,
        isSymbol = pkg.env$isSymbol,
        min = pkg.env$min,
        isGeneEnrich = isGeneEnrich
      )
      
      return(enrichedGene)
    }
  }
