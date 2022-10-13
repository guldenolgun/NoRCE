#' Pearson correlation coefficient value of the miRNA genes between
#' miRNA:mRNA for a given correlation cut-off and cancer.
#'
#' @param mirnagene Data frame of the miRNA genes in mature format
#' @param cancer Name of the TCGA project code such as 'BRCA' that is analyzed
#'      for miRNA-mRNA correlation. Possible cancer types ACC, BLCA, BRCA,
#'      CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN,
#'      KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC,
#'      SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of
#'      the miRNA-mRNA
#' @param databaseFile Path of the miRcancer.db file
#'
#' @return Data frame of the miRNA-mRNA correlation result
#' 
#' @import dbplyr
#' @import stringr
#'
#' @export
corrbased <- function(mirnagene,
                      cancer,
                      minAbsCor,
                      databaseFile) {
  colnames(mirnagene) <- c('g')
  
  a <-
    as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "",
                       mirnagene$g))
  
  a <- data.frame(str_replace_all(a[,1], 'miR', 'mir'))
  colnames(a) <- 'g'
  a <- unique(rbind(a, mirnagene))
  
  conn <- DBI::dbConnect(RSQLite::SQLite(), databaseFile)
  
  dat <- conn %>%
    dplyr::tbl('cor_mir') %>%
    dplyr::select(mirna_base, feature, cancer) %>%
    dplyr::filter(mirna_base %in% local(a$g)) %>%
    dplyr::collect() %>%
    tidyr::gather(cancer, cor, -mirna_base, -feature) %>%
    dplyr::mutate(cor = cor / 100) %>% dplyr::filter(abs(cor) > minAbsCor) %>%
    dplyr::arrange(dplyr::desc(abs(cor))) %>% na.omit()

  return(dat)
}

#' Pearson correlation coefficient value of the mRNA genes between 
#' miRNA:mRNA for a given correlation cut-off and cancer.
#'
#' @param mRNAgene Data frame of the mRNA genes
#' @param cancer Name of the TCGA project code such as 'BRCA' that is analyzed
#'      for miRNA-mRNA correlation. Possible cancer types ACC, BLCA, BRCA,
#'      CESC, CHOL, COAD, COADREAD, DLBC, ESCA, GBMLGG, HNSC, KICH, KIPAN,
#'      KIRC, KIRP, LGG, LIHC, LUAD, LUSC, OV, PAAD, PCPG, PRAD, READ, SARC,
#'      SKCM, STAD, STES, TGCT, THCA, THYM, UCEC, UCS, UVM
#' @param minAbsCor Cut-off value for the Pearson correlation coefficient of
#'      the miRNA-mRNA
#' @param databaseFile Path of miRcancer.db file
#'
#' @return Data frame of the miRNA-mRNA correlation result
#' 
#' @importFrom dplyr %>%
#'
#' @export
corrbasedMrna <-
  function(mRNAgene, cancer, minAbsCor, databaseFile) {
    colnames(mRNAgene) <- c('g')
    
    conn <- DBI::dbConnect(RSQLite::SQLite(), databaseFile)
    
    dat <- conn %>%
      dplyr::tbl('cor_mir') %>%
      dplyr::select(mirna_base, feature, cancer) %>%
      dplyr::filter(feature %in% !!mRNAgene$g) %>%
      dplyr::collect() %>%
      tidyr::gather(cancer, cor,-mirna_base,-feature) %>%
      dplyr::mutate(cor = cor / 100) %>%
      dplyr::filter(abs(cor) > minAbsCor) %>%
      dplyr::arrange(dplyr::desc(abs(cor))) %>% na.omit()
    
    return(dat)
  }

#' Get TCGA miRNAseq expression of miRNA genes for the given cancer
#'
#' @param mirnagene Data frame of the mature format
#' @param cancer Name of the TCGA project code such as 'BRCA'
#' @param databaseFile Path of miRcancer.db file
#'
#' @return Data frame of the raw read count of the given miRNA genes
#'       for different patients
#'       
#' @import dbplyr
#' @import stringr
#'
#' @export
getmiRNACount <- function(mirnagene, cancer, databaseFile) {
  colnames(mirnagene) <- c('g')
  
  a <-
    as.data.frame(gsub(paste(c("-3p", "-5p"), collapse = "|"), "",
                       mirnagene$g))
  
  a <- data.frame(str_replace_all(a[,1], 'miR', 'mir'))
  colnames(a) <- 'g'
  a <- unique(rbind(a, mirnagene))
  
  conn <- DBI::dbConnect(RSQLite::SQLite(), databaseFile)
  
  dat <-
    conn %>%
    dplyr::tbl('profiles') %>%
    dplyr::select(study, mirna_base, count) %>%
    dplyr::filter(mirna_base %in% !!a$g) %>%
    dplyr::filter(study %in% cancer) %>%
    dplyr::collect() %>% na.omit()
  
  return(dat)
}

#' Calculates the correlation coefficient values between two custom
#' expression data.
#'
#' @param exp1 Custom expression data matrix or SummarizedExperiment data.
#'      Columns must be genes and rows must be patients.
#' @param exp2 Custom expression data matrix or SummarizedExperiment data.
#'      Columns must be genes and rows must be patients.
#' @param label1 Gene names of the custom exp1 expression data. If it is
#'      not provided, column name of the exp1 data will be taken.
#' @param label2 Gene names of the custom exp2 expression data. If it is
#'      not provided, column name of the exp2 data will be taken.
#' @param corrMethod Correlation coeffient method that will be used for
#'      evaluation. Possible values are "pearson", "kendall", "spearman"
#' @param varCutoff Variance cut off that genes have less variance than
#'      this value will be trimmed
#' @param corCutoff Correlation cut off values for the given correlation
#'      method
#' @param pcut P-value cut off for the correlation values
#' @param alternate Holds the alternative hypothesis and "two.sided",
#'      "greater" or "less" are the possible values.
#' @param conf Confidence level for the returned confidence interval. It is
#'      only used for the Pearson correlation coefficient if there are at least
#'      4 complete pairs of observations.
#'
#' @return Pairwise relations between gene-gene with corresponding correlation
#'       value and pvalue
#'
#' @examples
#'
#' #Assume that mirnanorce and mrnanorce are custom patient by gene data
#' a<-calculateCorr(exp1 = mirna, exp2 = mrna )
#'
#' @export
calculateCorr <-
  function(exp1,
           exp2,
           label1 = '',
           label2 = '',
           corrMethod = "pearson",
           varCutoff = 0.0025,
           corCutoff = 0.3,
           pcut = 0.05,
           alternate = 'greater',
           conf = 0.95) {
    if (class(exp1)[[1]] == "RangedSummarizedExperiment") {
      tmp1 <- SummarizedExperiment::assay(exp1)
      exp1 <- t(tmp1)
    }
    if (class(exp2)[[1]] == "RangedSummarizedExperiment") {
      tmp2 <- SummarizedExperiment::assay(exp2)
      exp2 <- t(tmp2)
    }
    
    ifelse(label1 == '',
           label1 <- data.frame(colnames(exp1)),
           colnames(exp1) <- t(label1))
    
    ifelse(label2 == '',
           label2 <- data.frame(colnames(exp2)),
           colnames(exp2) <- t(label2))
    
    
    var1 <- apply(exp1, 2, var)
    var2 <- apply(exp2, 2, var)
    
    exp1 <- exp1[, which(var1 > varCutoff)]
    exp2 <- exp2[, which(var2 > varCutoff)]
    label1 <- data.frame(label1[which(var1 > varCutoff),])
    label2 <- data.frame(label2[which(var2 > varCutoff),])
    extractData <- data.frame()
    for (i in seq_len(ncol(exp2))) {
      #for (i in 1:dim(exp2)[2]) {
      tmp <-
        apply(exp1, 2, function(x)
          cor.test(
            x,
            exp2[[i]],
            alternative = alternate,
            conf.level = conf,
            method = corrMethod
          ))
      out <- lapply(tmp, function(x)
        c(x$estimate, x$p.value))
      corpValue <- data.frame(do.call(rbind, out))
      index <- which(abs(corpValue$cor) > corCutoff &
                       corpValue$V2 < pcut)
      if (!IRanges::isEmpty(index)) {
        extractData <-
          rbind(extractData, data.frame(label1[index,],
                                        label2[i,], corpValue[index, ]))
      }
    }
    colnames(extractData) <-
      c('firstExp', 'SecondExp', 'Cor', 'Pvalue')
    return(extractData)
  }
