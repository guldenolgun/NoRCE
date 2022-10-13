#' Draw dot plot of the enrichment object
#'
#' @param mrnaObject Object of the enrichment result
#' @param type Draw the dot plot according to the p-value or adjusted p-value
#'      ("pvalue", "pAdjust")
#' @param n Number of GO terms or pathways, that ordered by type and has
#'      least number of top p-value
#'
#' @return Dot plot of the top n enrichment results
#' @importFrom AnnotationDbi Term
#' @importFrom ggplot2 aes element_text geom_point ggplot labs theme theme_bw
#'
#' @export
drawDotPlot <- function(mrnaObject, type = "pAdjust", n) {
  if (missing(mrnaObject)) {
    message("Expected input is NoRCE object.")
  }
  
  if (missing(n)) {
    n <- length(mrnaObject@ID)
  }
  if (n > length(mrnaObject@ID))
    n <- length(mrnaObject@ID)
  tmp <- topEnrichment(mrnaObject, type, n)
  if (type == "pvalue") {
    p <-
      ggplot(tmp, aes(x = Pvalue, y = reorder(Term, Pvalue))) +
      geom_point(aes(size = EnrichGeneNumber, color =
                       Pvalue), position = "dodge") + theme_bw() + theme(
                         axis.text.x = element_text(
                           angle = 85,
                           hjust = 1,
                           size = 15
                         ),
                         axis.text.y = element_text(size = 15),
                         legend.text = element_text(size = 13),
                         legend.title = element_text(size = 13)
                       ) + labs(
                         x =
                           "p-value",
                         y = "Terms",
                         color = "p-value",
                         size = "Gene Count"
                       )
  }
  else{
    p <-
      ggplot(tmp, aes(y = reorder(Term, PAdjust), x = PAdjust)) +
      geom_point(aes(size = EnrichGeneNumber, color =
                       PAdjust), position = "dodge") + theme_bw() + theme(
                         axis.text.x = element_text(
                           angle = 85,
                           hjust = 1,
                           size = 15
                         ),
                         axis.text.y = element_text(size = 15),
                         legend.text = element_text(size = 13),
                         legend.title = element_text(size = 13)
                       ) + labs(
                         x =
                           "pAdjust-value",
                         y = "Terms",
                         color = "p-Adj",
                         size = "Gene Count"
                       )
  }
  return(p)
}


#' Write the tabular form of the pathway or GO term enrichment results
#'
#' @param mrnaObject Object of the enrichment result
#' @param fileName File name of the txt file
#' @param sept File separator, by default, it is tab('\\t')
#' @param type Draw the dot plot according to the p-value or adjusted
#'     p-value ("pvalue", "pAdjust"). Default value is "pAdjust".
#' @param n Number of GO terms or pathways, that ordered by type and has least
#'      number of top p-value
#'
#' @return Text file of the enrichment results in a tabular format
#' 
#' @examples 
#' ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, org_assembly='hg19',
#'    near=TRUE, genetype = 'Ensembl_gene')
#'    
#' writeEnrichment(mrnaObject = ncGO,fileName = "a.txt",sept = '\t')
#'
#' @export
writeEnrichment <-
  function(mrnaObject,
           fileName,
           sept = "\t",
           type = "pAdjust",
           n) {
    if (missing(n)) {
      n <- length(mrnaObject@ID)
    }
    if (n > length(mrnaObject@ID))
      n <- length(mrnaObject@ID)
    if (missing(mrnaObject)) {
      message("Expected input is NoRCE object.")
    }
    if (missing(fileName)) {
      message("Please specify the name of the output file.")
    }
    dd <- topEnrichment(mrnaObject, type, n)
    write.table(dd,
                file = fileName,
                sep = sept,
                row.names = FALSE)
  }

#' Number of top enrichment results of the pathway or GO terms for the given
#' object and the order type - p-value or adjusted p-value.
#'
#' @param mrnaObject Object of the enrichment result
#' @param type Draw the dot plot according to the p-value or adjusted p-value
#'     ("pvalue", "pAdjust")
#' @param n Number of GO terms or pathways, that ordered by type and has
#'     least number of top p-value
#'
#' @return Give top n enrichment results
#' 
#' @importFrom dplyr %>%
#' 
#' @examples 
#' ncGO<-geneGOEnricher(gene = brain_disorder_ncRNA, org_assembly='hg19',
#'    near=TRUE, genetype = 'Ensembl_gene')
#'    
#' result = topEnrichment(mrnaObject = ncGO, type = "pvalue", n = 10)
#'
#' @export
topEnrichment <- function(mrnaObject, type, n) {
  if (missing(mrnaObject)) {
    message("Expected input is NoRCE object.")
  }
  if (missing(type)) {
    message("Type of the sorting method is missing.")
  }
  if (missing(n)) {
    message("Number of top enrichment is missing.")
  }
  if (n > length(mrnaObject@ID))
    n <- length(mrnaObject@ID)
  
  table <- data.frame(gene = rep(
    names(mrnaObject@geneList),
    lapply(mrnaObject@geneList, length)
  ),
  go = unlist(mrnaObject@geneList))
  
  tmp <- mrnaObject@ncGeneList
  is.na(tmp) <- lengths(tmp) == 0
  table1 <- data.frame(gene = rep(
    names(mrnaObject@geneList),
    lapply(tmp, length)
  ),
  go = unlist(tmp))
  
  table <- table[!duplicated(table), ]
  
  table1 <- table1[!duplicated(table1), ]
  
  #split_tibble function is copied from https://stackoverflow.com/a/39638933
  split_tibble <- function(tibble, column = 'col') {
    tibble %>% split(., .[, column]) %>%
      lapply(., function(x)
        x[, setdiff(names(x), column)])
  }
  
  #xy.list <- split(table$go, table$gene)
  #xy.list1 <- split(table1$go, table1$gene)
  
  xy.list <- split_tibble(table, 'gene')
  xy.list1 <- split_tibble(table1, 'gene')
  
  ft <- lapply(xy.list, paste0, collapse = " ")
  ft1 <- lapply(xy.list1, paste0, collapse = " ")
  #x_nc <- as.data.frame(lengths(xy.list1))
  x.1 <- as.data.frame(lengths(xy.list))
  
  rt <- row.names(x.1)
  # rt1 <- row.names(x_nc)
  x.1 <- data.frame(x.1, rt)
  #x_nc <- data.frame(x_nc, rt1)
  x.2 <-
    data.frame(
      go = rep(names(ft), lapply(ft, length)),
      gene = unlist(ft),
      ncgene = unlist(ft1)
    )
  
  dd <- merge(x.1, x.2, by.x = "rt", by.y = "go")
  
  newf <-
    data.frame(
      as.data.frame(mrnaObject@ID),
      as.data.frame(mrnaObject@Term),
      as.double(as.character(mrnaObject@pvalue)),
      as.data.frame(mrnaObject@pAdj),
      as.data.frame(mrnaObject@GeneRatio),
      as.data.frame(mrnaObject@BckRatio)
    )
  
  dd <- merge(newf, dd, by.x = "mrnaObject.ID", by.y = "rt")
  colnames(dd) <-
    c(
      "ID",
      "Term",
      "Pvalue",
      "PAdjust",
      "GeneRatio",
      "BackGroundRatio",
      "EnrichGeneNumber",
      "GeneList",
      "ncGeneList"
    )
  
  ifelse(type == "pvalue", dd <-
           dd[order(dd$Pvalue), ], dd <- dd[order(dd$PAdjust), ])
  
  return(dd[seq_len(n), ])
}

#' Create interaction network for top n enriched GO term:coding RNA or GO-term:noncoding RNA interaction.
#' Nodes are GO term and RNA, edges are interactions between them. Each
#' GO-term is annotated and enriched with the mRNAs provided from the input
#' list.
#'
#' @param mrnaObject Output of enrichment results
#' @param type Sort in terms of p-values or FDR. Possible values
#'     "pvalue", "padjust"
#' @param n Number of top enrichments
#' @param isNonCode Boolean value that checks whether node of the network is
#'      GO-term\& coding or GO-term\& noncoding genes. By default, it is FALSE
#'      so node of the network is GO-term\& coding gene. Otherwise, nodes are
#'      GO-term\& noncoding genes.
#' @param takeID Boolean value that checks the name decision of the GO/pathway
#'      node, GO-term/pathway-term or GO ID-pathway ID. If it is true, name of
#'      the GO/pathway node will be GO ID/pathway ID will be used, otherwise,
#'      name of the GO/pathway node is GO-term. By default, it is FALSE. It is
#'      suggested to used when the GO-term is two long or the GO-term is
#'      missing for the custom enrichment database.
#'
#'
#' @return Network
#' 
#' @importFrom igraph cluster_optimal degree graph_from_data_frame 
#' @importFrom igraph layout_with_fr norm_coords V E V<- E<-
#' @importFrom ggplot2 aes element_text geom_point ggplot labs theme theme_bw
#'
#' @export
createNetwork <-
  function(mrnaObject,
           type = 'pvalue',
           n,
           isNonCode = FALSE,
           takeID = FALSE) {
    if (missing(mrnaObject)) {
      message("Expected input is NoRCE object.")
    }
    if (missing(n)) {
      message(
        "Number of top enrichment is missing.
        Please specify maximum number of enrichment of interest"
      )
    }
    if (n > length(mrnaObject@ID))
      n <- length(mrnaObject@ID)
    tf <- topEnrichment(mrnaObject = mrnaObject,
                        type = type,
                        n = n)
    
    grap <- data.frame()
    if (isNonCode)
    {
      for (i in seq_len(n)) {
        if (!takeID) {
          foo <-
            data.frame(do.call('cbind', strsplit(
              as.character(tf$ncGeneList[i]), ' ', fixed = FALSE
            )))
          ng <- data.frame(go = rep(tf[i, 2], ncol(foo)),
                           gene = foo)
          grap <- rbind(grap, ng)
        }
        else{
          foo <-
            data.frame(do.call('cbind', strsplit(
              as.character(tf$ncGeneList[i]), ' ', fixed = FALSE
            )))
          ng <- data.frame(go = rep(tf[i, 1], ncol(foo)),
                           gene = foo)
          grap <- rbind(grap, ng)
        }
      }
    }
    else
    {
      for (i in seq_len(n)) {
        if (!takeID) {
          foo <-
            data.frame(do.call('cbind', strsplit(
              as.character(tf$GeneList[i]), ' ', fixed = FALSE
            )))
          ng <- data.frame(go = rep(tf[i, 2], ncol(foo)),
                           gene = foo)
          grap <- rbind(grap, ng)
        }
        else{
          foo <-
            data.frame(do.call('cbind', strsplit(
              as.character(tf$GeneList[i]), ' ', fixed = FALSE
            )))
          ng <- data.frame(go = rep(tf[i, 1], ncol(foo)),
                           gene = foo)
          grap <- rbind(grap, ng)
        }
      }
    }
    
    g <- graph_from_data_frame(grap, directed = FALSE)
    deg <- degree(g, mode = "all")
    V(g)$size <- deg * 1.5
    V(g)$frame.color <- "white"
    V(g)$color <- "orange"
    E(g)$arrow.mode <- 0
    E(g)$width <- E(g)$weight * 2
    E(g)$edge.color <- "grey20"
    graphics::plot(g, vertex.label.color = "black")
    
    clp <- cluster_optimal(g)
    V(g)$community <- clp$membership
    colrs <-
      grDevices::adjustcolor(
        c(
          "gray50",
          "tomato",
          "gold",
          "yellowgreen",
          "dark orange",
          "red",
          "blue",
          "green"
        ),
        alpha.f = .6
      )
    l <-  layout_with_fr(g)
    l <-  norm_coords(
      l,
      ymin = -1,
      ymax = 1,
      xmin = -1,
      xmax = 1
    )
    p <-
      graphics::plot(
        g,
        vertex.color = colrs[V(g)$community],
        vertex.label.color = "black",
        vertex.label.cex = .75,
        rescale = FALSE,
        layout = l * 1.0
      )
    return(p)
  }

#' Plot and save the GO term DAG of the top n enrichments in terms of 
#' p-values or adjusted p-values with an user provided format
#'
#' @param mrnaObject Output of enrichment results
#' @param type Sort in terms of p-values or FDR. possible values "pvalue",
#'     "padjust"
#' @param n Number of top enrichments
#' @param filename Name of the DAG file
#' @param imageFormat Image format of the DAG. possible values "png" or "svg"
#' @param p_range Break points for the p-values or FDR. By default
#'     [0.05, 0.001, 0.0005, 0.0001, 0.00005,0.00001,0] is used
#'
#' @return Saves image file in a given format
#'
#'
#' @examples
#' ncRNAPathway<-mirnaPathwayEnricher(gene = brain_mirna,
#'                                    org_assembly = 'hg19',near = TRUE)
#'
#' @export
getGoDag <-
  function(mrnaObject,
           type,
           n,
           filename,
           imageFormat,
           p_range = seq(0, 0.05, by = 0.001))
  {
    if (missing(mrnaObject)) {
      message("Expected input is NoRCE object.")
    }
    if (missing(type)) {
      message("Type of the sorting method is missing.")
    }
    if (missing(n))
      message("Number of top enrichment is missing. ")
    
    if (missing(imageFormat))
      message("Image format is missing. ")
    
    if (n > length(mrnaObject@ID))
      n <- length(mrnaObject@ID)
    if (missing(filename)) {
      message("Name of the output file is missing. ")
    }
    if (missing(imageFormat)) {
      message("Format of the output image is missing.")
    }
    
    dt <- topEnrichment(mrnaObject = mrnaObject,
                        type = type,
                        n = n)
    
    node_color <-
      grDevices::colorRampPalette(c("lightgoldenrodyellow", "orangered"),
                                  bias = 0.5)(length(p_range))
    
    color <- seq_along(dt$Pvalue)#seq_len(2)
    if (type == 'pvalue') {
      for (i in seq_along(dt$Pvalue)) {
        for (j in seq_along(p_range)) {
          if (dt$Pvalue[i] <= p_range[j]) {
            color[i] <- node_color[j]
          }
        }
      }
    }
    else{
      for (i in seq_along(color)) {
        for (j in seq_along(p_range)) {
          if (dt$PAdjust[i] <= p_range[j]) {
            color[i] <- node_color[j]
          }
        }
      }
    }
    gojson <-
      paste0("\"", dt$ID, "\":{\"fill\": \"", color, "\"},", collapse = "")
    
    gojson1 <-
      paste0("{", substr(gojson, 1, nchar(gojson) - 1), "}")
    
    goGraph <-
      RCurl::postForm(
        "http://amigo.geneontology.org/visualize",
        term_data = gojson1,
        inline = "false",
        format = imageFormat,
        mode = 'amigo',
        term_data_type = "json"
      )
    
    if (imageFormat == 'png') {
      png::writePNG(png::readPNG(goGraph), filename)
    }
    else{
      writeLines(goGraph, filename)
    }
  }

#' Display the enriched KEGG diagram of the KEGG pathway. This function is
#' specific to only one KEGG pathway id and identifies the enriched genes
#' in the diagram.
#'
#' @param mrnaObject Output of enrichment results
#' @param pathway Kegg pathway term such as 'hsa04010'
#' @param org_assembly Genome assembly of interest for the analysis.
#'     Possible assemblies are "mm10" for mouse, "dre10" for zebrafish,
#'     "rn6" for rat, "dm6" for fruit fly, "ce11" for worm, "sc3" for yeast,
#'     "hg19" and "hg38" for human
#'
#' @return Shows kegg diagram marked with an enriched genes in a browser
#' 
#' @importFrom utils browseURL read.table write.table
#' 
#' @examples
#' ncRNAPathway<-mirnaPathwayEnricher(gene = brain_mirna,
#'                                    org_assembly = 'hg19',near = TRUE)
#'
#' getKeggDiagram(mrnaObject = ncRNAPathway, org_assembly ='hg19',
#'                pathway = ncRNAPathway@ID[1])
#' @export
#'
getKeggDiagram <-
  function(mrnaObject,
           pathway,
           org_assembly = c("hg19",
                            "hg38",
                            "mm10",
                            "dre10",
                            "rn6",
                            "dm6",
                            "ce11",
                            "sc3")) {
    if (missing(mrnaObject)) {
      message("Expected input is NoRCE object.")
    }
    if (missing(pathway))
      message("Expected pathway ID is missing. Please specify pathway ID")
    
    assembly(org_assembly)
    
    path_index <- which(names(mrnaObject@geneList) == pathway)
    ns <- unique(
      getBM(
        attributes = c("hgnc_symbol", "entrezgene_id"),
        filters = "hgnc_symbol",
        values = mrnaObject@geneList[path_index],
        mart = pkg.env$mart
      )
    )
    n <- paste(ns$entrezgene_id, collapse = '/')
    if(identical(interactive(), TRUE)){
      browseURL(
        paste0(
          "http://www.kegg.jp/kegg-bin/show_pathway?",
          pathway,
          '/',
          n,
          collapse = ''
        )
      )
    }
  }

#' Display the enriched Reactome diagram of the given Reactome pathway id.
#' This function is specific to only one pathway id and identifies the 
#' enriched genes in the diagram.
#'
#' @param mrnaObject Output of enrichment results
#' @param pathway Reactome pathway term
#' @param imageFormat Image format of the diagram. Possible image formats are
#'      'png', 'svg'
#'
#' @return Shows reactome diagram marked with an enriched genes in a browser
#'
#' @examples
#'
#' br_enr<-reactomeEnrichment(genes = breastmRNA,org_assembly='hg19')
#'
#' getReactomeDiagram(mrnaObject = br_enr,pathway = br_enr@ID[1],
#'                    imageFormat = 'png')
#'
#'@export
getReactomeDiagram <- function(mrnaObject, pathway, imageFormat) {
  if (missing(mrnaObject)) {
    message("Expected input is NoRCE object.")
  }
  if (missing(pathway))
    message("Expected pathway ID is missing. Please specify pathway ID")
  if (missing(imageFormat))
    message("Image format is missing. The format of the image should be ")
  
  path_index <- which(names(mrnaObject@geneList) == pathway)
  n <-
    paste(mrnaObject@geneList[[path_index]], collapse = ',', sep = '')
  a <-
    paste0(
      "https://reactome.org/ContentService/exporter/diagram/",
      mrnaObject@ID[path_index],
      ".",
      imageFormat,
      "?flg=",
      n
    )
  if(identical(interactive(), TRUE))
    browseURL(a)
}
