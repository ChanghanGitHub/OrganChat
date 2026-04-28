#' @title Calculate the mean expression level from single-cell metabolite data (for two conditions), for multi-organ OrganChat (MOOC) analysis
#' @description
#' If single-cell metabolite data is available, then this function can be applied to calculate the
#' average metabolite amount in each cell cluster in both conditions. The input must be the normalized amount matrix.
#' It should be applied to each organ if the data is available for multiple organs.
#'
#' @param object an OrganChat object
#' @param scMeta.data normalized metabolite amount matrix (metabolite-by-cell).
#' @param organ.name a char contains the organ names to be calculated.
#' @param idents a char of cluster labels to group the cells.
#' @param conditions.label a char of condition labels.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#' @param type for metabolite data, OrganChat accepts two types of measurements, "amount" and "flux". The default type is "amount",
#' which means the data should be non-negative values. The "flux" data may have both negative (uptake) and positive (release) values.
#' @param scale.max The maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#'
#' @export
#'
MO_cal_scMeta.bugroup <- function(object,
                                  scMeta.data,
                                  organ.name = NULL,
                                  idents = NULL,
                                  conditions.label = NULL,
                                  mean_method = NULL,
                                  type = "amount",
                                  scale.max = 10){

  object <- CA_cal_scMeta.bugroup(object = object,
                                  scMeta.data = scMeta.data,
                                  organ.name = organ.name,
                                  idents = idents,
                                  conditions.label = conditions.label,
                                  mean_method = mean_method,
                                  type = type,
                                  scale.max = scale.max)

  return(object)
}


#' @title Calculate the mean expression level from single-cell transcriptomics data (for two conditions), for multi-organ OrganChat (MOOC) analysis
#' @description
#' Single-cell transcriptomics data of two organs are requested to perform OOC. This function can be applied to
#' calculate the average gene expression of the long-range signals (LS) and downstream genes (DG) in each cell cluster, respectively (for both organs under BOTH conditions).
#' The LS and DG will be inferred using the corresponding OrganChatDB.
#' The outputs will be stored in "Object - LSdata/DGdata - organ".
#'
#' @param object an OrganChat object
#' @param DB the corresponding OrganChatDB used in the analysis.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#' @param scale.max the maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#' @export
MO_cal_scRNA.bygroup <- function(object,
                                 DB = OrganChatDB,
                                 mean_method = NULL,
                                 scale.max = 10){

  # get the LS names from the DB
  LS.list = unique(DB[[1]]$from)

  # for each condition
  for (i in 1:2) {
    for (j in 1:length(object@organs)) {

      organ.select = object@organs[j]

      # select the cells for organ j in condition i
      cell.select = names(object@conditions.list[[j]])[ object@conditions.list[[j]] == object@conditions[i] ]
      idents.select = object@idents.list[[j]][cell.select]
      data.select = object@data.list[[j]][ , cell.select]

      # find the cell barcodes of organ1 and organ2
      celluse <- names(idents.select)
      # find all the genes that will be used
      geneuse <- rownames(data.select)

      meanExpr.df = cal_meanExpr(data = data.select,
                                 celluse = celluse,
                                 geneuse = geneuse,
                                 idents = idents.select,
                                 mean_method = mean_method)

      object@DGdata[[i]][[organ.select]] = meanExpr.df

      if(!is.null(scale.max)){
        meanExpr.df = meanExpr.df/max(meanExpr.df)*scale.max
      }

      # find the matched LS (it may be empty)
      meanExpr.LS = meanExpr.df[rownames(meanExpr.df) %in% LS.list, ]

      if( nrow(meanExpr.LS)>0 ){

        # merge
        data.new = meanExpr.LS
        data.new$Type = "amount"

        if(nrow(object@LSdata[[i]][[organ.select]])==0){
          object@LSdata[[i]][[organ.select]] = data.new
        }else if( length( setdiff(colnames(object@LSdata[[i]][[j]]), colnames(data.new)) )==0 ){
          data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[i]][[j]])) , colnames(object@LSdata[[i]][[j]])]
          object@LSdata[[i]][[j]] = rbind(object@LSdata[[i]][[j]], data.new)
        }else{
          warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
        }
      }

    }

  }
  return(object)
}


#' @title Merge existing cluster-level metabolite amount/flux data (for two conditions), for multi-organ OrganChat (MOOC) analysis
#' @description
#' If the metabolite amount/flux data for both organs are ready, and are not calculated using th cal_scMeta.bygroup function,
#' they can still be integrated into the OrganChat object using this function for comparison analysis of two conditions.
#' This is usually because when single-cell metabolite data is unavailable, then the metabolite flux data can be inferred via METAFlux using OrganChat built-in modules.
#' Note the metabolites are only considered as the LS, the outputs will be stored in "Object - LSdata - organ".
#'
#' @param object an OrganChat object
#' @param metabolitedata.list.condition1 a list of metabolite data (average amount/flux for each cluster) for two organs in condition 1.
#' @param metabolitedata.list.condition2 a list of metabolite data (average amount/flux for each cluster) for two organs in condition 2.
#' @param type for metabolite data, OrganChat accepts two types of measurements, "amount" and "flux". The default type is "amount",
#' which means the data should be non-negative values. The "flux" data may have both negative (uptake) and positive (release) values.
#' @param scale.max the maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#'
#' @return an OrganChat object
#' @export
MO_merge_metabolite.bygroup <- function(object,
                                        metabolitedata.list.condition1 = NULL,
                                        metabolitedata.list.condition2 = NULL,
                                        type = "amount",
                                        scale.max = 10){

  if(length(metabolitedata.list.condition1)!=length(object@organs) | length(metabolitedata.list.condition2)!=length(object@organs)){
    stop("The input metabolitedata.list must have the data for every organ.")
  }

  # check the input 'type'
  if(is.null(type)){ type = "amount" }
  if( !type %in% c("amount", "flux") ){
    stop("The input 'type' must be either 'amount' or 'flux'.")
  }

  # scale the data
  if(!is.null(scale.max)){
    m1 = c()
    m2 = c()
    for (i in 1:length(object@organs)) {
      m1 = c( m1, max(abs(metabolitedata.list.condition1[[i]])) )
      m2 = c( m2, max(abs(metabolitedata.list.condition2[[i]])) )
    }
    m1.max = max( m1 )
    m2.max = max( m2 )

    for (i in 1:length(object@organs)) {
      metabolitedata.list.condition1[[i]] = metabolitedata.list.condition1[[i]]/m1.max*scale.max
      metabolitedata.list.condition2[[i]] = metabolitedata.list.condition2[[i]]/m2.max*scale.max
    }
  }

  for (i in 1:2) { # for condition i
    for (j in 1:length(object@organs)) { # for organ j

      if(i == 1){
        data.new = metabolitedata.list.condition1[[j]]
      }else{
        data.new = metabolitedata.list.condition2[[j]]
      }

      data.new$Type = type

      if(nrow(object@LSdata[[i]][[j]])==0){
        object@LSdata[[i]][[j]] = data.new
      }else if( length( setdiff(colnames(object@LSdata[[i]][[j]]), colnames(data.new)) )==0 ){
        data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[i]][[j]])) , colnames(object@LSdata[[i]][[j]])]
        object@LSdata[[i]][[j]] = rbind(object@LSdata[[i]][[j]], data.new)
      }else{
        warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
      }
    }
  }

  return(object)
}


#' @title Multi-organ OrganChat (MOOC) analysis of two conditions
#' @description
#' Calculate the communication probability, perform bootstrap tests, and filter the LS-T pathways for each organ and condition. Store the results as a list.
#'
#' @param object an OrganChat object
#' @param DB the corresponding OrganChatDB used in the analysis.
#' @param LS the list of LS used for inferring the pathways, the input is NULL by default, which means all will be used.
#' @param receptor the list of Receptor used for inferring the pathways, the input is NULL by default, which means all will be used.
#' @param tf the list of Signaling Effector used for inferring the pathways, the input is NULL by default, which means all will be used.
#' @param target  the list of Target used for inferring the pathways, the input is NULL by default, which means all will be used.
#' @param K half-saturation constant of the Hill function, the default value is 0.5.
#' @param N Hill coefficient, the default value is 2.
#' @param R the number of bootstrap testing to perfomr, the default value is 5.
#' @param CI value of confidence interval, CI = 0.95 by default.
#' @param CA_cutoff_value a vector of 4 cutoff values for LS, R, SE, T, respectively. For each LS-T pathway, each component must have a higher
#' expression than the cutoff value in at least one condition. The default value is NULL.
#' @param CA_cutoff_chatP LS-T pathways with a communication probability lower than this cutoff value under both conditions will be removed.
#' The default value is NULL.
#' @param CA_cutoff_sd LS-T pathways with a standard deviation (calculated based on bootstrap results) higher than this cutoff value under at least one condition will be removed.
#' The default value is NULL.
#' @param CA_cutoff_CI S-T pathways with a total number of components lying within the confidence interval smaller than this value under at least one condition will be removed.
#' The default value is 0. The input must be 0,1,2, or 3.
#'
#' @importFrom utils combn
#' @importFrom forcats fct_drop
#' @importFrom boot boot
#'
#' @return a list
#' @export
MO_pairwise_analysis <- function(object,
                                 DB = OrganChatDB,
                                 LS = NULL,
                                 receptor = NULL,
                                 tf = NULL,
                                 target = NULL,
                                 K = 0.5,
                                 N = 2,
                                 R = 5,
                                 CI = 0.95,
                                 CA_cutoff_value = NULL, # "&" logic between cutoff conditions
                                 CA_cutoff_chatP = NULL, # "|" logic between cutoff conditions
                                 CA_cutoff_sd = NULL,
                                 CA_cutoff_CI = 0){

  MO.net.list <- list()

  organ.combn <- combn(object@organs, 2)
  n.combn <- ncol(combn(object@organs, 2))

  for (i in 1:n.combn) {
    organ1 = organ.combn[1, i]
    organ2 = organ.combn[2, i]

    print(paste("Analyzing the ", i, " of ", n.combn, " combination."))

    CAOC.object = create_OrganChat(object = object@data.list[c(organ1, organ2)],
                                   organs = c(organ1, organ2),
                                   idents.list = object@idents.list[c(organ1, organ2)],
                                   ComparisonAnalysis = TRUE,
                                   conditions = object@conditions,
                                   conditions.list = object@conditions.list[c(organ1, organ2)],
                                   assay = NULL,
                                   do.sparse = T)

    CAOC.object@LSdata[[1]] = object@LSdata[[1]][c(organ1, organ2)]
    CAOC.object@LSdata[[2]] = object@LSdata[[2]][c(organ1, organ2)]

    CAOC.object@DGdata[[1]] = object@DGdata[[1]][c(organ1, organ2)]
    CAOC.object@DGdata[[2]] = object@DGdata[[2]][c(organ1, organ2)]

    CAOC.object = CA_net0_infer(object = CAOC.object,
                                DB = DB,
                                LS = LS,
                                receptor = receptor,
                                tf = tf,
                                target = target)

    CAOC.object = CA_net_analysis(object = CAOC.object,
                                  K = K,
                                  N = N,
                                  R = R,
                                  CI = CI,
                                  CA_cutoff_value = CA_cutoff_value,
                                  CA_cutoff_chatP = CA_cutoff_chatP,
                                  CA_cutoff_sd = CA_cutoff_sd,
                                  CA_cutoff_CI = CA_cutoff_CI)

    MO.net.list[[i]] = CAOC.object@net
  }

  MO.net.list[[(n.combn+1)]] = organ.combn

  return(MO.net.list)
}


#' @title Construct a list for the directed & weighted GRN from the output of MO_pairwise_analysis function.
#' @description
#' This function generates a list of directed & weighted GRNs from MO_pairwise_analysis function outputs. Each will have a data frame for edges and nodes, respectively.
#'
#' @param object an OrganChat object, normally the output of MO_merge_metabolite.bygroup function.
#' @param MO.net.list output (a list) of the MO_pairwise_analysis function.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#'
#' @importFrom data.table rbindlist
#' @import dplyr
#'
#' @return a list
#'
#' @export
#'
MO_networkconstruct <- function(object,
                                MO.net.list,
                                category = "All",
                                chatP_cutoff = NULL){

  organ.combn <- MO.net.list[[length(MO.net.list)]]
  n.combn <- ncol(organ.combn)

  EdgeNode.list = list()
  Edge.list = list(Condition1 = list(), Condition2 = list())
  Node.list = list(Condition1 = list(), Condition2 = list())

  for (i in 1:n.combn){
    organ1 = organ.combn[1, i]
    organ2 = organ.combn[2, i]

    CAOC.object = create_OrganChat(object = object@data.list[c(organ1, organ2)],
                                   organs = c(organ1, organ2),
                                   idents.list = object@idents.list[c(organ1, organ2)],
                                   ComparisonAnalysis = TRUE,
                                   conditions = object@conditions,
                                   conditions.list = object@conditions.list[c(organ1, organ2)],
                                   assay = NULL,
                                   do.sparse = T)

    CAOC.object@net <- MO.net.list[[i]]

    # generate the edge-node table for both conditions
    EdgeNode.list[[1]] <- construct_EdgeNode(CAOC.object,
                                             category = category,
                                             chatP_cutoff = chatP_cutoff,
                                             single.condition = TRUE)

    EdgeNode.list[[2]] <- construct_EdgeNode(CAOC.object,
                                             category = category,
                                             chatP_cutoff = chatP_cutoff,
                                             single.condition = FALSE)

    # here we need to add the cell population for each organ in each condition manually
    Population_organ1 <- unname(table(CAOC.object@conditions.list[[1]]))
    Population_organ2 <- unname(table(CAOC.object@conditions.list[[2]]))
    # condition 1
    EdgeNode.list[[1]]$nodes$Population <- c(Population_organ1[1], Population_organ2[1])
    # condition 2
    EdgeNode.list[[2]]$nodes$Population <- c(Population_organ1[2], Population_organ2[2])

    Edge.list[[1]][[i]] = EdgeNode.list[[1]][[1]]
    Node.list[[1]][[i]] = EdgeNode.list[[1]][[2]]

    Edge.list[[2]][[i]] = EdgeNode.list[[2]][[1]]
    Node.list[[2]][[i]] = EdgeNode.list[[2]][[2]]
  }

  output = list(Condition1 = list(), Condition2 = list())

  output[[1]][[1]] = rbindlist(Edge.list[[1]], use.names = TRUE)
  output[[1]][[2]] = rbindlist(Node.list[[1]], use.names = TRUE)
  output[[1]][[1]] = output[[1]][[1]][!duplicated(output[[1]][[1]]), ]
  output[[1]][[2]] = output[[1]][[2]][!duplicated(output[[1]][[2]]), ]
  names(output[[1]]) = c("edge", "node")

  output[[2]][[1]] = rbindlist(Edge.list[[2]], use.names = TRUE)
  output[[2]][[2]] = rbindlist(Node.list[[2]], use.names = TRUE)
  output[[2]][[1]] = output[[2]][[1]][!duplicated(output[[2]][[1]]), ]
  output[[2]][[2]] = output[[2]][[2]][!duplicated(output[[2]][[2]]), ]
  names(output[[2]]) = c("edge", "node")

  for (k in 1:length(object@organs)) {
    max1 = max(output[[1]]$edge[ output[[1]]$edge$from == object@organs[k] & output[[1]]$edge$to == object@organs[k] ]$ave.weight)
    max2 = max(output[[2]]$edge[ output[[2]]$edge$from == object@organs[k] & output[[2]]$edge$to == object@organs[k] ]$ave.weight)

    output[[1]]$edge = output[[1]]$edge[ !(output[[1]]$edge$from == object@organs[k] &
                                             output[[1]]$edge$to == object@organs[k] &
                                             output[[1]]$edge$ave.weight < max1) , ]

    output[[2]]$edge = output[[2]]$edge[ !(output[[2]]$edge$from == object@organs[k] &
                                             output[[2]]$edge$to == object@organs[k] &
                                             output[[2]]$edge$ave.weight < max2) , ]
  }

  names(output) = object@conditions

  return(output)
}


#' @title Chord diagrams of the organ-level communications of two conditions
#' @description
#' Chord diagrams of the organ-level communications of two conditions.
#'
#' @param network.input output (a list) of the MO_networkconstruct function.
#' @param style style select which metric to showcase using the edge width, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", "ave.sd", or NULL (default, meaning "total.counts" will be used).
#' @param color.platte select the colors, the default selection is brewer.pal(8, "Set2").
#' @param scale.edge.width a multiplier to scale the edge width, the default value is 5.
#' @param scale.vertex.size a multiplier to scale the vertex size, the default value is 30.
#' @param vertex.frame.color argument of the plot function, the default value is "white".
#' @param vertex.label.color argument of the plot function, the default value is "black".
#' @param vertex.label.dist argument of the plot function, the default value is 3.
#' @param edge.arrow.size argument of the plot function, the default value is 0.6.
#' @param edge.curved argument of the plot function, the default value is 0.2.
#' @param margin argument of the plot function, the default value is 0.1.
#' @param vertex.max scale the vertex size using the maximum in "condition1" or "condition2" or NULL (default, meaning using the maximum of both conditions).
#' @param edge.max scale the edge width using the maximum in "condition1" or "condition2" or NULL (default, meaning using the maximum of both conditions).
#'
#' @importFrom graphics legend par strwidth title
#' @import RColorBrewer
#' @importFrom scales rescale
#' @importFrom igraph graph_from_data_frame V E ends
#' @import dplyr
#' @import ggplot2
#' @import stringr
#'
#' @export
#'
MO_networkplot <- function(network.input,
                           style = "total.counts",
                           color.platte = NULL,
                           scale.edge.width = 5,
                           scale.vertex.size = 30,
                           vertex.frame.color = "white",
                           vertex.label.color="black",
                           vertex.label.dist = 3,
                           edge.arrow.size = 0.6,
                           edge.curved = 0.2,
                           margin = 0.1,
                           vertex.max = NULL,
                           edge.max = NULL){

  if( is.null(style)){
    style = "total.counts"
    print("The default 'style' is 'total.counts'.")
  }else if( !style %in% c("ave.weight", "total.counts", "total.chatP") ){
    stop("The input 'style' should be 'total.counts', 'ave.weight', or 'total.chatP'.")
  }

  if(is.null(color.platte)){
    color.platte = RColorBrewer::brewer.pal(8, "Set2")
  }

  # function for label location
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x+start) %% (2*pi)*direction
    c.rotate(scales::rescale(x, c(0, 2*pi), range(x)))
  }

  par(mfrow=c(1,2))
  Vmax = c()
  Emax = c()

  for (i in 1:2) {
    df.plot = network.input[[i]]
    # generate links data frame
    links <- data.frame(source = df.plot$edge$from, target = df.plot$edge$to, value = df.plot$edge[[style]])
    links <- links[links$value > 0, ]
    # generate nodes data frame
    nodes <- data.frame(name = df.plot$node$Organ, population = df.plot$node$Population)
    nodes <- nodes[!duplicated(nodes), ]
    # constrict network object
    network <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)

    Vmax = c(Vmax, max(igraph::V(network)$population))
    Emax = c(Emax, max(igraph::E(network)$value))
  }

  if(!is.null(vertex.max)){
    if(vertex.max == "condition1"){
      Vmax = Vmax[1]
    }else if(vertex.max == "condition2"){
      Vmax = Vmax[2]
    }else{
      Vmax = max(Vmax)
    }
  }

  if(!is.null(edge.max)){
    if(edge.max == "condition1"){
      Emax = Emax[1]
    }else if(edge.max == "condition2"){
      Emax = Emax[2]
    }else{
      Emax = max(Emax)
    }
  }

  for (i in 1:2) {

    df.plot = network.input[[i]]
    # generate links data frame
    links <- data.frame(source = df.plot$edge$from, target = df.plot$edge$to, value = df.plot$edge[[style]])
    links <- links[links$value > 0, ]
    # generate nodes data frame
    nodes <- data.frame(name = df.plot$node$Organ, population = df.plot$node$Population)
    nodes <- nodes[!duplicated(nodes), ]
    # constrict network object
    network <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)

    vcolor <- color.platte[1:nrow(nodes)]
    igraph::V(network)$color <- vcolor[factor(igraph::V(network)$name)]
    igraph::V(network)$size <- igraph::V(network)$population/max(Vmax)
    igraph::E(network)$width <- igraph::E(network)$value/max(Emax)
    # for edge color
    edge.start <- igraph::ends(network, es=igraph::E(network), names=F)[,1]
    edge.col <- igraph::V(network)$color[edge.start]
    # for label location
    lab.locs <- radian.rescale(x=1:nrow(nodes), direction=-1, start=0)

    # for self-loop angle
    myMatrix = igraph::ends(network, es=igraph::E(network), names=F)
    selfloopAngles <- numeric(nrow(myMatrix))
    M <- nrow(myMatrix)
    for(j in 1:nrow(myMatrix)) {
      if(myMatrix[j,1]==myMatrix[j,2]){
        selfloopAngles[j] <- (2*pi*(M - j)/M)
      }
    }

    title = paste0(names(network.input)[i]," - ", style)

    plot(network,
         edge.width = igraph::E(network)$width*scale.edge.width,
         vertex.size = igraph::V(network)$size*scale.vertex.size,
         vertex.frame.color = vertex.frame.color,
         vertex.label.color = vertex.label.color,
         vertex.label.dist = vertex.label.dist,
         vertex.label.degree = lab.locs,
         edge.color = edge.col,
         edge.arrow.size = edge.arrow.size,
         edge.loop.angle = selfloopAngles,
         layout = layout.circle,
         main = title,
         edge.curved = edge.curved,
         margin = margin)
  }

}
