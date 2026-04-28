#' @title Select LS-T pathways for visualization
#' @description
#' Select the preferred pathways from a data.frame, based on user-select preference and LS type.
#' @param df a data frame of LS-T pathways
#' @param pref a char indicating the preference when selecting the pathways, the input should be one of the following: "chatP", "LS", "Receptor", "TF", "Target"
#' The default value is "LS", meaning OrganChat will try to cover as many LS as possible, by selecting one LS-T pathway of each LS (ranked by the communication probability) until reaching the limit n.path.
#' If the input is "chatP", then OrganChat will select the top n.path LS-T pathways based on the communication probability.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param n.path total number pathways to plot, the default value is 20.
#'
#' @import dplyr
#'
#' @return a data frame
#'
#' @export
#'
select_paths <- function(df,
                         pref = "LS",
                         category = "All",
                         n.path = 20){

  # check input 'pref'
  if( !(pref %in% c("chatP", "LS", "Receptor", "TF", "Target")) == TRUE ){
    message("The input 'pref' is invalid, set it to 'LS'.")
    pref = "LS"
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    pref = "All"
  }

  if( category != "All"){
    ind = which(grepl(category, df$category, fixed=TRUE))
    df = df[ind, ]
  }

  # check input 'n.path'
  if(is.null(n.path) & nrow(df)>=n.path){
    n.path = 20
  }else if( nrow(df)<n.path ){
    message("The input 'n.path' is larger than the nrow(df), output the 'df'.")
    return(df)
  }

  if( pref=="chatP" ){
    # select the top n.path rows
    df %>%
      dplyr::arrange(-chatP) %>%
      slice_head(n = n.path) -> df.s
  }else{

    df$select <- unlist(as.vector(df[ , pref]))

    if( length(unique(df$select))>=n.path ){
      message("The number of components involoved are larger then n.path, select the row with highest chatP.")
      # for each component, select the pathway with the highest chatP
      df %>%
        group_by(select) %>%
        dplyr::arrange(-chatP) %>%
        slice_head(n = 1) %>%
        ungroup() -> df.s

      # select the top n.path rows
      df.s %>%
        dplyr::arrange(-chatP) %>%
        slice_head(n = n.path) -> df.s
    }else{

      df %>%
        group_by(select) %>%
        dplyr::arrange(-chatP) %>%
        slice_head(n = 1) %>%
        ungroup() -> df.0
      df.1 <- anti_join(df, df.0, by = join_by(Path))

      q = nrow(df.0)
      p = n.path - q
      r = length(unique(df.1$select))

      if(p<0){
        df.1 %>%
          dplyr::arrange(-chatP) -> df.1
        df.s = df.1[ c(1:n.path), ]
      }else{
        while (p>=r & p>0) {
          df.1 %>%
            group_by(select) %>%
            dplyr::arrange(-chatP) %>%
            slice_head(n = 1) %>%
            ungroup() -> df.new

          df.0 <- rbind(df.0, df.new)
          df.1 <- anti_join(df.1, df.new, by = join_by(Path))

          q = nrow(df.0)
          p = n.path - q
          r = length(unique(df.1$select))
        }

        df.1 %>%
          group_by(select) %>%
          dplyr::arrange(-chatP) %>%
          slice_head(n = 1) %>%
          ungroup() -> df.new
        # select the top n.path rows
        df.new %>%
          dplyr::arrange(-chatP) %>%
          slice_head(n = p) -> df.new

        df.s <- rbind(df.0, df.new)
      }

    }
  }
  return(df.s[ , 1:ncol(df.s)-1])
}


#' @title River plot of selected LS-T pathways
#' @description
#' River plot of the selected pathways from a data.frame (using ggplot), based on user-select preference and LS type.
#'
#' @param df a data frame of LS-T pathways from OrganChat object
#' @param pref a char indicating the preference when selecting the pathways, the input should be one of the following: "chatP", "LS", "Receptor", "TF", "Target"
#' The default value is "LS", meaning OrganChat will try to cover as many LS as possible, by selecting one LS-T pathway of each LS (ranked by the communication probability) until reaching the limit n.path.
#' If the input is "chatP", then OrganChat will select the top n.path LS-T pathways based on the communication probability.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param n.path total number pathways to plot, the default value is 20.
#' @param text.color the default value is TRUE.
#' @param LS_name a boolean value to specify if the metabolite name is used, the default value is TRUE. Otherwise the HMDB ID will be presented.
#' @param width argument of the geom_stratum function, the default value is 1/3.
#' @param size argument of the geom_stratum function, the default value is 0.1.
#' @param expand argument of the scale_x_discrete function, the default value is c(.01, .05).
#' @param type argument of the scale_fill_continuous function, the default value is "viridis".
#' @param palette argument of the scale_color_distiller function, the default value is "RdPu".
#'
#' @importFrom reshape2 melt
#' @import ggalluvial
#' @import ggplot2
#'
#' @return a ggplot2 object
#'
#' @export
#'
OC_riverplot <- function(df,
                         pref = "LS",
                         category = "All",
                         n.path = 20,
                         text.color = TRUE,
                         LS_name = TRUE,
                         width = 1/3,
                         size = 0.1,
                         expand = c(.01, .05),
                         type = "viridis",
                         palette = "RdPu"){

  df.plot <- select_paths(df, pref, category, n.path)

  df_long_1 <- melt(df.plot, id = c("Path", "LS_name",
                                    "Sender.group", "Receiver.group",
                                    "LS.value", "Receptor.value", "TF.value", "Target.value", "chatP"),
                    measure.vars = c("LS", "Receptor", "TF", "Target"),
                    variable.name = "Role", value.name = "Name")

  df_long_2 <- melt(df.plot, id = c("Path", "LS_name", "Sender.group", "chatP"),
                    measure.vars = c( "LS.value", "Receptor.value", "TF.value", "Target.value"),
                    variable.name = "Role.value", value.name = "Expression")

  df_long <- df_long_1[ , c("Path", "LS_name", "Sender.group", "Receiver.group", "chatP", "Role", "Name")]
  df_long$Expression <- df_long_2$Expression

  # check the input of 'text.color'
  if(is.logical(text.color)==FALSE){
    text.color = TRUE
  }

  # check the input of 'LS_name'
  if(is.logical(LS_name)==FALSE){
    LS_name = TRUE
  }

  if(LS_name == TRUE){

    df_long[df_long$Role == "LS", ]$Name <- df_long[df_long$Role == "LS", ]$LS_name
    df_long$Name <- as.factor(df_long$Name)
  }else{ df_long$Name <- as.factor(df_long$Name) }

  # df_long %>%
  #   arrange(desc(chatP) )

  df_long = df_long[order(-df_long$chatP), ]

  # adapted ggplot function
  if(text.color == TRUE){
    ggplot(df_long,
           aes(x = Role, stratum = Name, alluvium = Path, color = Expression)) +
      geom_flow(aes(fill = chatP)) +
      geom_stratum(width = width, size = size) +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("LS", "Receptor", "TF", "Target"), expand = expand) +
      scale_fill_continuous(type = type) +
      scale_color_distiller(palette = palette) +
      theme_minimal()
  }else{
    ggplot(df_long,
           aes(x = Role, stratum = Name, alluvium = Path)) +
      scale_fill_continuous(type = type) +
      geom_flow(stat = "alluvium", aes(fill = chatP)) +
      geom_stratum(width = width, size = size) +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("LS", "Receptor", "TF", "Target"), expand = expand) +
      theme_minimal()
  }
}


#' @title Select LS-T pathways from one direction table for plot
#' @description
#' Select LS-T pathways from one direction, calculate the total count (total.counts), average communication probability (ave.weight), and total communication probability (total.chatP).
#'
#' @param object an OrganChat object
#' @param net a list of data frames under "object - net - Condition - direction" layer
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#'
#' @import dplyr
#'
#' @return a data frame
#'
#' @export
#'
OC_singledir_netanalysis <- function(object,
                                     net,
                                     dir,
                                     chatP_cutoff = NULL,
                                     category = "All"){
  # check the input 'dir'
  if( !dir %in% c("forward", "backward", "self_organ1", "self_organ2") | is.null(dir)){
    stop("The input 'dir' should be one of the following:
         forward, backward, self_organ1, or self_organ2.")
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    category = "All"
  }

  ave.weight <- c()
  total.counts <- c()
  total.chatP <- c()
  ave.sd <- c()
  sender.group <- c()
  receiver.group <- c()
  sender.organ <- c()
  receiver.organ <- c()

  # assign organs & cluster labels based on the direction
  if( dir=="forward" ){
    S.organ <- object@organs[1] # from organ1 (Metabolite) to organ2 (RNA)
    R.organ <- object@organs[2]
  } else if( dir=="backward" ){
    S.organ <- object@organs[2] # from organ2 (Metabolite) to organ1 (RNA)
    R.organ <- object@organs[1]
  } else if( dir=="self_organ1" ){
    S.organ <- object@organs[1]
    R.organ <- object@organs[1]
  } else if( dir=="self_organ2" ){
    S.organ <- object@organs[2]
    R.organ <- object@organs[2]
  }
  S.group <- levels(object@idents.list[[S.organ]])
  R.group <- levels(object@idents.list[[R.organ]])

  for (j in 1:length(net)) { # check for each cluster pair

    df.path <- net[[j]] # for specific cluster-cluster network table

    if( nrow(df.path)==0){
      ave.weight <- c(ave.weight, 0)
      total.counts <- c(total.counts, 0)
      total.chatP <- c(total.chatP, 0)
      ave.sd <- c(ave.sd, 0)
    }else{
      # apply filter "chetP_cutoff"
      if(!is.null(chatP_cutoff)){
        df.path = df.path[ df.path$chatP >= chatP_cutoff, ]
      }
      # apply filter "category"
      if( category != "All"){
        ind = which(grepl(category, df.path$category, fixed=TRUE))
        df.path = df.path[ind, ]
      }
      # ave.weight & total.counts
      if( nrow(df.path)>0 ){
        ave.weight <- c(ave.weight, sum(df.path$chatP)/nrow(df.path))
        total.counts <- c(total.counts, nrow(df.path))
        total.chatP <- c(total.chatP, sum(df.path$chatP))
        ave.sd <- c(ave.sd, sum(df.path$bootstrap_sd)/nrow(df.path))
      }else{
        ave.weight <- c(ave.weight, 0)
        total.counts <- c(total.counts, 0)
        total.chatP <- c(total.chatP, 0)
        ave.sd <- c(ave.sd, 0)
      }
    }
  }


  # sender.group $ receiver.group
  for (l in 1:length(S.group)) {
    for (m in 1:length(R.group)) {
      sender.group <- c(sender.group, S.group[l])
      receiver.group <- c(receiver.group, R.group[m])

      sender.organ <- c(sender.organ, S.organ)
      receiver.organ <- c(receiver.organ, R.organ)
    }
  }

  df.output = data.frame(from = paste0(sender.group, "(s)"),
                         to = paste0(receiver.group, "(r)"),
                         ave.weight = ave.weight, total.counts = total.counts,
                         total.chatP = total.chatP, ave.sd = ave.sd,
                         sender.organ = sender.organ, receiver.organ = receiver.organ)

  return(df.output)
}


#' @title Chord diagram of cluster-level OOC
#' @description
#' Chord diagram of cluster-level OOC under the selected direction (using chordDiagram function).
#'
#' @param object an OrganChat object
#' @param net a list of data frames under "object - net - Condition - direction" layer
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param style select which metric to showcase using the edge width, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", "ave.sd", or NULL (default, meaning "total.counts" will be used).
#' @param big.gap argument of the chordDiagram function, the default value is 20.
#' @param directional argument of the chordDiagram function, the default value is 1.
#' @param grid.col argument of the chordDiagram function, the default value is NULL.
#' @param direction.type argument of the chordDiagram function, the default value is c("diffHeight", "arrows").
#' @param link.arr.type argument of the chordDiagram function, the default value is  "big.arrow".
#' @param link.target.prop argument of the chordDiagram function, the default value is TRUE.
#' @param small.gap argument of the chordDiagram function, the default value is 1.
#' @param annotationTrack argument of the chordDiagram function, the default value is "grid".
#' @param scale argument of the chordDiagram function, the default value is FALSE.
#' @param title.font argument of the circos.trackPlotRegion function, the default value is 2.
#' @param title.cex argument of the circos.trackPlotRegion function, the default value is 0.2.
#' @param title.line argument of the circos.trackPlotRegion function, the default value is 0.
#' @param title.adj argument of the circos.trackPlotRegion function, the default value is 0.5.
#' @param label.cex argument of the circos.trackPlotRegion function, the default value is 0.75.
#' @param circle.margin argument of the circos.par function, the default value is 0.5.
#'
#' @importFrom graphics legend par strwidth title
#' @import circlize
#'
#' @export
#'
OC_chorddiagram <- function(object,
                            net, # a list of sig network of one direction
                            dir = NULL,
                            category = "All",
                            chatP_cutoff = NULL,
                            style = NULL,
                            big.gap = 20,
                            directional = 1,
                            grid.col = NULL,
                            direction.type = c("diffHeight", "arrows"),
                            link.arr.type = "big.arrow",
                            link.target.prop = TRUE,
                            small.gap = 1,
                            annotationTrack = "grid",
                            scale = FALSE,
                            title.font = 2,
                            title.cex = 0.2,
                            title.line = 0,
                            title.adj = 0.5,
                            label.cex = 0.75,
                            circle.margin = 0.5){


  # check the input 'dir'
  if( !dir %in% c("forward", "backward", "self_organ1", "self_organ2") | is.null(dir)){
    stop("The input 'dir' should be one of the following:
         forward, backward, self_organ1, or self_organ2.")
  }

  if( length(net)>0 ){ # the results for this direction is not empty

    df.output <- OC_singledir_netanalysis(object = object,
                                          net = net,
                                          dir = dir,
                                          chatP_cutoff = chatP_cutoff,
                                          category = category)

  }else{
    stop("No signaling pathway found!")
  }

  # check the input 'style'
  if( is.null(style)){
    style = "total.counts"
    print("The default 'style' is 'total.counts'.")
  }else if( !style %in% c("ave.weight", "total.counts", "total.chatP", "ave.sd") ){
    stop("The input 'style' should be 'total.counts', 'ave.weight', 'total.chatP', or 'ave.sd'.")
  }

  if(style == "total.counts"){
    df.plot = data.frame(from = df.output$from,
                         to = df.output$to,
                         value = df.output$total.counts)
    df.plot = df.plot[df.plot$value >0, ]
  }else if(style == "ave.weight"){
    df.plot = data.frame(from = df.output$from,
                         to = df.output$to,
                         value = df.output$ave.weight)
    df.plot = df.plot[df.plot$value >0, ]
  }else if(style == "total.chatP"){
    df.plot = data.frame(from = df.output$from,
                         to = df.output$to,
                         value = df.output$total.chatP)
    df.plot = df.plot[df.plot$value >0, ]
  }else{
    df.plot = data.frame(from = df.output$from,
                         to = df.output$to,
                         value = df.output$ave.sd)
    df.plot = df.plot[df.plot$value >0, ]
  }

  if(nrow(df.plot)==0){
    stop("No pathway found.")
  }

  circos.par(start.degree = 180, circle.margin = circle.margin)
  chordDiagram(df.plot,
               big.gap = big.gap,
               directional = directional,
               grid.col = grid.col,
               direction.type = direction.type,
               link.arr.type = link.arr.type,
               link.target.prop = link.target.prop,
               small.gap = small.gap,
               annotationTrack = annotationTrack,
               preAllocateTracks = list( track.height = max(strwidth(unlist(dimnames(df.plot)))) ),
               scale = scale)

  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = label.cex)
    # circos.axis(h = "top", labels.cex = 0.5, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)

  # legend("right",
  #        pch = pch,
  #        legend = c( unique(df.plot$from), unique(df.plot$to) ),
  #        bty = bty,
  #        col = col,
  #        cex = cex,
  #        pt.cex = pt.cex,
  #        border = border,
  #        xjust = xjust,
  #        text.width = text.width,
  #        text.font = text.font)

  title(paste0("From top (", unique(df.output$sender.organ), ") to bottom (", unique(df.output$receiver.organ), ")") ,
        font = title.font, cex = title.cex, line = title.line, adj = title.adj)

  circos.clear()

  # return(p)
}


#' @title Bubble-heatmap showing two metrics of the cluster-level OOC
#' @description
#' Bubble-heatmap showing two metrics by bubble size and color of the cluster-level OOC (using ggplot function).
#'
#' @param object an OrganChat object
#' @param net a list of data frames under "object - net - Condition - direction" layer
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param size select which metric to showcase using the bubble size, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", or "ave.sd". The default value is "ave.weight".
#' @param color select which metric to showcase using the bubble color, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", or "ave.sd". The default value is "total.chatP".
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return a ggplot2 object
#'
#' @export
#'
OC_bubbleheatmap <- function(object,
                             net, # a list of sig network of one direction
                             dir = NULL,
                             chatP_cutoff = NULL,
                             category = "All",
                             size = "ave.weight",
                             color = "total.chatP"){

  # check the input 'dir'
  if( !dir %in% c("forward", "backward", "self_organ1", "self_organ2") | is.null(dir)){
    stop("The input 'dir' should be one of the following:
         forward, backward, self_organ1, or self_organ2.")
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    pref = "All"
  }

  if( length(net)>0 ){ # the results for this direction is not empty

    df.output <- OC_singledir_netanalysis(object = object,
                                          net = net,
                                          dir = dir,
                                          chatP_cutoff = chatP_cutoff,
                                          category = category)

  }else{
    stop("No signaling pathway found!")
  }

  if( !size %in% c("ave.weight", "total.counts", "total.chatP", "ave.sd") |
      !color %in% c("ave.weight", "total.counts", "total.chatP", "ave.sd") ){
    stop("The input 'size' and 'color' should be 'total.counts', 'ave.weight', 'total.chatP', or 'ave.sd'.")
  }

  indexname <- c("ave.weight", "total.counts", "total.chatP", "ave.sd")
  legendtitle <- c("Ave. sig. prob.", "Total path", "Total sig. prob.", "Ave. bootstrap sd")
  size.name <- legendtitle[which(indexname == size)]
  color.name <- legendtitle[which(indexname == color)]

  df.plot = data.frame(from = df.output$from,
                       to = df.output$to,
                       size = df.output[, size],
                       color = df.output[, color])
  df.plot = df.plot[df.plot$size>0 | df.plot$color>0, ]

  p <- ggplot(df.plot, aes(to, from, size = size, color = color)) +
    geom_point() +
    xlab("Receiver group") + ylab("Sender group") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_viridis_c(option = "C", name = color.name) +
    guides(size = guide_legend(title = size.name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}


#' @title Construct a list of data frames recording the edges & nodes from a LS-T table
#' @description
#' From any data frame of the LS-T pathways in the OrganChat object, construct a list of data frames recording the directed & weighted network.
#' The direction indicates the activation along the LS-R-SE-T path, the weight is predicted based on the Hill function model.
#'
#' @param df a data frame of the LS-T pathways in the OrganChat object
#'
#' @import dplyr
#'
#' @return a list
#'
#' @export
#'
constrcuct_EdgeNode_df <- function(df){

  # find the names corresponding to the Metabolite ID
  df$LS <- df$LS_name

  # check if there are duplicated genes
  R.list = unique(df$Receptor)
  TF.list = unique(df$TF)
  T.list =unique(df$Target)

  All.list = unique( c(R.list, TF.list, T.list) )
  A.matrix = matrix(0, length(All.list), 3)
  rownames(A.matrix) = All.list

  A.matrix[ rownames(A.matrix) %in% R.list, 1] = 1
  A.matrix[ rownames(A.matrix) %in% TF.list, 2] = 1
  A.matrix[ rownames(A.matrix) %in% T.list, 3] = 1

  for (i in 1:length(All.list)) {
    if( sum(A.matrix[i, ])>1 ){
      # check receptor
      if( A.matrix[i, 1]==1 ){
        df$Receptor[df$Receptor == rownames(A.matrix)[i]] = paste0(df$Receptor[df$Receptor == rownames(A.matrix)[i]], "*")
      }
      # check TF
      if( A.matrix[i, 2]==1 ){
        df$TF[df$TF == rownames(A.matrix)[i]] = paste0(df$TF[df$TF == rownames(A.matrix)[i]], "**")
      }
      # check target
      if( A.matrix[i, 3]==1 ){
        df$Target[df$Target == rownames(A.matrix)[i]] = paste0(df$Target[df$Target == rownames(A.matrix)[i]], "***")
      }
    }
  }

  # construct node data frame
  df.node1 = df[ , c("LS", "Sender.group", "LS.value")]
  df.node1$nodeID = df.node1$LS
  df.node1$role = "LS"
  df.node1$size = df.node1$LS.value
  df.node1$group = df.node1$Sender.group
  df.node1 = df.node1[ , c("nodeID", "role", "size", "group")]

  df.node2 = df[ , c("Receptor", "Receiver.group", "Receptor.value")]
  colnames(df.node2) = c("nodeID", "group", "size")
  df.node2$role = "Receptor"
  df.node2 = df.node2[ , c("nodeID", "role", "size", "group")]

  df.node3 = df[ , c("TF", "Receiver.group", "TF.value")]
  colnames(df.node3) = c("nodeID", "group", "size")
  df.node3$role = "TF"
  df.node3 = df.node3[ , c("nodeID", "role", "size", "group")]

  df.node4 = df[ , c("Target", "Receiver.group", "Target.value")]
  colnames(df.node4) = c("nodeID", "group", "size")
  df.node4$role = "Target"
  df.node4 = df.node4[ , c("nodeID", "role", "size", "group")]

  df.node = rbind(df.node1, df.node2, df.node3, df.node4)
  df.node = df.node[!duplicated(df.node), ]
  df.node$node = c(0:(nrow(df.node)-1))

  # construct edge data frame
  df.edge1 = df[ , c("LS", "Receptor", "chatP_layer1")]
  colnames(df.edge1) = c("from", "to", "value")
  df.edge2 = df[ , c("Receptor", "TF", "chatP_layer2")]
  colnames(df.edge2) = c("from", "to", "value")
  df.edge3 = df[ , c("TF", "Target", "chatP_layer3")]
  colnames(df.edge3) = c("from", "to", "value")

  df.edge = rbind(df.edge1, df.edge2, df.edge3)
  df.edge = df.edge[!duplicated(df.edge), ]

  df.edge = left_join(df.edge, df.node[, c("nodeID", "node")], join_by(from == nodeID), relationship = "many-to-many")
  colnames(df.edge)[4] = "source"
  df.edge = left_join(df.edge, df.node[, c("nodeID", "node")], join_by(to == nodeID), relationship = "many-to-many")
  colnames(df.edge)[5] = "target"
  df.edge$total.counts = nrow(df)
  df.edge$ave.sd = sum(df$bootstrap_sd)/nrow(df)

  output <- list()
  # save the results
  output[[1]] = df.edge
  output[[2]] = df.node
  names(output) = c("edge", "node")

  return(output)
}


#' @title Construct a list for the directed & weighted GRN with user-selected filters.
#' @description
#' This function is basically constrcuct_EdgeNode_df function with additional options, generating a list of outputs from constrcuct_EdgeNode_df function.
#'
#' @param df a data frame of the LS-T pathways in the OrganChat object
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param restrict.by use this parameter to indicate whether each GRN has a single LS, a single Target, or both. The input should be one of the following:
#' "LS", "Target", "Both". The default value is "Both"
#'
#' @import dplyr
#'
#' @return a list
#'
#' @export
#'
construct_Network_df <- function(df,
                                 category = "All",
                                 restrict.by = "Both"){

  if(nrow(df)==0){ stop('No inferred pathway found!') }

  if( is.null(restrict.by)){
    restrict.by = "Both"
  }else if( !restrict.by %in% c("LS", "Target", "Both")){
    stop("Input 'restrict.by' should be 'LS', 'Target', or 'Both'.")
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    pref = "All"
  }

  if( category != "All"){
    ind = which(grepl(category, df$category, fixed=TRUE))
    df = df[ind, ]
  }

  df.output <- list()
  MT.names <- c()

  M.list <- unique(df$LS)
  T.list <- unique(df$Target)

  if( restrict.by == "Both" ){
    k = 1
    for (i in 1:length(M.list)) {
      for (j in 1:length(T.list)) {

        MT.network <- list()

        df.MT = df[df$LS == M.list[i] & df$Target == T.list[j], ]
        if(nrow(df.MT)==0){
          # df.output[[k]] = data.frame()
          # MT.names = c( MT.names, paste0(M.list[i], "->", T.list[j]) )

          # k = k+1
          next
        }

        MT.network <- constrcuct_EdgeNode_df(df = df.MT)
        df.output[[k]] = MT.network
        MT.names = c( MT.names, paste0(M.list[i], "->", T.list[j]) )
        k = k+1
      }
    }
  }else if( restrict.by == "LS"){
    k = 1
    for (i in 1:length(M.list)){
      MT.network <- list()

      df.MT = df[df$LS == M.list[i], ]
      if(nrow(df.MT)==0){
        # df.output[[k]] = data.frame()
        # MT.names = c( MT.names, paste0(M.list[i], "->All") )

        # k = k+1
        next
      }

      MT.network <- constrcuct_EdgeNode_df(df = df.MT)
      df.output[[k]] = MT.network
      MT.names = c( MT.names, paste0(M.list[i], "->All") )
      k = k+1
    }
  }else if( restrict.by == "Target" ){
    k = 1
    for (j in 1:length(T.list)){
      MT.network <- list()

      df.MT = df[df$Target == T.list[j], ]
      if(nrow(df.MT)==0){
        # df.output[[k]] = data.frame()
        # MT.names = c( MT.names, paste0("All->", T.list[j]) )

        # k = k+1
        next
      }

      MT.network <- constrcuct_EdgeNode_df(df = df.MT)
      df.output[[k]] = MT.network
      MT.names = c( MT.names, paste0("All->", T.list[j]) )
      k = k+1
    }
  }

  names(df.output) <- MT.names
  return(df.output)
}


#' @title Calculate the signal inflow & outflow of a GRN
#' @description
#' For a cluster-cluster communication table of one direction, calculate the following three metrics for each "LS & T" pair: (1) total counts -- total number of pathways involved;
#' (2) ave. sd -- average bootstrap sd of these pathways; (3) max flow --  maximum flow from the LS to the T (using the graph_from_data_frame function).
#'
#' @param df a data frame of the LS-T pathways in the OrganChat object
#'
#' @import dplyr
#' @importFrom igraph graph_from_data_frame
#'
#' @return a data frame
#'
#' @export
#'
construct_MTmaxflow_df <- function(df){
  # construct the network data frame
  df.list <- construct_Network_df(df, category = "All", restrict.by = "Both")

  names <- c()
  total.counts <- c()
  ave.sd <- c()
  max.flow <- c()

  for (i in 1:length(df.list)) {
    # x <- str_split( names(df.list[i]), fixed("->"))
    x.from <- df.list[[i]]$node$nodeID[df.list[[i]]$node$role == "LS"]
    x.to <- df.list[[i]]$node$nodeID[df.list[[i]]$node$role == "Target"]

    names <- c(names, names(df.list[i]))

    if(length(df.list[[i]])==0){
      total.counts <- c(total.counts, 0)
      ave.sd <- c(ave.sd, 0)
      max.flow <- c(max.flow, 0)
      next
    }

    total.counts <- c(total.counts, df.list[[i]]$edge$total.counts[1])
    ave.sd <- c(ave.sd, df.list[[i]]$edge$ave.sd[1])

    d = df.list[[i]]$edge[ , c("from", "to", "value")]
    g = graph_from_data_frame(d, directed = TRUE, vertices = NULL)
    max.flow <- c(max.flow, max_flow(g, x.from, x.to, E(g)$value)$value)

  }

  df.output = data.frame(ST.pair = names, total.counts = total.counts, max.flow = max.flow, ave.sd = ave.sd)

  return(df.output)
}


#' @title Plot the network structure
#' @description
#' Plot the network structure (using igraph package), the input should be one GRN of the output from 'construct_Network_df' function.
#'
#' @param df a data frame of the LS-T pathways in the OrganChat object.
#' @param edge.width.scalemax scale the edge width to a selected maximum value, the default value is 2.
#' @param node.size.scalemax scale the node size to a selected maximum value, the default value is 20.
#' @param edge.color argument of the plot function, the default value is "grey50".
#' @param edge.arrow.size argument of the plot function, the default value is 0.5.
#' @param edge.arrow.width argument of the plot function, the default value is 0.5.
#' @param vertex.frame.color argument of the plot function, the default value is "white".
#' @param vertex.label.color argument of the plot function, the default value is "black".
#' @param vertex.label.cex argument of the plot function, the default value is 0.65.
#' @param vertex.label.dist argument of the plot function, the default value is 0.
#' @param vertex.label.degree argument of the plot function, the default value is 0.
#'
#' @importFrom graphics legend par strwidth title
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @importFrom igraph graph_from_data_frame V E layout_as_tree
#'
#' @export
#'
network_treeplot <- function(df,
                             edge.width.scalemax = 2,
                             node.size.scalemax = 20,
                             edge.color = "grey50",
                             edge.arrow.size=0.5,
                             edge.arrow.width=0.5,
                             vertex.frame.color = "white",
                             vertex.label.color= "black",
                             vertex.label.cex=0.65,
                             vertex.label.dist=0,
                             vertex.label.degree=0){

  links <- data.frame( source=df[[1]]$edge$from,
                       target=df[[1]]$edge$to,
                       importance=df[[1]]$edge$value/max(df[[1]]$edge$value)*edge.width.scalemax)

  nodes <- data.frame( name=df[[1]]$node$nodeID,
                       carac=df[[1]]$node$role,
                       size=df[[1]]$node$size/max(df[[1]]$node$size)*node.size.scalemax)

  net.plot <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
  net.plot$layout <- igraph::layout_as_tree
  igraph::V(net.plot)$carac <- factor(igraph::V(net.plot)$carac, levels = c("LS", "Receptor", "TF", "Target"))

  # Make a palette of 4 colors
  coul  <- brewer.pal(4, "Set1")

  # Create a vector of color
  my_color <- coul[as.numeric(as.factor(igraph::V(net.plot)$carac))]

  plot(net.plot, vertex.color = my_color,
       edge.width = igraph::E(net.plot)$importance,
       edge.color = edge.color,
       edge.arrow.size = edge.arrow.size,                           # Arrow size, defaults to 1
       edge.arrow.width = edge.arrow.width,
       vertex.size = igraph::V(net.plot)$size,
       vertex.frame.color = vertex.frame.color,
       vertex.label.color = vertex.label.color,
       vertex.label.cex = vertex.label.cex,
       vertex.label.dist = vertex.label.dist,                          # Distance between the label and the vertex
       vertex.label.degree = vertex.label.degree) # A simple plot of the network - we'll talk more about plots later

  # Add a legend
  legend("topright", legend=levels(igraph::V(net.plot)$carac)  , col = coul , bty = "n", pch=20 , pt.cex = 1.5, cex = 0.75, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
}


#' @title Rose diagram to show the signaling inflow/outflow
#' @description
#' Rose diagram to show the components of a cluster-cluster communication with the most inflow/outflow (using ggplot2).
#'
#' @param df a data frame of the LS-T pathways in the OrganChat object.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param n the number of components to show in each layer (LS, R, SE, T), the default value is 3.
#' @param outflow.colors argument of the scale_fill_manual function, defining the colors for LS, R, and SE. The default value is c("#e41a1c", "#377eb8", "#4daf4a").
#' @param inflow.colors argument of the scale_fill_manual function, defining the colors for R, SE, and T. The default value is c("#377eb8", "#4daf4a", "#984ea3").
#' @param start.pos argument of the coord_radial function, the default value is 0.25*pi.
#' @param rotate.angle argument of the coord_radial function, the default value is TRUE.
#' @param expand argument of the coord_radial function, the default value is FALSE.
#' @param inner.radius argument of the coord_radial function, the default value is 0.1.
#' @param r.axis.inside argument of the coord_radial function, the default value is TRUE.
#' @param legend.position argument of the theme function, the default value is "left".
#' @param grid.width argument of the grid.arrange function, the default value is 0.6.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom forcats fct_drop
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
components_rosediagram <- function(df,
                                   category = "All",
                                   n = 3,
                                   outflow.colors = c("#e41a1c", "#377eb8", "#4daf4a"),
                                   inflow.colors = c("#377eb8", "#4daf4a", "#984ea3"),
                                   start.pos = 0.25*pi,
                                   rotate.angle = TRUE,
                                   expand = FALSE,
                                   inner.radius = 0.1,
                                   r.axis.inside = TRUE,
                                   legend.position="left",
                                   grid.width = 0.6){

  if( is.null(df) | nrow(df)==0){
    stop("The input data frame cannot be empty!")
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    pref = "All"
  }

  if( category != "All"){
    ind = which(grepl(category, df$category, fixed=TRUE))
    df = df[ind, ]
  }

  # find the names corresponding to the LS ID
  df$LS <- df$LS_name

  Role = c("LS", "Receptor", "TF", "Target")

  outflow_df = data.frame(Name = as.character(), chatP_sum = as.numeric(), Role = as.character())
  inflow_df = data.frame(Name = as.character(), chatP_sum = as.numeric(), Role = as.character())

  for (i in 1:3) {

    df_layer = df[, c(Role[i], Role[i+1], paste0("chatP_layer", i))]
    df_layer = df_layer[!duplicated(df_layer, by=c(Role[i], Role[i+1])), ]

    df.outflow = df_layer[, c(Role[i], paste0("chatP_layer", i))]
    colnames(df.outflow) = c("Name", "chatP_layer")
    df.inflow = df_layer[, c(Role[i+1], paste0("chatP_layer", i))]
    colnames(df.inflow) = c("Name", "chatP_layer")

    df.outflow = aggregate(chatP_layer ~ Name, df.outflow, sum)
    names(df.outflow) = c("Name", "chatP_sum")
    df.inflow = aggregate(chatP_layer ~ Name, df.inflow, sum)
    names(df.inflow) = c("Name", "chatP_sum")

    df.outflow <- df.outflow[ , c("Name", "chatP_sum")]
    df.outflow <- df.outflow[!duplicated(df.outflow), ]
    df.outflow$Role = Role[i]
    df.outflow$chatP_sum = df.outflow$chatP_sum/max(df.outflow$chatP_sum)

    df.inflow <- df.inflow[ , c("Name", "chatP_sum")]
    df.inflow <- df.inflow[!duplicated(df.inflow), ]
    df.inflow$Role = Role[i+1]
    df.inflow$chatP_sum = df.inflow$chatP_sum/max(df.inflow$chatP_sum)

    outflow_df <- rbind(outflow_df, df.outflow)
    inflow_df <- rbind(inflow_df, df.inflow)

  }

  outflow_df$Role = factor(outflow_df$Role, levels = c("LS", "Receptor", "TF", "Target"))
  inflow_df$Role = factor(inflow_df$Role, levels = c("LS", "Receptor", "TF", "Target"))

  outflow_df <- outflow_df %>%
    group_by(Role) %>%
    dplyr::arrange(desc(chatP_sum)) %>%
    slice_head(n = n)

  inflow_df <- inflow_df %>%
    group_by(Role) %>%
    dplyr::arrange(desc(chatP_sum)) %>%
    slice_head(n = n)

  outflow_df$Name = factor(outflow_df$Name, levels = unique(outflow_df$Name) )
  inflow_df$Name = factor(inflow_df$Name, levels = unique(inflow_df$Name) )

  p1 <- ggplot(outflow_df, aes(x=Name, y=chatP_sum, fill=Role) ) + geom_bar(stat='identity') +
    scale_fill_manual(values = outflow.colors) +
    coord_radial(start = start.pos, rotate.angle = rotate.angle,
                 expand = expand, inner.radius = inner.radius, r.axis.inside = r.axis.inside) +
    labs(title = "Outflow of the components",
         y = "Flow",
         fill = "Role") +
    theme(legend.position=legend.position)

  p2 <- ggplot(inflow_df, aes(x=Name, y=chatP_sum, fill=Role) ) + geom_bar(stat='identity') +
    scale_fill_manual(values = inflow.colors) +
    coord_radial(start = start.pos, rotate.angle = rotate.angle,
                 expand = expand, inner.radius = inner.radius, r.axis.inside = r.axis.inside) +
    labs(title = "Inflow of the components",
         y = "Flow",
         fill = "Role") +
    theme(legend.position=legend.position)

  grid.arrange(p1, p2, nrow = 1, widths=c(grid.width, grid.width))
}


#' @title Scatter plot of the top associated GO terms
#' @description
#' Scatter plot via ggplot2 to show the top-related GO terms based on the OC_GOanalysis output.
#'
#' @param df a data frame of the OC_GOanalysis output.
#' @param n total number of GO terms to show, the default value is 10.
#' @param size size of the dots, the default value is 6.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
GOanalysis_scatterplot <- function(df,
                                   n = 10,
                                   size=6){

  df.plot = df
  df.plot %>%
    slice_head(n = n) -> df.plot

  df.plot$significance = -log10(df.plot$p_value)
  df.plot$term_name = factor(df.plot$term_name, levels = df.plot$term_name)

  ggplot(df.plot, aes(x=df.plot$significance, y=df.plot$precision, color = df.plot$term_name)) +
    geom_point(size=size) +
    theme_bw()
}



