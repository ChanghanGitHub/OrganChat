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








