#' @title Plot the differential index (DI) of LS-T pathways
#' @description
#' Bar plot of the differential index (DI) of LS-T pathways in each cluster-cluster pair. "DI>0: up-regulated in the first condition."
#' @param df a data frame of differentially expressed LS-T pathways "object - DEpath".
#' @param n number of pathways to show in each cluster-cluster pair, the default value is 5.
#' @param axis.text.size argument of the theme function, the default value is 20.
#'
#' @import dplyr
#' @import ggplot2
#' @import stringr
#'
#' @export
#'
CAOC_DI_barplot <- function(df,
                            n = 5,
                            axis.text.size = 20){

  if(is.null(df)|nrow(df)==0){stop("The input table cannot be empty!")}
  # switch the metabolite ID by its name
  df$pair = paste0(df$Sender.group, " -> ", df$Receiver.group)
  df$pair = factor(df$pair)

  df1 = df[df$DI <=0, ]
  df1 %>%
    group_by(pair) %>%
    arrange( DI ) %>%
    slice_head(n = n) %>%
    ungroup() -> df1

  df2 = df[df$DI >=0, ]
  df2 %>%
    group_by(pair) %>%
    arrange( desc(DI) ) %>%
    slice_head(n = n) %>%
    ungroup() -> df2

  df <- rbind(df1, df2)
  df <- df[!duplicated(df), ]

  for (i in 1:nrow(df)) {
    if(df[i, "id"]=="HMDB"){
      str_sub(df[i , "Path"], 1, 11) <- df$LS_name[i]
    }
  }
  df$ID = paste0(df$Path, "_", df$Sender.group, "_", df$Receiver.group)

  df %>%
    arrange(DI) %>%
    group_by(pair) -> df

  df$ID = factor(df$ID, levels = unique(df$ID) )

  # make.title = paste0("DI>0: up-regulated in the first condition.")

  ggplot(df, aes(x=ID, y=DI, fill=pair)) +
    geom_bar(stat="identity", colour="white") +
    guides(fill=guide_legend(reverse=TRUE)) +
    scale_x_discrete(label = df$Path) +
    labs(y = "Differential index (DI)",
         x = "Pathway") +
    theme(axis.text=element_text(size=axis.text.size)) +
    theme_classic() +
    coord_flip()
}


#' @title Chord diagrams of the cluster-level communications of two conditions
#' @description
#' Chord diagrams of the cluster-level communications of two conditions.
#' @param object an CAOC object (OrganChat object of comparison analysis).
#' @param dir specify which direction to select. The default value is NULL.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param style style select which metric to showcase using the edge width, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", "ave.sd", or NULL (default, meaning "total.counts" will be used).
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
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
#' @import dplyr
#' @import ggplot2
#' @import stringr
#'
#' @export
#'
CAOC_chorddiagram <- function(object,
                              dir = NULL,
                              category = "All",
                              style = NULL,
                              chatP_cutoff = NULL,
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

  OC.list = list()
  net.list = list()
  section.name = c()
  df.plot.list = list()
  for(i in 1:2){
    # create OC objects for each condition
    OC.list[[i]] <- create_OrganChat(object = object@data.list,
                                     organs = object@organs,
                                     idents.list = object@idents.list,
                                     ComparisonAnalysis = FALSE,
                                     assay = NULL,
                                     do.sparse = T)
    # merge LS data
    OC.list[[i]]@LSdata[[1]] = object@LSdata[[i]]
    # merge scRNA data
    OC.list[[i]]@DGdata[[1]] = object@DGdata[[i]]
    # merge ident info
    # OC.list[[i]]@idents = object@idents.bycondition[[i]]
    # merge organ info
    # OC.list[[i]]@organs.metabolite = object@organs.metabolite
    # OC.list[[i]]@organs.rna = object@organs.rna
    # merge net
    OC.list[[i]]@net[[1]] = object@net[[i]]
    # data frame for the plot
    net.list[[i]] = OC.list[[i]]@net[[1]][[dir]]

    if( length(net.list[[i]])>0 ){ # the results for this direction is not empty

      df.output <- OC_singledir_netanalysis(object = OC.list[[i]],
                                            net = net.list[[i]],
                                            dir = dir,
                                            chatP_cutoff = chatP_cutoff,
                                            category = category)
    }else{
      stop(paste0("No signaling pathway found for cindition ", object@conditions[i]))
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
    }else{
      df.plot.list[[i]] = df.plot
    }

    section.name <- unique( c(section.name, df.plot$from, df.plot$to) )
  }

  # set the color code
  color.list = rand_color(length(section.name), transparency = 0.4)
  names(color.list) = section.name

  # plot the chord diagram
  for(j in 1:length(object@conditions)){

    df.plot = df.plot.list[[j]]
    section.plot = unique(c(df.plot$from, df.plot$to))

    if(is.null(grid.col)){
      grid.col = color.list[section.plot]
    }

    circos.clear()
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

    # title(paste0("From top (", unique(df.output$sender.organ), ") to bottom (", unique(df.output$receiver.organ), ")") ,
    #       sub = paste0("Condition: ", object@conditions[j]),
    #       font = title.font, cex = title.cex, line = title.line, adj = title.adj)

    circos.clear()
  }
}


#' @title Bubble-heatmap showing two metrics of the cluster-level OOC of two conditions
#' @description
#' Bubble-heatmap showing two metrics by bubble size and color of the cluster-level OOC of two conditions (using ggplot function).
#'
#' @param object an OrganChat object
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param size select which metric to showcase using the bubble size, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", or "ave.sd". The default value is "ave.weight".
#' @param color select which metric to showcase using the bubble color, the input should be one of the following:
#' "ave.weight", "total.counts", "total.chatP", or "ave.sd". The default value is "total.chatP".
#' @param ncol argument of the wrap_plots function, the default value is 2, meaning placing two plots in 2 columns, respectively.
#'
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @import dplyr
#'
#' @return a ggplot2 object
#'
#' @export
#'
CAOC_bubbleheatmap <- function(object,
                               dir = NULL,
                               chatP_cutoff = NULL,
                               category = "All",
                               size = "ave.weight",
                               color = "total.chatP",
                               ncol = 2){

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

  OC.list = list()
  net.list = list()
  p.list = list()
  df.plot.list = c()
  for(i in 1:length(object@conditions)){
    # create OC objects for each condition
    OC.list[[i]] <- create_OrganChat(object = object@data.list,
                                     organs = object@organs,
                                     idents.list = object@idents.list,
                                     ComparisonAnalysis = FALSE,
                                     assay = NULL,
                                     do.sparse = T)
    # merge LS data
    OC.list[[i]]@LSdata[[1]] = object@LSdata[[i]]
    # merge scRNA data
    OC.list[[i]]@DGdata[[1]] = object@DGdata[[i]]
    # merge ident info
    # OC.list[[i]]@idents = object@idents.bycondition[[i]]
    # merge organ info
    # OC.list[[i]]@organs.metabolite = object@organs.metabolite
    # OC.list[[i]]@organs.rna = object@organs.rna
    # merge net
    OC.list[[i]]@net[[1]] = object@net[[i]]
    # data frame for the plot
    net.list[[i]] = OC.list[[i]]@net[[1]][[dir]]

    if( length(net.list[[i]])>0 ){ # the results for this direction is not empty

      df.output <- OC_singledir_netanalysis(object = OC.list[[i]],
                                            net = net.list[[i]],
                                            dir = dir,
                                            chatP_cutoff = chatP_cutoff,
                                            category = category)
    }else{
      stop(paste0("No signaling pathway found for cindition ", object@conditions[i]))
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

    df.plot.list[[i]] = df.plot

    p.list[[i]] <- ggplot(df.plot, aes(to, from, size = size)) +
      geom_point(aes(color = color)) +
      xlab("Receiver group") + ylab("Sender group") +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_color_viridis_c(option = "C", name = color.name) +
      guides(size = guide_legend(title = size.name)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  }

  wrap_plots(p.list, ncol = ncol)
  # grid.arrange(grobs = p.list, nrow = 1, widths=c( rep(grid.width, length(object@conditions)) ) )

  # names(df.plot.list) = object@conditions
  # return(df.plot.list)
}




