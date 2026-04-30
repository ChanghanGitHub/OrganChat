#' @title Calculate the mean expression level from a given expression matrix
#' @description
#' Calculate the mean expression level in each cluster based on the normalized expression matrix (gene-by-cell), for one condition.
#' The output is a data frame (gene-by-cluster).
#'
#' @param data normalized expression matrix (gene-by-cell).
#' @param celluse a char contains cell barcodes to be used to select cells. It is NULL by default, meaning all the cells will be used.
#' @param geneuse a char contains gene symbols to be used to select genes. It is NULL by default, meaning all the genes will be calculated.
#' @param idents a char of cluster labels to group the cells.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#'
#' @importFrom forcats fct_drop
#' @importFrom stats quantile
#'
#' @return a data frame
#' @export
cal_meanExpr <- function(data,
                         celluse = NULL,
                         geneuse = NULL,
                         idents = NULL,
                         mean_method = NULL){

  # expression data of organ1 (row: cell; column: gene; the last column: cell group label)
  M = as.matrix(data[ geneuse, celluse])

  if(length(geneuse)==1){
    df1 <- as.data.frame(M)
    colnames(df1) = geneuse
  }else if(length(geneuse)>1){
    df1 <- as.data.frame(t(M))
  }

  df1$group_label <- fct_drop(idents[names(idents) %in% celluse])
  df1$group_label <- as.character(df1$group_label)

  q = c(.25, .5, .75)

  if(is.null(mean_method)){

    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    q2 <- setDT(df1)[ , lapply(.SD, quantile,q[2]), keyby = group_label]
    q3 <- setDT(df1)[ , lapply(.SD, quantile,q[3]), keyby = group_label]

    yy <- 0.25*q1[,geneuse,with=FALSE]+0.5*q2[,geneuse,with=FALSE]+0.25*q3[,geneuse,with=FALSE]
    # rownames(yy) <- q1$group_label
    meanExpr.df = as.data.frame(t(yy))
  } else if(mean_method=="mean"){
    q1 <- setDT(df1)[ , lapply(.SD, quantile,q[1]), keyby = group_label]
    yy <- setDT(df1)[, lapply(.SD, mean), keyby = group_label]
    # rownames(yy) <- q1$group_label
    meanExpr.df = as.data.frame(t(yy[, -1]))
  }

  # meanExpr.df = as.data.frame(t(yy[, -1]))
  colnames(meanExpr.df) <- q1$group_label

  return(meanExpr.df)
}


#' @title Calculate the mean expression level from single-cell metabolite data
#' @description
#' If single-cell metabolite data is available, then this function can be applied to calculate the
#' average metabolite amount in each cell cluster. The input must be the normalized amount matrix, for one organ under one condition.
#' It should be applied to each organ if the data is available for multiple organs.
#'
#' @param object an OrganChat object
#' @param scMeta.data normalized metabolite amount matrix (metabolite-by-cell).
#' @param organ.name a char contains the organ names to be calculated.
#' @param idents a char of cluster labels to group the cells.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#' @param scale.max the maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#' @param type for metabolite data, OrganChat accepts two types of measurements, "amount" and "flux". The default type is "amount",
#' which means the data should be non-negative values. The "flux" data may have both negative (uptake) and positive (release) values.
#' @param single.condition the input is TRUE by default.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#' @export
cal_scMeta.bygroup <- function(object,
                               scMeta.data,
                               organ.name = NULL,
                               idents = NULL,
                               mean_method = NULL,
                               scale.max = 10,
                               type = "amount",
                               single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check the input 'type'
  if(is.null(type)){ type = "amount" }
  if( !type %in% c("amount", "flux") ){
    stop("The input 'type' must be either 'amount' or 'flux'.")
  }

  if ( inherits(x = scMeta.data, what = c("matrix", "Matrix", "dgCMatrix")) ) {
    data <- scMeta.data
  }

  if( is.null(organ.name) ){
    stop("Must specify the organ name!")
  }else if( !(organ.name %in% object@organs) ){
    stop("The organ name must be one of the two selected organs.")
  }

  idents = factor(idents)
  idents <- fct_drop(idents)

  x = setdiff( levels(idents), levels(object@idents.list[[organ.name]]) )
  if( length(x)>0 ){
    stop("The cell clusters in the scMETA.data and the object@data.list$organ.name must be the same.")
  }else{
    idents = factor(idents, levels = levels(object@idents.list[[organ.name]]) )
  }

  if (!identical(names(idents), colnames(scMeta.data))) {
    warning("The cell barcodes in 'idents' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the names of 'idents'!")
    names(idents) <- colnames(scMeta.data)
  }

  ##############################################################################

  meanExpr.df = cal_meanExpr(data = data,
                             celluse = names(idents),
                             geneuse = rownames(data),
                             idents = idents,
                             mean_method = mean_method)

  if(!is.null(scale.max)){
    meanExpr.df = meanExpr.df/max(meanExpr.df)*scale.max
  }

  data.new = meanExpr.df
  data.new$Type = type

  if(nrow(object@LSdata[[con.layer]][[organ.name]])==0){
    object@LSdata[[con.layer]][[organ.name]] = data.new
  }else if( length( setdiff(colnames(object@LSdata[[con.layer]][[i]]), colnames(data.new)) )==0 ){
    data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[con.layer]][[i]])) , colnames(object@LSdata[[con.layer]][[i]])]
    object@LSdata[[con.layer]][[i]] = rbind(object@LSdata[[con.layer]][[i]], data.new)
  }else{
    warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
  }

  return(object)
}


#' @title Calculate the mean expression level from single-cell transcriptomics data
#' @description
#' Single-cell transcriptomics data of two organs are requested to perform OOC. This function can be applied to
#' calculate the average gene expression of the long-range signals (LS) and downstream genes (DG) in each cell cluster, respectively (for both organs under one condition).
#' The LS and DG will be inferred using the corresponding OrganChatDB.
#' The outputs will be stored in "Object - LSdata/DGdata - organ".
#'
#' @param object an OrganChat object
#' @param DB the corresponding OrganChatDB used in the analysis.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#' @param scale.max the maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#' @param single.condition the input is TRUE by default.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#' @export
cal_scRNA.bygroup <- function(object,
                              DB = OrganChatDB,
                              mean_method = NULL,
                              scale.max = 10,
                              single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # get the LS names from the DB
  LS.list = unique(DB[[1]]$from)

  for (i in 1:2) {

    organ.select = object@organs[i]
    idents.select = object@idents.list[[i]]
    data.select = object@data.list[[i]]

    # find the cell barcodes of organ1 and organ2
    celluse <- names(idents.select)
    # find all the genes that will be used
    geneuse <- rownames(data.select)

    meanExpr.df = cal_meanExpr(data = data.select,
                               celluse = celluse,
                               geneuse = geneuse,
                               idents = idents.select,
                               mean_method = mean_method)

    object@DGdata[[1]][[organ.select]] = meanExpr.df

    if(!is.null(scale.max)){
      meanExpr.df = meanExpr.df/max(meanExpr.df)*scale.max
    }

    # find the matched LS (it may be empty)
    meanExpr.LS = meanExpr.df[rownames(meanExpr.df) %in% LS.list, ]

    if( nrow(meanExpr.LS)>0 ){

      # merge
      data.new = meanExpr.LS
      data.new$Type = "amount"

      if(nrow(object@LSdata[[con.layer]][[organ.select]])==0){
        object@LSdata[[con.layer]][[organ.select]] = data.new
      }else if( length( setdiff(colnames(object@LSdata[[con.layer]][[i]]), colnames(data.new)) )==0 ){
        data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[con.layer]][[i]])) , colnames(object@LSdata[[con.layer]][[i]])]
        object@LSdata[[con.layer]][[i]] = rbind(object@LSdata[[con.layer]][[i]], data.new)
      }else{
        warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
      }
    }

  }

  return(object)
}


#' @title Merge existing cluster-level metabolite amount/flux data
#' @description
#' If the metabolite amount/flux data for both organs are ready, and are not calculated using th cal_scMeta.bygroup function,
#' they can still be integrated into the OrganChat object using this function.
#' This is usually because when single-cell metabolite data is unavailable, then the metabolite flux data can be inferred via METAFlux using OrganChat built-in modules.
#' Note the metabolites are only considered as the LS, the outputs will be stored in "Object - LSdata - organ".
#'
#' @param object an OrganChat object
#' @param metabolitedata.list a list of metabolite data (average amount/flux for each cluster) for two organs
#' @param scale.max the maximum value used to scale the data, which is 10 by default. No scaling will be performed if the input is NULL.
#' @param type for metabolite data, OrganChat accepts two types of measurements, "amount" and "flux". The default type is "amount",
#' which means the data should be non-negative values. The "flux" data may have both negative (uptake) and positive (release) values.
#' @param single.condition the input is TRUE by default.
#'
#' @return an OrganChat object
#' @export
merge_metabolite.bygroup <- function(object,
                                     metabolitedata.list = NULL,
                                     scale.max = 10,
                                     type = "flux",
                                     single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check the input 'type'
  if(is.null(type)){ type = "amount" }
  if( !type %in% c("amount", "flux") ){
    stop("The input 'type' must be either 'amount' or 'flux'.")
  }

  if( length(metabolitedata.list)!=2 ){
    stop("The input 'metabolitedata.list' should have two data frames representing the flux data of the two organs.")
  }

  # scale the data
  if(!is.null(scale.max)){
    m1 = max( abs(metabolitedata.list[[1]]) )
    m2 = max( abs(metabolitedata.list[[2]]) )
    metabolitedata.list[[1]] = metabolitedata.list[[1]]/max(m1, m2)*scale.max
    metabolitedata.list[[2]] = metabolitedata.list[[2]]/max(m1, m2)*scale.max
  }

  # merge the data
  for (i in 1:2) {

    data.new = metabolitedata.list[[i]]
    data.new$Type = type

    if(nrow(object@LSdata[[con.layer]][[i]])==0){
      object@LSdata[[con.layer]][[i]] = data.new
    }else if( length( setdiff(colnames(object@LSdata[[con.layer]][[i]]), colnames(data.new)) )==0 ){
      data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[con.layer]][[i]])) , colnames(object@LSdata[[con.layer]][[i]])]
      object@LSdata[[con.layer]][[i]] = rbind(object@LSdata[[con.layer]][[i]], data.new)
    }else{
      warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
    }
  }

  return(object)
}


#' @title LS-R-SE-T pathway inference
#' @description
#' Infer LS-R-SE-T pathways using OrganChatDB, pathways with duplicated genes will be removed.
#' The output will be a data frame and stored in "Object - net0".
#'
#' @param object an OrganChat object
#' @param DB the corresponding OrganChatDB used in the analysis.
#' @param LS the list of LS used for inferring the pathways, the input is NULL by default, which means all will be used
#' @param receptor the list of Receptor used for inferring the pathways, the input is NULL by default, which means all will be used
#' @param tf the list of Signaling Effector used for inferring the pathways, the input is NULL by default, which means all will be used
#' @param target  the list of Target used for inferring the pathways, the input is NULL by default, which means all will be used
#' @param single.condition the input is TRUE by default.
#' @param hmdb_dictionary OrganChat-provided database matching "metabolite name" and the HMDB ID.
#'
#' @importFrom forcats fct_drop
#'
#' @return a data frame
#' @export
OCnet0_infer <- function(object,
                         DB = OrganChatDB,
                         LS = NULL,
                         receptor = NULL,
                         tf = NULL,
                         target = NULL,
                         single.condition = TRUE,
                         hmdb_dictionary = hmdb_dictionary){

  # "FALSE" is designed mainly for the CAOC
  if( isTRUE(single.condition) ){
    con.layer = 1
    LS.list <- intersect( rownames(object@LSdata[[con.layer]][[1]]), rownames(object@LSdata[[con.layer]][[2]]) )
  }else{
    LS.list <- intersect( rownames(object@LSdata[[1]][[1]]),
                          rownames(object@LSdata[[1]][[2]]),
                          rownames(object@LSdata[[2]][[1]]),
                          rownames(object@LSdata[[2]][[2]]),)
  }


  # get the DS list from the data
  if( isTRUE(single.condition) ){
    con.layer = 1
    DG.list <- intersect( rownames(object@DGdata[[con.layer]][[1]]), rownames(object@DGdata[[con.layer]][[2]]) )
  }else{
    DG.list <- intersect( rownames(object@DGdata[[1]][[1]]),
                          rownames(object@DGdata[[1]][[2]]),
                          rownames(object@DGdata[[2]][[1]]),
                          rownames(object@DGdata[[2]][[2]]),)
  }

  if( is.null(LS.list) ){
    stop("No shared long-range signal (LS) was found")
  }

  if( is.null(DG.list) ){
    stop("No down-stream gene (DG) was found")
  }

  name.list <- unique(c(LS.list, DG.list))

  # filter the database
  Role <- c("LS", "Receptor", "TF", "Target")

  DB[[1]] = DB[[1]][ (DB[[1]]$from %in% name.list)&(DB[[1]]$to %in% name.list) , c("from", "to", "category", "id")]
  colnames(DB[[1]])[1:2] <- c(Role[1], Role[2])

  for (i in 2:3) {
    DB[[i]] = DB[[i]][ (DB[[i]]$from %in% name.list)&(DB[[i]]$to %in% name.list) , 1:2]
    colnames(DB[[i]]) <- c(Role[i], Role[i+1])
  }

  SigPath <- inner_join(DB[[1]], DB[[2]], by = join_by(Receptor), multiple = "all", relationship = "many-to-many")
  SigPath <- inner_join(SigPath, DB[[3]], by = join_by(TF), multiple = "all", relationship = "many-to-many")

  if(!is.null(receptor)){
    SigPath <- SigPath[ (SigPath$Receptor %in% receptor) , ]
  }

  if(!is.null(tf)){
    SigPath <- SigPath[ (SigPath$TF %in% tf) , ]
  }

  if(!is.null(target)){
    SigPath <- SigPath[ (SigPath$Target %in% target) , ]
  }

  if(nrow(SigPath)==0){
    stop("No signal network inferred")
  }

  SigPath = SigPath[ , c("LS", "Receptor", "TF", "Target", "category", "id")]

  # remove pathways with duplicated genes
  boolean_dup <- apply(SigPath[ , 1:4], 1, duplicated)
  index = apply(boolean_dup, 2, sum)
  SigPath$index = index
  SigPath = SigPath[SigPath$index == 0, c("LS", "Receptor", "TF", "Target", "category", "id")]
  SigPath = SigPath[!duplicated(SigPath), ]

  SigPath$Path <- paste(SigPath$LS, SigPath$Receptor, SigPath$TF, SigPath$Target, sep="*")
  SigPath <- SigPath[!duplicated(SigPath$Path), ]

  # find the metabolite name
  M.name <- hmdb_dictionary[hmdb_dictionary$ID %in% SigPath$LS, ]
  M.name %>%
    group_by(ID) %>%
    slice_head(n = 1) %>%
    ungroup() -> M.name
  colnames(M.name) = c("LS", "LS_name")

  SigPath <- left_join(SigPath, M.name, join_by(LS))
  for (i in 1:nrow(SigPath)) {
    if(is.na(SigPath$LS_name[i])){
      SigPath$LS_name[i] = SigPath$LS[i]
    }
  }

  object@net0 <- SigPath

  return(object)
}


#' @title LS-R-SE-T pathway filtering (for one direction)
#' @description
#' Filter the 'net0' by removing lowly/zero expressed LS & DG for one direction, one condition.
#' The output will be a list of data frames for cluster-cluster communications.
#'
#' @param object an OrganChat object
#' @param cutoff_value a vector of 4 cutoff values for LS, R, SE, T, respectively.
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param single.condition the input is TRUE by default.
#'
#' @importFrom forcats fct_drop
#'
#' @return a list
#' @export
filter_OCnet0_onedir <- function(object, cutoff_value = NULL, dir = NULL, single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check input 'dir'
  if( dir!="forward" & dir!="backward" & dir!="self_organ1" & dir!="self_organ2"){
    stop("The input 'dir' should be 'forward', 'backward', 'self_organ1', or 'self_organ2'.")
  }

  # check input 'cutoff_value'
  if(is.null(cutoff_value)){ cutoff_value = c(0,0,0,0) }
  if( length(cutoff_value)!=4 ){
    stop("The input 'cutoff_value' should be a 1x4 vector.")
  }

  # save the results
  net.list <- list()

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

  k = 1
  name = c()

  if( dir %in% c("forward", "backward") ){
    for (i in 1:length(S.group)) {
      for (j in 1:length(R.group)) {
        # sender and receiver groups
        sender = S.group[i]
        receiver = R.group[j]

        # LS data for the sender and receiver groups
        LSdata.s = data.frame(LS = rownames(object@LSdata[[con.layer]][[S.organ]]),
                              LS.value_s = object@LSdata[[con.layer]][[S.organ]][ , sender],
                              Type_s =  object@LSdata[[con.layer]][[S.organ]][ , "Type"])

        LSdata.r = data.frame(LS = rownames(object@LSdata[[con.layer]][[R.organ]]),
                              LS.value_r = object@LSdata[[con.layer]][[R.organ]][ , receiver],
                              Type_r =  object@LSdata[[con.layer]][[R.organ]][ , "Type"])

        LSdata.df <- inner_join(LSdata.s, LSdata.r, by = join_by(LS), relationship = "many-to-many")

        # for flux data, sender--release, receiver--uptake
        LSdata.flux <- LSdata.df[ LSdata.df$Type_s == "flux" & LSdata.df$Type_r == "flux" &
                                    LSdata.df$LS.value_s > 0 & LSdata.df$LS.value_r < 0 , ]
        LSdata.amount <- LSdata.df[ LSdata.df$Type_s == "amount" & LSdata.df$Type_r == "amount", ]

        if( nrow(LSdata.flux)==0 & nrow(LSdata.amount)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(sender,"_to_", receiver))
          k = k+1
          next
        }

        if( nrow(LSdata.flux)>0 & nrow(LSdata.amount)>0 ){
          LSdata.df = rbind(LSdata.flux, LSdata.amount)
        }else if( nrow(LSdata.flux)>0 ){
          LSdata.df = LSdata.flux
        }else if( nrow(LSdata.amount)>0 ){
          LSdata.df = LSdata.amount
        }

        LSdata.df$LS.value = LSdata.df$LS.value_s - LSdata.df$LS.value_r
        LSdata.df = LSdata.df[ , c("LS", "LS.value")]

        LSdata.df <- LSdata.df[ LSdata.df$LS.value >= cutoff_value[1] , ] # nrow can be 0

        if( nrow(LSdata.df)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(sender,"_to_", receiver))
          k = k+1
          next
        }

        # expression data for the R, TF, T in the receiver group
        DGdata <- data.frame(DG = rownames(object@DGdata[[con.layer]][[R.organ]]),
                             Value = object@DGdata[[con.layer]][[R.organ]][ , receiver])
        DGdata.r <- DGdata[ DGdata$Value >= cutoff_value[2] , ]
        DGdata.tf <- DGdata[ DGdata$Value >= cutoff_value[3] , ]
        DGdata.t <- DGdata[ DGdata$Value >= cutoff_value[4] , ]

        if( nrow(DGdata.r)==0 | nrow(DGdata.tf)==0 | nrow(DGdata.t)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(sender,"_to_", receiver))
          k = k+1
          next
        }

        colnames(DGdata.r) = c("Receptor", "Receptor.value")
        colnames(DGdata.tf) = c("TF", "TF.value")
        colnames(DGdata.t) = c("Target", "Target.value")

        # filter the net0
        net.filtered <- object@net0
        net.filtered$Sender.group <- sender
        net.filtered$Receiver.group <- receiver
        net.filtered <- net.filtered[ net.filtered$LS %in% LSdata.df$LS &
                                        net.filtered$Receptor %in% DGdata.r$Receptor &
                                        net.filtered$TF %in% DGdata.tf$TF &
                                        net.filtered$Target %in% DGdata.t$Target , ]
        # combine net.filtered and the data
        net.filtered <- left_join(net.filtered, LSdata.df, by = join_by(LS), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.r, by = join_by(Receptor), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.tf, by = join_by(TF), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.t, by = join_by(Target), relationship = "many-to-many")

        # add the data frame to the object
        net.filtered = net.filtered[!duplicated(net.filtered), ]
        net.list[[k]] = net.filtered
        name = c(name, paste0(sender,"_to_", receiver))
        k = k+1
      }
    }
  }else if( dir %in% c("self_organ1", "self_organ2") ){
    for (i in 1:length(S.group)) {
      for (j in 1:length(R.group)) {
        # sender and receiver groups
        sender = S.group[i]
        receiver = R.group[j]

        # LS data for the sender and receiver groups
        LSdata.s = data.frame(LS = rownames(object@LSdata[[con.layer]][[S.organ]]),
                              LS.value_s = object@LSdata[[con.layer]][[S.organ]][ , sender],
                              Type_s =  object@LSdata[[con.layer]][[S.organ]][ , "Type"])
        LSdata.df = LSdata.s
        # for flux data, sender--release, receiver--0
        LSdata.flux <- LSdata.df[ LSdata.df$Type_s == "flux" & LSdata.df$LS.value_s > 0, ]
        LSdata.amount <- LSdata.df[ LSdata.df$Type_s == "amount", ]

        if( nrow(LSdata.flux)==0 & nrow(LSdata.amount)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(sender,"_to_", receiver))
          k = k+1
          next
        }

        if( nrow(LSdata.flux)>0 & nrow(LSdata.amount)>0 ){
          LSdata.df = rbind(LSdata.flux, LSdata.amount)
        }else if( nrow(LSdata.flux)>0 ){
          LSdata.df = LSdata.flux
        }else if( nrow(LSdata.amount)>0 ){
          LSdata.df = LSdata.amount
        }

        LSdata.df$LS.value = LSdata.df$LS.value_s - 0
        LSdata.df = LSdata.df[ , c("LS", "LS.value")]

        LSdata.df <- LSdata.df[ LSdata.df$LS.value >= cutoff_value[1] , ] # nrow can be 0

        if( nrow(LSdata.df)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(S.organ, "_", sender,"_to_", receiver))
          k = k+1
          next
        }

        # expression data for the R, TF, T in the receiver group
        DGdata <- data.frame(DG = rownames(object@DGdata[[con.layer]][[R.organ]]),
                             Value = object@DGdata[[con.layer]][[R.organ]][ , receiver])
        DGdata.r <- DGdata[ DGdata$Value >= cutoff_value[2] , ]
        DGdata.tf <- DGdata[ DGdata$Value >= cutoff_value[3] , ]
        DGdata.t <- DGdata[ DGdata$Value >= cutoff_value[4] , ]

        if( nrow(DGdata.r)==0 | nrow(DGdata.tf)==0 | nrow(DGdata.t)==0 ){
          net.list[[k]] = list()
          name = c(name, paste0(S.organ, "_", sender,"_to_", receiver))
          k = k+1
          next
        }

        colnames(DGdata.r) = c("Receptor", "Receptor.value")
        colnames(DGdata.tf) = c("TF", "TF.value")
        colnames(DGdata.t) = c("Target", "Target.value")

        # filter the net0
        net.filtered <- object@net0
        net.filtered$Sender.group <- sender
        net.filtered$Receiver.group <- receiver
        net.filtered <- net.filtered[ net.filtered$LS %in% LSdata.df$LS &
                                        net.filtered$Receptor %in% DGdata.r$Receptor &
                                        net.filtered$TF %in% DGdata.tf$TF &
                                        net.filtered$Target %in% DGdata.t$Target , ]
        # combine net.filtered and the data
        net.filtered <- left_join(net.filtered, LSdata.df, by = join_by(LS), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.r, by = join_by(Receptor), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.tf, by = join_by(TF), relationship = "many-to-many")
        net.filtered <- left_join(net.filtered, DGdata.t, by = join_by(Target), relationship = "many-to-many")

        # add the data frame to the object
        net.filtered = net.filtered[!duplicated(net.filtered), ]
        net.list[[k]] = net.filtered
        name = c(name, paste0(S.organ, "_", sender,"_to_", receiver))
        k = k+1
      }
    }
  }

  names(net.list) <- name
  return(net.list)
}


#' @title LS-R-SE-T pathway filtering for all directions
#' @description
#' Filter the 'net0' by removing lowly/zero expressed LS & DG for all four directions, one condition.
#' The output will be stored under "object - net".
#'
#' @param object an OrganChat object
#' @param cutoff_value a vector of 4 cutoff values for LS, R, SE, T, respectively.
#' @param single.condition the input is TRUE by default.
#'
#' @return an OrganChat object
#' @export
filter_OCnet0 <- function(object, cutoff_value = NULL, single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  object@net[[con.layer]][[1]] = filter_OCnet0_onedir(object, cutoff_value = cutoff_value,
                                                      dir = "forward", single.condition = single.condition)
  object@net[[con.layer]][[2]] = filter_OCnet0_onedir(object, cutoff_value = cutoff_value,
                                                      dir = "backward", single.condition = single.condition)
  object@net[[con.layer]][[3]] = filter_OCnet0_onedir(object, cutoff_value = cutoff_value,
                                                      dir = "self_organ1", single.condition = single.condition)
  object@net[[con.layer]][[4]] = filter_OCnet0_onedir(object, cutoff_value = cutoff_value,
                                                      dir = "self_organ2", single.condition = single.condition)
  names(object@net[[con.layer]]) <- c("forward", "backward", "self_organ1", "self_organ2")
  return(object)
}


#' @title Remove the receptors activated by non-LS signals
#' @description
#' (Optional step) Remove the receptors (and the associated LS-T pathways) that: (1) they can be activated by the ligands outside
#' the OrganChatDB, based on other ligand-receptor sources such as CellChatDB; (2) these ligands are recorded in the data and highly expressed.
#' The additional ligand-receptor database is need (LRDB) to perform this step.
#'
#' @param object an OrganChat object
#' @param DB the corresponding OrganChatDB used in the analysis.
#' @param LRDB an additional ligand-receptor database, such as CellChatDB, should be a data frame contains ligand-receptor pairs.
#' @param mean_method it is NULL by default and a weighted mean will be calculated (which is statisticallt more robust against noise).
#' If the input is "mean", then the arithmetic mean will be calculated.
#' @param exp_cutoff a scalar used to identify highly expressed ligands.
#' @param single.condition the input is TRUE by default.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#' @export
filter_OCreceptor <- function(object,
                              DB = OrganChatDB,
                              LRDB = LRDB,
                              mean_method = NULL,
                              exp_cutoff = NULL,
                              single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check the inputs
  if(is.null(LRDB)){
    stop("The input 'LRDB' must be a data frame of the ligand-receptor pairs.")
  }

  # for organ1
  data.select = object@data.list[[1]]
  idents.select = object@idents.list[[1]]

  celluse <- colnames(data.select)

  # get the ligand list
  geneuse <- unique(LRDB$from[LRDB$to %in% rownames(data.select)])
  geneuse <- setdiff(geneuse, unique(DB[[1]]$from))
  geneuse <- geneuse[geneuse %in% rownames(data.select)]

  meanExpr.organ1 = cal_meanExpr(data = data.select,
                                 celluse = celluse,
                                 geneuse = geneuse,
                                 idents = idents.select,
                                 mean_method = mean_method)
  meanExpr.organ1$Max = apply(meanExpr.organ1, 1, max)
  meanExpr.organ1$Gene = rownames(meanExpr.organ1)

  # select the acceptable receptors based on the cutoff value
  if(is.null(exp_cutoff)){
    R.organ1 <- rownames(object@data.list[[1]])
  }else{
    L.organ1 <- meanExpr.organ1$Gene[meanExpr.organ1$Max > exp_cutoff]
    R.organ1 <- unique(LRDB$to[ (LRDB$to %in% rownames(data.select)) & (LRDB$from %in% L.organ1) ])
  }

  # for organ 2
  data.select = object@data.list[[2]]
  idents.select = object@idents.list[[2]]

  celluse <- colnames(data.select)

  # get the ligand list
  geneuse <- unique(LRDB$from[LRDB$to %in% rownames(data.select)])
  geneuse <- setdiff(geneuse, unique(DB[[1]]$from))
  geneuse <- geneuse[geneuse %in% rownames(data.select)]

  meanExpr.organ2 = cal_meanExpr(data = data.select,
                                 celluse = celluse,
                                 geneuse = geneuse,
                                 idents = idents.select,
                                 mean_method = mean_method)
  meanExpr.organ2$Max = apply(meanExpr.organ2, 1, max)
  meanExpr.organ2$Gene = rownames(meanExpr.organ2)

  # select the acceptable receptors based on the cutoff value
  if(is.null(exp_cutoff)){
    R.organ2 <- rownames(object@data.list[[1]])
  }else{
    L.organ2 <- meanExpr.organ2$Gene[meanExpr.organ2$Max > exp_cutoff]
    R.organ2 <- unique(LRDB$to[ (LRDB$to %in% rownames(data.select)) & (LRDB$from %in% L.organ2) ])
  }

  ##############################################################################
  if(length(object@net[[con.layer]]$forward)>0){
    for (i in 1:length(object@net[[con.layer]]$forward)) {
      if( !is.null(dim(object@net[[con.layer]]$forward[[i]])) & nrow(object@net[[con.layer]]$forward[[i]])>0 ){
        object@net[[con.layer]]$forward[[i]] <- object@net[[con.layer]]$forward[[i]][object@net[[con.layer]]$forward[[i]]$Receptor %in% R.organ2, ]
      }
    }
  }

  if(length(object@net[[con.layer]]$backward)>0){
    for (i in 1:length(object@net[[con.layer]]$backward)) {
      if( !is.null(dim(object@net[[con.layer]]$backward[[i]])) & nrow(object@net[[con.layer]]$backward[[i]])>0 ){
        object@net[[con.layer]]$backward[[i]] <- object@net[[con.layer]]$backward[[i]][object@net[[con.layer]]$backward[[i]]$Receptor %in% R.organ1, ]
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ1)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ1)) {
      if( !is.null(dim(object@net[[con.layer]]$self_organ1[[i]])) & nrow(object@net[[con.layer]]$self_organ1[[i]])>0 ){
        object@net[[con.layer]]$self_organ1[[i]] <- object@net[[con.layer]]$self_organ1[[i]][object@net[[con.layer]]$self_organ1[[i]]$Receptor %in% R.organ1, ]
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ2)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ2)) {
      if( !is.null(dim(object@net[[con.layer]]$self_organ2[[i]])) & nrow(object@net[[con.layer]]$self_organ2[[i]])>0  ){
        object@net[[con.layer]]$self_organ2[[i]] <- object@net[[con.layer]]$self_organ2[[i]][object@net[[con.layer]]$self_organ2[[i]]$Receptor %in% R.organ2, ]
      }
    }
  }

  return(object)
}


#' @title Calculate the communication probability (individual LS-T pathway)
#' @description
#' Calculate the communication probability of an individual LS-T pathway, the signaling probability between two components are
#' modeled using Hill function x^N/(x^N+K^N).
#'
#' @param object an OrganChat object
#' @param df a data frame of LS-T pathways.
#' @param K half-saturation constant of the Hill function, the default value is 0.5.
#' @param N Hill coefficient, the default value is 2.
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#'
#' @return a data frame
#' @export
cal_singlechat_prob <- function(object, df, K = 0.5, N = 2,
                                dir = NULL){

  # check input 'dir'
  if( dir!="forward" & dir!="backward" & dir!="self_organ1" & dir!="self_organ2"){
    stop("The input 'dir' should be 'forward', 'backward', 'self_organ1', or 'self_organ2'.")
  }

  # define hill function
  hill <- function(x,K,N){
    return(x^N/(x^N+K^N))
  }

  # calculate the LS-receptor signaling probability
  x = df$LS.value
  x[x<=0] = 0

  P1 <- hill(x*df$Receptor.value, K, N)
  P2 <- hill(df$Receptor.value*df$TF.value, K, N)
  P3 <- hill(df$TF.value*df$Target.value, K, N)

  df$chatP_layer1 <- P1
  df$chatP_layer2 <- P2
  df$chatP_layer3 <- P3

  df$chatP <- P1*P2*P3

  return(df)
}


#' @title Calculate the communication probability (all LS-T pathways in the object)
#' @description
#' Calculate the communication probability of all LS-T pathways in the object, the signaling probability between two components are
#' modeled using Hill function x^N/(x^N+K^N).
#'
#' @param object an OrganChat object
#' @param K half-saturation constant of the Hill function, the default value is 0.5.
#' @param N Hill coefficient, the default value is 2.
#' @param single.condition the input is TRUE by default.
#'
#' @return an OrganChat object
#' @export
cal_netchat_prob <- function(object, K = 0.5, N = 2,
                             single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  if(length(object@net[[con.layer]]$forward)>0){
    for (i in 1:length(object@net[[con.layer]]$forward)) {
      if( !is.null(dim(object@net[[con.layer]]$forward[[i]])) & nrow(object@net[[con.layer]]$forward[[i]])>0 ){
        object@net[[con.layer]]$forward[[i]] <- cal_singlechat_prob(object, object@net[[con.layer]]$forward[[i]],K,N,dir="forward")
      }
    }
  }

  if(length(object@net[[con.layer]]$backward)>0){
    for (i in 1:length(object@net[[con.layer]]$backward)) {
      if( !is.null(dim(object@net[[con.layer]]$backward[[i]])) & nrow(object@net[[con.layer]]$backward[[i]])>0 ){
        object@net[[con.layer]]$backward[[i]] <- cal_singlechat_prob(object, object@net[[con.layer]]$backward[[i]],K,N,dir="backward")
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ1)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ1)) {
      if( !is.null(dim(object@net[[con.layer]]$self_organ1[[i]])) & nrow(object@net[[con.layer]]$self_organ1[[i]])>0 ){
        object@net[[con.layer]]$self_organ1[[i]] <- cal_singlechat_prob(object, object@net[[con.layer]]$self_organ1[[i]],K,N,dir="self_organ1")
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ2)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ2)) {
      if( !is.null(dim(object@net[[con.layer]]$self_organ2[[i]])) & nrow(object@net[[con.layer]]$self_organ2[[i]])>0 ){
        object@net[[con.layer]]$self_organ2[[i]] <- cal_singlechat_prob(object, object@net[[con.layer]]$self_organ2[[i]],K,N,dir="self_organ2")
      }
    }
  }

  return(object)
}


#' @title Perform the bootstrap test ((individual LS-T pathway))
#' @description
#' Perform the bootstrap test and calculated the standard deviation for an individual LS-T pathway stored under "object - net".
#'
#' @param object an OrganChat object
#' @param R the number of bootstrap testing to perfomr, the default value is 5.
#' @param dir the direction of comunication to analyze, The input should be one of the following: 'forward', 'backward', 'self_organ1', or 'self_organ2'.
#' @param CI value of confidence interval, CI = 0.95 by default.
#'
#' @importFrom forcats fct_drop
#' @importFrom boot boot
#' @importFrom stats sd quantile
#'
#' @return a list
#' @export
cal_singlechat_bootstrap <- function(object, R = 5, dir = NULL, CI = 0.95){

  # check input 'dir'
  if( dir!="forward" & dir!="backward" & dir!="self_organ1" & dir!="self_organ2"){
    stop("The input 'dir' should be 'forward', 'backward', 'self_organ1', or 'self_organ2'.")
  }

  # check the R value
  if(is.null(R)){
    R = 5
    message("The default value of 'R' (number of bootstrap) is 5.")
  }

  # check the CI value
  if(is.null(CI)){
    CI = 0.95
  }

  geneuse <- unique(c(object@net0$Receptor, object@net0$TF, object@net0$Target))

  # define a function to calculate the mean for each column
  model_mean <- function(data, index){apply(data[index, ], 2, mean)}

  # generate a data.frame to store the results
  bootstrap_sd.df <- data.frame(Gene = geneuse)
  bootstrap_CI.df <- data.frame(Gene = geneuse)

  if(dir=="backward"|dir=="self_organ1"){ # "backward"--organ 1 is the receiver

    # for organ1
    data.select = object@data.list[[1]]
    idents.select = object@idents.list[[1]]
    celluse <- colnames(data.select)

    # expression data of organ1 (row: cell; column: gene; the last column: cell group label)
    M1 = as.matrix(data.select[ geneuse, celluse])
    if(length(geneuse)==1){
      df1 <- as.data.frame(M1)
      colnames(df1) = geneuse
    }else if(length(geneuse)>1){
      df1 <- as.data.frame(t(M1))
    }

    df1$group_label <- fct_drop(idents.select[names(idents.select) %in% celluse])
    df1$group_label <- as.factor(df1$group_label)

    # cell groups in organ 1
    group1 <- levels(df1$group_label)

    for (i in 1:length(group1)) {
      df.bs <- df1[df1$group_label == levels(df1$group_label)[i], -ncol(df1)]
      bs.results <- boot(data = df.bs, statistic = model_mean, R)
      x <- apply(bs.results$t, 2, sd)
      bs.results_df <- data.frame(sd = x)
      colnames(bs.results_df) <- paste0(levels(df1$group_label)[i], "_sd")

      # calculate the CI
      q = c(0+(1-CI)/2 , 1-(1-CI)/2)
      CI_min = apply(bs.results$t, 2, quantile, q[1])
      CI_max = apply(bs.results$t, 2, quantile, q[2])
      bs.results_CI = (bs.results$t0 - CI_min)*(CI_max - bs.results$t0)
      bs.results_CI[ bs.results_CI>0 ] = TRUE
      bs.results_CI[ bs.results_CI<=0 ] = FALSE
      bs.results_CI <- data.frame(CI = bs.results_CI)
      colnames(bs.results_CI) <- paste0(levels(df1$group_label)[i], "_CI")

      bootstrap_sd.df <- cbind(bootstrap_sd.df, bs.results_df)
      bootstrap_CI.df <- cbind(bootstrap_CI.df, bs.results_CI)
    }

  }else if(dir=="forward"|dir=="self_organ2"){ # organ 2 is the receiver

    # for organ2
    data.select = object@data.list[[2]]
    idents.select = object@idents.list[[2]]
    celluse <- colnames(data.select)

    # expression data of organ2 (row: cell; column: gene; the last column: cell group label)
    M2 = as.matrix(data.select[ geneuse, celluse])
    if(length(geneuse)==1){
      df2 <- as.data.frame(M2)
      colnames(df2) = geneuse
    }else if(length(geneuse)>1){
      df2 <- as.data.frame(t(M2))
    }
    df2$group_label <- fct_drop(idents.select[names(idents.select) %in% celluse])
    df2$group_label <- as.factor(df2$group_label)

    # cell groups in organ 2
    group2 <- levels(df2$group_label)

    for (i in 1:length(group2)) {
      df.bs <- df2[df2$group_label == levels(df2$group_label)[i], -ncol(df2)]
      bs.results <- boot(data = df.bs, statistic = model_mean, R)
      # calculate the sd
      x <- apply(bs.results$t, 2, sd)
      bs.results_df <- data.frame(sd = x)
      colnames(bs.results_df) <- paste0(levels(df2$group_label)[i], "_sd")

      # calculate the CI
      q = c(0+(1-CI)/2 , 1-(1-CI)/2)
      CI_min = apply(bs.results$t, 2, quantile,q[1])
      CI_max = apply(bs.results$t, 2, quantile,q[2])
      bs.results_CI = (bs.results$t0 - CI_min)*(CI_max - bs.results$t0)
      bs.results_CI[ bs.results_CI>0 ] = TRUE
      bs.results_CI[ bs.results_CI<=0 ] = FALSE
      bs.results_CI <- data.frame(CI = bs.results_CI)
      colnames(bs.results_CI) <- paste0(levels(df2$group_label)[i], "_CI")

      bootstrap_sd.df <- cbind(bootstrap_sd.df, bs.results_df)
      bootstrap_CI.df <- cbind(bootstrap_CI.df, bs.results_CI)
    }

  }

  bootstrap_output <- list(bootstrap_sd = bootstrap_sd.df, bootstrap_CI = bootstrap_CI.df)
  # net.df <- df
  #
  # net.df <- left_join(net.df, bootstrap_sd.df1, join_by("Receptor"=="Gene"))
  # net.df <- left_join(net.df, bootstrap_sd.df1, join_by("TF"=="Gene"))
  # net.df <- left_join(net.df, bootstrap_sd.df1, join_by("Target"=="Gene"))
  #
  # df$bootstrap_sd <- apply(net.df[, c(paste0(receiver, "_sd"), paste0(receiver, "_sd.x"), paste0(receiver, "_sd.y"))], 1, sum)
  # df$bootstrap_CI <- apply(net.df[, c(paste0(receiver, "_CI"), paste0(receiver, "_CI.y"), paste0(receiver, "_CI.y"))], 1, sum)

  return(bootstrap_output)
}


#' @title Perform the bootstrap test (all LS-T pathways in the object)
#' @description
#' Perform the bootstrap test and calculated the standard deviation for all LS-T pathways stored under "object - net",
#' and filter them based on selected cutoff values.
#'
#' @param object an OrganChat object
#' @param R the number of bootstrap testing to perfomr, the default value is 5.
#' @param sd_cutoff LS-T pathways with a standard deviation (calculated based on bootstrap results) higher than this cutoff value will be removed. The default value is 0.05.
#' @param p_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param CI_cutoff LS-T pathways with a total number of components lying within the confidence interval smaller than this value will be removed. The default value is 0.
#' The input must be 0,1,2, or 3.
#' @param single.condition the input is TRUE by default.
#' @param CI the value of confidence interval, CI = 0.95 by default.
#'
#' @importFrom forcats fct_drop
#' @importFrom boot boot
#'
#' @return an OrganChat object
#' @export
cal_netchat_bootstrap <- function(object, R = 5, sd_cutoff = 0.05, p_cutoff = NULL, CI_cutoff = 0,
                                  single.condition = TRUE, CI = 0.95){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check the sd_cutoff value
  if(is.null(sd_cutoff)){
    sd_cutoff = 0.05
    print("The default value of 'sd_cutoff' is 0.05.")
  }

  # check the p_cutoff value
  if(is.null(p_cutoff)){
    p_cutoff = 0
    print("The default value of 'p_cutoff' is 0.")
  }

  if( !CI_cutoff %in% c(0,1,2,3) ){
    CI_cutoff = 0
    print("Ignore the 'CI_cutoff'.")
  }

  ##############################################################################
  print("Perform bootstrap test (1/2)")

  bootstrap_output_organ2 = cal_singlechat_bootstrap(object = object, R = R, dir = "forward", CI = CI)
  bootstrap_output_organ1 = cal_singlechat_bootstrap(object = object, R = R, dir = "backward", CI = CI)

  ##############################################################################
  print("Integrate results (2/2)")

  # direction -- "forward"
  S.organ <- object@organs[1] # from organ1 (Metabolite) to organ2 (RNA)
  R.organ <- object@organs[2]
  S.group <- levels(object@idents.list[[S.organ]])
  R.group <- levels(object@idents.list[[R.organ]])

  if( length(S.group)>0 & length(R.group)>0 ){
    k = 1
    for (i in 1:length(S.group)) {
      for (j in 1:length(R.group)) {

        if( !is.null(dim(object@net[[con.layer]]$forward[[k]])) & nrow(object@net[[con.layer]]$forward[[k]])>0  ){
          df = object@net[[con.layer]]$forward[[k]]

          net.df <- df
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Target"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Target"=="Gene"))

          df$bootstrap_sd <- apply(net.df[, c(paste0(R.group[j], "_sd"),
                                              paste0(R.group[j], "_sd.x"),
                                              paste0(R.group[j], "_sd.y"))], 1, sum)
          df$bootstrap_CI <- apply(net.df[, c(paste0(R.group[j], "_CI"),
                                              paste0(R.group[j], "_CI.y"),
                                              paste0(R.group[j], "_CI.y"))], 1, sum)

          object@net[[con.layer]]$forward[[k]] <- df[df$bootstrap_sd<sd_cutoff & df$bootstrap_CI>=CI_cutoff , ]
          object@net[[con.layer]]$forward[[k]] <- object@net[[con.layer]]$forward[[k]][object@net[[con.layer]]$forward[[k]]$chatP>=p_cutoff, ]
        }
        k = k+1
      }
    }
  }
  ##############################################################################
  # direction -- "backward"
  S.organ <- object@organs[2] # from organ1 (Metabolite) to organ2 (RNA)
  R.organ <- object@organs[1]
  S.group <- levels(object@idents.list[[S.organ]])
  R.group <- levels(object@idents.list[[R.organ]])

  if( length(S.group)>0 & length(R.group)>0 ){
    k = 1
    for (i in 1:length(S.group)) {
      for (j in 1:length(R.group)) {

        if( !is.null(dim(object@net[[con.layer]]$backward[[k]])) & nrow(object@net[[con.layer]]$backward[[k]])>0 ){
          df = object@net[[con.layer]]$backward[[k]]

          net.df <- df
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Target"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Target"=="Gene"))

          df$bootstrap_sd <- apply(net.df[, c(paste0(R.group[j], "_sd"),
                                              paste0(R.group[j], "_sd.x"),
                                              paste0(R.group[j], "_sd.y"))], 1, sum)
          df$bootstrap_CI <- apply(net.df[, c(paste0(R.group[j], "_CI"),
                                              paste0(R.group[j], "_CI.y"),
                                              paste0(R.group[j], "_CI.y"))], 1, sum)

          object@net[[con.layer]]$backward[[k]] <- df[df$bootstrap_sd<sd_cutoff & df$bootstrap_CI>=CI_cutoff & df$chatP>=p_cutoff, ]
        }

        k = k+1
      }
    }
  }
  ##############################################################################
  # direction -- "self_organ1"
  S.organ <- object@organs[1] # from organ1 (Metabolite) to organ2 (RNA)
  R.organ <- object@organs[1]
  S.group <- levels(object@idents.list[[S.organ]])
  R.group <- levels(object@idents.list[[R.organ]])

  if( length(S.group)>0 & length(R.group)>0 ){
    k = 1
    for (i in 1:length(S.group)) {
      for (j in 1:length(R.group)) {

        if( !is.null(dim(object@net[[con.layer]]$self_organ1[[k]])) & nrow(object@net[[con.layer]]$self_organ1[[k]])>0  ){
          df = object@net[[con.layer]]$self_organ1[[k]]

          net.df <- df
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Target"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ1[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Target"=="Gene"))

          df$bootstrap_sd <- apply(net.df[, c(paste0(R.group[j], "_sd"),
                                              paste0(R.group[j], "_sd.x"),
                                              paste0(R.group[j], "_sd.y"))], 1, sum)
          df$bootstrap_CI <- apply(net.df[, c(paste0(R.group[j], "_CI"),
                                              paste0(R.group[j], "_CI.y"),
                                              paste0(R.group[j], "_CI.y"))], 1, sum)

          object@net[[con.layer]]$self_organ1[[k]] <- df[df$bootstrap_sd<sd_cutoff & df$bootstrap_CI>=CI_cutoff & df$chatP>=p_cutoff, ]
        }

        k = k+1
      }
    }
  }
  ##############################################################################
  # direction -- "self_organ2"
  S.organ <- object@organs[2] # from organ1 (Metabolite) to organ2 (RNA)
  R.organ <- object@organs[2]
  S.group <- levels(object@idents.list[[S.organ]])
  R.group <- levels(object@idents.list[[R.organ]])

  if( length(S.group)>0 & length(R.group)>0 ){
    k = 1
    for (i in 1:length(S.group)){
      for (j in 1:length(R.group)){

        if( !is.null(dim(object@net[[con.layer]]$self_organ2[[k]])) & nrow(object@net[[con.layer]]$self_organ2[[k]])>0  ){
          df = object@net[[con.layer]]$self_organ2[[k]]

          net.df <- df
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[1]][ , c("Gene", paste0(R.group[j],"_sd"))], join_by("Target"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Receptor"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("TF"=="Gene"))
          net.df <- left_join(net.df, bootstrap_output_organ2[[2]][ , c("Gene", paste0(R.group[j],"_CI"))], join_by("Target"=="Gene"))

          df$bootstrap_sd <- apply(net.df[, c(paste0(R.group[j], "_sd"),
                                              paste0(R.group[j], "_sd.x"),
                                              paste0(R.group[j], "_sd.y"))], 1, sum)
          df$bootstrap_CI <- apply(net.df[, c(paste0(R.group[j], "_CI"),
                                              paste0(R.group[j], "_CI.y"),
                                              paste0(R.group[j], "_CI.y"))], 1, sum)

          object@net[[con.layer]]$self_organ2[[k]] <- df[df$bootstrap_sd<sd_cutoff & df$bootstrap_CI>=CI_cutoff & df$chatP>=p_cutoff, ]
        }
        k = k+1
      }
    }
  }

  return(object)
}


#' @title Construct the edge and node data.frames
#' @description
#' Construct the edge and node data.frames from all data frames stored under "object - net".
#'
#' @param object an OrganChat object
#' @param category use this parameter to select the types of LS. The input should be one of the following: "All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite".
#' The default value is "All", meaning all the LS will be included.
#' @param chatP_cutoff LS-T pathways with a communication probability lower than this cutoff value will be removed. The default value is NULL.
#' @param single.condition the input is TRUE by default.
#'
#' @importFrom forcats fct_drop
#'
#' @return a list
#' @export
construct_EdgeNode <- function(object,
                               category = "All",
                               chatP_cutoff = NULL,
                               single.condition = TRUE){

  # "FALSE" is designed mainly for the CAOC
  if(isTRUE(single.condition)){
    con.layer = 1
  }else{
    con.layer = 2
  }

  # check input 'category'
  if( !(category %in% c("All", "hormone_peptide", "cytokines", "hormone_nonpeptide", "metabolite")) == TRUE ){
    message("The input 'category' is invalid, set it to 'All'.")
    pref = "All"
  }

  x1f <- c()
  x2f <- c()
  x1b <- c()
  x2b <- c()
  x1s1 <- c()
  x2s1 <- c()
  x1s2 <- c()
  x2s2 <- c()

  if(length(object@net[[con.layer]]$forward)>0){
    for (i in 1:length(object@net[[con.layer]]$forward)) {
      if( nrow(object@net[[con.layer]]$forward[[i]])>0 ){
        df = object@net[[con.layer]]$forward[[i]]
        # check category
        if( category != "All"){
          ind = which(grepl(category, df$category, fixed=TRUE))
          df = df[ind, ]
        }
        # check chatP_cutoff
        if(!is.null(chatP_cutoff)){
          df = df[df$chatP > chatP_cutoff, ]
        }

        x1f <- c(x1f, df$Path)
        x2f <- c(x2f, df$chatP)
      }
    }
  }

  if(length(object@net[[con.layer]]$backward)>0){
    for (i in 1:length(object@net[[con.layer]]$backward)) {
      if( nrow(object@net[[con.layer]]$backward[[i]])>0 ){
        df = object@net[[con.layer]]$backward[[i]]
        # check category
        if( category != "All"){
          ind = which(grepl(category, df$category, fixed=TRUE))
          df = df[ind, ]
        }
        # check chatP_cutoff
        if(!is.null(chatP_cutoff)){
          df = df[df$chatP > chatP_cutoff, ]
        }

        x1b <- c(x1b, df$Path)
        x2b <- c(x2b, df$chatP)
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ1)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ1)) {
      if( nrow(object@net[[con.layer]]$self_organ1[[i]])>0 ){
        df = object@net[[con.layer]]$self_organ1[[i]]
        # check category
        if( category != "All"){
          ind = which(grepl(category, df$category, fixed=TRUE))
          df = df[ind, ]
        }
        # check chatP_cutoff
        if(!is.null(chatP_cutoff)){
          df = df[df$chatP > chatP_cutoff, ]
        }

        x1s1 <- c(x1s1, df$Path)
        x2s1 <- c(x2s1, df$chatP)
      }
    }
  }

  if(length(object@net[[con.layer]]$self_organ2)>0){
    for (i in 1:length(object@net[[con.layer]]$self_organ2)) {
      if( nrow(object@net[[con.layer]]$self_organ2[[i]])>0 ){
        df = object@net[[con.layer]]$self_organ2[[i]]
        # check category
        if( category != "All"){
          ind = which(grepl(category, df$category, fixed=TRUE))
          df = df[ind, ]
        }
        # check chatP_cutoff
        if(!is.null(chatP_cutoff)){
          df = df[df$chatP > chatP_cutoff, ]
        }

        x1s2 <- c(x1s2, df$Path)
        x2s2 <- c(x2s2, df$chatP)
      }
    }
  }

  # "forward"
  if( length(x1f)>0 ){
    df.aggregate_forward <- data.frame(from = object@organs[1], to = object@organs[2], Path = x1f, chatP = x2f)
    df.edge_forward <- data.frame(from = object@organs[1], to = object@organs[2],
                                  total.counts = length(df.aggregate_forward$Path),
                                  total.chatP = sum(df.aggregate_forward$chatP),
                                  ave.weight = sum(df.aggregate_forward$chatP)/length(df.aggregate_forward$Path))
  }else{
    df.edge_forward = data.frame()
  }

  # "backward"
  if( length(x1b)>0  ){
    df.aggregate_backward <- data.frame(from = object@organs[2], to = object@organs[1], Path = x1b, chatP = x2b)
    df.edge_backward <- data.frame(from = object@organs[2], to = object@organs[1],
                                   total.counts = length(df.aggregate_backward$Path),
                                   total.chatP = sum(df.aggregate_backward$chatP),
                                   ave.weight = sum(df.aggregate_backward$chatP)/length(df.aggregate_backward$Path))
  }else{
    df.edge_backward = data.frame()
  }

  # "self_organ1"
  if( length(x1s1)>0  ){
    df.aggregate_self_organ1 <- data.frame(from = object@organs[1], to = object@organs[1], Path = x1s1, chatP = x2s1)
    df.edge_self_organ1 <- data.frame(from = object@organs[1], to = object@organs[1],
                                      total.counts = length(df.aggregate_self_organ1$Path),
                                      total.chatP = sum(df.aggregate_self_organ1$chatP),
                                      ave.weight = sum(df.aggregate_self_organ1$chatP)/length(df.aggregate_self_organ1$Path))
  }else{
    df.edge_self_organ1 = data.frame()
  }

  # "self_organ2"
  if( length(x1s2)>0 ){
    df.aggregate_self_organ2 <- data.frame(from = object@organs[2], to = object@organs[2], Path = x1s2, chatP = x2s2)
    df.edge_self_organ2 <- data.frame(from = object@organs[2], to = object@organs[2],
                                      total.counts = length(df.aggregate_self_organ2$Path),
                                      total.chatP = sum(df.aggregate_self_organ2$chatP),
                                      ave.weight = sum(df.aggregate_self_organ2$chatP)/length(df.aggregate_self_organ2$Path))
  }else{
    df.edge_self_organ2 = data.frame()
  }

  # merge
  df.edge <- rbind(df.edge_forward, df.edge_backward, df.edge_self_organ1, df.edge_self_organ2)

  # the cell populations for each organ
  Population <- c(length(object@idents.list[[1]]), length(object@idents.list[[2]]))
  df.node <- data.frame(Organ = object@organs, Population = Population)

  EdgeNode <- list(edges = df.edge, nodes = df.node)
  return(EdgeNode)
}


#' @title GO enrichment analysis
#' @description
#' Perform GO enrichment analysis of a data frame of LS-T pathways using the package "gprofiler2".
#'
#' @param df a data frame of LS-T pathways
#' @param LS selected LS to be included in the GO analysis. The default value is NULL, meaning all will be used.
#' @param layer specify which of these will be included: "Receptor", "TF", or "Target". The default value is NULL, meaning all will be used.
#' @param organism the species used in the analysis. The default value is "mmusculus".
#' @param ordered_query argument of the gost function in gprofiler2, the default value is TRUE.
#' @param multi_query argument of the gost function in gprofiler2, the default value is FALSE.
#' @param significant argument of the gost function in gprofiler2, the default value is TRUE.
#' @param exclude_iea argument of the gost function in gprofiler2, the default value is FALSE.
#' @param measure_underrepresentation argument of the gost function in gprofiler2, the default value is FALSE.
#' @param evcodes argument of the gost function in gprofiler2, the default value is FALSE.
#' @param user_threshold argument of the gost function in gprofiler2, the default value is 0.05.
#' @param correction_method argument of the gost function in gprofiler2, the default value is "g_SCS".
#' @param domain_scope argument of the gost function in gprofiler2, the default value is "annotated".
#' @param custom_bg argument of the gost function in gprofiler2, the default value is NULL.
#' @param numeric_ns argument of the gost function in gprofiler2, the default value is "".
#' @param sources argument of the gost function in gprofiler2, the default value is NULL.
#' @param as_short_link argument of the gost function in gprofiler2, the default value is FALSE.
#' @param highlight argument of the gost function in gprofiler2, the default value is FALSE.
#'
#' @importFrom data.table setorder
#' @importFrom gprofiler2 gost
#'
#' @return a list
#' @export
OC_GOanalysis <- function(df,
                          LS = NULL,
                          layer = NULL,
                          organism = "mmusculus",
                          ordered_query = TRUE,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "g_SCS",
                          domain_scope = "annotated", custom_bg = NULL,
                          numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE){

  # select the long-range signals
  if(!is.null(LS)){
    LS = intersect(LS, unique(df$LS))
  }else{
    LS = unique(df$LS)
  }

  # check the input 'layer'
  if( !is.null(layer) & length(setdiff(layer, c("Receptor", "TF", "Target")))>0 ){
    layer = NULL
  }

  # obtain the downstream genes and their expression values
  df = df[df$LS %in% LS, c("Receptor", "TF", "Target", "Receptor.value", "TF.value", "Target.value")]
  df.Receptor = df[, c("Receptor", "Receptor.value")]
  colnames(df.Receptor) = c("Gene", "Value")
  df.TF = df[, c("TF", "TF.value")]
  colnames(df.TF) = c("Gene", "Value")
  df.Target = df[, c("Target", "Target.value")]
  colnames(df.Target) = c("Gene", "Value")
  df.list = list(Receptor = df.Receptor, TF = df.TF, Target = df.Target)

  if( is.null(layer) ){
    df.GO <- rbind(df.Receptor, df.TF, df.Target)
  }else{
    layer.list = list()
    for (i in 1:length(layer)) {
      layer.list[[i]] = df.list[[layer[i]]]
    }
    df.GO <- rbindlist(layer.list)
  }

  # sort by the expression values
  df.GO <- df.GO[!duplicated(df.GO), ]
  df.GO <- setorder(df.GO, -Value)

  # perform GO enrichment analysis
  gostres <- gost(query = df.GO$Gene,
                  organism = organism, ordered_query = ordered_query,
                  multi_query = multi_query, significant = significant, exclude_iea = exclude_iea,
                  measure_underrepresentation = measure_underrepresentation, evcodes = evcodes,
                  user_threshold = user_threshold, correction_method = correction_method,
                  domain_scope = domain_scope, custom_bg = custom_bg,
                  numeric_ns = numeric_ns, sources = sources, as_short_link = as_short_link, highlight = highlight)

  return(gostres)
}


