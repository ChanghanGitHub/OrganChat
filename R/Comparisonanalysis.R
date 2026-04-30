#' @title Calculate the mean expression level from single-cell metabolite data (for two conditions)
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
#' @export
CA_cal_scMeta.bygroup <- function(object,
                                  scMeta.data,
                                  organ.name = NULL,
                                  idents = NULL,
                                  conditions.label = NULL,
                                  mean_method = NULL,
                                  type = "amount",
                                  scale.max = 10){

  conditions.label = factor(conditions.label)
  conditions.label <- fct_drop(conditions.label)

  x = setdiff( object@conditions, levels(conditions.label) )
  if( length(x)>0 ){
    stop("The conditions in the scMETA.data and the object@conditions must be the same.")
  }else{
    conditions.label = factor(conditions.label, levels = object@conditions )
  }

  if (!identical(names(conditions.label), colnames(scMeta.data))) {
    warning("The cell barcodes in 'conditions.label' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the names of 'conditions.label'!")
    names(conditions.label) <- colnames(scMeta.data)
  }

  # select the data for each condition
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

  cell.condition1 = names(conditions.label)[ conditions.label == object@conditions[1] ]
  cell.condition2 = names(conditions.label)[ conditions.label == object@conditions[2] ]

  idents.condition1 = idents[cell.condition1]
  idents.condition2 = idents[cell.condition2]

  scMeta.data.condition1 = scMeta.data[ , cell.condition1]
  scMeta.data.condition2 = scMeta.data[ , cell.condition2]

  # for condition 1
  object = cal_scMeta.bygroup(object = object,
                              scMeta.data = scMeta.data.condition1,
                              organ.name = organ.name,
                              idents = idents.condition1,
                              mean_method = mean_method,
                              scale.max = scale.max,
                              type = type,
                              single.condition = TRUE)

  # for condition 2
  object = cal_scMeta.bygroup(object = object,
                              scMeta.data = scMeta.data.condition2,
                              organ.name = organ.name,
                              idents = idents.condition2,
                              mean_method = mean_method,
                              scale.max = scale.max,
                              type = type,
                              single.condition = FALSE)

  return(object)
}


#' @title Calculate the mean expression level from single-cell transcriptomics data (for two conditions)
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
CA_cal_scRNA.bygroup <- function(object,
                                 DB = OrganChatDB,
                                 mean_method = NULL,
                                 scale.max = 10){
  # for each condition
  for (i in 1:2) {

    # in condition i, select the cells for each organ
    cell.organ1 = names(object@conditions.list[[1]])[ object@conditions.list[[1]] == object@conditions[i] ]
    cell.organ2 = names(object@conditions.list[[2]])[ object@conditions.list[[2]] == object@conditions[i] ]

    idents.organ1 = object@idents.list[[1]][cell.organ1]
    idents.organ2 = object@idents.list[[2]][cell.organ2]

    data.organ1 = object@data.list[[1]][ , cell.organ1]
    data.organ2 = object@data.list[[2]][ , cell.organ2]

    data.singlecondition = list(organ1 = data.organ1, organ2 = data.organ2)
    idents.singlecondition = list(organ1 = idents.organ1, organ2 = idents.organ2)

    object.singlecondition = suppressMessages( create_OrganChat(object =data.singlecondition,
                                                                organs = object@organs,
                                                                idents.list = idents.singlecondition,
                                                                ComparisonAnalysis = FALSE,
                                                                assay = NULL,
                                                                do.sparse = T) )

    object.singlecondition = cal_scRNA.bygroup(object = object.singlecondition,
                                               DB = DB,
                                               mean_method = mean_method,
                                               scale.max = scale.max,
                                               single.condition = TRUE)

    for (j in 1:2) {
      data.new = object.singlecondition@LSdata[[1]][[j]]
      if(nrow(object@LSdata[[i]][[j]])==0){
        object@LSdata[[i]][[j]] = data.new
      }else if( length( setdiff(colnames(object@LSdata[[i]][[j]]), colnames(data.new)) )==0 ){
        data.new = data.new[ setdiff(rownames(data.new), rownames(object@LSdata[[i]][[j]])) , colnames(object@LSdata[[i]][[j]])]
        object@LSdata[[i]][[j]] = rbind(object@LSdata[[i]][[j]], data.new)
      }else{
        warning("The colnames of the existing and new data frames are not identical, skip the new data frame.")
      }
    }
    object@DGdata[[i]] = object.singlecondition@DGdata[[1]]
  }
  return(object)
}


#' @title Merge existing cluster-level metabolite amount/flux data (for two conditions)
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
CA_merge_metabolite.bygroup <- function(object,
                                        metabolitedata.list.condition1 = NULL,
                                        metabolitedata.list.condition2 = NULL,
                                        type = "amount",
                                        scale.max = 10){

  if( length(metabolitedata.list.condition1)!=2 | length(metabolitedata.list.condition2)!=2 ){
    stop("Each 'metabolitedata.list.condition' input should have two data frames representing the flux data of the two organs.")
  }

  for (i in 1:2) {
    metabolitedata.list.condition1[[i]] = metabolitedata.list.condition1[[i]][ , levels(object@idents.list[[i]])]
    metabolitedata.list.condition2[[i]] = metabolitedata.list.condition2[[i]][ , levels(object@idents.list[[i]])]
  }

  object = merge_metabolite.bygroup(object = object,
                                    metabolitedata.list = metabolitedata.list.condition1,
                                    scale.max = scale.max,
                                    type = type,
                                    single.condition = TRUE)

  object = merge_metabolite.bygroup(object = object,
                                    metabolitedata.list = metabolitedata.list.condition2,
                                    scale.max = scale.max,
                                    type = type,
                                    single.condition = FALSE)

  return(object)
}


#' @title LS-R-SE-T pathway inference (for two conditions)
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
#' @param hmdb_dictionary OrganChat-provided database matching "metabolite name" and the HMDB ID.
#'
#' @importFrom forcats fct_drop
#'
#' @return an OrganChat object
#' @export
CA_net0_infer <- function(object,
                          DB = OrganChatDB,
                          LS = NULL,
                          receptor = NULL,
                          tf = NULL,
                          target = NULL,
                          hmdb_dictionary = OrganChatDB_hmdb_dictionary){

  net0.list = list()
  # for each condition
  for (i in 1:2) {

    # in condition i, select the cells for each organ
    cell.organ1 = names(object@conditions.list[[1]])[ object@conditions.list[[1]] == object@conditions[i] ]
    cell.organ2 = names(object@conditions.list[[2]])[ object@conditions.list[[2]] == object@conditions[i] ]

    idents.organ1 = object@idents.list[[1]][cell.organ1]
    idents.organ2 = object@idents.list[[2]][cell.organ2]

    data.organ1 = object@data.list[[1]][ , cell.organ1]
    data.organ2 = object@data.list[[2]][ , cell.organ2]

    data.singlecondition = list(organ1 = data.organ1, organ2 = data.organ2)
    idents.singlecondition = list(organ1 = idents.organ1, organ2 = idents.organ2)

    object.singlecondition = suppressMessages( create_OrganChat(object =data.singlecondition,
                                                                organs = object@organs,
                                                                idents.list = idents.singlecondition,
                                                                ComparisonAnalysis = FALSE,
                                                                assay = NULL,
                                                                do.sparse = T) )

    object.singlecondition@LSdata[[1]] = object@LSdata[[i]]
    object.singlecondition@DGdata[[1]] = object@DGdata[[i]]

    object.singlecondition = OCnet0_infer(object = object.singlecondition,
                                          DB = DB,
                                          LS = LS,
                                          receptor = receptor,
                                          tf = tf,
                                          target = target,
                                          single.condition = TRUE,
                                          hmdb_dictionary = hmdb_dictionary)

    net0.list[[i]] = object.singlecondition@net0
  }

  net0 = net0.list[[1]][ net0.list[[1]]$Path %in% intersect(net0.list[[1]]$Path, net0.list[[2]]$Path) , ]
  object@net0 = net0

  return(object)
}


#' @title Comparison analysis of two conditions
#' @description
#' Calculate the communication probability, perform bootstrap tests, and filter the LS-T pathways for each organ and condition.
#' The output will be a data frame and stored in "Object - net".
#'
#' @param object an OrganChat object
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
#' @importFrom forcats fct_drop
#' @importFrom boot boot
#'
#' @return an OrganChat object
#' @export
CA_net_analysis <- function(object,
                            K = 0.5,
                            N = 2,
                            R = 5,
                            CI = 0.95,
                            CA_cutoff_value = NULL, # "&" logic between cutoff conditions
                            CA_cutoff_chatP = NULL, # "|" logic between cutoff conditions
                            CA_cutoff_sd = NULL,   # "&" logic between cutoff conditions
                            CA_cutoff_CI = 0){ # "&" logic between cutoff conditions

  # check the cutoff values
  if( !is.null(CA_cutoff_value) & length(CA_cutoff_value)!=4 ){
    CA_cutoff_value = NULL
    warning("The input 'CA_cutoff_value' has been ignored, or change it to a 1x4 numeric vector.")
  }

  if( !is.null(CA_cutoff_chatP) & length(CA_cutoff_chatP)!=1 ){
    CA_cutoff_chatP = NULL
    warning("The input 'CA_cutoff_chatP' has been ignored, or change it to a numeric scale.")
  }

  if( !is.null(CA_cutoff_sd) & length(CA_cutoff_sd)!=1 ){
    CA_cutoff_sd = NULL
    warning("The input 'CA_cutoff_sd' has been ignored, or change it to a numeric scale.")
  }

  if( !CA_cutoff_CI %in% c(0,1,2,3) ){
    CA_cutoff_sd = 0
    warning("The input 'CA_cutoff_CI' has been ignored, or select from {0, 1, 2, 3}.")
  }

  # for each condition
  for (i in 1:2) {

    # in condition i, select the cells for each organ
    cell.organ1 = names(object@conditions.list[[1]])[ object@conditions.list[[1]] == object@conditions[i] ]
    cell.organ2 = names(object@conditions.list[[2]])[ object@conditions.list[[2]] == object@conditions[i] ]

    idents.organ1 = object@idents.list[[1]][cell.organ1]
    idents.organ2 = object@idents.list[[2]][cell.organ2]

    data.organ1 = object@data.list[[1]][ , cell.organ1]
    data.organ2 = object@data.list[[2]][ , cell.organ2]

    data.singlecondition = list(organ1 = data.organ1, organ2 = data.organ2)
    idents.singlecondition = list(organ1 = idents.organ1, organ2 = idents.organ2)

    object.singlecondition = suppressMessages( create_OrganChat(object =data.singlecondition,
                                                                organs = object@organs,
                                                                idents.list = idents.singlecondition,
                                                                ComparisonAnalysis = FALSE,
                                                                assay = NULL,
                                                                do.sparse = T) )

    object.singlecondition@LSdata[[1]] = object@LSdata[[i]]
    object.singlecondition@DGdata[[1]] = object@DGdata[[i]]
    object.singlecondition@net0 = object@net0

    object.singlecondition = filter_OCnet0(object = object.singlecondition,
                                           cutoff_value = c(0, 0, 0, 0),
                                           single.condition = TRUE)

    object.singlecondition = cal_netchat_prob(object = object.singlecondition,
                                              K = K,
                                              N = N,
                                              single.condition = TRUE)

    object@net[[i]] = object.singlecondition@net[[1]]
  }

  # filter the data frames
  for (j in 1:4) { # for four directions
    for (k in 1:length(object@net[[1]][[j]])) {
      ncol.max = ncol(object@net[[1]][[j]][[k]])
      # condition 1
      df1 = object@net[[1]][[j]][[k]]
      # condition 2
      df2 = object@net[[2]][[j]][[k]]

      df.merge = inner_join(df1, df2, join_by(Path))

      if( !is.null(CA_cutoff_value) ){
        df.merge = df.merge[ df.merge$LS.value.x >CA_cutoff_value[1]  |
                               df.merge$LS.value.y >CA_cutoff_value[1] , ]

        df.merge = df.merge[ df.merge$Receptor.value.x > CA_cutoff_value[2] |
                               df.merge$Receptor.value.y > CA_cutoff_value[2] , ]

        df.merge = df.merge[  df.merge$TF.value.x > CA_cutoff_value[3] |
                                df.merge$TF.value.y > CA_cutoff_value[3] , ]


        df.merge = df.merge[ df.merge$Target.value.x > CA_cutoff_value[4] |
                               df.merge$Target.value.y > CA_cutoff_value[4] , ]
      }

      if( !is.null(CA_cutoff_chatP) ){
        df.merge = df.merge[ df.merge$chatP.x > CA_cutoff_chatP |
                               df.merge$chatP.y > CA_cutoff_chatP , ]
      }

      object@net[[1]][[j]][[k]] = object@net[[1]][[j]][[k]][ object@net[[1]][[j]][[k]]$Path %in% df.merge$Path,  ]
      object@net[[2]][[j]][[k]] = object@net[[2]][[j]][[k]][ object@net[[2]][[j]][[k]]$Path %in% df.merge$Path,  ]
    }
  }

  # bootstrap test
  # for each condition
  for (i in 1:2) {

    # in condition i, select the cells for each organ
    cell.organ1 = names(object@conditions.list[[1]])[ object@conditions.list[[1]] == object@conditions[i] ]
    cell.organ2 = names(object@conditions.list[[2]])[ object@conditions.list[[2]] == object@conditions[i] ]

    idents.organ1 = object@idents.list[[1]][cell.organ1]
    idents.organ2 = object@idents.list[[2]][cell.organ2]

    data.organ1 = object@data.list[[1]][ , cell.organ1]
    data.organ2 = object@data.list[[2]][ , cell.organ2]

    data.singlecondition = list(organ1 = data.organ1, organ2 = data.organ2)
    idents.singlecondition = list(organ1 = idents.organ1, organ2 = idents.organ2)

    object.singlecondition = suppressMessages( create_OrganChat(object =data.singlecondition,
                                                                organs = object@organs,
                                                                idents.list = idents.singlecondition,
                                                                ComparisonAnalysis = FALSE,
                                                                assay = NULL,
                                                                do.sparse = T) )

    object.singlecondition@LSdata[[1]] = object@LSdata[[i]]
    object.singlecondition@DGdata[[1]] = object@DGdata[[i]]
    object.singlecondition@net0 = object@net0

    object.singlecondition@net[[1]] = object@net[[i]]

    object.singlecondition = cal_netchat_bootstrap(object = object.singlecondition,
                                                   R = R,
                                                   sd_cutoff = Inf,
                                                   p_cutoff = -Inf,
                                                   CI_cutoff = 0,
                                                   single.condition = TRUE,
                                                   CI = CI)

    object@net[[i]] = object.singlecondition@net[[1]]
  }

  # filter the data frames
  for (j in 1:4) { # for four directions
    for (k in 1:length(object@net[[1]][[j]])) {
      ncol.max = ncol(object@net[[1]][[j]][[k]])
      # condition 1
      df1 = object@net[[1]][[j]][[k]]
      # condition 2
      df2 = object@net[[2]][[j]][[k]]

      if( nrow(df1)>0 & nrow(df2)>0 ){
        df.merge = inner_join(df1, df2, join_by(Path))
      }else{
        next
      }

      df.merge = inner_join(df1, df2, join_by(Path))

      if( !is.null(CA_cutoff_sd) ){
        df.merge = df.merge[ df.merge$bootstrap_sd.x < CA_cutoff_sd &
                               df.merge$bootstrap_sd.y < CA_cutoff_sd , ]
      }

      df.merge = df.merge[ df.merge$bootstrap_CI.x >= CA_cutoff_CI &
                             df.merge$bootstrap_CI.y >= CA_cutoff_CI , ]


      object@net[[1]][[j]][[k]] = object@net[[1]][[j]][[k]][ object@net[[1]][[j]][[k]]$Path %in% df.merge$Path,  ]
      object@net[[2]][[j]][[k]] = object@net[[2]][[j]][[k]][ object@net[[2]][[j]][[k]]$Path %in% df.merge$Path,  ]
    }
  }

  return(object)
}


#' @title Export pathways
#' @description
#' Export all inferred pathways in one table
#'
#' @param object an OrganChat object
#'
#' @return a data frame
#' @export
Export_path_table <- function(object){

  path_table = list()
  # 'length(object@net[object@conditions[i]][[1]])' is just to collect all four directions
  for (i in 1:length(object@conditions)) { # controls the conditions
    df.list = list()
    ind = 1

    for (j in 1:4) { # controls the directions
      # all data frames for condition i direction j
      net.list <- object@net[[i]][[j]]

      for (k in 1:length(net.list)){
        if( !is.null(dim(net.list[[k]])) & nrow(net.list[[k]])>0 ){
          df.list[[ind]] = net.list[[k]]
          df.list[[ind]]$Condition = object@conditions[i]
          df.list[[ind]]$Direction = names(object@net[[i]][j])
          ind = ind + 1
        }
      }

      path_table[[i]] = rbindlist(df.list, use.names = TRUE)
      path_table[[i]] = path_table[[i]][!duplicated(path_table[[i]]), ]
    }

  }
  output = path_table
  names(output) = object@conditions

  return(output)
}


#' @title Calculate the quantile position
#' @description
#' Calculate the quantile position for every number in vector x
#'
#' @param x a vector
#'
#' @return a vector
#' @export
quant_vector <- function(x){

  if(class(x)[1]=="data.table"){
    x = as.numeric(unlist(x))
  }

  y = x
  for (i in 1:length(x)) {
    y[i] = length(which(x<x[i]))/length(x)
  }

  return(y)
}


#' @title Infer differential expressed pathways
#' @description
#' Generate a data frame for the LS-T pathways for comparison between conditions
#' The output will be a data frame and stored in "Object - DEpath".
#'
#' @param object an OrganChat object
#' @param Sender.group specify which Sender group to select. The default value is NULL.
#' @param Receiver.group specify which Receiver group to select. The default value is NULL.
#' @param Direction specify which direction to select. The default value is NULL.
#' @param p_cutoff LS-T pathways with a communication probability lower than this cutoff value under both conditions will be removed.
#' The default value is NULL.
#' @param sd_cutoff LS-T pathways with a standard deviation (calculated based on bootstrap results) higher than this cutoff value under at least one condition will be removed.
#' The default value is NULL.
#' @param DI_cutoff keep the LS-T pathways which the absolute RI value is higher then this cutoff value. The default value is NULL.
#' @param top keep the top pathways based on the positive and negative RI, respectively. The default value is NULL.
#'
#' @importFrom forcats fct_drop
#' @importFrom dplyr arrange desc slice_head
#'
#' @return an OrganChat object
#' @export
DEpath_table <- function(object,
                         Sender.group = NULL,
                         Receiver.group = NULL,
                         Direction = NULL,
                         p_cutoff = NULL,
                         sd_cutoff = NULL,
                         DI_cutoff = NULL,
                         top = NULL){

  # get the list of two tables for two conditions
  path.table = Export_path_table(object = object)

  if( nrow(path.table[[1]])==0 | nrow(path.table[[2]])==0 ){
    stop("No DE path inferred.")
  }

  # generate a Path_ID
  for (i in 1:2) {
    path.table[[i]]$Path_ID = paste0(path.table[[i]]$Path, "_",
                                     path.table[[i]]$Sender.group, "_",
                                     path.table[[i]]$Receiver.group, "_",
                                     path.table[[i]]$Direction)
    path.table[[i]] = path.table[[i]][ , c("Path_ID", "LS", "Receptor", "TF", "Target", "Path",
                                           "category", "id", "LS_name",
                                           "Sender.group", "Receiver.group", "Direction",
                                           "chatP", "bootstrap_sd")]
    colnames(path.table[[i]])[13:14] = c(paste0("chatP", "_", object@conditions[i]),
                                         paste0("bootstrap_sd", "_", object@conditions[i]))
  }

  DEpath.table = inner_join(path.table[[1]], path.table[[2]][ , c(1,13,14)], join_by(Path_ID))

  df <- data.frame(chatP_1 = DEpath.table[, 13],
                   chatP_2 = DEpath.table[, 15],
                   sd_1 = DEpath.table[, 14],
                   sd_2 = DEpath.table[, 16])
  colnames(df) = c("chatP_1", "chatP_2", "sd_1", "sd_2")

  # log2FC of the two chatP
  x = c(df$chatP_1, df$chatP_2)
  x = x[x>0]
  DEpath.table$log2FC = log2((df$chatP_1 + min(x))/(df$chatP_2 + min(x)))
  # adjust the log2FC
  Fx = 1.35*apply(df[, 1:2], 1, max)^2/(apply(df[, 1:2], 1, max)^2 + 0.125)
  Fy = -0.75*( (quant_vector(df[ ,3]) + quant_vector(df[ ,4]))/2 )^2 +1
  adjlog2FC = log2((df$chatP_1 + min(x))/(df$chatP_2 + min(x)))*Fx*Fy
  # calculate the RI
  DEpath.table$DI = 2/(1+exp(-2*adjlog2FC))-1

  if(!is.null(Sender.group)){
    DEpath.table = DEpath.table[DEpath.table$Sender.group %in% intersect(unique(DEpath.table$Sender.group), Sender.group), ]
  }

  if(!is.null(Receiver.group)){
    DEpath.table = DEpath.table[DEpath.table$Receiver.group %in% intersect(unique(DEpath.table$Receiver.group), Receiver.group), ]
  }

  if(!is.null(Direction)){
    DEpath.table = DEpath.table[DEpath.table$Direction %in% intersect(unique(DEpath.table$Direction), Direction), ]
  }

  if(!is.null(p_cutoff)){
    DEpath.table = DEpath.table[unlist(DEpath.table[ , 13]) >= p_cutoff | unlist(DEpath.table[ , 15]) >= p_cutoff, ]
  }

  if(!is.null(sd_cutoff)){
    DEpath.table = DEpath.table[unlist(DEpath.table[ , 14]) <= sd_cutoff & unlist(DEpath.table[ , 16]) <= sd_cutoff, ]
  }

  if(!is.null(DI_cutoff)){
    DEpath.table = DEpath.table[ abs(DEpath.table$RI) >= DI_cutoff, ]
  }

  if(!is.null(top)){
    DEpath.table %>%
      arrange(desc(DEpath.table$DI)) %>%
      slice_head(n = top) -> DEpath.1

    DEpath.table %>%
      arrange(DEpath.table$DI) %>%
      slice_head(n = top) -> DEpath.2

    DEpath.table <- rbind(DEpath.1, DEpath.2)
    DEpath.table <- DEpath.table[!duplicated(DEpath.table), ]
  }

  object@DEpath = DEpath.table

  return(object)
}



















