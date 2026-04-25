#' @title The OrganChat class
#' @description
#' The OrganChat object is created using the following inputs: (1) single-cell transcriptomics data matrix (normalized);
#' (2) organ names; (3) cell cluster and condition labels (named lists with cell barcodes).
#'
#' @import Matrix
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom data.table data.table
setClassUnion(name = 'AnyMatrix', members = c("matrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))
setClassUnion(name = 'AnyDF', members = c("data.frame"))

#' The key slots used in the OrganChat object are described below.
#'
#' @slot organs a char contains organ names (should be the same names recorded in the meta data)
#' @slot conditions a char contains condition labels (should be the same names recorded in the meta data)
#' @slot data.list a list of single-cell data matrices for each organ
#' @slot idents.list a list of cluster labels (named with cell barcodes) for each organ
#' @slot conditions.list a list of condition labels (named with cell barcodes) for each organ.
#' The input 'conditions.list' must have the same levels listed in the input "conditions".
#' @slot LSdata a list of data frame storing the cluster-level average expression of the long-range signal (LS) for each organ
#' @slot DGdata a list of data frame storing the cluster-level average expression of the downstream gene (DG, including R, SE, T) for each organ
#' @slot net0 a data frame storing all possible LS-T pathways inferred based on the OrganChatDB, without any filtering
#' @slot net a list of data frame with multiple layers, storing the filtered LS-T pathways for each cluster-cluster pair, in each condition and organ
#' @slot DEpath a data frame storing the comparison results of the LS-T pathways between conditions
#' @slot options a list of miscellaneous data, such as parameters used throughout analysis
#' #'
#' @exportClass OrganChat
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
object <- methods::setClass("OrganChat",
                            slots = c(organs = 'character',
                                      conditions = 'character',
                                      data.list = 'list',
                                      idents.list = 'list',
                                      conditions.list = 'list',
                                      LSdata = 'list',
                                      DGdata = 'list',
                                      net0 = "data.frame",
                                      net = 'list',
                                      DEpath = "data.frame",
                                      options = "list")
)


#' @title Create an OrganChat object
#' @description
#' Create an OrganChat object using single-cell transcriptomics data matrix (normalized), meta data of organ, condition, and cluster.
#'
#' @param object a normalized (NOT count) data matrix (gene-by-cell)
#' @param organs a char contains organ names (should be the same names recorded in the meta data)
#' @param idents.list a list of cluster labels (named with cell barcodes) for each organ
#' @param ComparisonAnalysis a boolean value to determine if the comparison between condition is performing,
#' comparison analysis will be used for multi-organ (more than 2) scenario. Only the clusters shared by both conditions will be used for comparison analysis.
#' @param conditions a char contains condition labels (should be the same names recorded in the meta data)
#' @param conditions.list a list of condition labels (named with cell barcodes) for each organ
#' @param assay an assay to use when the input is a Seurat object
#' @param do.sparse whether use sparse format
#'
#' @return an OrganChat object
#' @export
#' @importFrom methods as new
#' @importFrom Matrix t
create_OrganChat <- function(object,
                             organs = NULL,
                             idents.list = NULL,
                             ComparisonAnalysis = FALSE,
                             conditions = NULL,
                             conditions.list = NULL,
                             assay = NULL,
                             do.sparse = T){

  # data matrix as input
  if(class(object) == "list"){
    for (i in 1:length(object)) {
      if(inherits(x = object[[i]], what = c("matrix", "Matrix", "dgCMatrix"))){
        next
      }else{ stop("The input data should be a list of normalized expression matrix of the organs.") }
    }
    data.list <- object
  }

  if(class(idents.list) != "list" | length(idents.list) != length(data.list)){
    stop("The input 'idents.list' must be the list of the cluster labels for each organ.")
  }

  if(is.null(organs) | length(organs)!=length(data.list)){
    stop("The input 'organs' must be the organ names.")
  }else{
    names(data.list) = organs
    names(idents.list) = organs
  }

  for (i in 1:length(data.list)) {
    idents.list[[i]] = factor(idents.list[[i]])
    idents.list[[i]] <- fct_drop(idents.list[[i]])

    if (!identical(names(idents.list[[i]]), colnames(data.list[[i]]))) {
      warning("The cell barcodes in 'idents.list' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the names of 'idents.list'!")
      names(idents.list[[i]]) <- colnames(data.list[[i]])
    }
  }

  LSdata = list( Condition1 = list(),  Condition2 = list() )
  DGdata = list( Condition1 = list(),  Condition2 = list() )
  net = list( Condition1 = list(),  Condition2 = list() )

  for (i in 1:length(data.list)) {
    LSdata[[1]][[i]] = data.frame()
    LSdata[[2]][[i]] = data.frame()
    DGdata[[1]][[i]] = data.frame()
    DGdata[[2]][[i]] = data.frame()
  }

  names(LSdata[[1]]) = organs
  names(LSdata[[2]]) = organs
  names(DGdata[[1]]) = organs
  names(DGdata[[2]]) = organs

  ##############################################################################
  # for ComparisonAnalysis
  if(is.null(ComparisonAnalysis) & length(organs)==2){
    ComparisonAnalysis = FALSE
  }

  if(!isTRUE(ComparisonAnalysis) & length(organs)>2){
    ComparisonAnalysis = TRUE
    warning("For multi-organ (more than 2) scenario, only the comparison analysis is allowed.")
  }

  if(isFALSE(ComparisonAnalysis)){
    conditions = character()
    conditions.list = list()
  }else if( is.null(conditions.list) | length(conditions.list)!=length(data.list) | is.null(conditions) | length(conditions)!=2 ){
    stop("For comparison analysis,the name of two conditions must be specified in the input 'conditions', and 'conditions.list' must be the list of condition labels.")
  }

  if( isTRUE(ComparisonAnalysis) ){
    # name the list of condition labels
    names(conditions.list) = organs
    # check the levels of the factors
    for (i in 1:length(organs)) {
      if( length(setdiff(conditions, unique(conditions.list[[i]])))==0){
        conditions.list[[i]] = conditions.list[[i]][ conditions.list[[i]] %in% conditions ]
        conditions.list[[i]] = factor(conditions.list[[i]], levels = conditions)
        conditions.list[[i]] <- fct_drop(conditions.list[[i]])
      }else{ stop("The input 'conditions.list' must have the same levels.") }
    }

    # add cell parcodes
    for (i in 1:length(data.list)) {
      if (!identical(names(conditions.list[[i]]), colnames(data.list[[i]]))) {
        warning("The cell barcodes in 'conditions.list' is different from those in the used data matrix.")
        # names(conditions.list[[i]]) <- colnames(data.list[[i]])
        data.list[[i]] = data.list[[i]][ , names(conditions.list[[i]])]
      }
    }

    # for each organ, check the cell types in different conditions
    for (i in 1:length(organs)) {
      barcodes1 = names(conditions.list[[i]][ conditions.list[[i]] == conditions[1] ])
      barcodes2 = names(conditions.list[[i]][ conditions.list[[i]] == conditions[2] ])
      ident1 = idents.list[[i]][ barcodes1 ]
      ident2 = idents.list[[i]][ barcodes2 ]
      ident1 <- fct_drop(ident1)
      ident2 <- fct_drop(ident2)

      if( identical(levels(ident1), levels(ident2))==FALSE ){
        # select the shared cell types
        ident.shared = intersect(levels(ident1), levels(ident2))

        if(length(ident.shared)==0){
          message(paste0("There is no shared cell types in organ ", organs[i]))
          stop("Cannot perform comparison analysis!")
        }else{
          message("The following cell types in ", organs[i], " will be used: ")
          print(ident.shared)

          # update the idents.list, conditions.list, data.list
          idents.list[[i]] = idents.list[[i]][idents.list[[i]] %in% ident.shared]
          barcodes.shared = intersect(names(idents.list[[i]]), c(barcodes1, barcodes2))
          conditions.list[[i]] = conditions.list[[i]][barcodes.shared]
          data.list[[i]] = data.list[[i]][ , barcodes.shared]
          idents.list[[i]] = idents.list[[i]][barcodes.shared]

          idents.list[[i]] <- fct_drop(idents.list[[i]])
          conditions.list[[i]] <- fct_drop(conditions.list[[i]])
        }
      }

    }


    # change the names for LSdata, DGdata, net
    names(LSdata) = conditions
    names(DGdata) = conditions
    names(net) = conditions
  }


  object <- methods::new(Class = "OrganChat",
                         organs = organs,
                         conditions = conditions,
                         data.list = data.list,
                         idents.list = idents.list,
                         conditions.list = conditions.list,
                         LSdata = LSdata,
                         DGdata = DGdata,
                         net = net)
  object@options$mode <- "single"
  return(object)

}





