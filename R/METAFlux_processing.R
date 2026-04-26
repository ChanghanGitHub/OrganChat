#' @title Updated METAFlux generate_boots function
#' @description
#' Updated METAFlux generate_boots function to be compatible with Seurat V5
#' generate_boots: Generate bootstrap index data
#' @param celltype colnames of single cell data.Colnames should be labeled as cell type or cluster.
#' @param n  number of bootstrap
#'
#'
#' @return a data frame
#'
generate_boots <- function(celltype, n) {
dt <- data.frame(cluster = celltype, id = 1:length(celltype))
index <- do.call(cbind, sapply(1:n, function(x) {
  splits <- dt %>%
    group_by(cluster) %>%
    sample_n(dplyr::n(), replace = TRUE) %>%
    ungroup() %>%
    dplyr::select("id")
}))
return(index)
}


#' @title Updated METAFlux get_ave_exp function
#' @description
#' Updated METAFlux get_ave_exp function to be compatible with Seurat V5
#' get_ave_exp: Calculate mean expression for one bootstrap
#'
#' @param i index
#' @param myseurat single cell Seurat object.METAFlux will calculate on "data" slot
#' @param samples generated bootstrap index data
#' @param myident Seurat idents.This will be a character string indicating the grouping of the seurat object
#'
#' @return a data frame
#'
#' @import Seurat
#' @import dplyr
#' @importFrom SeuratObject CreateAssay5Object
#'
get_ave_exp_Updated <- function(i, myseurat, samples, myident) {
  meta.data = myseurat@meta.data[samples[,i],]
  sample <- myseurat@assays$RNA$counts[,samples[,i]]
  name <- colnames(sample)
  for (j in 1:length(name)) {
    name[j] <- paste0(name[j], "_", j)
  }
  colnames(sample) <- name
  rownames(meta.data) <- name
  sample <- CreateAssay5Object(sample)
  options(Seurat.object.assay.version = "v5")
  SeuratObject <- suppressWarnings(
    CreateSeuratObject(counts=sample, meta.data = meta.data))
  SeuratObject <- NormalizeData(SeuratObject, verbose = FALSE)
  ave <- GetAssayData(AverageExpression(SeuratObject, group.by = myident,return.seurat = T), assay = "RNA") %>% as.data.frame()
  return(ave)
}


#' @title Updated METAFlux calculate_avg_exp function
#' @description
#' Updated METAFlux calculate_avg_exp function to be compatible with Seurat V5
#' calculate_avg_exp: Calculate bootstrap mean expression for single cell data
#'
#' @param n_bootstrap number of bootstrap
#' @param seed random seed
#' @param myseurat Seurat object. METAFlux will calculate on "data" slot
#' @param myident Seurat idents for grouping.This will be a character string indicating the grouping of the seurat object
#'
#' @return mean expression data
#'
#' @import Seurat
#'
#' @export
#'
calculate_avg_exp_Updated <- function(myseurat, myident, n_bootstrap, seed) {
  set.seed(seed)
  samples = generate_boots(myseurat@meta.data[,myident], n_bootstrap)
  exp <- lapply(1:n_bootstrap, get_ave_exp_Updated, myseurat, samples, myident)
  exp <- do.call(cbind, exp)
  return(exp)
}


#' @title Perform METAFlux analysis in OrganChat
#' @description
#' Run METAFlux in OrganChat with updated functions
#'
#' @param Data.input Seurat object
#' @param myident Seurat idents for grouping.This will be a character string indicating the grouping of the seurat object
#' @param n_bootstrap number of the bootstrap test, the default value is 5.
#' @param seed random seed, the default value is 1.
#' @param medium = human_blood
#'
#' @return mean expression data
#'
#' @import Seurat
#'
#'
#'
#'
#'
#' @export
#'
run_METAFlux <- function(Data.input,
                         myident,
                         n_bootstrap=5,
                         seed=1,
                         medium = human_blood){

  # library(METAFlux)

  if(class(Data.input) != "Seurat"){
    stop("The input data should be a Seurat object.")
  }

  if(is.null(myident)){
    stop("The input 'myident' is missing.")
  }

  mean_exp = calculate_avg_exp_Updated(myseurat = Data.input, myident = myident, n_bootstrap = n_bootstrap, seed = seed)

  #calculate metabolic reaction scores
  scores <- calculate_reaction_score(data = mean_exp)

  #calculate the fractions of celltype/clusters
  x = round(table(Data.input@meta.data[ , myident])/nrow(Data.input@meta.data), 3)
  x[1] = 1- sum(x[2:length(x)])

  flux = compute_sc_flux(num_cell = length(x), fraction = unname(x), fluxscore = scores, medium = medium)

  return(flux)
}


#' @title Process and reformat the METAFlux output
#' @description
#' Calculate the average flux value based on the bootstrap results; scale the influx and outflux, respectively; adding cluster labels, reaction IDs, etc.
#'
#' @param flux the METAFlux output table
#' @param cell.types a char of cell cluster levels
#' @param exclude specify if neglect the the reactions driven by external_medium, internal_medium, or both.
#' The input should be one of the following: "None", "external_medium", "internal_medium", "Both". The default value is "None".
#' @param METAFlux_lookup use the METAFlux_lookup database to match the reaction ID.
#' @param scale.byflux a boolean value to determine if scale the average influx and outflux data to 0 and 1 or -1 and 0, respectively.
#'
#' @return a data frame
#'
#' @importFrom stringr str_split
#' @importFrom forcats fct_drop
#'
#'
#'
#'
#'
#' @import dplyr
#'
#' @export
#'
METAFlux_reformat <- function(flux,
                              cell.types,
                              exclude = "None",
                              METAFlux_lookup = METAFlux_lookup,
                              scale.byflux = TRUE){

  if(is.null(cell.types)){ stop("The file 'cell.types' is missing.") }

  if(is.null(exclude)){
    exclude = "None"
  }else if( !exclude %in% c("None", "external_medium", "internal_medium", "Both")){
    stop("The input 'exclude' must be one of the following: 'None', 'external_medium', 'internal_medium', or 'Both'.")
  }

  if(is.null(METAFlux_lookup)){ stop("The file 'METAFlux_lookup' is missing.") }
  if(is.null(scale.byflux)){ scale.byflux = TRUE }

  # separate it into two parts: celltype & external_medium
  flux_celltype <- flux[ (startsWith(rownames(flux), "celltype") ) ==TRUE, ]
  flux_exmedium <- flux[ (startsWith(rownames(flux), "external_medium") ) ==TRUE, ]

  # for celltype part, separate it by determining if its an internal_medium or not
  flux_celltype_normal <- flux_celltype[ grepl("internal_medium", rownames(flux_celltype), fixed=TRUE) ==FALSE,  ]
  flux_celltype_inmedium <- flux_celltype[ grepl("internal_medium", rownames(flux_celltype), fixed=TRUE) ==TRUE,  ]

  ##############################################################################
  # for "flux_celltype_normal"
  df <- flux_celltype_normal
  x <- str_split(rownames(df), fixed(" "), n=3)
  cluster <- c()
  ReactionID <- c()

  for (i in 1:length(x)) {
    cluster <- c(cluster, x[[i]][2])
    ReactionID <- c(ReactionID, x[[i]][3])
  }

  df$cluster <- as.factor(cluster)
  df$ReactionID <- ReactionID
  df$celltype_internal_medium <- "FALSE"
  rownames(df) <- 1:length(x)
  # output
  df.normal <- df

  ##############################################################################
  # for "flux_celltype_inmedium"
  df <- flux_celltype_inmedium
  x <- str_split(rownames(df), fixed(" "), n=4)
  cluster <- c()
  ReactionID <- c()

  for (i in 1:length(x)) {
    cluster <- c(cluster, x[[i]][2])
    ReactionID <- c(ReactionID, x[[i]][4])
  }

  df$cluster <- as.factor(cluster)
  df$ReactionID <- ReactionID
  df$celltype_internal_medium <- "TRUE"
  rownames(df) <- 1:length(x)
  # output
  df.inmedium <- df

  ##############################################################################
  # for "flux_exmedium"
  df <- flux_exmedium
  x <- str_split(rownames(df), fixed(" "), n=2)
  cluster <- c()
  ReactionID <- c()

  for (i in 1:length(x)) {
    cluster <- c(cluster, x[[i]][1])
    ReactionID <- c(ReactionID, x[[i]][2])
  }

  df$cluster <- as.factor(cluster)
  df$ReactionID <- ReactionID
  df$celltype_internal_medium <- "FALSE"
  rownames(df) <- 1:length(x)
  # output
  df.exmedium <- df

  ##############################################################################

  if( exclude == "None" ){
    flux_reformatted <- rbind(df.normal, df.inmedium, df.exmedium)
    flux_reformatted$cluster <- as.factor(flux_reformatted$cluster)
    levels(flux_reformatted$cluster) <- c(cell.types, "external_medium")
  }else if( exclude == "Both" ){
    flux_reformatted <- df.normal
    flux_reformatted$cluster <- as.factor(flux_reformatted$cluster)
    levels(flux_reformatted$cluster) <- c(cell.types)
  }else if( exclude == "internal_medium" ){
    flux_reformatted <- rbind(df.normal, df.exmedium)
    flux_reformatted$cluster <- as.factor(flux_reformatted$cluster)
    levels(flux_reformatted$cluster) <- c(cell.types, "external_medium")
  }else if( exclude == "external_medium" ){
    flux_reformatted <- rbind(df.normal, df.inmedium)
    flux_reformatted$cluster <- as.factor(flux_reformatted$cluster)
    levels(flux_reformatted$cluster) <- c(cell.types)
  }

  # sd & mean
  n_bootstrap = ncol(flux)
  flux_reformatted$Mean <- apply(flux_reformatted[, 1:n_bootstrap], 1, mean)
  flux_reformatted$Std <- apply(flux_reformatted[, 1:n_bootstrap], 1, sd)

  # find the rows that recorded in the "METAFlux_lookup" table
  matched_ID <- METAFlux_lookup$ReactionID[METAFlux_lookup$Status == 1]
  flux_reformatted <- flux_reformatted[ (flux_reformatted$ReactionID %in% matched_ID), ]

  if( scale.byflux == TRUE ){
    flux_r <- flux_reformatted[flux_reformatted$Mean > 0, ] # release
    flux_u <- flux_reformatted[flux_reformatted$Mean < 0, ] # uptake
    # scale the values
    flux_r$Mean <- (flux_r$Mean - min(flux_r$Mean)) / (max(flux_r$Mean) - min(flux_r$Mean))
    flux_u$Mean <- -(flux_u$Mean - min(flux_u$Mean)) / (max(flux_u$Mean) - min(flux_u$Mean))
    # merge
    flux_output <- rbind(flux_r, flux_u)
  }else if( scale.byflux == FALSE ){
    flux_output <- flux_reformatted
  }

  return(flux_output)
}


#' @title Convert the processed METAFlux output for OrganChat analysis
#' @description
#' Match the ReactionID (Metabolite name) to the HMDB ID, construct the data frame (row--metabolites (LS), col--clusters) that is ready
#' to be stored under "OrganChat object - LSdata - organ".
#'
#' @param cell.types a char of cell cluster levels
#' @param flux the processed and reformated METAFlux output table
#' @param METAFlux_lookup use the METAFlux_lookup database to match the reaction ID.
#' @param hmdb_dictionary OrganChat-provided database matching "metabolite name" and the HMDB ID.
#'
#' @return a data frame
#'
#' @importFrom stats aggregate
#' @importFrom stringr str_split
#' @importFrom forcats fct_drop
#' @import dplyr
#'
#' @export
#'
METAFlux_convert_OrganChat <- function(cell.types,
                                       flux,
                                       METAFlux_lookup = METAFlux_lookup,
                                       hmdb_dictionary = hmdb_dictionary){

  if(is.null(cell.types)){ stop("The file 'cell.types' is missing.") }
  if(is.null(METAFlux_lookup)){ stop("The file 'METAFlux_lookup' is missing.") }
  if(is.null(hmdb_dictionary)){ stop("The file 'hmdb_dictionary' is missing.") }

  flux <- flux[ flux$cluster %in% cell.types , c("cluster", "ReactionID", "Mean")]
  flux <- left_join(flux, METAFlux_lookup, by = join_by(ReactionID), relationship = "many-to-many")
  # fix the upper/lower case issue
  flux$Name <- tolower(flux$Name)
  hmdb_dictionary$Name <- tolower(hmdb_dictionary$Name)

  # fine the matched meta IDs
  flux <- left_join(flux, hmdb_dictionary, by = join_by(Name), multiple = "all", relationship = "many-to-many")
  flux_output <- data.frame(Metabolite = flux$ID, Value = flux$Mean, Cluster = flux$cluster)
  flux_output <- aggregate(flux_output$Value, by = list(flux_output$Metabolite, flux_output$Cluster), sum)
  colnames(flux_output) <- c("Metabolite", "Cluster", "Value")

  # check if the Number of metabolites between clusters are equal
  if( length(unique(flux_output$Cluster))!=length(cell.types) ){
    stop("The number of clusters in the flux data is not equal to the number of cell.types.")
  }

  flux_output$Cluster <- factor(flux_output$Cluster, levels = cell.types)
  if(length(unique(table(flux_output$Cluster)))!=1){
    warning("Number of metabolites in each cluster is not equal.")
    n.keep = min(table(flux_output$Cluster))
    flux_output %>%
      group_by(flux_output$Cluster) %>%
      slice_head(n = n.keep) %>%
      ungroup() -> flux_output
  }

  if(length(unique(table(flux_output$Metabolite)))!=1){
    warning("Number of clusters in each metabolite is not equal.")
    n.keep = min(table(flux_output$Metabolite))
    flux_output %>%
      group_by(flux_output$Metabolite) %>%
      slice_head(n = n.keep) %>%
      ungroup() -> flux_output
  }

  LS.df = data.frame(LS = unique(flux_output$Metabolite))
  for (i in 1:length(cell.types)) {
    LS.new = data.frame( LS = flux_output$Metabolite[flux_output$Cluster == cell.types[i]],
                         Value = flux_output$Value[flux_output$Cluster == cell.types[i]])
    colnames(LS.new)[2] = cell.types[i]
    LS.df = left_join(LS.df, LS.new, join_by(LS), relationship = "many-to-many")
  }

  rownames(LS.df) = LS.df$LS
  LS.df = LS.df[ , 2:ncol(LS.df)]

  return(LS.df)
}



