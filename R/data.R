#' Human gene ontology categories
#'
#' A list of functional categories from gene ontology and other sources, and the genes in each category.
#' This had been parsed from a gmt file.
#' The source file is updated monthly so this should be updated too.
#'
#' @format A named list of vectors
#' \describe{
#'   \item{name}{name of functional category }
#'   \item{genes}{vector of genes in the category}
#'   ...
#' }
#' @source \url{http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Human_GO_AllPathways_no_GO_iea_October_01_2019_symbol.gmt}
#' @examples
#' head(human_categories, n = 2)
"human_categories"



#' Mouse gene ontology categories
#'
#' A list of functional categories from gene ontology and other sources, and the genes in each category.
#' This had been parsed from a gmt file.
#' The source file is updated monthly so this should be updated too.
#'
#' @format A named list of vectors
#' \describe{
#'   \item{name}{name of functional category }
#'   \item{genes}{vector of genes in the category}
#'   ...
#' }
#' @source \url{http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/Mouse_GO_AllPathways_no_GO_iea_October_01_2019_symbol.gmt}
#' @examples
#' head(mouse_categories, n=3)
"mouse_categories"



#' Set of mouse genes
#'
#' A dataset containing a set of 100 mouse genes.
#'
#' @format A character vector containing 100 gene names
#'
#' @examples
#' x <- head(genes_1)
#' x
"genes_1"


#' Set of background genes
#'
#' A dataset containing a set of background mouse genes.
#' All protein coding genes in a certain release?
#'
#' @format A character vector containing 29125 gene names
#' @examples
#' head(bg_genes)
"bg_genes"


