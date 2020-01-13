# This is the code used to process the gene set files and save them
# as usable data files for the package

process_GMT <- function(file, min_genes = 3, max_genes = 30000) {

  gmt_file <- data.table::fread(
    file,
    sep = "\n",
    header = FALSE,
    data.table = FALSE
  )[, 1]

  categories <- strsplit(gmt_file, "\t")

  names(categories) <- sapply(categories, `[[`, 1)

  genes_in_categories <- lapply(categories, `[`, -(1:2))

  genes_in_categories <- sapply(genes_in_categories, toupper)

  category_lengths <- sapply(genes_in_categories, length)

  genes_in_categories[category_lengths >= min_genes & category_lengths <= max_genes]
}

human_categories <- process_GMT(file = "Human_GO_AllPathways_no_GO_iea_January_01_2020_symbol.gmt.txt", min_genes = 5)
mouse_categories <- process_GMT(file = "Mouse_GO_AllPathways_no_GO_iea_January_01_2020_symbol.gmt.txt", min_genes = 5)


# usethis::use_data(human_go, overwrite = TRUE) # 3.3MB
usethis::use_data(human_categories, compress = "xz", overwrite = TRUE) # 1.8MB but a lot slower to save
usethis::use_data(mouse_categories, compress = "xz", overwrite = TRUE) # 2.1MB but a lot slower to save
# usethis::use_data(mouse_categories, overwrite = TRUE) # 3.7MB

# using example genes from my GOcategoryStats package
devtools::load_all("M:/GOcategoryStats/")
usethis::use_data(genes_1)
usethis::use_data(bg_genes)
