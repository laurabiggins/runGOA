#' Clean text
#'
#' Remove erroneous text characters.
#'
#' @param text character vector
#' @param chars character vector of characters to remove from text
#' @param remove_empty remove empty elements from the vector.
#' @param remove_dup remove duplicated elements from the vector.
#' @param upper_case convert to upper case
#' @return A character vector
#' @examples
#' clean_text("\thello_world")
#' clean_text("\thello_world", upper_case = FALSE)
#' @export

clean_text <- function(text,
                       chars = c("\t", "\r", "\n", "\\", "\ ", "\""),
                       remove_empty = TRUE,
                       remove_dup = TRUE,
                       upper_case = TRUE
) {
  char_pattern <- paste(chars, collapse = "|")
  char_pattern <- paste0("[", char_pattern, "]")
  cleaned <- gsub(pattern = char_pattern, replacement = "", x = as.character(text))
  if (remove_dup) cleaned <- remove_duplicates(cleaned)
  if (upper_case) cleaned <- toupper(cleaned)
  ifelse(remove_empty, return(cleaned[sapply(cleaned, nchar) > 0]), return(cleaned))
}


#' remove_duplicates
#'
#' remove duplicates from vector or data frame
#'
#' @param x vector or data frame to remove duplicates from
#' @param column If a data frame is passed in, the column that is used
#' for deduplicating (default = 1)
#'
#' @return A character vector or data frame.
#' @examples
#' hello_string <- c("h", "e", "l", "l", "o")
#' hello_string
#' remove_duplicates(hello_string)
#'
#' df <- data.frame(hello_string, 1:5)
#' df
#' remove_duplicates(df)
#' @export

remove_duplicates <- function(x, column = 1) {
  if (is.data.frame(x) | is.matrix(x)) {
    return(x[!duplicated(x[,column]), ])
  } else if (is.list(x)) {
    stop("x must be a vector, data frame or matrix")
  } else if (is.vector(x)) {
    return(unique(x))
  } else{
    warning("data must be a vector, data frame or matrix")
  }
}
