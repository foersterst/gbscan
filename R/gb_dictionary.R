# Exported functions -----------------------------------------------------------

#' Get the marker dictionary
#'
#' @description
#' `gb_dictionary()` reads the marker dictionary. The dictionary tells gbscan
#' how to search for each marker and how to identify the marker in a record.
#'
#' The package supplies a default dictionary. It contains common Sanger
#' markers for invertebrate systematics. To change it, do these steps:
#'
#' 1. Write the default dictionary to a CSV file with
#'    `gb_dictionary(path = "markers.csv", write = TRUE)`.
#' 2. Edit the CSV file. You can add, remove, or change rows.
#' 3. Read your file with `gb_dictionary("markers.csv")`.
#' 4. Give the result to the `dictionary` argument of the other functions.
#'
#' You can also edit the tibble in R, for example with [dplyr::mutate()].
#'
#' @param path Path to a CSV file. If `NULL` (the default), the function
#'   reads the dictionary that the package supplies.
#' @param write If `TRUE`, the function writes the default dictionary to
#'   `path`. Use this to get a copy that you can edit. The function does not
#'   replace a file that exists.
#'
#' @returns A tibble with one row for each marker and these columns:
#'   * `marker`: Short marker name. gbscan uses it for folder and file names.
#'   * `genome`: `"mitochondrial"` or `"nuclear"`. gbscan uses only
#'     mitochondrial markers for features in mitochondrial records.
#'   * `description`: Full marker name.
#'   * `query`: Entrez search terms for GenBank.
#'   * `regex`: Regular expression (case-insensitive). gbscan matches it
#'     against the gene, product, note, and definition text of a feature.
#'     If `NA`, gbscan does not assign the marker to features.
#'   * `min_len`, `max_len`: Length limits for the GenBank search.
#'
#'   With `write = TRUE`, the function returns `path` invisibly.
#'
#' @details
#' The order of the rows is important. If more than one regular expression
#' matches a feature, gbscan uses the first matching row.
#'
#' The `mitogenome` row is special. It finds complete mitochondrial genomes.
#' It has no regular expression, so no feature gets the `mitogenome` marker.
#' The parsers cut the separate genes (for example COI) from the genomes.
#'
#' @export
#' @examples
#' dict <- gb_dictionary()
#' dict
#'
#' # Add a marker
#' dict <- dplyr::bind_rows(
#'   dict,
#'   tibble::tibble(
#'     marker = "ND2",
#'     genome = "mitochondrial",
#'     description = "NADH dehydrogenase subunit 2",
#'     query = 'ND2[Gene Name] OR "NADH dehydrogenase subunit 2"[Title]',
#'     regex = "^(nd2|nad2)$|nadh dehydrogenase subunit (2|ii)\\b",
#'     min_len = 100,
#'     max_len = 5000
#'   )
#' )
gb_dictionary <- function(path = NULL, write = FALSE) {
  default <- system.file("extdata", "markers.csv", package = "gbscan")

  if (isTRUE(write)) {
    if (is.null(path)) {
      cli::cli_abort("{.arg path} must be a file path when {.code write = TRUE}.")
    }
    if (file.exists(path)) {
      cli::cli_abort(c(
        "Can't write the dictionary.",
        "x" = "The file {.file {path}} exists.",
        "i" = "Use a different path or remove the file."
      ))
    }
    file.copy(default, path)
    cli::cli_inform(c("v" = "Wrote the marker dictionary to {.file {path}}."))
    return(invisible(path))
  }

  path <- path %||% default
  if (!file.exists(path)) {
    cli::cli_abort("Can't find the file {.file {path}}.")
  }

  dict <- utils::read.csv(
    path,
    na.strings = "",
    stringsAsFactors = FALSE,
    colClasses = c(min_len = "numeric", max_len = "numeric")
  )
  .check_dictionary(tibble::as_tibble(dict))
}

# Helpers ----------------------------------------------------------------------

#' Check a marker dictionary
#'
#' Makes sure that a dictionary has the necessary columns and unique marker
#' names. Stops with an error if it does not.
#'
#' @param dictionary A data frame, usually from [gb_dictionary()].
#' @param arg Argument name to show in error messages.
#' @param call Environment to show in error messages.
#' @returns The dictionary as a tibble.
#' @noRd
.check_dictionary <- function(
  dictionary,
  arg = rlang::caller_arg(dictionary),
  call = rlang::caller_env()
) {
  needed <- c(
    "marker", "genome", "query", "regex", "min_len", "max_len"
  )
  if (!is.data.frame(dictionary)) {
    cli::cli_abort(
      "{.arg {arg}} must be a data frame, not {.obj_type_friendly {dictionary}}.",
      call = call
    )
  }

  missing <- setdiff(needed, names(dictionary))
  if (length(missing) > 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must have the column{?s} {.field {missing}}.",
        "i" = "Use {.fn gb_dictionary} to see the necessary columns."
      ),
      call = call
    )
  }

  dup <- unique(dictionary$marker[duplicated(dictionary$marker)])
  if (length(dup) > 0) {
    cli::cli_abort(
      "{.arg {arg}} has duplicated marker name{?s}: {.val {dup}}.",
      call = call
    )
  }

  bad_genome <- setdiff(unique(dictionary$genome), c("mitochondrial", "nuclear"))
  if (length(bad_genome) > 0) {
    cli::cli_abort(
      c(
        "Column {.field genome} of {.arg {arg}} has incorrect values.",
        "i" = "Use {.val mitochondrial} or {.val nuclear}.",
        "x" = "Incorrect value{?s}: {.val {bad_genome}}."
      ),
      call = call
    )
  }

  tibble::as_tibble(dictionary)
}

#' Select markers from a dictionary
#'
#' Returns the dictionary rows for the given markers, in the given order.
#' Stops with an error if a marker is not in the dictionary.
#'
#' @param markers Character vector of marker names.
#' @param dictionary A marker dictionary.
#' @param call Environment to show in error messages.
#' @returns A tibble with one row for each marker.
#' @noRd
.select_markers <- function(markers, dictionary, call = rlang::caller_env()) {
  if (!is.character(markers) || length(markers) == 0) {
    cli::cli_abort(
      "{.arg markers} must be a character vector with one or more names.",
      call = call
    )
  }

  unknown <- setdiff(markers, dictionary$marker)
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        "{.val {unknown}} {?is/are} not in the marker dictionary.",
        "i" = "Available markers: {.val {dictionary$marker}}.",
        "i" = "Add markers to the dictionary with {.fn gb_dictionary}."
      ),
      call = call
    )
  }

  dictionary[match(unique(markers), dictionary$marker), ]
}
