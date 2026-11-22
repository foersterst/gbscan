# Exported functions -----------------------------------------------------------

#' Count records for taxa and markers
#'
#' @description
#' `gb_count()` counts the records in GenBank for each combination of taxon
#' and marker. It does not download the records. Use it before
#' [gb_download()] to see how large a download will be.
#'
#' The function prints the result with shortened queries. It returns the
#' full result invisibly. To see the full queries, save the result in an
#' object.
#'
#' @param taxa Character vector of taxon names, for example
#'   `c("Tityus", "Centruroides")`. You can also use NCBI taxonomy IDs in the
#'   format `"txid6887"`.
#' @param markers Character vector of marker names from the dictionary, for
#'   example `c("COI", "16S")`. If `NULL` (the default), the function uses
#'   all markers in the dictionary.
#' @param dictionary Marker dictionary. See [gb_dictionary()].
#'
#' @returns A tibble with one row for each taxon and marker, invisibly. The
#'   columns are `taxon`, `marker`, `query` (the full Entrez query), and
#'   `count`.
#'
#' @details
#' The function makes an Entrez query from the `query`, `min_len`, and
#' `max_len` columns of the dictionary. The `query` column of the result
#' has the full query. You can copy it into the GenBank web page to examine
#' the results.
#'
#' The NCBI rate limit is 3 requests per second. With an API key, the limit
#' is 10 requests per second. To use a key, add the line
#' `ENTREZ_KEY=your_key` to your `.Renviron` file. Then restart R.
#'
#' @family pipeline
#' @export
#' @examplesIf interactive()
#' gb_count(taxa = "Tityus", markers = c("COI", "16S", "28S"))
#'
#' # Save the result to see the full queries
#' counts <- gb_count(taxa = c("Tityus", "Centruroides"), markers = "COI")
#' counts$query
gb_count <- function(
  taxa,
  markers = NULL,
  dictionary = gb_dictionary()
) {
  dictionary <- .check_dictionary(dictionary)
  .check_taxa(taxa)
  dict <- .select_markers(markers %||% dictionary$marker, dictionary)

  grid <- expand.grid(
    taxon = taxa,
    marker = dict$marker,
    stringsAsFactors = FALSE
  )
  grid$query <- .build_query(grid$taxon, grid$marker, dict)
  grid$count <- NA_integer_

  cli::cli_progress_bar("Counting records", total = nrow(grid))
  for (i in seq_len(nrow(grid))) {
    grid$count[[i]] <- .ncbi_count(grid$query[[i]])
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  out <- tibble::as_tibble(grid)

  out |>
    dplyr::mutate(query = stringr::str_trunc(.data$query, width = 25)) |>
    print()

  invisible(out)
}

# Helpers ----------------------------------------------------------------------

#' Build search queries
#'
#' Makes one query for each pair of taxon and marker. The query has the form
#' `"taxon"[Organism] AND (marker terms) AND min:max[SLEN]`.
#'
#' @param taxon Character vector of taxon names.
#' @param marker Character vector of marker names (same length as `taxon`).
#' @param dict Dictionary rows for the markers.
#' @returns Character vector of queries.
#' @noRd
.build_query <- function(taxon, marker, dict) {
  row <- dict[match(marker, dict$marker), ]

  organism <- ifelse(
    grepl("^txid[0-9]+$", taxon),
    paste0(taxon, "[Organism:exp]"),
    paste0("\"", taxon, "\"[Organism]")
  )
  slen <- ifelse(
    is.na(row$min_len) | is.na(row$max_len),
    "",
    sprintf(" AND %d:%d[SLEN]", as.integer(row$min_len), as.integer(row$max_len))
  )
  sprintf("%s AND (%s)%s", organism, row$query, slen)
}

#' Check the taxa argument
#'
#' @param taxa Value to check.
#' @param call Environment to show in error messages.
#' @returns `taxa`, invisibly.
#' @noRd
.check_taxa <- function(taxa, call = rlang::caller_env()) {
  if (!is.character(taxa) || length(taxa) == 0 || anyNA(taxa) || any(!nzchar(taxa))) {
    cli::cli_abort(
      "{.arg taxa} must be a character vector of taxon names.",
      call = call
    )
  }
  invisible(taxa)
}
