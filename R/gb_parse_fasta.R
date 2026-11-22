# Exported functions -----------------------------------------------------------

#' Write the sequences of downloaded records to FASTA files
#'
#' @description
#' `gb_parse_fasta()` reads the files that [gb_download()] saved in `dir`.
#' It writes the sequences of each marker to a FASTA file. Each marker gets
#' its own folder: `out_dir/<marker>/<marker>.fasta`.
#'
#' @inheritParams gb_parse_meta
#' @param out_dir Folder for the FASTA files. The function makes the folder
#'   if it does not exist.
#' @param markers Character vector of marker names. If `NULL` (the default),
#'   the function writes all markers that it finds.
#' @param keep Optional. A data frame with a `feature_id` column, usually the
#'   `features` table from [gb_parse_meta()] after you filtered it. If you
#'   supply it, the function writes only these features. The function can
#'   also use the other columns of `keep` in the FASTA headers.
#' @param region `"feature"` (the default) writes only the sequence of the
#'   feature. For example, it cuts the COI gene from a mitogenome.
#'   `"record"` writes the full sequence of the record.
#' @param header Character vector of column names for the FASTA headers.
#'   The function joins the values with `"|"` and replaces spaces with
#'   underscores. The columns can come from `keep`, from the records table,
#'   or from the features table of [gb_parse_meta()]. The default is
#'   `c("accession", "species")`.
#' @param overwrite If `FALSE` (the default), the function stops if a FASTA
#'   file exists.
#'
#' @returns A tibble with one row for each marker, invisibly. The columns are
#'   `marker`, `n_sequences`, and `file`.
#'
#' @details
#' The function skips features without a marker and features without a
#' sequence. A feature has no sequence if its location points to a
#' different record.
#'
#' @family pipeline
#' @export
#' @examplesIf interactive()
#' # All sequences of all markers
#' gb_parse_fasta("tityus_data", out_dir = "fasta")
#'
#' # Only the longest COI sequence of each species
#' meta <- gb_parse_meta("tityus_data")
#' best <- meta$features |>
#'   dplyr::filter(marker == "COI", feature_length >= 500) |>
#'   dplyr::left_join(meta$records, by = "accession") |>
#'   dplyr::filter(!is.na(species)) |>
#'   dplyr::arrange(dplyr::desc(feature_length)) |>
#'   dplyr::distinct(species, .keep_all = TRUE)
#'
#' gb_parse_fasta(
#'   "tityus_data",
#'   out_dir = "fasta",
#'   keep = best,
#'   header = c("accession", "species", "country"),
#'   overwrite = TRUE
#' )
gb_parse_fasta <- function(
  dir,
  out_dir = "fasta",
  markers = NULL,
  keep = NULL,
  region = c("feature", "record"),
  header = c("accession", "species"),
  overwrite = FALSE,
  dictionary = gb_dictionary()
) {
  region <- rlang::arg_match(region)
  dictionary <- .check_dictionary(dictionary)
  if (!is.null(keep) && !("feature_id" %in% names(keep))) {
    cli::cli_abort(c(
      "{.arg keep} must have a {.field feature_id} column.",
      "i" = "Use the {.field features} table of {.fn gb_parse_meta}."
    ))
  }

  parsed <- .parse_dir(dir, dictionary, region = region)
  feats <- parsed$features
  feats <- feats[!is.na(feats$marker) & !is.na(feats$sequence) & nzchar(feats$sequence), ]
  if (!is.null(markers)) {
    feats <- feats[feats$marker %in% markers, ]
  }
  if (!is.null(keep)) {
    feats <- feats[feats$feature_id %in% keep$feature_id, ]
  }
  if (nrow(feats) == 0) {
    cli::cli_abort("No sequences to write. Check {.arg markers} and {.arg keep}.")
  }

  labels <- .fasta_headers(feats, parsed$records, keep, header)

  found <- unique(feats$marker)
  files <- file.path(out_dir, .slug(found), paste0(.slug(found), ".fasta"))
  exists <- file.exists(files)
  if (!overwrite && any(exists)) {
    cli::cli_abort(c(
      "{sum(exists)} FASTA file{?s} exist{?s/}: {.file {files[exists]}}.",
      "i" = "Use {.code overwrite = TRUE} or a different {.arg out_dir}."
    ))
  }

  n <- integer(length(found))
  for (i in seq_along(found)) {
    rows <- feats$marker == found[[i]]
    dir.create(dirname(files[[i]]), recursive = TRUE, showWarnings = FALSE)
    .write_fasta(labels[rows], feats$sequence[rows], files[[i]])
    n[[i]] <- sum(rows)
  }

  out <- tibble::tibble(marker = found, n_sequences = n, file = files)
  cli::cli_inform(c("v" = "Wrote {sum(n)} sequence{?s} for {length(found)} marker{?s} to {.path {out_dir}}."))
  invisible(out)
}

# Helpers ----------------------------------------------------------------------

#' Make FASTA headers
#'
#' Gets the values of the header columns for each feature. The helper looks
#' for each column first in `keep`, then in the features table, and then in
#' the records table.
#'
#' @param feats Features tibble (filtered).
#' @param records Records tibble.
#' @param keep Data frame with `feature_id`, or `NULL`.
#' @param header Character vector of column names.
#' @param call Environment to show in error messages.
#' @returns Character vector of headers, without `>`.
#' @noRd
.fasta_headers <- function(feats, records, keep, header, call = rlang::caller_env()) {
  rec_idx <- match(feats$accession, records$accession)
  keep_idx <- if (is.null(keep)) NULL else match(feats$feature_id, keep$feature_id)

  values <- lapply(header, function(col) {
    if (!is.null(keep) && col %in% names(keep)) {
      as.character(keep[[col]][keep_idx])
    } else if (col %in% names(feats)) {
      as.character(feats[[col]])
    } else if (col %in% names(records)) {
      as.character(records[[col]][rec_idx])
    } else {
      cli::cli_abort(
        c(
          "Can't find the header column {.field {col}}.",
          "i" = "Use a column of {.arg keep} or of the tables from {.fn gb_parse_meta}."
        ),
        call = call
      )
    }
  })

  labels <- do.call(paste, c(values, sep = "|"))
  gsub("\\s+", "_", labels)
}

#' Write sequences to a FASTA file
#'
#' Writes each sequence on one line.
#'
#' @param labels Character vector of headers, without `>`.
#' @param sequences Character vector of sequences.
#' @param file Path of the file.
#' @returns `file`, invisibly.
#' @noRd
.write_fasta <- function(labels, sequences, file) {
  lines <- as.vector(rbind(paste0(">", labels), sequences))
  writeLines(lines, file)
  invisible(file)
}
