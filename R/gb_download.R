# Exported functions -----------------------------------------------------------

#' Download records for taxa and markers
#'
#' @description
#' `gb_download()` downloads the records for each combination of taxon and
#' marker and saves them in `dir`. Use [gb_parse_meta()] and
#' [gb_parse_fasta()] to read the saved files.
#'
#' The function is made for large downloads (more than 50,000 records for
#' one taxon and marker). It downloads the records in batches. If the
#' download stops (for example, because of a network problem), run the same
#' command again. The function continues from the last complete batch.
#'
#' @inheritParams gb_count
#' @param markers Character vector of marker names from the dictionary.
#' @param dir Folder for the downloaded files. The function makes the folder
#'   if it does not exist.
#' @param batch Number of records in one download request. The default is
#'   500. Use a smaller value if NCBI often stops a request. The maximum
#'   value is 10,000.
#' @param taxonomy If `TRUE` (the default), the function also downloads the
#'   NCBI Taxonomy data for the records. [gb_parse_meta()] uses these data
#'   to find the taxonomic ranks (order, family, subfamily, tribe,
#'   genus, and species).
#'
#' @returns A tibble with one row for each taxon and marker, invisibly. The
#'   columns are `taxon`, `marker`, `n_records` (number of records
#'   that the search found), and `path` (the folder with the files).
#'
#' @details
#' ## Folder structure
#'
#' The function saves the files in this structure:
#'
#' ```
#' dir/
#'   genbank/
#'     <marker>/<taxon>/
#'       query.txt            the query, the date, and the number of records
#'       ids.txt              all record IDs that the search found
#'       done.txt             the IDs that are downloaded
#'       records_00001.xml.gz the records (GBSeq XML), one file per batch
#'     taxonomy/
#'       taxonomy_00001.xml.gz
#' ```
#'
#' The function compresses the GenBank files with gzip to save disk space.
#'
#' ## Resume and update
#'
#' The function saves the record IDs when it first runs the search. When you
#' run it again, it downloads only the IDs that are not in `done.txt`. If
#' the query changed (for example, because you edited the dictionary), the
#' function runs the search again and downloads only the new IDs.
#'
#' To get records that were added to GenBank after the first search, remove
#' the file `ids.txt`. Then run the function again.
#'
#' A record can be in more than one folder. For example, a mitogenome can be
#' in the `COI` folder and in the `mitogenome` folder. The parsers remove
#' these duplicates.
#'
#' @family pipeline
#' @export
#' @examplesIf interactive()
#' # COI and complete mitogenomes of the scorpion genus Tityus
#' gb_download(
#'   taxa = "Tityus",
#'   markers = c("COI", "mitogenome"),
#'   dir = "tityus_data"
#' )
gb_download <- function(
  taxa,
  markers,
  dir,
  batch = 500,
  taxonomy = TRUE,
  dictionary = gb_dictionary()
) {
  dictionary <- .check_dictionary(dictionary)
  .check_taxa(taxa)
  dict <- .select_markers(markers, dictionary)
  if (!rlang::is_scalar_integerish(batch) || batch < 1 || batch > 10000) {
    cli::cli_abort("{.arg batch} must be a whole number from 1 to 10000.")
  }
  if (!rlang::is_string(dir)) {
    cli::cli_abort("{.arg dir} must be a folder path (one string).")
  }

  grid <- expand.grid(taxon = taxa, marker = dict$marker, stringsAsFactors = FALSE)
  grid$query <- .build_query(grid$taxon, grid$marker, dict)
  grid$path <- file.path(dir, "genbank", .slug(grid$marker), .slug(grid$taxon))
  grid$n_records <- NA_integer_

  for (i in seq_len(nrow(grid))) {
    cli::cli_h2("{grid$taxon[[i]]} | {grid$marker[[i]]}")
    dir.create(grid$path[[i]], recursive = TRUE, showWarnings = FALSE)
    grid$n_records[[i]] <- .download_genbank(grid$query[[i]], grid$path[[i]], batch, dir)
  }

  if (isTRUE(taxonomy)) {
    cli::cli_h2("NCBI Taxonomy")
    .download_taxonomy(dir)
  }

  invisible(tibble::as_tibble(grid[c("taxon", "marker", "n_records", "path")]))
}

# Helpers ----------------------------------------------------------------------

#' Download GenBank records for one query
#'
#' Gets the record IDs (or reads the saved IDs), and then downloads the IDs
#' that are not in `done.txt`, in batches. After each batch, the helper adds
#' the IDs to `done.txt` and the NCBI taxonomy IDs to the taxonomy list.
#'
#' @param query Entrez query.
#' @param path Folder for the files of this query.
#' @param batch Number of records in one request.
#' @param dir Root download folder. The taxonomy list is in
#'   `dir/genbank/taxonomy`.
#' @returns Number of records that the search found.
#' @noRd
.download_genbank <- function(query, path, batch, dir) {
  ids_file <- file.path(path, "ids.txt")
  query_file <- file.path(path, "query.txt")
  done_file <- file.path(path, "done.txt")
  tax_dir <- file.path(dir, "genbank", "taxonomy")
  dir.create(tax_dir, recursive = TRUE, showWarnings = FALSE)

  old_query <- .read_lines_safe(query_file)[1]
  if (file.exists(ids_file) && identical(old_query, query)) {
    ids <- .read_lines_safe(ids_file)
    cli::cli_inform(c("i" = "Using {length(ids)} saved record ID{?s}."))
  } else {
    cli::cli_progress_step("Searching GenBank")
    ids <- .ncbi_search_ids(query)
    writeLines(ids, ids_file)
    writeLines(c(query, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(ids)), query_file)
    cli::cli_progress_done()
    cli::cli_inform(c("i" = "Found {length(ids)} record{?s}."))
  }

  todo <- setdiff(ids, .read_lines_safe(done_file))
  if (length(todo) == 0) {
    cli::cli_inform(c("v" = "All records are downloaded."))
    return(length(ids))
  }

  chunks <- split(todo, ceiling(seq_along(todo) / batch))
  n_old <- length(list.files(path, pattern = "^records_[0-9]+\\.xml\\.gz$"))
  n_got <- 0L

  cli::cli_progress_bar(
    format = "Downloading {cli::pb_bar} {cli::pb_current}/{cli::pb_total} batches | ETA {cli::pb_eta}",
    total = length(chunks)
  )
  for (j in seq_along(chunks)) {
    txt <- .ncbi_fetch(chunks[[j]], db = "nuccore")
    file <- file.path(path, sprintf("records_%05d.xml.gz", n_old + j))
    .write_gz(txt, file)

    taxids <- unique(sub("taxon:", "", regmatches(txt, gregexpr("taxon:[0-9]+", txt))[[1]]))
    cat(taxids, file = file.path(tax_dir, "taxids.txt"), sep = "\n", append = TRUE)
    cat(chunks[[j]], file = done_file, sep = "\n", append = TRUE)

    n_got <- n_got + lengths(regmatches(txt, gregexpr("<GBSeq>", txt, fixed = TRUE)))
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  if (n_got < length(todo)) {
    cli::cli_inform(c(
      "!" = "NCBI sent {n_got} of {length(todo)} requested record{?s}.",
      "i" = "NCBI does not send records that are removed or suppressed."
    ))
  }
  cli::cli_inform(c("v" = "Saved the records in {.path {path}}."))
  length(ids)
}

#' Download NCBI Taxonomy records
#'
#' Downloads the taxonomy records for the taxonomy IDs in
#' `dir/genbank/taxonomy/taxids.txt` that are not downloaded yet.
#'
#' @param dir Root download folder.
#' @param batch Number of taxonomy IDs in one request.
#' @returns Number of new taxonomy IDs, invisibly.
#' @noRd
.download_taxonomy <- function(dir, batch = 200) {
  tax_dir <- file.path(dir, "genbank", "taxonomy")
  done_file <- file.path(tax_dir, "done.txt")
  wanted <- unique(.read_lines_safe(file.path(tax_dir, "taxids.txt")))
  todo <- setdiff(wanted[nzchar(wanted)], .read_lines_safe(done_file))

  if (length(todo) == 0) {
    cli::cli_inform(c("v" = "Taxonomy data are complete."))
    return(invisible(0L))
  }

  chunks <- split(todo, ceiling(seq_along(todo) / batch))
  n_old <- length(list.files(tax_dir, pattern = "^taxonomy_[0-9]+\\.xml\\.gz$"))
  cli::cli_progress_bar("Downloading taxonomy", total = length(chunks))
  for (j in seq_along(chunks)) {
    txt <- .ncbi_fetch(chunks[[j]], db = "taxonomy")
    .write_gz(txt, file.path(tax_dir, sprintf("taxonomy_%05d.xml.gz", n_old + j)))
    cat(chunks[[j]], file = done_file, sep = "\n", append = TRUE)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  cli::cli_inform(c("v" = "Downloaded taxonomy data for {length(todo)} tax{?on/a}."))
  invisible(length(todo))
}
