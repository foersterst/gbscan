#' Get a classification from the Catalogue of Life
#'
#' Uses the name match service of ChecklistBank for the latest release of
#' the Catalogue of Life (dataset `3LR`).
#'
#' @param genus Genus name (string).
#' @returns A data frame with the columns `name` and `rank`, or `NULL`.
#' @noRd
.clife_classification <- function(genus) {
  resp <- httr2::request("https://api.checklistbank.org") |>
    httr2::req_url_path_append("dataset", "3LR", "match", "nameusage") |>
    httr2::req_url_query(q = genus, rank = "genus") |>
    httr2::req_user_agent("gbscan R package") |>
    httr2::req_throttle(rate = 5, realm = "checklistbank") |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()
  
  body <- httr2::resp_body_json(resp)
  usage <- body$usage
  if (is.null(usage)) {
    return(NULL)
  }
  # The classification can leave out the matched name, so add it.
  cl <- c(usage$classification, list(usage))
  out <- data.frame(
    name = vapply(cl, \(x) as.character(x$name %||% NA_character_)[1], character(1)),
    rank = vapply(cl, \(x) as.character(x$rank %||% NA_character_)[1], character(1))
  )
  out[!duplicated(out$rank), ]
}

#' Cut a feature sequence from a record sequence
#'
#' @param sequence The full record sequence (string).
#' @param location The feature location string.
#' @returns The feature sequence (string), in the 5' to 3' direction of the
#'   feature. `NA` if the location is not supported.
#' @noRd
.extract_feature <- function(sequence, location) {
  m <- .parse_location(location)
  if (is.null(m) || is.na(sequence)) {
    return(NA_character_)
  }
  pieces <- substring(sequence, m[, "start"], m[, "end"])
  minus <- m[, "strand"] < 0
  pieces[minus] <- vapply(pieces[minus], .revcomp, character(1))
  paste(pieces, collapse = "")
}

#' Ranks that the second database can fill
#'
#' @returns Character vector of rank names.
#' @noRd
.fill_rank_columns <- function() {
  c("tribe", "subfamily", "family", "order", "class", "phylum")
}

#' Fill empty ranks from a second database
#'
#' Looks up each genus that has empty ranks, compares the result with the
#' ranks that the record already has, and fills only the empty ranks. Adds
#' the column `filled_ranks`.
#'
#' @param records Records tibble from `.parse_dir()`.
#' @param db Database code: `"gbif"`, `"clife"`, `"itis"`, or `"worms"`.
#' @param dir Root download folder. The lookup cache goes here.
#' @returns The records tibble with filled ranks.
#' @noRd
.fill_ranks <- function(records, db, dir) {
  targets <- .fill_rank_columns()
  check <- c("family", "order", "class", "phylum")
  records$filled_ranks <- NA_character_
  
  empty <- is.na(as.matrix(records[targets]))
  need <- !is.na(records$genus) & rowSums(empty) > 0
  if (!any(need)) {
    cli::cli_inform(c("i" = "No empty ranks to fill."))
    return(records)
  }
  
  lookup <- .lookup_genera(unique(records$genus[need]), db, dir)
  idx <- match(records$genus, lookup$genus)
  use <- need & !is.na(idx) & lookup$found[idx] %in% TRUE
  
  # Do not change a record if a known rank is different in the lookup.
  # A difference means that the genus name possibly is a homonym.
  conflict <- rep(FALSE, nrow(records))
  for (r in check) {
    known <- records[[r]]
    other <- lookup[[r]][idx]
    conflict <- conflict | (!is.na(known) & !is.na(other) & known != other)
  }
  n_conflict <- sum(use & conflict)
  use <- use & !conflict
  
  filled <- matrix(FALSE, nrow(records), length(targets), dimnames = list(NULL, targets))
  for (r in targets) {
    rows <- use & is.na(records[[r]]) & !is.na(lookup[[r]][idx])
    records[[r]][rows] <- lookup[[r]][idx[rows]]
    filled[rows, r] <- TRUE
  }
  any_filled <- rowSums(filled) > 0
  records$filled_ranks[any_filled] <- apply(
    filled[any_filled, , drop = FALSE],
    1,
    \(x) paste(targets[x], collapse = ";")
  )
  
  cli::cli_inform(c(
    "v" = "Filled ranks of {sum(any_filled)} record{?s} from {.val {db}}.",
    "!" = if (n_conflict > 0) {
      "Did not change {n_conflict} record{?s}: the {.val {db}} classification of the genus is different."
    }
  ))
  records
}

#' Summarize feature locations
#'
#' Gets the start, end, strand, partial status, and length of each feature
#' from its location string.
#'
#' @param location Character vector of location strings.
#' @returns A data frame with the columns `start`, `end`, `strand`,
#'   `partial`, and `feature_length`. The length is the sum of the interval
#'   lengths.
#' @noRd
.location_summary <- function(location) {
  parsed <- lapply(location, .parse_location)
  ok <- !vapply(parsed, is.null, logical(1))
  n <- length(location)
  out <- data.frame(
    start = rep(NA_integer_, n),
    end = rep(NA_integer_, n),
    strand = rep(NA_integer_, n),
    partial = grepl("[<>]", location),
    feature_length = rep(NA_integer_, n)
  )
  
  if (any(ok)) {
    p <- parsed[ok]
    out$start[ok] <- vapply(p, \(m) min(m[, c("start", "end")]), integer(1))
    out$end[ok] <- vapply(p, \(m) max(m[, c("start", "end")]), integer(1))
    out$strand[ok] <- vapply(p, \(m) m[1, "strand"], integer(1))
    out$feature_length[ok] <- vapply(
      p,
      \(m) sum(abs(m[, "end"] - m[, "start"]) + 1L),
      integer(1)
    )
  }
  out
}

#' Look up the ranks of genera, with a cache
#'
#' Reads the cache file `dir/ranks_<db>.csv`, looks up the genera that are
#' not in it, and writes the new results to it. A genus that the database
#' does not have is also saved (with `found = FALSE`). Thus the helper does
#' not look it up again.
#'
#' @param genera Character vector of genus names.
#' @param db Database code.
#' @param dir Root download folder.
#' @returns A data frame with the columns `genus`, `found`, and one column
#'   for each rank in `.fill_rank_columns()`.
#' @noRd
.lookup_genera <- function(genera, db, dir) {
  cols <- c("genus", "found", .fill_rank_columns())
  cache_file <- file.path(dir, paste0("ranks_", db, ".csv"))
  cache <- if (file.exists(cache_file)) {
    utils::read.csv(cache_file, colClasses = "character", na.strings = "")
  } else {
    as.data.frame(stats::setNames(replicate(length(cols), character(), simplify = FALSE), cols))
  }
  cache$found <- as.logical(cache$found)
  
  todo <- setdiff(genera, cache$genus)
  if (length(todo) > 0) {
    if (db != "clife") {
      rlang::check_installed("taxize", reason = paste0("to get ranks from ", db, "."))
    }
    new <- vector("list", length(todo))
    cli::cli_progress_bar("Looking up genera in {db}", total = length(todo))
    for (i in seq_along(todo)) {
      new[[i]] <- .lookup_genus(todo[[i]], db)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()
    
    new <- do.call(rbind, new)
    cache <- rbind(cache[cols], new[cols])
    utils::write.csv(cache, cache_file, row.names = FALSE, na = "")
    cli::cli_inform(c("i" = "Saved the {db} lookups in {.file {cache_file}}."))
  }
  cache[cache$genus %in% genera, cols]
}

#' Look up the ranks of one genus
#'
#' Gets the classification of the genus from the database. Uses the
#' Catalogue of Life API for `"clife"` and the taxize package for the other
#' databases. If the lookup fails or gives a different genus name, the
#' result has `found = FALSE`.
#'
#' @param genus Genus name (string).
#' @param db Database code.
#' @returns A one-row data frame with the columns `genus`, `found`, and one
#'   column for each rank in `.fill_rank_columns()`.
#' @noRd
.lookup_genus <- function(genus, db) {
  cl <- tryCatch(
    if (db == "clife") .clife_classification(genus) else .taxize_classification(genus, db),
    error = \(e) NULL
  )
  
  out <- data.frame(genus = genus, found = FALSE)
  for (r in .fill_rank_columns()) {
    out[[r]] <- NA_character_
  }
  if (is.null(cl) || !is.data.frame(cl) || nrow(cl) == 0) {
    return(out)
  }
  
  cl$rank <- tolower(cl$rank)
  if (!identical(cl$name[match("genus", cl$rank)], genus)) {
    return(out)
  }
  out$found <- TRUE
  for (r in .fill_rank_columns()) {
    out[[r]] <- cl$name[match(r, cl$rank)]
  }
  out
}

#' Count GenBank records for a query
#'
#' @param query Entrez query string.
#' @param db NCBI database.
#' @returns Number of records (integer).
#' @noRd
.ncbi_count <- function(query, db = "nuccore") {
  resp <- .ncbi_request(
    "esearch.fcgi",
    list(db = db, term = query, retmax = 0, retmode = "json")
  )
  as.integer(httr2::resp_body_json(resp)$esearchresult$count)
}

#' Download records as XML
#'
#' Gets the records for a set of IDs with EFetch. The helper checks that the
#' XML is complete. NCBI sometimes sends a truncated response with a status
#' of 200. In that case the helper tries again.
#'
#' @param ids Character vector of IDs.
#' @param db NCBI database (`"nuccore"` or `"taxonomy"`).
#' @param tries Number of tries for a truncated response.
#' @returns The XML text.
#' @noRd
.ncbi_fetch <- function(ids, db = "nuccore", tries = 3) {
  params <- list(db = db, id = paste(ids, collapse = ","), retmode = "xml")
  if (db == "nuccore") {
    params$rettype <- "gb"
  }
  end_tag <- if (db == "nuccore") "</GBSet>" else "</TaxaSet>"
  
  for (i in seq_len(tries)) {
    txt <- httr2::resp_body_string(.ncbi_request("efetch.fcgi", params))
    if (grepl(end_tag, substr(txt, max(1, nchar(txt) - 200), nchar(txt)), fixed = TRUE)) {
      return(txt)
    }
    Sys.sleep(10 * i)
  }
  cli::cli_abort(c(
    "NCBI sent incomplete data {tries} time{?s}.",
    "i" = "Run the same command again. The download continues from the last complete batch."
  ))
}

#' Get the NCBI API key
#'
#' Reads the key from the environment variable `ENTREZ_KEY` (the name that
#' the rentrez package also uses) or `NCBI_API_KEY`. Put the key in your
#' `.Renviron` file. Do not put it in a script.
#'
#' @returns The key as a string, or `NULL` if no key is set.
#' @noRd
.ncbi_key <- function() {
  key <- Sys.getenv("ENTREZ_KEY", unset = Sys.getenv("NCBI_API_KEY"))
  if (nzchar(key)) key else NULL
}

#' Send one request to NCBI E-utilities
#'
#' Sends a POST request. POST lets the request hold many record IDs. The
#' request obeys the NCBI rate limit (3 requests per second, or 10 with an
#' API key). It tries again after a temporary failure, with a longer wait
#' each time.
#'
#' @param endpoint E-utility name, for example `"esearch.fcgi"`.
#' @param params Named list of request parameters.
#' @param timeout Maximum time for one try, in seconds.
#' @returns An httr2 response.
#' @noRd
.ncbi_request <- function(endpoint, params, timeout = 300) {
  key <- .ncbi_key()
  params <- c(params, list(tool = "gbscan", api_key = key))
  
  httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils") |>
    httr2::req_url_path_append(endpoint) |>
    httr2::req_body_form(!!!params) |>
    httr2::req_user_agent("gbscan R package") |>
    httr2::req_throttle(rate = if (is.null(key)) 3 else 10, realm = "ncbi") |>
    httr2::req_timeout(timeout) |>
    httr2::req_retry(
      max_tries = 6,
      retry_on_failure = TRUE,
      is_transient = \(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504),
      backoff = \(i) min(5 * 2^(i - 1), 120)
    ) |>
    httr2::req_perform()
}

#' Get all record IDs for a query
#'
#' ESearch gives a maximum of 10,000 IDs in one response. This helper sends
#' as many requests as necessary and increases `retstart` each time. Thus it
#' gets all the IDs, also for queries with more than 50,000 records.
#'
#' @param query Entrez query string.
#' @param db NCBI database.
#' @param page Number of IDs in one request (maximum 10,000).
#' @returns Character vector of unique IDs (NCBI UIDs).
#' @noRd
.ncbi_search_ids <- function(query, db = "nuccore", page = 10000) {
  n <- .ncbi_count(query, db)
  if (n == 0) {
    return(character())
  }
  
  starts <- seq(0, n - 1, by = page)
  ids <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    resp <- .ncbi_request(
      "esearch.fcgi",
      list(
        db = db, term = query, retstart = starts[[i]], retmax = page,
        retmode = "json"
      )
    )
    ids[[i]] <- unlist(httr2::resp_body_json(resp)$esearchresult$idlist)
  }
  unique(as.character(unlist(ids)))
}

#' Parse a GenBank feature location
#'
#' Converts a location string to a matrix of intervals. The helper supports
#' these forms and all combinations of them:
#'
#' * `123..456` and `123` (one base)
#' * `<123..>456` (partial ends)
#' * `complement(...)`
#' * `join(...)` and `order(...)`
#'
#' For `complement(join(a, b))`, the helper reverses the order of the
#' intervals and sets the strand of each one to -1. Thus the result is
#' correct for features on the minus strand that have more than one part.
#'
#' @param location A single location string.
#' @returns An integer matrix with the columns `start`, `end`, and `strand`
#'   (1 or -1), one row for each interval in reading order. `NULL` if the
#'   location points to a different record (`J00194.1:100..202`) or is
#'   between two bases (`123^124`).
#' @noRd
.parse_location <- function(location) {
  if (is.na(location) || grepl("[:^]", location)) {
    return(NULL)
  }
  .parse_location_part(gsub("[<> ]", "", location))
}

#' Parse one part of a location
#'
#' Recursive worker for `.parse_location()`.
#'
#' @param x Location string without `<`, `>`, and spaces.
#' @returns See `.parse_location()`.
#' @noRd
.parse_location_part <- function(x) {
  if (startsWith(x, "complement(")) {
    inner <- .parse_location_part(substr(x, 12, nchar(x) - 1))
    if (is.null(inner)) {
      return(NULL)
    }
    inner <- inner[rev(seq_len(nrow(inner))), , drop = FALSE]
    inner[, "strand"] <- -inner[, "strand"]
    return(inner)
  }
  
  if (grepl("^(join|order)\\(", x)) {
    inner <- sub("^(join|order)\\(", "", substr(x, 1, nchar(x) - 1))
    parts <- lapply(.split_top_level(inner), .parse_location_part)
    if (any(vapply(parts, is.null, logical(1)))) {
      return(NULL)
    }
    return(do.call(rbind, parts))
  }
  
  pos <- regmatches(x, regexec("^([0-9]+)(\\.\\.([0-9]+))?$", x))[[1]]
  if (length(pos) == 0) {
    return(NULL)
  }
  start <- as.integer(pos[[2]])
  end <- if (nzchar(pos[[4]])) as.integer(pos[[4]]) else start
  cbind(start = start, end = end, strand = 1L)
}

#' Read lines from a file if it exists
#'
#' @param file Path to a text file.
#' @returns Character vector. Empty if the file does not exist.
#' @noRd
.read_lines_safe <- function(file) {
  if (file.exists(file)) readLines(file, warn = FALSE) else character()
}

#' Reverse complement a DNA sequence
#'
#' Supports the IUPAC ambiguity codes.
#'
#' @param x A single DNA sequence (string).
#' @returns The reverse complement (string).
#' @noRd
.revcomp <- function(x) {
  x <- chartr("ACGTURYKMBDHVNacgturykmbdhvn", "TGCAAYRMKVHDBNtgcaayrmkvhdbn", x)
  intToUtf8(rev(utf8ToInt(x)))
}

#' Make a name safe for folder names
#'
#' Replaces each run of characters that are not letters or digits with one
#' underscore. For example, `"Tityus serrulatus"` becomes
#' `"Tityus_serrulatus"`.
#'
#' @param x Character vector.
#' @returns Character vector.
#' @noRd
.slug <- function(x) {
  x <- gsub("[^A-Za-z0-9.]+", "_", x)
  gsub("^_+|_+$", "", x)
}

#' Split a string at commas that are not inside parentheses
#'
#' @param x A single string, for example `"complement(1..5),10..20"`.
#' @returns Character vector of the parts.
#' @noRd
.split_top_level <- function(x) {
  chars <- strsplit(x, "")[[1]]
  depth <- cumsum((chars == "(") - (chars == ")"))
  cut <- which(chars == "," & depth == 0)
  starts <- c(1, cut + 1)
  ends <- c(cut - 1, length(chars))
  substring(x, starts, ends)
}

#' Get a classification with taxize
#'
#' @param genus Genus name (string).
#' @param db taxize database name (`"gbif"`, `"itis"`, or `"worms"`).
#' @returns A data frame with the columns `name` and `rank`, or `NULL`.
#' @noRd
.taxize_classification <- function(genus, db) {
  res <- suppressMessages(suppressWarnings(
    taxize::classification(genus, db = db, rows = 1, messages = FALSE)
  ))
  cl <- res[[1]]
  if (!is.data.frame(cl)) NULL else cl[c("name", "rank")]
}

#' Write text to a gzip file safely
#'
#' Writes to a temporary file first, and then renames it. Thus the final
#' file exists only when it is complete.
#'
#' @param text Character string.
#' @param file Path of the final file.
#' @returns `file`, invisibly.
#' @noRd
.write_gz <- function(text, file) {
  tmp <- paste0(file, ".part")
  con <- gzfile(tmp, open = "wb")
  writeLines(text, con, useBytes = TRUE)
  close(con)
  file.rename(tmp, file)
  invisible(file)
}

#' Get the text of the first node that matches an XPath
#'
#' Works on a node set. Gives one value for each node, with `NA` where no
#' node matches.
#'
#' @param nodes An xml2 node set.
#' @param xpath XPath expression, relative to each node.
#' @returns Character vector with the same length as `nodes`.
#' @noRd
.xml_value <- function(nodes, xpath) {
  if (length(nodes) == 0) {
    return(character())
  }
  xml2::xml_text(xml2::xml_find_first(nodes, xpath))
}
