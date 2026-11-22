# Exported functions -----------------------------------------------------------

#' Read the metadata of downloaded records
#'
#' @description
#' `gb_parse_meta()` reads the files that [gb_download()] saved in `dir`. It
#' returns the metadata in two tibbles: one for records and one for
#' features.
#'
#' @param dir Folder that [gb_download()] used.
#' @param dictionary Marker dictionary. See [gb_dictionary()]. The function
#'   uses it to assign a marker to each feature. If you change the
#'   dictionary, run the function again. You do not have to download again.
#' @param ranks Source of the taxonomic ranks. `"ncbi"` (the default) uses
#'   only the NCBI Taxonomy data. The other values also use a second database
#'   to fill the empty ranks (`tribe` to `phylum`):
#'   * `"gbif"`: GBIF Backbone Taxonomy.
#'   * `"clife"`: Catalogue of Life.
#'   * `"itis"`: Integrated Taxonomic Information System.
#'   * `"worms"`: World Register of Marine Species.
#'
#'   All values except `"ncbi"` and `"clife"` need the taxize package. See
#'   the section "Fill empty ranks".
#'
#' @returns A list with two tibbles:
#'
#' `records` has one row for each record (accession) and these columns:
#'   * `accession`: GenBank accession with version.
#'   * `organism`: Organism name in the record.
#'   * `species`, `genus`, `tribe`, `subfamily`, `family`, `order`,
#'     `class`, `phylum`: Taxonomic ranks. `NA` if the rank is not
#'     available.
#'   * `taxid`: NCBI taxonomy ID.
#'   * `lineage`: Full lineage text.
#'   * `definition`, `length`, `mol_type`, `organelle`.
#'   * `specimen_voucher`, `isolate`, `country`, `lat_lon`,
#'     `collection_date`.
#'   * `authors`, `ref_title`, `journal`: First reference of the record.
#'   * `create_date`, `update_date`.
#'   * `filled_ranks`: The ranks that came from the `ranks` database,
#'     separated by `;`. `NA` if no rank came from it. This column is only
#'     present if `ranks` is not `"ncbi"`.
#'
#' `features` has one row for each gene-level feature and these columns:
#'   * `feature_id`: Unique feature identifier. Use it to filter features
#'     before you use [gb_parse_fasta()].
#'   * `accession`.
#'   * `marker`: Marker name from the dictionary. `NA` if no marker
#'     matches.
#'   * `matched_field`: The field that matched the marker (`"gene"`,
#'     `"product"`, `"note"`, or `"definition"`).
#'   * `feature_key`: GenBank feature type (`"CDS"`, `"rRNA"`, and so on).
#'   * `gene`, `product`, `note`: Feature qualifiers.
#'   * `location`, `start`, `end`, `strand`, `partial`, `feature_length`.
#'
#' @details
#' ## Features
#'
#' The function reads features of these types: `CDS`, `rRNA`, `misc_RNA`,
#' `mRNA`, `misc_feature`, and `gene`. If a record has none of them, the
#' function makes one feature for the full record. If a record has more than
#' one feature for the same marker (for example, a `gene` and a `CDS`), the
#' function keeps the most specific one.
#'
#' To find features without a marker, use `dplyr::filter(features,
#' is.na(marker))`. Then add a row or change a regular expression in the
#' dictionary.
#'
#' ## Taxonomic ranks
#'
#' The ranks come from the NCBI Taxonomy data that [gb_download()]
#' saved. The `species` column has the species-rank name of the NCBI
#' taxon. Thus a subspecies gets the name of its species, and a record that
#' is identified only to genus gets `NA`. Names in open nomenclature, for
#' example "Tityus sp. 1" and "Tityus sp. 2", stay different species. The
#' `organism` column always has the name as it is in the record.
#'
#' If the NCBI Taxonomy data are not available, the function gets the genus and species from the organism name. It also gets family
#' (`-idae`), subfamily (`-inae`), and tribe (`-ini`) from the lineage text.
#' This method is correct for most animal names.
#'
#' ## Fill empty ranks
#'
#' NCBI Taxonomy often has no tribe or subfamily for invertebrates. Use the
#' `ranks` argument to get these ranks from a second database. The function
#' does these steps:
#'
#' 1. It finds the records that have a genus and one or more empty ranks.
#' 2. It looks up each genus one time in the database.
#' 3. It compares the lookup result with the ranks that the record already
#'    has (family, order, class, and phylum). If one of them is different,
#'    the genus name possibly belongs to a different group (a homonym). Then
#'    the function does not change the record.
#' 4. It fills only the empty ranks. It does not replace a rank that
#'    exists.
#'
#' The function saves the lookup results in `dir`, in the file
#' `ranks_<database>.csv`. When you run the function again, it uses the saved
#' results and looks up only new genera. To look up all genera again, remove
#' the file.
#'
#' Different databases can use different classifications. Use the
#' `filled_ranks` column to find the ranks that came from the second
#' database.
#'
#' ## Duplicates
#'
#' A record can be in more than one download folder. The function keeps one
#' copy of each record and each feature.
#'
#' @family pipeline
#' @export
#' @examplesIf interactive()
#' meta <- gb_parse_meta("tityus_data")
#' meta$records
#' meta$features
#'
#' # Fill empty tribe and subfamily ranks from GBIF
#' meta <- gb_parse_meta("tityus_data", ranks = "gbif")
#' dplyr::count(meta$records, subfamily, filled_ranks)
#'
#' # Join the two tables
#' dplyr::left_join(meta$features, meta$records, by = "accession")
gb_parse_meta <- function(
  dir,
  ranks = c("ncbi", "gbif", "clife", "itis", "worms"),
  dictionary = gb_dictionary()
) {
  ranks <- rlang::arg_match(ranks)
  dictionary <- .check_dictionary(dictionary)
  out <- .parse_dir(dir, dictionary, region = NULL)
  if (ranks != "ncbi") {
    out$records <- .fill_ranks(out$records, ranks, dir)
  }

  n_rec <- nrow(out$records)
  n_feat <- nrow(out$features)
  n_none <- sum(is.na(out$features$marker))
  cli::cli_inform(c(
    "v" = "Read {n_rec} record{?s} and {n_feat} feature{?s}.",
    "i" = if (n_none > 0) "{n_none} feature{?s} {?has/have} no marker."
  ))
  out
}

# Helpers ----------------------------------------------------------------------

#' Add taxonomic ranks to the records table
#'
#' Uses the NCBI Taxonomy data. For records without taxonomy data, gets the
#' ranks from the organism name and the lineage text. Then puts the columns
#' in a fixed order.
#'
#' @param records Records tibble.
#' @param taxonomy Tibble from `.read_taxonomy()`, or `NULL`.
#' @returns The records tibble with the rank columns.
#' @noRd
.add_ranks <- function(records, taxonomy) {
  ranks <- .rank_columns()
  for (r in ranks) {
    records[[r]] <- NA_character_
  }

  idx <- if (is.null(taxonomy)) {
    rep(NA_integer_, nrow(records))
  } else {
    match(records$taxid, taxonomy$taxid)
  }
  has_tax <- !is.na(idx)

  for (r in ranks) {
    records[[r]][has_tax] <- taxonomy[[r]][idx[has_tax]]
  }

  if (any(!has_tax)) {
    guess <- .guess_ranks(records$organism[!has_tax], records$lineage[!has_tax])
    for (r in names(guess)) {
      records[[r]][!has_tax] <- guess[[r]]
    }
  }

  # Some taxa have no genus rank in NCBI Taxonomy. Use the first word of the
  # organism name for these.
  miss_genus <- is.na(records$genus)
  records$genus[miss_genus] <- sub(" .*", "", records$organism[miss_genus])

  first <- c(
    "accession", "organism", ranks, "taxid", "lineage", "definition",
    "length", "mol_type", "organelle", "specimen_voucher", "isolate",
    "country", "lat_lon", "collection_date", "authors", "ref_title",
    "journal", "create_date", "update_date"
  )
  records[c(intersect(first, names(records)), setdiff(names(records), first))]
}

#' Assign markers to features
#'
#' Matches the regular expressions of the dictionary against the text of
#' each feature. The helper examines the fields in this order: gene,
#' product, note, and definition. In each field, it examines the markers in
#' dictionary order. The first match is the result.
#'
#' Mitochondrial records (organelle "mitochondrion") can get only
#' mitochondrial markers. Other records get nuclear markers first. Thus
#' "large subunit ribosomal RNA" becomes 16S in a mitochondrial record and
#' 28S in a nuclear record.
#'
#' @param gene,product,note,definition Character vectors with the same
#'   length. Use `NA` where a value is missing.
#' @param organelle Character vector (same length).
#' @param dictionary Marker dictionary.
#' @returns A list with the character vectors `marker` and `field`.
#' @noRd
.assign_markers <- function(gene, product, note, definition, organelle, dictionary) {
  dict <- dictionary[!is.na(dictionary$regex), ]
  n <- length(gene)
  marker <- rep(NA_character_, n)
  field <- rep(NA_character_, n)
  is_mt_record <- grepl("mitochondri", organelle, ignore.case = TRUE)
  is_mt_marker <- dict$genome == "mitochondrial"
  marker_order <- c(which(!is_mt_marker), which(is_mt_marker))

  fields <- list(gene = gene, product = product, note = note, definition = definition)
  for (f in names(fields)) {
    txt <- tolower(trimws(fields[[f]]))
    for (k in marker_order) {
      todo <- is.na(marker) & !is.na(txt)
      if (!is_mt_marker[[k]]) {
        todo <- todo & !is_mt_record
      }
      if (!any(todo)) {
        next
      }
      hit <- todo
      hit[todo] <- grepl(dict$regex[[k]], txt[todo], perl = TRUE)
      marker[hit] <- dict$marker[[k]]
      field[hit] <- f
    }
  }
  list(marker = marker, field = field)
}

#' Remove redundant features
#'
#' A record often has a `gene` feature and a `CDS` feature for the same
#' gene. This helper keeps the most specific feature type for each record
#' and marker. The order is `CDS`, `rRNA`, `misc_RNA`, `mRNA`,
#' `misc_feature`, `gene`, `source`. For features without a marker, it
#' removes `gene` features that have the same location as another feature.
#'
#' @param features Features tibble with `accession`, `marker`,
#'   `feature_key`, and `location`.
#' @returns The filtered tibble.
#' @noRd
.drop_redundant_features <- function(features) {
  priority <- c("CDS", "rRNA", "misc_RNA", "mRNA", "misc_feature", "gene", "source")
  rank <- match(features$feature_key, priority)

  group <- paste(features$accession, features$marker)
  best <- stats::ave(rank, group, FUN = min)
  keep_assigned <- !is.na(features$marker) & rank == best

  dup_loc <- duplicated(paste(features$accession, features$location)) |
    duplicated(paste(features$accession, features$location), fromLast = TRUE)
  keep_unassigned <- is.na(features$marker) & !(features$feature_key == "gene" & dup_loc)

  features[keep_assigned | keep_unassigned, ]
}

#' Find the downloaded files in a folder
#'
#' @param dir Root download folder.
#' @param call Environment to show in error messages.
#' @returns A list with the elements `genbank` and `taxonomy`. Each
#'   element is a character vector of file paths.
#' @noRd
.find_files <- function(dir, call = rlang::caller_env()) {
  if (!dir.exists(dir)) {
    cli::cli_abort("Can't find the folder {.path {dir}}.", call = call)
  }
  gb_dir <- file.path(dir, "genbank")
  all_gb <- list.files(gb_dir, pattern = "\\.xml\\.gz$", recursive = TRUE, full.names = TRUE)
  is_tax <- grepl("^taxonomy_", basename(all_gb))

  files <- list(genbank = all_gb[!is_tax], taxonomy = all_gb[is_tax])
  if (length(files$genbank) == 0) {
    cli::cli_abort(
      c(
        "Can't find downloaded records in {.path {dir}}.",
        "i" = "Use {.fn gb_download} to download records first."
      ),
      call = call
    )
  }
  files
}

#' Complete a features table
#'
#' Adds location data, assigns markers, removes redundant features, and
#' adds the `feature_id` column.
#'
#' @param features Tibble with `accession`, `feature_key`, `gene`, `product`,
#'   `note`, and `location`.
#' @param records Records tibble (for `organelle` and `definition`).
#' @param dictionary Marker dictionary.
#' @returns The completed features tibble.
#' @noRd
.finish_features <- function(features, records, dictionary) {
  idx <- match(features$accession, records$accession)
  n_per_record <- table(features$accession)[features$accession]
  definition <- ifelse(n_per_record == 1, records$definition[idx], NA_character_)

  assigned <- .assign_markers(
    gene = features$gene,
    product = features$product,
    note = features$note,
    definition = definition,
    organelle = records$organelle[idx],
    dictionary = dictionary
  )
  features$marker <- assigned$marker
  features$matched_field <- assigned$field
  features <- .drop_redundant_features(features)

  features <- dplyr::bind_cols(features, .location_summary(features$location))
  features$feature_id <- paste(features$accession, features$feature_key, features$location, sep = "|")

  features[c(
    "feature_id", "accession", "marker", "matched_field",
    "feature_key", "gene", "product", "note", "location", "start", "end",
    "strand", "partial", "feature_length"
  )]
}

#' Get the authors of the first reference of each record
#'
#' @param recs Node set of GBSeq records.
#' @param accession Accessions of the records, in the same order.
#' @returns Character vector of author lists (comma-separated).
#' @noRd
.first_ref_authors <- function(recs, accession) {
  au <- xml2::xml_find_all(recs, "GBSeq_references/GBReference[1]/GBReference_authors/GBAuthor")
  if (length(au) == 0) {
    return(rep(NA_character_, length(accession)))
  }
  owner <- .xml_value(au, "ancestor::GBSeq/GBSeq_accession-version")
  joined <- tapply(xml2::xml_text(au), owner, paste, collapse = ", ")
  unname(joined[accession])
}

#' Guess taxonomic ranks without taxonomy data
#'
#' Gets the genus (first word) and species (first two words) from the
#' organism name. Names in open nomenclature (`sp.`, `cf.`, `aff.`, `nr.`,
#' `gr.`, with or without the dot) keep the full organism name. Thus
#' "Tityus sp. 1" and "Tityus sp. 2" stay different species. Gets family,
#' subfamily, and tribe from the lineage text by the suffixes `-idae`,
#' `-inae`, and `-ini`. The suffixes are correct for
#' animals. They are not correct for all plants, fungi, or bacteria.
#'
#' @param organism Character vector of organism names.
#' @param lineage Character vector of lineage texts (`;`-separated).
#' @returns A list of character vectors: `genus`, `species`, `tribe`,
#'   `subfamily`, and `family`.
#' @noRd
.guess_ranks <- function(organism, lineage) {
  words <- strsplit(organism, " ", fixed = TRUE)
  open_name <- grepl(" (sp|spp|cf|aff|nr|gr)\\.?( |$)", organism)
  species <- vapply(words, \(w) if (length(w) >= 2) paste(w[1:2], collapse = " ") else NA_character_, character(1))
  species[open_name] <- organism[open_name]

  parts <- strsplit(ifelse(is.na(lineage), "", lineage), "; ?")
  by_suffix <- function(suffix) {
    vapply(parts, function(p) {
      hit <- grep(paste0(suffix, "$"), p, value = TRUE)
      if (length(hit) > 0) hit[[length(hit)]] else NA_character_
    }, character(1))
  }

  list(
    genus = vapply(words, `[`, character(1), 1),
    species = species,
    tribe = by_suffix("ini"),
    subfamily = by_suffix("inae"),
    family = by_suffix("idae")
  )
}

#' Read all downloaded files in a folder
#'
#' Worker for [gb_parse_meta()] and [gb_parse_fasta()]. Reads every file,
#' combines the results, removes duplicates, and adds taxonomic ranks.
#'
#' @param dir Root download folder.
#' @param dictionary Marker dictionary.
#' @param region `NULL` for no sequences. `"feature"` or `"record"` to add a
#'   `sequence` column to the features table.
#' @param call Environment to show in error messages.
#' @returns A list with the tibbles `records` and `features`.
#' @noRd
.parse_dir <- function(dir, dictionary, region = NULL, call = rlang::caller_env()) {
  files <- .find_files(dir, call = call)
  parts <- vector("list", length(files$genbank))

  cli::cli_progress_bar("Reading files", total = length(files$genbank))
  for (i in seq_along(files$genbank)) {
    parts[[i]] <- .parse_genbank_file(files$genbank[[i]], dictionary, region)
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  records <- dplyr::bind_rows(lapply(parts, `[[`, "records"))
  features <- dplyr::bind_rows(lapply(parts, `[[`, "features"))
  records <- records[!duplicated(records$accession), ]
  features <- features[!duplicated(features$feature_id), ]

  taxonomy <- .read_taxonomy(files$taxonomy)
  records <- .add_ranks(records, taxonomy)

  list(records = records, features = features)
}

#' Read one GenBank XML file
#'
#' Reads a file of GBSeq records. The helper uses vectorized XPath queries,
#' so it is fast also for files with many records.
#'
#' @param file Path to a GBSeq XML file (can be gzip-compressed).
#' @param dictionary Marker dictionary.
#' @param region `NULL`, `"feature"`, or `"record"`. See `.parse_dir()`.
#' @returns A list with the tibbles `records` and `features`.
#' @noRd
.parse_genbank_file <- function(file, dictionary, region = NULL) {
  doc <- xml2::read_xml(file, options = c("HUGE", "NOBLANKS"))
  recs <- xml2::xml_find_all(doc, "/GBSet/GBSeq")
  if (length(recs) == 0) {
    return(list(records = NULL, features = NULL))
  }

  src_q <- paste0(
    "GBSeq_feature-table/GBFeature[GBFeature_key='source']/GBFeature_quals/",
    "GBQualifier[GBQualifier_name='%s']/GBQualifier_value"
  )
  src <- \(name) .xml_value(recs, sprintf(src_q, name))
  ref <- \(name) .xml_value(recs, sprintf("GBSeq_references/GBReference[1]/%s", name))

  accession <- .xml_value(recs, "GBSeq_accession-version")
  country <- src("geo_loc_name")
  country[is.na(country)] <- src("country")[is.na(country)]

  records <- tibble::tibble(
    accession = accession,
    organism = .xml_value(recs, "GBSeq_organism"),
    taxid = sub("^taxon:", "", .xml_value(
      recs,
      paste0(
        "GBSeq_feature-table/GBFeature[GBFeature_key='source']/GBFeature_quals/",
        "GBQualifier[GBQualifier_name='db_xref' and starts-with(GBQualifier_value, 'taxon:')]/",
        "GBQualifier_value"
      )
    )),
    lineage = .xml_value(recs, "GBSeq_taxonomy"),
    definition = .xml_value(recs, "GBSeq_definition"),
    length = as.integer(.xml_value(recs, "GBSeq_length")),
    mol_type = src("mol_type"),
    organelle = src("organelle"),
    specimen_voucher = src("specimen_voucher"),
    isolate = src("isolate"),
    country = country,
    lat_lon = src("lat_lon"),
    collection_date = src("collection_date"),
    authors = .first_ref_authors(recs, accession),
    ref_title = ref("GBReference_title"),
    journal = ref("GBReference_journal"),
    create_date = .xml_value(recs, "GBSeq_create-date"),
    update_date = .xml_value(recs, "GBSeq_update-date")
  )

  # Gene-level features ---------------------------------------------------
  keys <- c("CDS", "rRNA", "misc_RNA", "mRNA", "misc_feature", "gene")
  feat_q <- sprintf(
    "GBSeq_feature-table/GBFeature[%s]",
    paste0("GBFeature_key='", keys, "'", collapse = " or ")
  )
  feats <- xml2::xml_find_all(recs, feat_q)
  qual <- \(name) .xml_value(
    feats,
    sprintf("GBFeature_quals/GBQualifier[GBQualifier_name='%s']/GBQualifier_value", name)
  )

  features <- tibble::tibble(
    accession = .xml_value(feats, "../../GBSeq_accession-version"),
    feature_key = .xml_value(feats, "GBFeature_key"),
    gene = qual("gene"),
    product = qual("product"),
    note = qual("note"),
    location = .xml_value(feats, "GBFeature_location")
  )

  # Records without gene-level features get one feature for the full record
  bare <- setdiff(accession, features$accession)
  if (length(bare) > 0) {
    bare_len <- records$length[match(bare, accession)]
    features <- dplyr::bind_rows(
      features,
      tibble::tibble(
        accession = bare,
        feature_key = "source",
        location = paste0("1..", bare_len)
      )
    )
  }

  features <- .finish_features(features, records, dictionary)

  if (!is.null(region)) {
    seqs <- toupper(.xml_value(recs, "GBSeq_sequence"))
    full <- seqs[match(features$accession, accession)]
    features$sequence <- if (region == "record") {
      full
    } else {
      mapply(.extract_feature, full, features$location, USE.NAMES = FALSE)
    }
  }

  list(records = records, features = features)
}

#' Names of the rank columns
#'
#' @returns Character vector of rank names, from low to high.
#' @noRd
.rank_columns <- function() {
  c("species", "genus", "tribe", "subfamily", "family", "order", "class", "phylum")
}

#' Read NCBI Taxonomy files
#'
#' @param files Paths to taxonomy XML files.
#' @returns A tibble with one row for each taxonomy ID and the columns
#'   `taxid`, `phylum`, `class`, `order`, `family`, `subfamily`, `tribe`,
#'   `genus`, and `species`. Merged (old) IDs also have rows.
#' @noRd
.read_taxonomy <- function(files) {
  ranks <- .rank_columns()
  if (length(files) == 0) {
    return(NULL)
  }

  long <- lapply(files, function(file) {
    doc <- xml2::read_xml(file, options = c("HUGE", "NOBLANKS"))
    top <- xml2::xml_find_all(doc, "/TaxaSet/Taxon")
    lin <- xml2::xml_find_all(top, "LineageEx/Taxon")
    aka <- xml2::xml_find_all(top, "AkaTaxIds/TaxId")
    list(
      nodes = data.frame(
        taxid = c(.xml_value(top, "TaxId"), .xml_value(lin, "../../TaxId")),
        rank = c(.xml_value(top, "Rank"), .xml_value(lin, "Rank")),
        name = c(.xml_value(top, "ScientificName"), .xml_value(lin, "ScientificName"))
      ),
      aka = data.frame(alias = xml2::xml_text(aka), taxid = .xml_value(aka, "../../TaxId"))
    )
  })
  nodes <- do.call(rbind, lapply(long, `[[`, "nodes"))
  aka <- do.call(rbind, lapply(long, `[[`, "aka"))

  out <- data.frame(taxid = unique(nodes$taxid))
  for (r in ranks) {
    sub <- nodes[nodes$rank == r, ]
    out[[r]] <- sub$name[match(out$taxid, sub$taxid)]
  }
  if (nrow(aka) > 0) {
    alias_rows <- out[match(aka$taxid, out$taxid), ]
    alias_rows$taxid <- aka$alias
    out <- rbind(out, alias_rows[!alias_rows$taxid %in% out$taxid, ])
  }
  tibble::as_tibble(out)
}
