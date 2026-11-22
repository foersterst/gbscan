fake_gbset <- function(ids) {
  recs <- paste0(
    "<GBSeq><GBSeq_accession-version>", ids, "</GBSeq_accession-version>",
    "<GBQualifier_value>taxon:", ids, "</GBQualifier_value></GBSeq>",
    collapse = ""
  )
  paste0("<GBSet>", recs, "</GBSet>")
}

test_that("gb_download downloads in batches and resumes", {
  dir <- withr::local_tempdir()
  calls <- 0
  local_mocked_bindings(
    .ncbi_search_ids = function(query, db = "nuccore") as.character(1:7),
    .ncbi_fetch = function(ids, db = "nuccore") {
      calls <<- calls + 1
      if (db == "taxonomy") "<TaxaSet></TaxaSet>" else fake_gbset(ids)
    }
  )

  res <- suppressMessages(gb_download("Tityus", "COI", dir, batch = 3))
  path <- file.path(dir, "genbank", "COI", "Tityus")
  expect_equal(res$n_records, 7)
  expect_length(list.files(path, pattern = "^records_"), 3)
  expect_setequal(readLines(file.path(path, "done.txt")), as.character(1:7))
  expect_length(list.files(file.path(dir, "genbank", "taxonomy"), pattern = "^taxonomy_"), 1)

  calls <- 0
  suppressMessages(gb_download("Tityus", "COI", dir, batch = 3))
  expect_equal(calls, 0)
})

test_that("gb_download downloads only missing IDs after an interruption", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "genbank", "COI", "Tityus")
  dir.create(path, recursive = TRUE)
  query <- .build_query("Tityus", "COI", gb_dictionary())
  writeLines(c(query, "date", "5"), file.path(path, "query.txt"))
  writeLines(as.character(1:5), file.path(path, "ids.txt"))
  writeLines(as.character(1:2), file.path(path, "done.txt"))

  got <- character()
  local_mocked_bindings(
    .ncbi_search_ids = function(...) stop("must not search"),
    .ncbi_fetch = function(ids, db = "nuccore") {
      got <<- c(got, ids)
      fake_gbset(ids)
    }
  )
  suppressMessages(gb_download("Tityus", "COI", dir, taxonomy = FALSE))
  expect_equal(got, as.character(3:5))
})

test_that("gb_download checks batch", {
  expect_snapshot(gb_download("Tityus", "COI", tempdir(), batch = 0), error = TRUE)
})
