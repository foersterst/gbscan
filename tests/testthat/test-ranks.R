fake_lookup <- function(genus, db) {
  out <- data.frame(
    genus = genus, found = TRUE, tribe = NA_character_, subfamily = NA_character_,
    family = "Buthidae", order = "Scorpiones", class = "Arachnida", phylum = "Arthropoda"
  )
  if (genus == "Tityus") {
    out$tribe <- "Tityini"
    out$subfamily <- "Tityinae"
  }
  out
}

copy_fixture <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  file.copy(list.files(fixture_dir(), full.names = TRUE), d, recursive = TRUE)
  d
}

test_that("ranks = 'ncbi' keeps the old behavior", {
  rec <- suppressMessages(gb_parse_meta(fixture_dir()))$records
  expect_false("filled_ranks" %in% names(rec))
  expect_true(all(is.na(rec$subfamily)))
})

test_that("empty ranks are filled and existing ranks are kept", {
  d <- copy_fixture()
  calls <- 0
  local_mocked_bindings(.lookup_genus = function(genus, db) {
    calls <<- calls + 1
    fake_lookup(genus, db)
  })
  rec <- suppressMessages(gb_parse_meta(d, ranks = "gbif"))$records
  expect_equal(unique(rec$subfamily), "Tityinae")
  expect_equal(unique(rec$tribe), "Tityini")
  expect_equal(rec$filled_ranks[rec$accession == "AB000001.1"], "tribe;subfamily")
  expect_equal(rec$filled_ranks[rec$accession == "KX000003.1"], "tribe;subfamily;order;class;phylum")
  expect_equal(calls, 1)
  expect_true(file.exists(file.path(d, "ranks_gbif.csv")))

  calls <- 0
  suppressMessages(gb_parse_meta(d, ranks = "gbif"))
  expect_equal(calls, 0)
})

test_that("a different classification blocks the fill", {
  d <- copy_fixture()
  local_mocked_bindings(.lookup_genus = function(genus, db) {
    out <- fake_lookup(genus, db)
    out$class <- "Insecta"
    out
  })
  rec <- suppressMessages(gb_parse_meta(d, ranks = "itis"))$records
  ncbi_rows <- rec$accession %in% c("AB000001.1", "MT000002.1")
  expect_true(all(is.na(rec$subfamily[ncbi_rows])))
  expect_true(all(is.na(rec$filled_ranks[ncbi_rows])))
})

test_that("lookup_genus checks the genus name and rank names", {
  local_mocked_bindings(.taxize_classification = function(genus, db) {
    data.frame(name = c("Arthropoda", "Buthidae", "Tityinae", "Tityus"), rank = c("Phylum", "Family", "Subfamily", "Genus"))
  })
  out <- .lookup_genus("Tityus", "gbif")
  expect_true(out$found)
  expect_equal(out$subfamily, "Tityinae")
  expect_false(.lookup_genus("Centruroides", "gbif")$found)
})

test_that("lookup failures give found = FALSE", {
  local_mocked_bindings(.clife_classification = function(genus) stop("no network"))
  expect_false(.lookup_genus("Tityus", "clife")$found)
})

test_that("gb_parse_meta checks ranks", {
  expect_snapshot(gb_parse_meta(fixture_dir(), ranks = "col"), error = TRUE)
})
