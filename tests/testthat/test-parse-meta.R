test_that("gb_parse_meta returns records and features", {
  meta <- suppressMessages(gb_parse_meta(fixture_dir()))
  expect_named(meta, c("records", "features"))
  expect_equal(nrow(meta$records), 4)
  expect_false("source" %in% names(meta$records))
  expect_equal(anyDuplicated(meta$features$feature_id), 0)
})

test_that("ranks come from taxonomy data, merged IDs, or guesses", {
  rec <- suppressMessages(gb_parse_meta(fixture_dir()))$records
  r <- \(acc, col) rec[[col]][rec$accession == acc]
  expect_equal(r("AB000001.1", "order"), "Scorpiones")
  expect_equal(r("MT000002.1", "species"), "Tityus obscurus")
  expect_equal(r("KX000003.1", "family"), "Buthidae")
  expect_equal(r("KX000003.1", "genus"), "Tityus")
})

test_that("markers are assigned and redundant features are removed", {
  f <- suppressMessages(gb_parse_meta(fixture_dir()))$features
  m <- \(acc) sort(stats::na.omit(f$marker[f$accession == acc]))
  expect_equal(m("AB000001.1"), "COI")
  expect_equal(f$feature_key[f$accession == "AB000001.1"], "CDS")
  expect_equal(m("MT000002.1"), c("16S", "COI", "COIII"))
  expect_equal(m("KX000003.1"), "28S")
  expect_equal(f$matched_field[f$accession == "KX000004.1"], "definition")
})

test_that("organelle decides between 16S and 28S", {
  dict <- gb_dictionary()
  out <- .assign_markers(
    gene = c(NA, NA),
    product = c("large subunit ribosomal RNA", "large subunit ribosomal RNA"),
    note = c(NA, NA),
    definition = c(NA, NA),
    organelle = c("mitochondrion", NA),
    dictionary = dict
  )
  expect_equal(out$marker, c("16S", "28S"))
})

test_that("a changed dictionary changes the assignment", {
  dict <- gb_dictionary()
  dict <- dict[dict$marker != "28S", ]
  f <- suppressMessages(gb_parse_meta(fixture_dir(), dictionary = dict))$features
  expect_equal(f$marker[f$accession == "KX000003.1"], "16S")
})

test_that("gb_parse_meta fails on an empty folder", {
  d <- withr::local_tempdir()
  expect_snapshot(gb_parse_meta(d), error = TRUE, transform = \(x) sub(" in .*\\.$", " in <dir>.", x))
})

test_that("open names stay different species without taxonomy data", {
  org <- c("Tityus sp. 1", "Tityus sp. 2", "Tityus sp 1", "Tityus cf bahiensis", "Tityus serrulatus", "Tityus")
  g <- .guess_ranks(org, rep(NA, length(org)))
  expect_equal(g$species, c(org[1:4], "Tityus serrulatus", NA))
})
