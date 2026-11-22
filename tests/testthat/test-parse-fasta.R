test_that("gb_parse_fasta writes one folder per marker", {
  out_dir <- withr::local_tempdir()
  res <- suppressMessages(gb_parse_fasta(fixture_dir(), out_dir = out_dir))
  expect_setequal(res$marker, c("COI", "16S", "COIII", "28S", "18S"))
  expect_true(all(file.exists(file.path(out_dir, res$marker, paste0(res$marker, ".fasta")))))

  coiii <- readLines(file.path(out_dir, "COIII", "COIII.fasta"))
  expect_equal(coiii, c(">MT000002.1|Tityus_obscurus", "CCCCGT"))
})

test_that("keep filters features and gives header columns", {
  out_dir <- withr::local_tempdir()
  meta <- suppressMessages(gb_parse_meta(fixture_dir()))
  keep <- meta$features[meta$features$accession == "AB000001.1", ]
  keep$label <- "my_label"
  suppressMessages(gb_parse_fasta(fixture_dir(), out_dir, keep = keep, header = c("label", "family")))
  expect_equal(readLines(file.path(out_dir, "COI", "COI.fasta"))[1], ">my_label|Buthidae")
  expect_equal(list.files(out_dir), "COI")
})

test_that("region = 'record' writes full sequences", {
  out_dir <- withr::local_tempdir()
  suppressMessages(gb_parse_fasta(fixture_dir(), out_dir, markers = "COIII", region = "record"))
  expect_equal(nchar(readLines(file.path(out_dir, "COIII", "COIII.fasta"))[2]), 56)
})

test_that("gb_parse_fasta does not overwrite files and checks headers", {
  out_dir <- withr::local_tempdir()
  suppressMessages(gb_parse_fasta(fixture_dir(), out_dir, markers = "28S"))
  expect_snapshot(gb_parse_fasta(fixture_dir(), out_dir, markers = "28S"), error = TRUE, transform = \(x) gsub("'[^']*'", "<file>", x))
  expect_snapshot(gb_parse_fasta(fixture_dir(), out_dir, markers = "18S", header = "nope"), error = TRUE)
})
