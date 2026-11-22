test_that("default dictionary has the expected markers and columns", {
  dict <- gb_dictionary()
  expect_s3_class(dict, "tbl_df")
  expect_contains(dict$marker, c("COI", "12S", "16S", "18S", "28S", "ITS1", "ITS2", "H3", "mitogenome"))
  expect_named(dict, c("marker", "genome", "description", "query", "regex", "min_len", "max_len"))
  expect_equal(anyDuplicated(dict$marker), 0)
})

test_that("all regular expressions compile", {
  dict <- gb_dictionary()
  rx <- stats::na.omit(dict$regex)
  expect_no_error(for (r in rx) grepl(r, "test", perl = TRUE))
})

test_that("write and read a dictionary copy", {
  path <- withr::local_tempfile(fileext = ".csv")
  expect_message(gb_dictionary(path, write = TRUE), "Wrote")
  expect_equal(gb_dictionary(path), gb_dictionary())
  expect_snapshot(gb_dictionary(path, write = TRUE), error = TRUE, transform = \(x) sub("file .* exists", "file <path> exists", x))
})

test_that("check_dictionary finds problems", {
  dict <- gb_dictionary()
  expect_snapshot(.check_dictionary(dict[, -1]), error = TRUE)
  expect_snapshot(.check_dictionary(rbind(dict, dict[1, ])), error = TRUE)
  bad <- dict
  bad$genome[1] <- "plastid"
  expect_snapshot(.check_dictionary(bad), error = TRUE)
})

test_that("select_markers keeps order and rejects unknown markers", {
  dict <- gb_dictionary()
  expect_equal(.select_markers(c("16S", "COI"), dict)$marker, c("16S", "COI"))
  expect_snapshot(.select_markers("XYZ", dict[1:3, ]), error = TRUE)
})
