test_that("simple and partial locations", {
  expect_equal(unname(.parse_location("10..20")[1, ]), c(10L, 20L, 1L))
  expect_equal(unname(.parse_location("<10..>20")[1, ]), c(10L, 20L, 1L))
  expect_equal(unname(.parse_location("7")[1, ]), c(7L, 7L, 1L))
})

test_that("complement and join combinations", {
  m <- .parse_location("complement(join(1..3,7..9))")
  expect_equal(unname(m[, "start"]), c(7L, 1L))
  expect_equal(unname(m[, "strand"]), c(-1L, -1L))

  m2 <- .parse_location("join(complement(7..9),complement(1..3))")
  expect_equal(m, m2)
})

test_that("unsupported locations give NULL", {
  expect_null(.parse_location("J00194.1:100..202"))
  expect_null(.parse_location("123^124"))
  expect_null(.parse_location(NA))
})

test_that("extract_feature cuts and reverse-complements", {
  s <- "AAACCCGGGTTTA"
  expect_equal(.extract_feature(s, "4..6"), "CCC")
  expect_equal(.extract_feature(s, "complement(1..4)"), "GTTT")
  expect_equal(.extract_feature(s, "join(1..2,12..13)"), "AATA")
  expect_equal(.extract_feature(s, "complement(join(1..2,12..13))"), "TATT")
})

test_that("location_summary computes length and partial status", {
  out <- .location_summary(c("<1..>20", "join(1..5,11..15)", "X:1..5"))
  expect_equal(out$feature_length, c(20L, 10L, NA))
  expect_equal(out$partial, c(TRUE, FALSE, FALSE))
})

test_that("revcomp handles ambiguity codes", {
  expect_equal(.revcomp("ACGTRYN"), "NRYACGT")
})
