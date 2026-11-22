test_that("build_query makes GenBank queries", {
  dict <- gb_dictionary()
  q <- .build_query(c("Tityus", "txid6887"), c("COI", "COI"), dict)
  expect_match(q[1], '^"Tityus"\\[Organism\\] AND \\(COI\\[Gene Name\\]', perl = TRUE)
  expect_match(q[1], "AND 100:5000\\[SLEN\\]$")
  expect_match(q[2], "^txid6887\\[Organism:exp\\]")
})

test_that("gb_count uses the counter for each taxon and marker", {
  local_mocked_bindings(.ncbi_count = function(query, db = "nuccore") nchar(query))
  out <- gb_count(c("Tityus", "Centruroides"), c("COI", "16S"))
  expect_equal(nrow(out), 4)
  expect_equal(out$count, nchar(out$query))
})

test_that("gb_count checks its arguments", {
  expect_snapshot(gb_count(character()), error = TRUE)
  expect_snapshot(gb_count(NA_character_), error = TRUE)
})
