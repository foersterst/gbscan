# gb_parse_fasta does not overwrite files and checks headers

    Code
      gb_parse_fasta(fixture_dir(), out_dir, markers = "28S")
    Condition
      Error in `gb_parse_fasta()`:
      ! 1 FASTA file exists: <file>.
      i Use `overwrite = TRUE` or a different `out_dir`.

---

    Code
      gb_parse_fasta(fixture_dir(), out_dir, markers = "18S", header = "nope")
    Condition
      Error in `gb_parse_fasta()`:
      ! Can't find the header column nope.
      i Use a column of `keep` or of the tables from `gb_parse_meta()`.

