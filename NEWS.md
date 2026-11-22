# gbscan 0.0.0.9000

* First version of the package.
* `gb_count()` counts GenBank records for taxa and markers. It prints the
  counts with shortened queries and returns the full result invisibly.
* `gb_dictionary()` reads and copies the marker dictionary.
* `gb_download()` downloads records in batches and continues a stopped download.
* `gb_parse_fasta()` writes the sequences of each marker to a FASTA file.
* `gb_parse_meta()` reads the record and feature metadata into two tibbles.
  Its `ranks` argument fills empty ranks from GBIF, Catalogue of Life, ITIS,
  or WoRMS.
