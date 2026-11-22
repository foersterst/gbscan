# gbscan

gbscan counts and downloads nucleotide records for genetic markers and taxa
from GenBank. It reads the downloaded records into two tibbles
(records and features) and writes the sequences of each marker to FASTA
files.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("<github-user>/gbscan")
```

## Functions

| Function | Use |
|---|---|
| `gb_dictionary()` | Get or edit the marker dictionary. |
| `gb_count()` | Count the records for taxa and markers. |
| `gb_download()` | Download the records in batches. A stopped download continues from the last batch. |
| `gb_parse_meta()` | Read the metadata into a records tibble and a features tibble. |
| `gb_parse_fasta()` | Write the sequences of each marker to a FASTA file. |

## Example

```r
library(gbscan)

gb_count(taxa = "Tityus", markers = c("COI", "16S"))

gb_download(taxa = "Tityus", markers = c("COI", "mitogenome"), dir = "tityus_data")

meta <- gb_parse_meta("tityus_data")
gb_parse_fasta("tityus_data", out_dir = "fasta")
```

The vignette `vignette("tityus-coi")` shows the full pipeline.

## NCBI API key

NCBI accepts 3 requests per second without an API key, and 10 with a key.
Add the line `ENTREZ_KEY=your_key` to your `.Renviron` file. Do not write
the key in a script.
