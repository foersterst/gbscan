# Build the marker dictionary shipped in inst/extdata/markers.csv.
# Run this script from the package root after you change a marker.
#
# Columns:
#   marker       Short marker name. gbscan uses it for folder and file names.
#   genome       "mitochondrial" or "nuclear". It controls marker assignment.
#   description  Full marker name.
#   query        Entrez search terms for GenBank, joined with OR.
#   regex        Case-insensitive pattern. gbscan matches it against the
#                gene, product, note, and definition fields of a feature.
#   min_len      Minimum record length in the GenBank search ([SLEN]).
#   max_len      Maximum record length in the GenBank search ([SLEN]).

gene <- function(...) paste0(c(...), "[Gene Name]")
title <- function(...) paste0("\"", c(...), "\"[Title]")
or <- function(...) paste(c(...), collapse = " OR ")

markers <- tibble::tribble(
  ~marker, ~genome, ~description, ~query, ~regex, ~min_len, ~max_len,

  # Mitochondrial protein-coding genes ------------------------------------
  "COI", "mitochondrial", "cytochrome c oxidase subunit I",
  or(
    gene("COI", "COX1", "CO1", "COXI"),
    title(
      "cytochrome oxidase subunit I", "cytochrome c oxidase subunit I",
      "cytochrome oxidase subunit 1", "cytochrome c oxidase subunit 1"
    )
  ),
  "^(coi|cox1|co1|coxi|cox-i|coi-5p)$|cytochrome (c )?oxidase subunit (i|1)\\b",
  100, 5000,
  "COII", "mitochondrial", "cytochrome c oxidase subunit II",
  or(
    gene("COII", "COX2", "CO2", "COXII"),
    title(
      "cytochrome oxidase subunit II", "cytochrome c oxidase subunit II",
      "cytochrome oxidase subunit 2", "cytochrome c oxidase subunit 2"
    )
  ),
  "^(coii|cox2|co2|coxii|cox-ii)$|cytochrome (c )?oxidase subunit (ii|2)\\b",
  100, 5000,
  "COIII", "mitochondrial", "cytochrome c oxidase subunit III",
  or(
    gene("COIII", "COX3", "CO3", "COXIII"),
    title(
      "cytochrome oxidase subunit III", "cytochrome c oxidase subunit III",
      "cytochrome oxidase subunit 3", "cytochrome c oxidase subunit 3"
    )
  ),
  "^(coiii|cox3|co3|coxiii|cox-iii)$|cytochrome (c )?oxidase subunit (iii|3)\\b",
  100, 5000,
  "CytB", "mitochondrial", "cytochrome b",
  or(gene("CYTB", "COB"), title("cytochrome b")),
  "^(cytb|cob|cyt b|cyt-b)$|cytochrome b\\b",
  100, 5000,
  "ND1", "mitochondrial", "NADH dehydrogenase subunit 1",
  or(gene("ND1", "NAD1"), title("NADH dehydrogenase subunit 1")),
  "^(nd1|nad1)$|nadh dehydrogenase subunit (1|i)\\b",
  100, 5000,
  "ND4", "mitochondrial", "NADH dehydrogenase subunit 4",
  or(gene("ND4", "NAD4"), title("NADH dehydrogenase subunit 4")),
  "^(nd4|nad4)$|nadh dehydrogenase subunit (4|iv)\\b",
  100, 5000,
  "ND5", "mitochondrial", "NADH dehydrogenase subunit 5",
  or(gene("ND5", "NAD5"), title("NADH dehydrogenase subunit 5")),
  "^(nd5|nad5)$|nadh dehydrogenase subunit (5|v)\\b",
  100, 5000,
  "ATP6", "mitochondrial", "ATP synthase F0 subunit 6",
  or(gene("ATP6", "ATPASE6"), title("ATP synthase F0 subunit 6", "ATPase 6")),
  "^(atp6|atpase6|atpase 6)$|atp synthase (f0 )?subunit (6|a)\\b",
  100, 5000,

  # Mitochondrial ribosomal RNA genes --------------------------------------
  "12S", "mitochondrial", "12S ribosomal RNA (small subunit)",
  or(
    gene("12S", "rrnS", "s-rRNA"),
    title("12S ribosomal RNA", "12S rRNA", "small subunit ribosomal RNA")
  ),
  "^(12s|rrns|s-rrna|rrn12)$|12s|small subunit ribosomal rna",
  100, 5000,
  "16S", "mitochondrial", "16S ribosomal RNA (large subunit)",
  or(
    gene("16S", "rrnL", "l-rRNA"),
    title("16S ribosomal RNA", "16S rRNA", "large subunit ribosomal RNA")
  ),
  "^(16s|rrnl|l-rrna|rrn16)$|16s|large subunit ribosomal rna",
  100, 5000,

  # Nuclear ribosomal RNA genes and spacers --------------------------------
  "18S", "nuclear", "18S ribosomal RNA (small subunit)",
  or(gene("18S"), title("18S ribosomal RNA", "18S rRNA", "small subunit ribosomal RNA")),
  "^(18s|ssu)$|18s|small subunit ribosomal rna",
  100, 5000,
  "ITS1", "nuclear", "internal transcribed spacer 1",
  or(gene("ITS1"), title("internal transcribed spacer 1", "ITS1")),
  "^its-?1$|internal transcribed spacer 1\\b",
  100, 5000,
  "5.8S", "nuclear", "5.8S ribosomal RNA",
  or(gene("5.8S"), title("5.8S ribosomal RNA", "5.8S rRNA")),
  "5\\.8s",
  100, 5000,
  "ITS2", "nuclear", "internal transcribed spacer 2",
  or(gene("ITS2"), title("internal transcribed spacer 2", "ITS2")),
  "^its-?2$|internal transcribed spacer 2\\b",
  100, 5000,
  "28S", "nuclear", "28S ribosomal RNA (large subunit)",
  or(gene("28S"), title("28S ribosomal RNA", "28S rRNA", "large subunit ribosomal RNA")),
  "^(28s|lsu)$|28s|large subunit ribosomal rna",
  100, 5000,

  # Nuclear protein-coding genes -------------------------------------------
  "H3", "nuclear", "histone H3",
  or(gene("H3", "His3", "H3A"), title("histone H3")),
  "^(h3|his3|h3a)$|histone h3\\b",
  100, 5000,
  "EF1a", "nuclear", "elongation factor 1 alpha",
  or(
    gene("EF1a", "EF-1a", "EF-1alpha", "EF1alpha"),
    title(
      "elongation factor 1 alpha", "elongation factor-1 alpha",
      "elongation factor 1-alpha"
    )
  ),
  "^ef-?1-?(a|alpha)$|elongation factor.?1.?(a|alpha)",
  100, 5000,
  "wingless", "nuclear", "wingless (Wnt-1)",
  or(gene("wingless", "wg"), title("wingless")),
  "^wg$|wingless",
  100, 5000,
  "ArgK", "nuclear", "arginine kinase",
  or(gene("ArgK"), title("arginine kinase")),
  "^argk$|arginine kinase",
  100, 5000,
  "CAD", "nuclear", "carbamoyl-phosphate synthetase (CAD, rudimentary)",
  or(
    gene("CAD", "rudimentary"),
    title(
      "carbamoylphosphate synthetase", "carbamoyl-phosphate synthetase",
      "carbamoyl phosphate synthetase", "carbamoylphosphate synthase"
    )
  ),
  "^cad$|rudimentary|carbamoyl.?phosphate synth",
  100, 5000,
  "MDH", "nuclear", "malate dehydrogenase",
  or(gene("MDH"), title("malate dehydrogenase")),
  "^c?mdh\\d?$|malate dehydrogenase",
  100, 5000,
  "RpS5", "nuclear", "ribosomal protein S5",
  or(gene("RpS5"), title("ribosomal protein S5")),
  "^rps5$|ribosomal protein s5\\b",
  100, 5000,
  "GAPDH", "nuclear", "glyceraldehyde-3-phosphate dehydrogenase",
  or(
    gene("GAPDH"),
    title(
      "glyceraldehyde-3-phosphate dehydrogenase",
      "glyceraldehyde 3-phosphate dehydrogenase"
    )
  ),
  "^gapdh$|glyceraldehyde[ -]3[ -]phosphate dehydrogenase",
  100, 5000,
  "IDH", "nuclear", "isocitrate dehydrogenase",
  or(gene("IDH"), title("isocitrate dehydrogenase")),
  "^idh\\d?$|isocitrate dehydrogenase",
  100, 5000,
  "CaATPase", "nuclear", "sarco/endoplasmic reticulum calcium ATPase",
  or(
    gene("Ca-ATPase", "SERCA"),
    title("Ca-ATPase", "calcium ATPase", "sarco/endoplasmic reticulum calcium ATPase")
  ),
  "ca-?atpase|calcium[- ]atpase|^serca|sarco.?endoplasmic reticulum ca",
  100, 5000,
  "Nex9", "nuclear", "sorting nexin 9",
  or(gene("Nex9", "SNX9"), title("sorting nexin 9", "sorting nexin-9")),
  "^nex9$|^snx9$|sorting nexin.?9",
  100, 5000,
  "DDC", "nuclear", "dopa decarboxylase",
  or(gene("DDC"), title("dopa decarboxylase", "DOPA decarboxylase")),
  "^ddc$|dopa decarboxylase",
  100, 5000,
  "PEPCK", "nuclear", "phosphoenolpyruvate carboxykinase",
  or(gene("PEPCK"), title("phosphoenolpyruvate carboxykinase")),
  "^pepck$|phosphoenolpyruvate carboxykinase",
  100, 5000,

  # Complete mitochondrial genomes -----------------------------------------
  # gbscan does not assign this marker to features (regex is NA). It uses
  # the search to get mitogenomes. The parsers then cut each gene from them.
  "mitogenome", "mitochondrial", "complete mitochondrial genome",
  "mitochondrion[Filter] AND \"complete genome\"[Title]",
  NA,
  10000, 30000
)

utils::write.csv(markers, "inst/extdata/markers.csv", row.names = FALSE, na = "")
