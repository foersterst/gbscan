# write and read a dictionary copy

    Code
      gb_dictionary(path, write = TRUE)
    Condition
      Error in `gb_dictionary()`:
      ! Can't write the dictionary.
      x The file <path> exists.
      i Use a different path or remove the file.

# check_dictionary finds problems

    Code
      .check_dictionary(dict[, -1])
    Condition
      Error:
      ! `dict[, -1]` must have the column marker.
      i Use `gb_dictionary()` to see the necessary columns.

---

    Code
      .check_dictionary(rbind(dict, dict[1, ]))
    Condition
      Error:
      ! `rbind(dict, dict[1, ])` has duplicated marker name: "COI".

---

    Code
      .check_dictionary(bad)
    Condition
      Error:
      ! Column genome of `bad` has incorrect values.
      i Use "mitochondrial" or "nuclear".
      x Incorrect value: "plastid".

# select_markers keeps order and rejects unknown markers

    Code
      .select_markers("XYZ", dict[1:3, ])
    Condition
      Error:
      ! "XYZ" is not in the marker dictionary.
      i Available markers: "COI", "COII", and "COIII".
      i Add markers to the dictionary with `gb_dictionary()`.

