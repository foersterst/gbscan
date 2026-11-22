# Cool loading message: works with both load_all() and library()
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    cli::format_inline(
      "\U0001F9EC {.pkg gbscan} {utils::packageVersion(pkgname)} \U0001F9EC"
    )
  )
}
