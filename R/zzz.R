# zzz.R — package lifecycle hooks
# Sourced last (alphabetically) so all package internals are available.

# .onAttach ---------------------------------------------------------------
# Fires when the user calls library(micromicon). Prints a compact startup
# screen. Suppressed by suppressPackageStartupMessages().

#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  if (!interactive()) return(invisible(NULL))

  ver <- tryCatch(
    as.character(utils::packageVersion("micromicon")),
    error = function(e) "?"
  )

  w    <- .micromicon_console_width()
  dash <- .micromicon_dash()

  # Title rule: micRomicon v0.x.x ────────
  title_text   <- paste0("micRomicon v", ver)
  title_styled <- .micromicon_bold(title_text)
  used <- nchar(title_text, type = "width") + 1L
  rem  <- max(0L, w - used)
  cat("\n\n", title_styled, " ", 
      paste(rep(dash, rem), collapse = ""), "\n", sep = "")

  # Welcome
  cat("Welcome to micRomicon!\n\n")

  # Getting-started hints
  cat("To get started, load a microbial genome with annotation:\n\n")
  code <- .micromicon_col256
  col  <- getOption("micromicon.color.code", 39L)
  cat("  ", code('read_genome("reference.gbk")', col),
      "                     # GenBank\n", sep = "")
  cat("  ", code('read_genome("reference.gff3", "reference.fasta")', col),
      "  # GFF3 + FASTA\n", sep = "")
  cat("\n")
  cat("Call ", .micromicon_bold("micromicon_functions()"),
      " for a full function index.\n", sep = "")

  # Closing rule
  cat(paste(rep(dash, w), collapse = ""), "\n", sep = "")
  invisible(NULL)
}
