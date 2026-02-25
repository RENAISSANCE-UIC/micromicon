# zzz.R — package lifecycle hooks
# Sourced last (alphabetically) so all package internals are available.

# Startup message ---------------------------------------------------------

#' @keywords internal
#' @noRd
.micromicon_print_welcome <- function() {
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

# .onAttach ---------------------------------------------------------------
# Fires when the user calls library(micromicon). Prints a compact startup
# screen. Suppressed by suppressPackageStartupMessages().

#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  # Suppress during devtools::load_all() — the package is already reachable
  # via :: in that context; the greeting would be noise on every reload.
  # pkgload::is_loading() is TRUE for the duration of load_all().
  if (isNamespaceLoaded("pkgload") &&
      isTRUE(tryCatch(pkgload::is_loading(pkgname), error = function(e) FALSE))) {
    return(invisible(NULL))
  }

  # Suppress in non-interactive sessions (Rscript, CI, R CMD INSTALL subprocesses).
  if (!interactive()) return(invisible(NULL))

  .micromicon_print_welcome()
}

# micromicon_welcome ------------------------------------------------------

#' Display the micromicon welcome screen
#'
#' Prints the same startup message shown when the package is first attached
#' with \code{library(micromicon)}. Useful for revisiting the quick-start
#' examples or for developers testing the message without restarting R.
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of printing.
#' @examples
#' micromicon_welcome()
#'
#' @export
micromicon_welcome <- function() {
  .micromicon_print_welcome()
}
