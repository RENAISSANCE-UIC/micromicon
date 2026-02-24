#' micromicon: Genomic Analysis for Microbial Genomes
#'
#' @description
#' Tools for loading, querying, and analyzing microbial reference genomes and
#' variant calls (breseq genome diffs). Built around two core S3 objects:
#' \code{genome_entity} (reference genome) and \code{genome_entity_gd}
#' (genome + variant calls).
#'
#' ## Typical workflow
#'
#' ```r
#' # 1. Load a reference genome
#' ref <- read_genome("reference.gbk")
#'
#' # 2. Explore features and sequences
#' get_features(ref, type = "CDS")
#' get_gene_dna(ref, "dnaA")
#'
#' # 3. Load variant calls
#' gd <- read_variants(ref, "annotated.gd")
#'
#' # 4. Analyse mutations
#' pm  <- predict_variants(gd)
#' eff <- annotate_variants(gd, pm)
#' ```
#'
#' Call \code{\link{micromicon_functions}()} at the console for a full
#' categorised function listing.
#'
#' @keywords internal
"_PACKAGE"

##' Print the micromicon function index in the console
#'
#' Renders a tidy, ANSI-aware, dependency-free index of the package’s
#' core capabilities. Prints to stdout (no tinted message background),
#' wraps descriptions to the current console width, and optionally
#' uses ANSI bold and a single function name color.
#'
#' @section Styling:
#' Set options to adjust styling at runtime:
#' \itemize{
#'   \item \code{options(micromicon.bold = TRUE)} to toggle bold headers.
#'   \item \code{options(micromicon.color.enabled = TRUE)} to enable colors.
#'   \item \code{options(micromicon.color.code = 39)} to set the xterm-256 color.
#'   \item \code{options(micromicon.ascii = FALSE)} for ASCII-only symbols.
#' }
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of printing.
#' @examples
#' micromicon_functions()
#'
#' @export
micromicon_functions <- function() {
  width <- .micromicon_console_width()
  
  .micromicon_rule_title("micromicon", width)
  cat("A concise index of the package’s core capabilities.\n\n", sep = "")
  
  .micromicon_section("Reading", list(
    c("read_genome()",   "import a reference genome (GenBank or GFF3 + FASTA)"),
    c("read_variants()", "load variant calls (GenomeDiff; VCF forthcoming)")
  ), width)
  
  .micromicon_section("Features", list(
    c("get_features()",     "retrieve all annotated features"),
    c("search_features()",  "filter by type, name, or coordinates"),
    c("get_roi_features()", "features intersecting a region of interest"),
    c("gene_info()",        "structured metadata for a single gene")
  ), width)
  
  .micromicon_section("Sequences", list(
    c("get_gene_dna()", "CDS nucleotide sequence"),
    c("get_gene_aa()",  "translated amino acid sequence"),
    c("get_roi_dna()",  "arbitrary region sequence"),
    c("sequences()",    "all sequences (chromosomes / contigs)")
  ), width)
  
  .micromicon_section("Coordinate mapping", list(
    c("map_cds_to_genomic()", "CDS coordinate \u2192 genomic coordinate"),
    c("map_genomic_to_cds()", "genomic coordinate \u2192 CDS coordinate"),
    c("extract_by_name()",    "extract sequence and features by gene name"),
    c("extract_by_coords()",  "extract sequence and features by coordinates")
  ), width)
  
  .micromicon_section("Metadata", list(
    c("get_genome_metadata()", "organism, topology, span, and provenance"),
    c("get_genome_seqnames()", "contig / chromosome names")
  ), width)
  
  .micromicon_section("Variants", list(
    c("predict_variants()",         "tidy variant table (GenomeDiff parser)"),
    c("annotate_variants()",        "molecular consequences of variants"),
    c("filter_variants()",          "filter by gene, frequency, type, or consequence"),
    c("validate_variant_in_gene()", "check whether a variant resides in a gene")
  ), width)
  
  .micromicon_section("Export", list(
    c("export_genome()", "write genome as GFF3, FASTA, or GenBank")
  ), width)
  
  .micromicon_rule_full(width)
  cat(.micromicon_bold("\u2139"), " For details: use ?<function> or see ?micromicon.\n", sep = "")
  invisible(NULL)
}

# Helpers ----

#' ---- ANSI helpers (no external deps) ----
#' @keywords internal
#' @noRd
.micromicon_console_width <- function() {
  w <- getOption("width", 80L)
  if (is.null(w) || !is.numeric(w) || length(w) != 1L || is.na(w)) 80L else as.integer(w)
}

#' @keywords internal
#' @noRd
.micromicon_ansi_on <- function() {
  opt <- getOption("micromicon.color.enabled", NULL)
  if (isTRUE(opt)) return(TRUE)
  if (identical(opt, FALSE)) return(FALSE)
  
  is_tty <- sink.number(type = "output") == 0L && isatty(stdout())
  is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")
  is_tty || is_rstudio
}

#' @keywords internal
#' @noRd
.micromicon_bold <- function(x) {
  if (.micromicon_ansi_on()) paste0("\033[1m", x, "\033[22m") else x
}

#' @keywords internal
#' @noRd
.micromicon_col256 <- function(x, code) {
  if (.micromicon_ansi_on()) paste0("\033[38;5;", code, "m", x, "\033[39m") else x
}

#' @keywords internal
#' @noRd
.micromicon_dash <- function() "\u2500"  # ─

# ---- Rules (headers) ----

#' @keywords internal
#' @noRd
.micromicon_rule_title <- function(text, width = .micromicon_console_width()) {
  dash <- .micromicon_dash()
  label <- .micromicon_bold(text)
  used <- nchar(text, type = "width") + 1L
  rem <- max(0L, width - used)
  cat(label, " ", paste(rep(dash, rem), collapse = ""), "\n", sep = "")
}

#' @keywords internal
#' @noRd
.micromicon_rule_section <- function(name, width = .micromicon_console_width()) {
  dash <- .micromicon_dash()
  left <- "─  "
  mid  <- .micromicon_bold(name)
  right_pad <- "  "
  used <- nchar("─  ", type = "width") +
    nchar(name, type = "width") +
    nchar("  ", type = "width")
  rem <- max(0L, width - used)
  cat(left, mid, right_pad, paste(rep(dash, rem), collapse = ""), "\n", sep = "")
}


#' @keywords internal
#' @noRd
.micromicon_rule_full <- function(width = .micromicon_console_width()) {
  cat(paste(rep(.micromicon_dash(), width), collapse = ""), "\n", sep = "")
}

# ---- Bullets with wrapping ----

#' @keywords internal
#' @noRd
.micromicon_bullet <- function(fun, desc, width = .micromicon_console_width(), color_code = 39) {
  bullet <- "\u2022 "  # •
  fun_col <- .micromicon_col256(fun, color_code)
  sep <- " — "
  
  prefix_visible <- nchar(bullet, type = "width") +
    nchar(fun,    type = "width") +
    nchar(sep,    type = "width")
  avail <- max(20L, width - prefix_visible)
  
  wrapped <- strwrap(desc, width = avail, exdent = 2L)
  if (length(wrapped) == 0) wrapped <- ""
  
  cat(bullet, fun_col, sep, wrapped[1L], "\n", sep = "")
  
  if (length(wrapped) > 1L) {
    cont_prefix <- paste0(
      paste(rep(" ", nchar(bullet, type = "width")), collapse = ""),
      paste(rep(" ", nchar(fun,    type = "width")), collapse = ""),
      paste(rep(" ", nchar(sep,    type = "width")), collapse = "")
    )
    for (i in 2:length(wrapped)) cat(cont_prefix, wrapped[i], "\n", sep = "")
  }
}

#' @keywords internal
#' @noRd
.micromicon_section <- function(title, entries, width = .micromicon_console_width(),
                                color_code = getOption("micromicon.color.code", 39)) {
  .micromicon_rule_section(title, width)
  if (length(entries)) {
    for (i in seq_along(entries)) {
      e <- entries[[i]]
      .micromicon_bullet(e[1], e[2], width = width, color_code = color_code)
    }
  }
}
