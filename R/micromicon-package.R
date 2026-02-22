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
#' pm  <- predict_mutations(gd)
#' eff <- compute_effects(gd, pm)
#' ```
#'
#' Call \code{\link{micromicon_functions}()} at the console for a full
#' categorised function listing.
#'
#' @keywords internal
"_PACKAGE"

# ── micromicon_functions ───────────────────────────────────────────────────────

# micromicon_functions <- function() {
#   
#   cli::cli_h1("micromicon")
#   
#   cli::cli_h2("Reading data")
#   cli::cli_ul(c(
#     "{.fn read_genome}    \u2014 load a reference genome (.gbk, .gff3 + .fasta)",
#     "{.fn read_variants}  \u2014 load variant calls against a reference (GD; VCF coming)"
#   ))
#   
#   cli::cli_h2("Genomic features")
#   cli::cli_ul(c(
#     "{.fn get_features}      \u2014 all annotated features as a data.frame",
#     "{.fn search_features}   \u2014 filter features by type, name pattern, or coordinates",
#     "{.fn get_roi_features}  \u2014 features overlapping a genomic region",
#     "{.fn gene_info}         \u2014 full details for a single gene"
#   ))
#   
#   cli::cli_h2("Sequences")
#   cli::cli_ul(c(
#     "{.fn get_gene_dna}  \u2014 CDS nucleotide sequence for a gene",
#     "{.fn get_gene_aa}   \u2014 amino acid sequence for a gene",
#     "{.fn get_roi_dna}   \u2014 nucleotide sequence for any genomic region",
#     "{.fn sequences}     \u2014 all contig/chromosome sequences"
#   ))
#   
#   cli::cli_h2("Coordinate mapping")
#   cli::cli_ul(c(
#     "{.fn map_cds_to_genomic}   \u2014 CDS position \u2192 genomic position",
#     "{.fn map_genomic_to_cds}   \u2014 genomic position \u2192 CDS position",
#     "{.fn extract_by_name}      \u2014 extract sequence and features by gene name",
#     "{.fn extract_by_coords}    \u2014 extract sequence and features by coordinates"
#   ))
#   
#   cli::cli_h2("Genome metadata")
#   cli::cli_ul(c(
#     "{.fn get_genome_metadata}  \u2014 organism name, length, topology, etc.",
#     "{.fn get_genome_seqnames}  \u2014 contig / chromosome names"
#   ))
#   
#   cli::cli_h2("Variants & mutations")
#   cli::cli_ul(c(
#     "{.fn predict_mutations}        \u2014 tidy mutation table from a genome diff",
#     "{.fn compute_effects}          \u2014 molecular consequences (aa change, consequence class)",
#     "{.fn validate_variant_in_gene} \u2014 check whether a variant falls within a gene"
#   ))
#   
#   cli::cli_h2("Export")
#   cli::cli_ul(c(
#     "{.fn export_genome}  \u2014 write genome as GFF3, FASTA, or GenBank"
#   ))
#   
#   cli::cli_inform(c(
#     "i" = "Use {.code ?<function>} for full documentation, e.g. {.code ?read_genome}.",
#     "i" = "Use {.code ?micromicon} for the package overview."
#   ))
#   
#   invisible(NULL)
# }


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

#' Print the micromicon function index in the console
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
    c("predict_mutations()",        "tidy mutation table (GenomeDiff parser)"),
    c("compute_effects()",          "molecular consequences of variants"),
    c("validate_variant_in_gene()", "check whether a variant resides in a gene")
  ), width)
  
  .micromicon_section("Export", list(
    c("export_genome()", "write genome as GFF3, FASTA, or GenBank")
  ), width)
  
  .micromicon_rule_full(width)
  cat(.micromicon_bold("\u2139"), " For details: use ?<function> or see ?micromicon.\n", sep = "")
  invisible(NULL)
}
