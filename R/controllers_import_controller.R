#' Read Genome Data from Various Formats
#'
#' @description
#' Unified entry point for reading genomic data from GenBank or GFF3+FASTA files.
#' Automatically detects format and returns a genome_entity object.
#'
#' This is a controller function that coordinates gateways and use cases.
#'
#' @details
#' ## Format Differences
#' - **GenBank**: Rich metadata (organism, taxonomy, references, comments)
#' - **GFF3+FASTA**: Minimal metadata (sequences and features only)
#'
#' ## Conversion Rules
#' - ALLOWED: GenBank to GFF3+FASTA (use `write_gff3()`, `write_fasta()`)
#' - FORBIDDEN: GFF3+FASTA to GenBank (no export function - intentional)
#'
#' Always keep original GenBank files. See package documentation for details.
#'
#' @param path Character string or named vector specifying input file(s):
#'   - Single file path with .gb/.gbk extension: GenBank format
#'   - Named vector c(gff = "path.gff3", fasta = "path.fasta"): GFF3+FASTA
#' @param gff Path to GFF3 file (alternative to using named vector in path)
#' @param fasta Path to FASTA file (alternative to using named vector in path)
#' @param format Character string to override auto-detection:
#'   "auto" (default), "genbank", or "gff3_fasta"
#' @param ... Additional arguments passed to use cases (e.g., auto_harmonize, verbose)
#'
#' @return
#' A \code{genome_entity} object. In interactive sessions, also prints a
#' formatted summary to the console:
#' \itemize{
#'   \item \strong{Organism} -- scientific name (GenBank sources only)
#'   \item \strong{Source} -- file format detected or specified
#'   \item \strong{Contigs} -- count, total base pairs, and topology
#'   \item \strong{Features} -- total count and breakdown by type (CDS, rRNA, tRNA, other)
#'   \item \strong{Next} -- three suggested follow-on function calls
#' }
#' Suppress the summary in scripts with \code{options(micromicon.quiet = TRUE)}
#' (not yet implemented) or by running in a non-interactive session.
#'
#' @seealso \code{\link{read_variants}()} to load variant calls against the
#'   returned genome; \code{\link{get_features}()} to browse the annotation table.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # GenBank format
#' genome1 <- read_genome("data/ecoli.gbk")
#'
#' # GFF3 + FASTA (named vector)
#' genome2 <- read_genome(c(gff = "data/ecoli.gff3", fasta = "data/ecoli.fna"))
#'
#' # GFF3 + FASTA (separate arguments)
#' genome3 <- read_genome(gff = "data/ecoli.gff3", fasta = "data/ecoli.fna")
#'
#' # Explicit format specification
#' genome4 <- read_genome("data/genome.txt", format = "genbank")
#' }
read_genome <- function(path = NULL, gff = NULL, fasta = NULL,
                        format = c("auto", "genbank", "gff3_fasta"), ...) {
  format <- match.arg(format)

  # Handle different input modes
  if (!is.null(gff) && !is.null(fasta)) {
    # Mode 1: gff= and fasta= arguments
    format <- "gff3_fasta"
    gff_path <- gff
    fasta_path <- fasta
  } else if (!is.null(path) && is.character(path) && length(path) == 2 &&
             !is.null(names(path)) && all(c("gff", "fasta") %in% names(path))) {
    # Mode 2: named vector c(gff = "...", fasta = "...")
    format <- "gff3_fasta"
    gff_path <- path["gff"]
    fasta_path <- path["fasta"]
  } else if (!is.null(path) && is.character(path) && length(path) == 1) {
    # Mode 3: single file path
    if (format == "auto") {
      # Auto-detect based on extension
      if (grepl("\\.(gb|gbk|genbank)$", path, ignore.case = TRUE)) {
        format <- "genbank"
      } else if (grepl("\\.gff3?$", path, ignore.case = TRUE)) {
        cli::cli_abort(c(
          "GFF3 file detected but FASTA file required.",
          "i" = "Provide both: read_genome(gff = '{path}', fasta = 'genome.fasta')"
        ))
      } else {
        cli::cli_abort(c(
          "Cannot auto-detect format from path: {path}",
          "i" = "Use format= argument or provide gff= and fasta= arguments"
        ))
      }
    }

    if (format == "genbank") {
      gbk_path <- path
    } else {
      cli::cli_abort("Single file path only works for GenBank format")
    }
  } else {
    cli::cli_abort(c(
      "Invalid arguments to read_genome().",
      "i" = "Provide either: path (GenBank file), or gff= and fasta= arguments, or c(gff = ..., fasta = ...)"
    ))
  }

  # Read based on detected format
  if (format == "genbank") {
    # Create GenBank gateway
    gateway <- create_genbank_gateway()

    # Call use case
    entity <- execute_import_genbank(gateway, gbk_path, options = list(...))
  } else if (format == "gff3_fasta") {
    # Create gateways
    gff_gateway <- create_gff_gateway(use_bioconductor = TRUE)
    fasta_gateway <- create_fasta_gateway(use_bioconductor = TRUE)

    # Call use case
    entity <- execute_import_gff_fasta(
      gff_gateway, fasta_gateway,
      gff_path, fasta_path,
      options = list(...)
    )
  }

  .read_genome_report(entity, format)
  entity
}

# Post-load summary -------------------------------------------------------

#' @keywords internal
#' @noRd
.read_genome_report <- function(entity, format) {
  if (!interactive()) return(invisible(NULL))

  w   <- .micromicon_console_width()
  col <- getOption("micromicon.color.code", 39L)

  # -- Compute stats --------------------------------------------------------

  n_contigs <- length(entity$indices$seqnames)
  total_bp  <- sum(nchar(entity$sequences$dna_raw))
  bp_fmt    <- formatC(total_bp, format = "d", big.mark = ",")

  feat_types <- if (nrow(entity$features) > 0) {
    table(entity$features$type)
  } else {
    table(character())
  }

  get_n <- function(key) {
    v <- feat_types[key]
    if (length(v) == 0L || is.na(v)) 0L else as.integer(v)
  }
  n_total <- sum(feat_types)
  n_cds   <- get_n("CDS")
  n_rrna  <- get_n("rRNA")
  n_trna  <- get_n("tRNA")
  n_other <- n_total - n_cds - n_rrna - n_trna

  # Topology
  tops     <- entity$metadata$topology
  tops     <- tops[!is.na(tops) & nzchar(tops)]
  topology <- if (length(tops) > 0L) tops[[1L]] else NA_character_

  # Contigs line
  contig_str <- paste0(
    formatC(n_contigs, format = "d", big.mark = ","),
    if (n_contigs == 1L) " contig" else " contigs",
    "  \u00b7  ", bp_fmt, " bp",
    if (!is.na(topology)) paste0("  \u00b7  ", topology) else ""
  )

  # Feature breakdown
  parts <- character()
  if (n_cds   > 0L) parts <- c(parts, paste0(formatC(n_cds,   format = "d", big.mark = ","), " CDS"))
  if (n_rrna  > 0L) parts <- c(parts, paste0(formatC(n_rrna,  format = "d", big.mark = ","), " rRNA"))
  if (n_trna  > 0L) parts <- c(parts, paste0(formatC(n_trna,  format = "d", big.mark = ","), " tRNA"))
  if (n_other > 0L) parts <- c(parts, paste0(formatC(n_other, format = "d", big.mark = ","), " other"))

  feat_str <- paste0(
    formatC(n_total, format = "d", big.mark = ","),
    "  \u2014  ",
    if (length(parts) > 0L) paste(parts, collapse = "  \u00b7  ") else "none"
  )

  # Organism (GenBank only)
  organism <- NULL
  if (format == "genbank") {
    orgs     <- entity$metadata$organism
    orgs     <- orgs[!is.na(orgs) & nzchar(orgs)]
    if (length(orgs) > 0L) organism <- orgs[[1L]]
  }

  fmt_label <- switch(format, genbank = "GenBank", gff3_fasta = "GFF3 + FASTA", format)

  # -- Print ----------------------------------------------------------------

  lbl <- function(x) sprintf("  %-11s", x)

  .micromicon_rule_title("Genome loaded", w)
  if (!is.null(organism)) cat(lbl("Organism"), organism, "\n", sep = "")
  cat(lbl("Source"),   fmt_label, "\n", sep = "")
  cat(lbl("Contigs"),  contig_str, "\n", sep = "")
  cat(lbl("Features"), feat_str, "\n", sep = "")
  cat("\n")

  fn1 <- sprintf("%-34s", "get_features(ref)")
  fn2 <- sprintf("%-34s", 'get_gene_dna(ref, "dnaA")')
  fn3 <- sprintf("%-34s", 'read_variants(ref, "out.gd")')
  cat(lbl("Next"), .micromicon_col256(fn1, col), "browse the annotation\n", sep = "")
  cat(lbl(""),     .micromicon_col256(fn2, col), "fetch a gene sequence\n",  sep = "")
  cat(lbl(""),     .micromicon_col256(fn3, col), "load variant calls\n",     sep = "")

  .micromicon_rule_full(w)
  invisible(NULL)
}

#' Read GenBank File
#'
#' @description
#' Read a GenBank (.gb, .gbk) file. By default returns a genome_entity object.
#' Set return_entity=FALSE for backward compatibility with old list format.
#'
#' @param path Path to GenBank file
#' @param return_entity Logical; if TRUE (default), return genome_entity,
#'   if FALSE, return legacy list format
#' @param ... Additional arguments (currently unused)
#'
#' @return A genome_entity object (if return_entity=TRUE) or list of records (if FALSE)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # New format (genome_entity)
#' genome <- read_gbk("data/ecoli.gbk")
#'
#' # Old format (for backward compatibility)
#' gbk_list <- read_gbk("data/ecoli.gbk", return_entity = FALSE)
#' }
read_gbk <- function(path, return_entity = TRUE, ...) {
  # Create gateway
  gateway <- create_genbank_gateway()

  if (return_entity) {
    # Return genome_entity
    entity <- execute_import_genbank(gateway, path, options = list(...))
    entity
  } else {
    # Return raw records (legacy format)
    records <- gateway$read(path)
    records
  }
}

#' Initialize Genome from GFF3 and FASTA Files
#'
#' @description
#' Initialize a genome object from GFF3 annotation and FASTA sequence files.
#' By default returns a genome_entity object. Set return_entity=FALSE for
#' backward compatibility with old format.
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param return_entity Logical; if TRUE (default), return genome_entity,
#'   if FALSE, return legacy genome_obj format (not yet implemented)
#' @param auto_harmonize Logical; attempt to harmonize seqname mismatches
#' @param verbose Logical; if TRUE, print progress messages
#' @param ... Additional arguments passed to use case
#'
#' @return A genome_entity object (if return_entity=TRUE)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # New format (genome_entity)
#' genome <- init_genome("data/ecoli.gff3", "data/ecoli.fna")
#'
#' # With options
#' genome <- init_genome(
#'   "data/ecoli.gff3",
#'   "data/ecoli.fna",
#'   auto_harmonize = TRUE,
#'   verbose = TRUE
#' )
#' }
init_genome <- function(gff_path, fasta_path,
                        return_entity = TRUE,
                        auto_harmonize = TRUE,
                        verbose = TRUE,
                        ...) {
  # Create gateways
  gff_gateway <- create_gff_gateway(use_bioconductor = TRUE)
  fasta_gateway <- create_fasta_gateway(use_bioconductor = TRUE)

  # Build options
  options <- list(
    auto_harmonize = auto_harmonize,
    verbose = verbose,
    ...
  )

  # Call use case
  entity <- execute_import_gff_fasta(
    gff_gateway, fasta_gateway,
    gff_path, fasta_path,
    options = options
  )

  if (return_entity) {
    entity
  } else {
    # Legacy format conversion (not yet implemented)
    cli::cli_warn("Legacy format (return_entity=FALSE) not yet implemented. Returning genome_entity.")
    entity
  }
}
