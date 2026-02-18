#' Convert genome_entity to legacy genome_obj format
#'
#' @description
#' Converts a genome_entity to the old init_genome() object format
#' for backward compatibility.
#'
#' @param entity A genome_entity object
#'
#' @return Legacy genome_obj list
#' @export
#'
#' @examples
#' \dontrun{
#' entity <- read_genome(gff = "annotation.gff3", fasta = "genome.fasta")
#' genome_obj <- entity_to_legacy_genome_obj(entity)
#' }
entity_to_legacy_genome_obj <- function(entity) {
  validate_genome_entity(entity)

  # Check for required Bioconductor packages
  if (!has_bioconductor()) {
    missing_pkgs <- character()
    required_pkgs <- c("GenomicRanges", "Biostrings", "rtracklayer", "Rsamtools", "IRanges", "S4Vectors")

    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_pkgs <- c(missing_pkgs, pkg)
      }
    }

    cli::cli_abort(c(
      "Converting to legacy genome_obj requires Bioconductor packages.",
      "x" = "Missing: {paste(missing_pkgs, collapse = ', ')}",
      "i" = "Install all required packages with:",
      " " = "if (!require('BiocManager')) install.packages('BiocManager')",
      " " = "BiocManager::install(c('GenomicRanges', 'Biostrings', 'rtracklayer', 'Rsamtools', 'IRanges', 'S4Vectors'))"
    ))
  }

  # Get GRanges (create if needed)
  gff <- features(entity, format = "GRanges")

  # Get DNAStringSet (create if needed)
  fasta <- sequences(entity, format = "DNAStringSet")

  # Create FaFile if we have a fasta_used path (only present for GFF+FASTA workflows)
  fa <- NULL
  if ("fasta_used" %in% names(entity$metadata)) {
    fasta_used <- entity$metadata$fasta_used[1]
    if (!is.na(fasta_used) && nzchar(fasta_used)) {
      fasta_path <- fasta_used
      if (file.exists(fasta_path)) {
        # Index if needed
        if (!file.exists(paste0(fasta_path, ".fai"))) {
          Rsamtools::indexFa(fasta_path)
        }
        fa <- Rsamtools::FaFile(fasta_path)
      }
    }
  }

  # For multi-row metadata, take first row's values if columns exist
  # (these columns only exist for GFF+FASTA workflows, not GenBank)
  gff_used_val <- ""
  if ("gff_used" %in% names(entity$metadata) && nrow(entity$metadata) > 0) {
    gff_used_val <- entity$metadata$gff_used[1]
    if (is.na(gff_used_val)) gff_used_val <- ""
  }

  import_errs_val <- ""
  if ("import_errors" %in% names(entity$metadata) && nrow(entity$metadata) > 0) {
    import_errs_val <- entity$metadata$import_errors[1]
    if (is.na(import_errs_val)) import_errs_val <- ""
  }

  list(
    gff = gff,
    fasta = fasta,
    fa = fa,
    seqnames = entity$indices$seqnames,
    gff_used = gff_used_val,
    import_errs = import_errs_val
  )
}
