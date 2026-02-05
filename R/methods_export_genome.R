#' Export Genome Methods
#'
#' @description
#' S3 method implementations for the export_genome() generic.
#' Provides polymorphic export functionality for different genome object types.
#'
#' @name export_genome-methods
NULL

#' @rdname export_genome-methods
#' @export
export_genome.genome_entity <- function(x, file, format = c("auto", "fasta", "gff3"),
                                       wrap_width = 80, source = "micromicon", ...) {
  format <- match.arg(format)

  # Auto-detect from file extension if needed
  if (format == "auto") {
    format <- detect_format_from_extension(file)
  }

  # Dispatch to appropriate use case
  switch(format,
    fasta = execute_export_fasta(x, file, wrap_width = wrap_width, ...),
    gff3 = execute_export_gff3(x, file, source = source, ...),
    cli::cli_abort("Unsupported export format: {format}")
  )
}

#' @rdname export_genome-methods
#' @export
export_genome.default <- function(x, file, format = c("auto", "fasta", "gff3"), ...) {
  cli::cli_abort("export_genome() not implemented for class {.cls {class(x)[1]}}")
}

#' Detect Export Format from File Extension
#'
#' @description
#' Helper function to automatically detect export format from file extension.
#'
#' @param file File path
#' @return Character string: "fasta" or "gff3"
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' detect_format_from_extension("output.fasta")  # "fasta"
#' detect_format_from_extension("output.gff3")   # "gff3"
#' detect_format_from_extension("output.gff")    # "gff3"
#' }
detect_format_from_extension <- function(file) {
  ext <- tolower(tools::file_ext(file))

  format <- switch(ext,
    "fasta" = "fasta",
    "fa" = "fasta",
    "fna" = "fasta",
    "faa" = "fasta",
    "gff3" = "gff3",
    "gff" = "gff3",
    {
      # Default to fasta if unrecognized
      cli::cli_warn(c(
        "Unrecognized file extension: {.file {ext}}",
        "i" = "Defaulting to FASTA format",
        "i" = "Specify format explicitly with format = 'fasta' or format = 'gff3'"
      ))
      "fasta"
    }
  )

  format
}
