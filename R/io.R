#' I/O Functions for Genome Data (Backward Compatibility)
#'
#' @description
#' This file provides backward-compatible I/O functions for reading genome data.
#' These functions are thin wrappers around the Clean Architecture controller implementations
#' located in R/controllers/import_controller.R.
#'
#' In Phase 6 of the Clean Architecture refactoring, these functions were moved to the
#' Interface Adapters layer (controllers). This file exists for documentation and to
#' maintain the historical API surface.
#'
#' @name io
#' @keywords internal
NULL

# Note: The actual implementations are in R/controllers/import_controller.R
#
# Exported functions:
# - read_genome(path, gff, fasta, format, ...) - Unified entry point for reading genomic data
# - read_gbk(path, return_entity, ...) - Read GenBank files
# - init_genome(gff_path, fasta_path, return_entity, ...) - Read GFF3 + FASTA files
#
# All functions are implemented in the controller layer and automatically
# exported via roxygen2 @export tags in those files.
