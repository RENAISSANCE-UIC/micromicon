#' Accessor Functions for genome_entity (Backward Compatibility)
#'
#' @description
#' This file provides backward-compatible accessor functions for genome_entity objects.
#' These functions are thin wrappers around the Clean Architecture controller implementations
#' located in R/controllers/accessor_controller.R.
#'
#' In Phase 6 of the Clean Architecture refactoring, these functions were moved to the
#' Interface Adapters layer (controllers). This file exists for documentation and to
#' maintain the historical API surface.
#'
#' @name accessors
#' @keywords internal
NULL

# Note: The actual implementations are in R/controllers/accessor_controller.R
# and R/frameworks/bioconductor_compat.R
#
# Exported functions:
# - has_bioconductor() - Check if Bioconductor packages are available (from bioconductor_compat.R)
# - sequences(x, format, ...) - Access sequences from genome_entity
# - features(x, format, type, ...) - Access features from genome_entity
# - metadata(x, ...) - Access metadata from genome_entity
# - seqnames(x, ...) - Get sequence names from genome_entity
#
# All functions are implemented in the controller layer and automatically
# exported via roxygen2 @export tags in those files.
