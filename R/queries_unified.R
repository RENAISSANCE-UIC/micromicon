#' Unified Query Functions (Backward Compatibility)
#'
#' @description
#' This file provides backward-compatible query functions for genome data.
#' These functions are thin wrappers around the Clean Architecture controller implementations
#' located in R/controllers/query_controller.R.
#'
#' In Phase 6 of the Clean Architecture refactoring, these functions were moved to the
#' Interface Adapters layer (controllers). This file exists for documentation and to
#' maintain the historical API surface.
#'
#' @name queries
#' @keywords internal
NULL

# Note: The actual implementations are in R/controllers/query_controller.R
#
# Exported functions:
# - extract_sequences_by_name(x, pattern, fields, type, ignore_case, translate, ...) - Extract sequences by name
# - extract_sequence_by_coords(x, seqname, start, end, strand, translate, names, ...) - Extract by coordinates
# - search_features(x, type, pattern, seqname, start, end, strand, ...) - Search features
# - analyze_roi(x, ...) - Analyze regions of interest (from R/queries.R if exists)
# - get_genomic_context(x, ...) - Get genomic context (from R/queries.R if exists)
#
# All functions are implemented in the controller layer and automatically
# exported via roxygen2 @export tags in those files.
