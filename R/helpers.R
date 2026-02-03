#' Extract Feature Attribute with NA Handling
#'
#' @param df Features data.frame
#' @param key Character; attribute name
#' @return Character vector (NA for missing)
#' @export
feat_attr <- function(df, key) {
  if (!key %in% names(df)) {
    return(rep(NA_character_, nrow(df)))
  }
  df[[key]]
}


#' Filter Features by Multiple Criteria
#'
#' @param df Features data.frame
#' @param type Feature type
#' @param name Gene name pattern (regex)
#' @param id Feature ID
#' @param attr Named list of attribute filters
#' @return Filtered data.frame
#' @export
feat_filter <- function(df, type = NULL, name = NULL, id = NULL, attr = list()) {
  if (!is.data.frame(df)) {
    cli::cli_abort("{.arg df} must be a data.frame")
  }

  # Filter by type
  if (!is.null(type) && "type" %in% names(df)) {
    df <- df[df$type == type, ]
  }

  # Filter by name pattern
  if (!is.null(name)) {
    matches <- rep(FALSE, nrow(df))
    for (field in c("gene", "locus_tag", "product")) {
      if (field %in% names(df)) {
        field_matches <- grepl(name, df[[field]], ignore.case = TRUE)
        field_matches[is.na(field_matches)] <- FALSE
        matches <- matches | field_matches
      }
    }
    df <- df[matches, ]
  }

  # Filter by ID
  if (!is.null(id) && "ID" %in% names(df)) {
    df <- df[!is.na(df$ID) & df$ID == id, ]
  }

  # Filter by attributes
  for (attr_name in names(attr)) {
    if (attr_name %in% names(df)) {
      df <- df[!is.na(df[[attr_name]]) & df[[attr_name]] == attr[[attr_name]], ]
    }
  }

  df
}


#' Create Standardized ROI Specification
#'
#' @param seqname Character; sequence name
#' @param start Integer; start position (1-based)
#' @param end Integer; end position (1-based, inclusive)
#' @param strand Character; strand (default "+")
#' @return roi_spec object with validated coordinates
#' @export
roi_coords <- function(seqname, start, end, strand = "+") {
  # Validate
  validate_genomic_coordinates(start, end, seqname)
  if (!is.null(strand)) validate_strand(strand)

  structure(
    list(
      seqname = seqname,
      start = as.integer(start),
      end = as.integer(end),
      strand = strand
    ),
    class = "roi_spec"
  )
}


#' Print ROI Specification
#'
#' @param x roi_spec object
#' @param ... Additional arguments (ignored)
#' @export
print.roi_spec <- function(x, ...) {
  cat("ROI: ", x$seqname, ":", x$start, "-", x$end,
      " (", x$strand, ")\n", sep = "")
  invisible(x)
}
