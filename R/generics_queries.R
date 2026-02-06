#' Search Features
#'
#' @param x Object to search
#' @param ... Additional arguments
#' @export
search_features <- function(x, ...) {
  UseMethod("search_features")
}

#' @export
search_features.genome_entity <- function(x, type = NULL, pattern = NULL,
                                         seqname = NULL, start = NULL,
                                         end = NULL, strand = NULL, ...) {
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (!is.null(pattern)) {
    matches <- rep(FALSE, nrow(feats))
    for (field in c("ID", "Name", "Alias", "gene", "locus_tag", "product")) {
      if (field %in% names(feats)) {
        field_matches <- grepl(pattern, feats[[field]], ignore.case = TRUE)
        matches <- matches | field_matches
      }
    }
    feats <- feats[matches, ]
  }

  if (!is.null(seqname) && "seqname" %in% names(feats)) {
    feats <- feats[feats$seqname == seqname, ]
  }

  if (!is.null(start) && "start" %in% names(feats)) {
    feats <- feats[feats$start >= start, ]
  }

  if (!is.null(end) && "end" %in% names(feats)) {
    feats <- feats[feats$end <= end, ]
  }

  if (!is.null(strand) && "strand" %in% names(feats)) {
    feats <- feats[feats$strand == strand, ]
  }

  feats
}

#' @export
search_features.default <- function(x, ...) {
  cli::cli_abort("search_features() not implemented for class {.cls {class(x)[1]}}")
}


#' Extract Sequences by Coordinates
#'
#' @param x Object to extract from
#' @param ... Additional arguments
#' @export
extract_by_coords <- function(x, ...) {
  UseMethod("extract_by_coords")
}

#' @export
extract_by_coords.genome_entity <- function(x, seqname, start, end,
                                           strand = "+",
                                           translate = FALSE,
                                           names = NULL, ...) {
  validate_genome_entity(x)
  options <- list(strand = strand, translate = translate, names = names)
  execute_extract_sequences_by_coords(x, seqname, start, end, options)
}

#' @export
extract_by_coords.default <- function(x, ...) {
  cli::cli_abort("extract_by_coords() not implemented for class {.cls {class(x)[1]}}")
}


#' Extract Sequences by Name
#'
#' @param x Object to extract from
#' @param ... Additional arguments
#' @export
extract_by_name <- function(x, ...) {
  UseMethod("extract_by_name")
}

#' @export
extract_by_name.genome_entity <- function(x, pattern,
                                         fields = c("gene", "locus_tag", "product"),
                                         type = NULL,
                                         ignore_case = TRUE,
                                         translate = FALSE, ...) {
  validate_genome_entity(x)
  options <- list(ignore_case = ignore_case, type_filter = type, translate = translate)
  execute_extract_sequences_by_name(x, pattern, fields, options)
}

#' @export
extract_by_name.default <- function(x, ...) {
  cli::cli_abort("extract_by_name() not implemented for class {.cls {class(x)[1]}}")
}


# Backward compatibility aliases
#' @rdname extract_by_coords
#' @export
extract_sequences_by_coords <- extract_by_coords

#' @rdname extract_by_name
#' @export
extract_sequences_by_name <- extract_by_name
