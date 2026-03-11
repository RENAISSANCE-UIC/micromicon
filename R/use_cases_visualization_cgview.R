# cgview_entity.R
# Adapters from micromicon genome_entity / genome_entity_gd objects to CGView.

#' Render a genome_entity or genome_entity_gd with CGView.js
#'
#' @description
#' Convenience wrapper that extracts the relevant data from a micromicon
#' `genome_entity` or `genome_entity_gd` object and passes it to
#' \code{build_cgview_json()} + \code{beta_cgview()}.
#'
#' For `genome_entity_gd` objects, breseq mutation events are automatically
#' overlaid as an innermost "Variants" track (SNP, DEL, INS, AMP, MOB, SUB).
#'
#' @param entity A `genome_entity` or `genome_entity_gd` object.
#' @param contig Character. Restrict to a single contig (by `seqname`).
#'   `NULL` (default) uses all contigs / first contig.
#' @param feature_types Character vector of feature types to display.
#'   Defaults to `c("CDS", "tRNA", "rRNA")`.
#' @param include_mutations Logical. For `genome_entity_gd` objects, whether
#'   to overlay mutation events as a Variants track. Default `TRUE`.
#' @param palette Named list of hex colours to override defaults.
#' @param width,height Widget dimensions passed to \code{beta_cgview()}.
#' @param elementId Optional fixed HTML id.
#' @param ... Additional arguments forwarded to \code{build_cgview_json()}.
#'
#' @return An htmlwidget.
#'
#' @seealso \code{beta_cgview()}, \code{build_cgview_json()},
#'   \code{beta_cgview_pair_from_entity()}
#'
#' @keywords internal
#' @noRd
beta_cgview_from_entity <- function(entity,
                                     contig            = NULL,
                                     feature_types     = c("CDS", "tRNA", "rRNA"),
                                     include_mutations = TRUE,
                                     palette           = NULL,
                                     width             = NULL,
                                     height            = 2000,
                                     elementId         = NULL,
                                     min_freq          = NULL,
                                     ...) {

  if (!inherits(entity, "genome_entity")) {
    stop("`entity` must be a genome_entity or genome_entity_gd.", call. = FALSE)
  }

  genome_length <- .cgview_genome_length(entity, contig)
  genome_name   <- .cgview_genome_name(entity, contig)

  features_df <- .cgview_features(entity, contig, feature_types)
  mutations_df <- NULL
  if (include_mutations && inherits(entity, "genome_entity_gd")) {
    mutations_df <- .cgview_mutations(entity, contig, min_freq = min_freq)
  }

  json <- build_cgview_json(
    genome_length = genome_length,
    genome_name   = genome_name,
    features_df   = features_df,
    mutations_df  = mutations_df,
    feature_types = feature_types,
    palette       = palette,
    ...
  )

  beta_cgview(json, width = width, height = height, elementId = elementId)
}

#' Side-by-side reference vs variants view from a genome_entity_gd
#'
#' @description
#' Builds a \code{beta_cgview_pair()} layout where the **left** panel shows
#' reference genome features (CDS, tRNA, rRNA) and the **right** panel shows
#' only the mutation overlay -- same genome coordinates, no feature rings.
#' Both panels share the same colour palette and genome length.
#'
#' @param entity A `genome_entity_gd` object (contains both reference
#'   annotation and breseq variant calls).
#' @param contig Character. Restrict to a single contig by `seqname`.
#' @param feature_types Feature types to show in the left (reference) panel.
#' @param palette Named list of hex colour overrides.
#' @param height Height in pixels for both panels.
#' @param left_title,right_title Optional panel titles (`NULL` = no title).
#' @param gap CSS gap between the two panels (default `"12px"`).
#' @param ... Additional arguments forwarded to \code{build_cgview_json()}.
#'
#' @return An `htmltools::browsable` tagList (two widgets side by side).
#'
#' @seealso \code{beta_cgview_pair()}, \code{build_cgview_json()},
#'   \code{beta_cgview_from_entity()}
#'
#' @keywords internal
#' @noRd
beta_cgview_pair_from_entity <- function(entity,
                                          contig        = NULL,
                                          feature_types = c("CDS", "tRNA", "rRNA"),
                                          palette       = NULL,
                                          height        = 700,
                                          left_title    = "Reference",
                                          right_title   = "Variants",
                                          gap           = "12px",
                                          min_freq      = NULL,
                                          ...) {

  if (!inherits(entity, "genome_entity_gd")) {
    stop("`entity` must be a genome_entity_gd.", call. = FALSE)
  }

  genome_length <- .cgview_genome_length(entity, contig)
  genome_name   <- .cgview_genome_name(entity, contig)
  features_df   <- .cgview_features(entity, contig, feature_types)
  mutations_df  <- .cgview_mutations(entity, contig, min_freq = min_freq)

  # Left panel: reference features only (no variants)
  left_json <- build_cgview_json(
    genome_length = genome_length,
    genome_name   = genome_name,
    features_df   = features_df,
    mutations_df  = NULL,
    feature_types = feature_types,
    palette       = palette,
    ...
  )

  # Right panel: variants only (no feature rings)
  right_json <- build_cgview_json(
    genome_length = genome_length,
    genome_name   = genome_name,
    features_df   = NULL,
    mutations_df  = mutations_df,
    feature_types = feature_types,
    palette       = palette,
    ...
  )

  beta_cgview_pair(
    left_json,  right_json,
    height      = height,
    left_title  = left_title,
    right_title = right_title,
    gap         = gap
  )
}

# ---- Internal helpers -------------------------------------------------------

#' @keywords internal
.cgview_features <- function(entity, contig, feature_types) {
  df <- entity$features
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

  if (!is.null(contig) && "seqname" %in% names(df)) {
    df <- df[df$seqname == as.character(contig), , drop = FALSE]
  }

  if (!is.null(feature_types)) {
    df <- df[df$type %in% feature_types, , drop = FALSE]
  }

  if (nrow(df) == 0) NULL else df
}

#' @keywords internal
.cgview_mutations <- function(entity, contig = NULL, min_freq = NULL) {
  events <- entity$events
  if (length(events) == 0) return(NULL)

  mutation_types <- c("SNP", "DEL", "INS", "MOB", "AMP", "SUB", "CON", "INV")

  mut_ev <- Filter(function(ev) {
    !isTRUE(ev$is_evidence) && isTRUE(ev$type %in% mutation_types)
  }, events)

  # Frequency filter via cached variants_predicted table
  if (!is.null(min_freq)) {
    vp       <- entity$variants_predicted
    freq_num <- suppressWarnings(as.numeric(sub("%$", "", vp$freq)) / 100)
    keep_pos <- suppressWarnings(
      as.integer(gsub(",", "", as.character(vp$position[!is.na(freq_num) & freq_num >= min_freq])))
    )
    mut_ev <- Filter(function(ev) as.integer(ev$position) %in% keep_pos, mut_ev)
  }

  if (!is.null(contig)) {
    mut_ev <- Filter(function(ev) {
      !is.null(ev$contig) && identical(as.character(ev$contig),
                                        as.character(contig))
    }, mut_ev)
  }

  if (length(mut_ev) == 0) return(NULL)

  rows <- lapply(mut_ev, function(ev) {
    size <- switch(ev$type,
      DEL = if (!is.null(ev$del_size) && !is.na(ev$del_size))
              as.integer(ev$del_size) else 1L,
      INS = if (!is.null(ev$ins_size) && !is.na(ev$ins_size))
              as.integer(ev$ins_size) else 1L,
      SUB = if (!is.null(ev$sub_size) && !is.na(ev$sub_size))
              as.integer(ev$sub_size) else 1L,
      1L
    )
    list(
      type   = as.character(ev$type),
      pos    = as.integer(ev$position),
      size   = size,
      contig = as.character(ev$contig %||% NA_character_)
    )
  })

  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

#' @keywords internal
.cgview_genome_length <- function(entity, contig = NULL) {
  md <- entity$metadata
  if (is.data.frame(md) && nrow(md) > 0) {
    rows <- if (!is.null(contig) && "seqname" %in% names(md))
              md[md$seqname == as.character(contig), , drop = FALSE]
            else
              md[1, , drop = FALSE]
    for (col in c("length_bp", "length", "seq_length")) {
      if (col %in% names(rows) && nrow(rows) > 0) {
        len <- suppressWarnings(as.integer(rows[[col]][1]))
        if (!is.na(len) && len > 0) return(len)
      }
    }
  }
  stop("Cannot determine genome length from entity metadata.", call. = FALSE)
}

#' @keywords internal
.cgview_genome_name <- function(entity, contig = NULL) {
  md <- entity$metadata
  if (is.data.frame(md) && nrow(md) > 0) {
    rows <- if (!is.null(contig) && "seqname" %in% names(md))
              md[md$seqname == as.character(contig), , drop = FALSE]
            else
              md[1, , drop = FALSE]
    if (nrow(rows) > 0) {
      for (col in c("organism", "accession", "definition", "seqname")) {
        if (col %in% names(rows)) {
          nm <- as.character(rows[[col]][1])
          if (!is.na(nm) && nzchar(nm)) return(nm)
        }
      }
    }
  }
  if (!is.null(contig)) paste0("Contig ", contig) else "Genome"
}
