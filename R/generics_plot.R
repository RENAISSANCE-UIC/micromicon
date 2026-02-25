#' Plot a circular genome map with CGView.js
#'
#' @description
#' S3 generic that dispatches on `genome_entity` or `genome_entity_gd`.
#'
#' @param entity A `genome_entity` or `genome_entity_gd`.
#' @param ... Passed to the class-specific method.
#'
#' @return An htmlwidget (single panel) or browsable tagList (paired panels).
#' @export
plot_cgview <- function(entity, ...) UseMethod("plot_cgview")

#' @param contig Character. Restrict to one contig by seqname. NULL = first contig.
#' @param feature_types Character vector of feature types to display.
#'   Default: `c("CDS", "tRNA", "rRNA")`.
#' @param palette Named list of hex colour overrides.
#' @param width,height Widget dimensions (px). height defaults to 2000.
#' @param elementId Optional fixed HTML element id.
#' @param viewer Where to display the plot. `"pane"` (default) opens in the
#'   RStudio Viewer pane; `"browser"` saves to a temporary HTML file and
#'   opens it in the system browser at full window width.
#'
#' @rdname plot_cgview
#' @export
plot_cgview.genome_entity <- function(entity,
                                       contig        = NULL,
                                       feature_types = c("CDS", "tRNA", "rRNA"),
                                       palette       = NULL,
                                       width         = 2000,
                                       height        = 2000,
                                       elementId     = NULL,
                                       viewer        = c("pane", "browser"),
                                       ...) {
  viewer <- match.arg(viewer)
  w <- beta_cgview_from_entity(
    entity,
    contig            = contig,
    feature_types     = feature_types,
    include_mutations = FALSE,
    palette           = palette,
    width             = width,
    height            = height,
    elementId         = elementId,
    ...
  )
  .cgview_open(w, viewer)
}

#' @param paired Logical. FALSE (default) = single map with mutations overlaid.
#'   TRUE = side-by-side reference vs variants panels.
#' @param left_title,right_title Panel titles for paired view. Ignored when
#'   `paired = FALSE`.
#' @param gap CSS gap between panels for paired view (default `"12px"`).
#'
#' @rdname plot_cgview
#' @export
plot_cgview.genome_entity_gd <- function(entity,
                                          contig        = NULL,
                                          feature_types = c("CDS", "tRNA", "rRNA"),
                                          paired        = FALSE,
                                          left_title    = "Reference",
                                          right_title   = "Variants",
                                          gap           = "12px",
                                          palette       = NULL,
                                          width         = 2000,
                                          height        = 2000,
                                          elementId     = NULL,
                                          viewer        = c("pane", "browser"),
                                          ...) {
  viewer <- match.arg(viewer)
  w <- if (paired) {
    beta_cgview_pair_from_entity(
      entity,
      contig        = contig,
      feature_types = feature_types,
      palette       = palette,
      height        = height,
      left_title    = left_title,
      right_title   = right_title,
      gap           = gap,
      ...
    )
  } else {
    beta_cgview_from_entity(
      entity,
      contig            = contig,
      feature_types     = feature_types,
      include_mutations = TRUE,
      palette           = palette,
      width             = width,
      height            = height,
      elementId         = elementId,
      ...
    )
  }
  .cgview_open(w, viewer)
}

#' @rdname plot_cgview
#' @export
plot_cgview.default <- function(entity, ...) {
  cli::cli_abort(
    "{.fn plot_cgview} is not implemented for class {.cls {class(entity)[1]}}."
  )
}
