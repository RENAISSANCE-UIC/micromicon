# controllers_plot_controller.R
# User-facing static plotting functions.

#' Plot a static genomic region of interest
#'
#' @description
#' Draws a linear map of a genomic region returned by [get_roi_features()].
#' Each feature is rendered as a directional arrow coloured uniquely per gene.
#' Gene labels are placed above the track with automatic 2-D collision
#' avoidance: overlapping labels are first nudged horizontally (within the
#' plot bounds) and, when horizontal room is exhausted, lifted to a higher
#' tier. Leader lines connect displaced labels back to their arrows.
#'
#' Requires the **ggplot2** package. Install it with
#' `install.packages("ggplot2")` if prompted.
#'
#' @param roi A `GRanges` or `data.frame` (or any object coercible to one)
#'   typically produced by [get_roi_features()]. Required columns:
#'   `seqnames`, `start`, `end`, `strand`. Optional columns `Name`, `gene`,
#'   `Alias`, and `type` are used when present to derive labels and feature
#'   colours.
#'
#'   If this argument is not the output of [get_roi_features()], you may see
#'   unexpected results. Run [get_roi_features()] first to obtain a compatible
#'   object.
#' @param title `character(1)` Plot title passed to [ggplot2::labs()].
#'   Default `NULL` (no title).
#' @param arrow_height `numeric(1)` Total vertical extent of each arrow in
#'   internal y-units. Increase for taller, more prominent arrows. Default
#'   `0.18`.
#' @param head_prop `numeric(1)` Fraction of each gene's length occupied by
#'   the arrowhead, in the range (0, 1]. Default `0.35`.
#' @param neck_prop `numeric(1)` Height of the arrow body as a fraction of
#'   the full arrowhead height, in the range (0, 1]. `1` gives a rectangular
#'   body; values below `1` create a narrowed neck. Default `0.60`.
#' @param label_size `numeric(1)` Text size for gene labels as accepted by
#'   [ggplot2::geom_text()]. Default `4.0`.
#' @param colors `character` or `NULL`. When `NULL` (default) colours are
#'   drawn from an evenly-spaced HCL palette matching ggplot2 defaults. Supply
#'   a named vector to map gene names to specific colours, or an unnamed
#'   vector that is recycled across genes in plotting order.
#'
#' @return A [`ggplot2::ggplot`] object. Further layers and theme adjustments
#'   can be added with `+` in the usual way.
#'
#' @seealso [get_roi_features()] to obtain a region of interest from a genome
#'   object. [plot_cgview()] for an interactive circular genome map.
#'
#' @export
#'
#' @examples
#' roi <- data.frame(
#'   seqnames = "U00096",
#'   start    = c(1971030, 1972836, 1973360, 1975329, 1976252,
#'                1977266, 1977847, 1978518, 1979753, 1980188, 1981587),
#'   end      = c(1972691, 1973339, 1975324, 1976255, 1977139,
#'                1977844, 1978197, 1978966, 1980181, 1981612, 1982387),
#'   strand   = c("-","-","-","-","-","-","-","-","+","-","-"),
#'   type     = "CDS",
#'   Name     = c("tar","cheW","cheA","motB","motA",
#'                "flhC","flhD","pgaptmp_001971","uspC","otsA","otsB"),
#'   stringsAsFactors = FALSE
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_roi(roi)
#' }
plot_roi <- function(roi,
                     title        = NULL,
                     arrow_height = 0.18,
                     head_prop    = 0.35,
                     neck_prop    = 0.60,
                     label_size   = 4.0,
                     colors       = NULL) {

  # ggplot2 availability check ------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_roi} requires the {.pkg ggplot2} package.",
      "i" = "Install it with: {.code install.packages(\"ggplot2\")}"
    ))
  }

  # Input validation ----------------------------------------------------------

  if (is.null(roi)) {
    cli::cli_abort(
      c("{.arg roi} is {.code NULL}.",
        "i" = "Run {.fn get_roi_features} first to obtain a compatible object.")
    )
  }

  df <- tryCatch(
    as.data.frame(roi, stringsAsFactors = FALSE),
    error = function(e)
      cli::cli_abort(
        c("Cannot coerce {.arg roi} to a data frame: {e$message}",
          "i" = "Run {.fn get_roi_features} first to obtain a compatible object.")
      )
  )

  if (nrow(df) == 0L) {
    cli::cli_abort(
      c("{.arg roi} contains no rows. Nothing to plot.",
        "i" = "Run {.fn get_roi_features} first to obtain a compatible object.")
    )
  }

  required <- c("seqnames", "start", "end", "strand")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c("{.arg roi} is missing required column{?s}: {.field {missing}}",
        "i" = "Run {.fn get_roi_features} first to obtain a compatible object.")
    )
  }

  if (!is.numeric(df$start) || !is.numeric(df$end)) {
    cli::cli_abort(
      "{.field start} and {.field end} must be numeric, not {.cls {class(df$start)}}."
    )
  }

  for (nm in c("arrow_height", "head_prop", "neck_prop", "label_size")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
      cli::cli_abort(
        "{.arg {nm}} must be a single positive finite number, not {.val {val}}."
      )
    }
  }

  if (head_prop > 1) {
    cli::cli_warn(
      "{.arg head_prop} is {head_prop}, which is > 1; arrowheads will exceed gene length."
    )
  }
  if (neck_prop > 1) {
    cli::cli_warn(
      "{.arg neck_prop} is {neck_prop}, which is > 1; body will be wider than head."
    )
  }

  if (!is.null(colors) && !is.character(colors)) {
    cli::cli_abort("{.arg colors} must be a character vector or {.code NULL}.")
  }

  # Informational messages ----------------------------------------------------

  n_feat  <- nrow(df)
  contig  <- unique(as.character(df$seqnames))[1L]
  span_kb <- round((max(df$end) - min(df$start)) / 1e3, 1)

  cli::cli_inform(c(
    "i" = "Plotting {n_feat} feature{?s} on {.val {contig}} \
           ({span_kb} kbp spanning {min(df$start)}\u2013{max(df$end)})"
  ))

  # Dispatch to use case -------------------------------------------------------

  plot_roi_impl(
    roi,
    title        = title,
    arrow_height = arrow_height,
    head_prop    = head_prop,
    neck_prop    = neck_prop,
    label_size   = label_size,
    colors       = colors
  )
}
