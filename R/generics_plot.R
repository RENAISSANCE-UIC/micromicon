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
                                       min_freq      = NULL,
                                       ...) {
  if (!is.null(min_freq)) {
    cli::cli_warn(c(
      "!" = "{.arg min_freq} was ignored.",
      "i" = "Frequency filtering is only available for {.cls genome_entity_gd} objects.",
      "i" = "Load variant calls with {.fn read_variants}, then call {.fn predict_variants}."
    ))
  }
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
#' @param min_freq Numeric (0–1) or NULL (default). When set, only variants
#'   with allele frequency >= this value are shown. Requires
#'   \code{gd$variants_predicted} to be populated first via
#'   \code{\link{predict_variants}()}.
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
                                          min_freq      = NULL,
                                          ...) {
  viewer <- match.arg(viewer)

  if (!is.null(min_freq)) {
    if (is.null(entity$variants_predicted)) {
      cli::cli_abort(c(
        "{.arg min_freq} requires {.code gd$variants_predicted} to be populated.",
        "i" = "Run {.fn predict_variants} on your {.cls genome_entity_gd} first."
      ))
    }
    stopifnot(is.numeric(min_freq), length(min_freq) == 1L,
              min_freq >= 0, min_freq <= 1)
  }

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
      min_freq      = min_freq,
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
      min_freq          = min_freq,
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

# =============================================================================
# plot_roi — genome_entity method
# =============================================================================

#' @rdname plot_roi
#'
#' @param entity A \code{genome_entity} (or \code{genome_entity_gd}) object.
#' @param gene `character(1)` Gene name (or pattern) to look up. Mutually
#'   exclusive with `contig`/`start`/`end`. The search is case-insensitive
#'   and matches `ID`, `Name`, `Alias`, `gene`, `locus_tag`, and `product`
#'   fields. Features of type `"gene"` are excluded from matching to avoid
#'   duplicate hits common in NCBI-style annotations. If exactly one feature
#'   matches, its coordinates plus `flank` define the plotted window.
#' @param contig `character(1)` Sequence name. Ignored when `gene` is supplied
#'   (the contig is inferred from the matched feature). In coordinate mode,
#'   when omitted the first contig returned by [get_contig_names()] is used
#'   with an informational message.
#' @param start,end `integer(1)` Region boundaries in coordinate mode. Both
#'   must be supplied together. Ignored when `gene` is supplied.
#' @param flank `integer(1)` Base pairs added on each side of the target
#'   region. Defaults to `5000` in gene mode and `0` in coordinate mode.
#'   Override explicitly to zoom in or out.
#' @param feature_type `character` Feature type(s) passed to
#'   [get_roi_features()] and drawn in the plot. Default `"CDS"`.
#'
#' @export
plot_roi.genome_entity <- function(entity,
                                   gene         = NULL,
                                   contig       = NULL,
                                   start        = NULL,
                                   end          = NULL,
                                   flank        = NULL,
                                   feature_type = "CDS",
                                   title        = NULL,
                                   arrow_height = 0.18,
                                   head_prop    = 0.35,
                                   neck_prop    = 0.60,
                                   label_size   = 4.0,
                                   colors       = NULL,
                                   ...) {

  using_gene   <- !is.null(gene)
  using_coords <- !is.null(start) || !is.null(end) || !is.null(contig)

  if (using_gene && using_coords) {
    cli::cli_abort(c(
      "Provide either {.arg gene} or coordinates \
       ({.arg contig}, {.arg start}, {.arg end}), not both."
    ))
  }

  if (!using_gene && !using_coords) {
    cli::cli_abort(c(
      "Provide either {.arg gene} or coordinates \
       ({.arg contig}, {.arg start}, {.arg end})."
    ))
  }

  # ---------------------------------------------------------------------------
  # Gene name mode
  # ---------------------------------------------------------------------------
  if (using_gene) {

    # Tiered gene lookup: broad search first, then narrow by priority.
    # Tier 1: exact case-insensitive match on gene/Name/ID/locus_tag.
    # Tier 2: word-boundary regex match across those fields + product.
    # Tier 3: original substring fallback (current behaviour).
    .narrow_hits <- function(feats, g) {
      feats <- feats[!tolower(feats$type) %in% "gene", , drop = FALSE]
      if (nrow(feats) <= 1L) return(feats)

      name_fields <- intersect(c("gene", "Name", "ID", "locus_tag"), names(feats))

      # Tier 1 — exact (case-insensitive)
      mask <- rep(FALSE, nrow(feats))
      for (f in name_fields) {
        v    <- as.character(feats[[f]])
        mask <- mask | (!is.na(v) & tolower(v) == tolower(g))
      }
      if (any(mask)) return(feats[mask, , drop = FALSE])

      # Tier 2 — word-boundary
      wb_pat  <- paste0("\\b", g, "\\b")
      wb_mask <- rep(FALSE, nrow(feats))
      for (f in intersect(c(name_fields, "product"), names(feats))) {
        v       <- as.character(feats[[f]])
        wb_mask <- wb_mask | (!is.na(v) & grepl(wb_pat, v, ignore.case = TRUE, perl = TRUE))
      }
      if (any(wb_mask)) return(feats[wb_mask, , drop = FALSE])

      # Tier 3 — original substring results
      feats
    }

    hits <- .narrow_hits(search_features(entity, pattern = gene), gene)

    if (nrow(hits) == 0L) {
      cli::cli_abort(c(
        "No features matching {.val {gene}} were found.",
        "i" = "Check spelling or try {.fn search_features} to explore available features."
      ))
    }

    if (nrow(hits) > 1L) {
      name_col <- if ("Name" %in% names(hits)) hits$Name else rep(NA_character_, nrow(hits))
      id_col   <- if ("ID"   %in% names(hits)) hits$ID   else rep(NA_character_, nrow(hits))
      labels   <- paste0(
        "[", hits$type, "] ",
        ifelse(!is.na(name_col) & nzchar(name_col), name_col,
               ifelse(!is.na(id_col) & nzchar(id_col), id_col, "<unnamed>")),
        " @ ", hits$seqname, ":", hits$start, "\u2013", hits$end
      )
      cli::cli_abort(c(
        "{nrow(hits)} features matched {.val {gene}}.",
        "i" = "Refine your query or use {.arg contig}, {.arg start}, \
               {.arg end} directly.",
        "i" = "Matches:",
        setNames(labels, rep("*", length(labels)))
      ))
    }

    feat_contig <- as.character(hits$seqname[1L])
    feat_start  <- hits$start[1L]
    feat_end    <- hits$end[1L]

    if (is.null(flank)) flank <- 5000L
    plot_start <- max(1L, feat_start - flank)
    plot_end   <- feat_end + flank

    if (is.null(title))
      title <- paste0("Region around ", gene)

    roi <- get_roi_features(entity,
                            contig       = feat_contig,
                            start        = plot_start,
                            end          = plot_end,
                            feature_type = feature_type)

    p <- plot_roi(roi,
                  title        = title,
                  arrow_height = arrow_height,
                  head_prop    = head_prop,
                  neck_prop    = neck_prop,
                  label_size   = label_size,
                  colors       = colors)

    cli::cli_inform(c(
      "i" = "Gene {.val {gene}} spans \
             {feat_contig}:{feat_start}\u2013{feat_end}.",
      "i" = "Plotted with {.arg flank} = {flank} bp on each side.",
      "i" = "To zoom out: \
             {.code plot_roi(entity, gene = \"{gene}\", flank = {flank * 2L})}"
    ))

    return(p)
  }

  # ---------------------------------------------------------------------------
  # Coordinate mode
  # ---------------------------------------------------------------------------
  if (is.null(start) || is.null(end)) {
    cli::cli_abort(
      "Both {.arg start} and {.arg end} must be provided in coordinate mode."
    )
  }

  if (is.null(contig)) {
    contig <- get_contig_names(entity)[1L]
    cli::cli_inform(
      "{.arg contig} not specified; using {.val {contig}} (first contig)."
    )
  }

  if (is.null(flank)) flank <- 0L

  if (is.null(title))
    title <- paste0("Region around ", contig, ":", start, "\u2013", end)

  roi <- get_roi_features(entity,
                          contig       = contig,
                          start        = start - flank,
                          end          = end   + flank,
                          feature_type = feature_type)

  plot_roi(roi,
           title        = title,
           arrow_height = arrow_height,
           head_prop    = head_prop,
           neck_prop    = neck_prop,
           label_size   = label_size,
           colors       = colors)
}

#' @rdname plot_roi
#' @export
plot_roi.default <- function(roi, ...) {
  cli::cli_abort(c(
    "{.fn plot_roi} does not know how to handle \
     {.cls {class(roi)[1]}} objects.",
    "i" = "Pass a {.cls data.frame} (or {.cls GRanges}) from \
           {.fn get_roi_features}, or a {.cls genome_entity}."
  ))
}
