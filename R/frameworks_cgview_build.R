# cgview_build.R
# Build a CGView-native JSON payload (as an R list) from tidy data frames.
#
# Design principles:
#   - R emits maximally neutral JSON: features, legend colors, track types,
#     and backbone/ruler layout. NO sequence data (suppresses GC auto-tracks).
#   - The JS binding (cgview.js) is the single source of truth for visual
#     styling such as per-track thicknessRatio.
#   - Feature color is driven by legend item swatchColor, which CGView uses
#     as the feature fill color when a feature's `legend` field matches.
#
# Track routing uses the `source` field on features:
#   CDS features              -> source = "CDS"   -> "CDS" track
#   tRNA/rRNA/tmRNA/repeat_*  -> source = "ncRNA" -> "ncRNA" track
#   Other feature types       -> source = type     -> their own track
#   Mutation events           -> source = "mutations" -> "Variants" track

# ---- Default colour palette ------------------------------------------------

.CGVIEW_PALETTE <- list(
  # Genomic features
  CDS           = "#4e79a7",   # Tableau blue
  tRNA          = "#8b6bd6",   # purple
  rRNA          = "#2ecc71",   # teal-green
  tmRNA         = "#e74c3c",   # red
  repeat_region = "#f39c12",   # amber
  gene          = "#95a5a6",   # grey
  other         = "#7f8c8d",   # dark grey
  # Mutation event types (warm palette, distinct from feature blues)
  SNP           = "#ff7f0e",   # orange
  DEL           = "#c0392b",   # crimson
  INS           = "#27ae60",   # green
  SUB           = "#e91e63",   # pink
  AMP           = "#1abc9c",   # teal
  MOB           = "#7d3c98"    # violet
)

# Feature types that get routed to the ncRNA track (source = "ncRNA")
.CGVIEW_NCRNA_TYPES <- c("tRNA", "rRNA", "tmRNA", "repeat_region")

# Feature types that should get arrow decoration
.CGVIEW_ARROW_TYPES <- c("CDS", "gene")

# ---- Public builder --------------------------------------------------------

#' Build a CGView JSON payload from tidy data frames
#'
#' @description
#' Produces the R list that CGView.js expects (top-level key `cgview`).
#' No GC tracks are emitted -- sequence length only, no raw DNA string --
#' so CGView will not auto-generate GC Content / GC Skew rings.
#'
#' **Track routing by `source` field:**
#' \describe{
#'   \item{`"CDS"` track}{Features whose `type` is `"CDS"` (or other types in
#'     `arrow_types`). Position `"both"`, strand-split.}
#'   \item{`"ncRNA"` track}{Features whose `type` is `"tRNA"`, `"rRNA"`,
#'     `"tmRNA"`, or `"repeat_region"`. Position `"both"`, strand-split.}
#'   \item{`"Variants"` track}{Mutation events. Position `"inside"`.}
#' }
#' Any feature type not in the above groups gets its own inside track.
#'
#' @param genome_length Integer. Genome length in base pairs.
#' @param genome_name Character. Sequence label shown in the viewer.
#' @param features_df Data frame of genomic features, or `NULL`. Expected
#'   columns: `start`, `end` (or `stop`), `type`, `strand` (`"+"` / `"-"`),
#'   and optionally `Name`, `Alias`, `ID`, `gene`, `locus_tag`, `product`.
#' @param mutations_df Data frame of mutation events, or `NULL`. Expected
#'   columns: `type` (e.g. `"SNP"`, `"DEL"`), `pos` (or `position`), and
#'   optionally `size` (bp span; defaults to 1).
#' @param feature_types Character vector of feature types to include.
#'   `NULL` includes everything present in `features_df`.
#' @param palette Named list mapping type names to hex colour strings.
#'   Merged over the built-in defaults.
#' @param mark_width_bp Integer width (bp) of mutation markers on the map.
#'   Defaults to `max(1, genome_length %/% 500)` (0.2\% of genome).
#' @param ruler_spacing Integer. Ruler tick spacing in bp.
#'   Defaults to a round number near 1/5 of genome length.
#'
#' @return A named R list with top-level key `cgview`, ready for
#'   [jsonlite::toJSON()] or \code{beta_cgview()}.
#'
#' @keywords internal
#' @noRd
build_cgview_json <- function(genome_length,
                               genome_name   = "Genome",
                               features_df   = NULL,
                               mutations_df  = NULL,
                               feature_types = NULL,
                               palette       = NULL,
                               mark_width_bp = NULL,
                               ruler_spacing = NULL) {

  pal <- .CGVIEW_PALETTE
  if (!is.null(palette)) pal[names(palette)] <- palette

  genome_length <- as.integer(genome_length)

  # Default mark width: 0.2% of genome (visible as a thin spike)
  mark_w <- if (is.null(mark_width_bp))
    max(1L, genome_length %/% 500L)
  else
    as.integer(mark_width_bp)

  # tickCount: CGView uses d3.ticks(start, stop, tickCount) which produces
  # "nice" round bp values automatically.  We target one tick per
  # ruler_spacing interval; for a 4.79 Mbp genome with spacing = 1 Mbp
  # that gives tickCount ~= 5, and d3 places ticks at 1, 2, 3, 4 Mbp.
  # NOTE: ruler$spacing is a *pixel* gap (default 2 px), NOT a bp interval --
  # do NOT pass ruler_sp as spacing or the ruler rings go millions of px
  # off-canvas.
  ruler_sp <- if (is.null(ruler_spacing))
    .cgview_nice_interval(genome_length)
  else
    as.integer(ruler_spacing)
  tick_count <- as.integer(max(4L, min(10L, round(genome_length / ruler_sp))))

  # ---- Filter features -------------------------------------------------------
  if (!is.null(features_df) && !is.null(feature_types)) {
    features_df <- features_df[features_df$type %in% feature_types, ,
                                drop = FALSE]
  }
  if (!is.null(features_df) && nrow(features_df) == 0) features_df <- NULL

  # ---- Normalise mutations_df column names -----------------------------------
  if (!is.null(mutations_df)) {
    if (!"pos" %in% names(mutations_df) && "position" %in% names(mutations_df))
      mutations_df$pos <- mutations_df$position
    if (!"size" %in% names(mutations_df))
      mutations_df$size <- 1L
  }
  if (!is.null(mutations_df) && nrow(mutations_df) == 0) mutations_df <- NULL

  # ---- Build feature objects -------------------------------------------------
  feature_list   <- list()
  legend_seen    <- character(0)  # track which types have legend entries

  if (!is.null(features_df)) {
    feature_list <- .cgview_build_feature_objects(features_df, pal)
    legend_seen  <- c(legend_seen, unique(as.character(features_df$type)))
  }

  if (!is.null(mutations_df)) {
    mut_feats    <- .cgview_build_mutation_objects(mutations_df, mark_w, pal)
    feature_list <- c(feature_list, mut_feats)
    legend_seen  <- c(legend_seen, unique(as.character(mutations_df$type)))
  }

  # ---- Build track list ------------------------------------------------------
  track_list <- .cgview_build_tracks(features_df, mutations_df)

  # ---- Legend items ----------------------------------------------------------
  # decoration on the legend item is the primary control for glyph shape in
  # CGView 1.8 -- it overrides any feature-level decoration field.
  legend_items <- lapply(unique(legend_seen), function(nm) {
    item <- list(
      name        = nm,
      swatchColor = if (nm %in% names(pal)) pal[[nm]] else pal$other
    )
    if (nm %in% .CGVIEW_ARROW_TYPES) item$decoration <- "arrow"
    item
  })

  # ---- Assemble CGView-native JSON -------------------------------------------
  list(
    cgview = list(
      version  = "1.0",
      settings = list(
        format          = "circular",
        backgroundColor = "#FFFFFF",
        arrowHeadLength = 0.3    # controls arrowhead prominence (0--1); default is near-invisible
      ),
      sequence = list(
        length = genome_length,
        name   = genome_name
      ),
      backbone = list(
        color     = "#999999",
        thickness = 6L
      ),
      ruler = list(
        visible    = TRUE,
        tickCount  = tick_count,   # d3.ticks() -> nice round Mbp/kbp labels
        color      = "#888888",
        tickLength = 10L,
        tickWidth  = 2L,
        font       = "sans-serif, plain, 18"
      ),
      features = feature_list,
      legend   = list(
        position = "top-right",
        items    = legend_items
      ),
      tracks   = track_list
    )
  )
}

# ---- Internal builders -----------------------------------------------------

#' @keywords internal
.cgview_build_feature_objects <- function(df, pal) {
  lapply(seq_len(nrow(df)), function(i) {
    row    <- df[i, , drop = FALSE]
    ft     <- as.character(row$type)
    nm     <- .cgview_feature_name(row, ft)
    end_v  <- if ("end"  %in% names(row)) row$end  else
              if ("stop" %in% names(row)) row$stop else row$start
    strand <- .cgview_strand_int(if ("strand" %in% names(row)) row$strand else NA)

    # Arrow decoration for CDS/gene; arc for everything else.
    # CGView 1.8 recognises "arrow" and auto-orients it from strand.
    deco <- if (ft %in% .CGVIEW_ARROW_TYPES) "arrow" else "arc"

    # Route to track via source field
    source_key <- if (ft %in% .CGVIEW_NCRNA_TYPES) "ncRNA" else ft

    list(
      name       = nm,
      type       = ft,
      start      = as.integer(row$start),
      stop       = as.integer(end_v),
      strand     = strand,
      source     = source_key,
      legend     = ft,
      decoration = deco,
      showLabel  = FALSE
    )
  })
}

#' @keywords internal
.cgview_build_mutation_objects <- function(df, mark_w, pal) {
  lapply(seq_len(nrow(df)), function(i) {
    row  <- df[i, , drop = FALSE]
    mt   <- as.character(row$type)
    pos  <- as.integer(row$pos)
    sz   <- if (!is.na(row$size) && as.integer(row$size) > 1L)
              as.integer(row$size)
            else
              mark_w
    list(
      name       = paste0(mt, "@", pos),
      type       = mt,
      start      = pos,
      stop       = pos + sz,
      strand     = 1L,
      source     = "mutations",
      legend     = mt,
      decoration = "arc",
      showLabel  = FALSE
    )
  })
}

#' @keywords internal
.cgview_build_tracks <- function(features_df, mutations_df) {
  tracks <- list()

  if (!is.null(features_df) && nrow(features_df) > 0) {
    types_present <- unique(as.character(features_df$type))

    # CDS track -- straddles the backbone, strand-split
    cds_types <- intersect(types_present, .CGVIEW_ARROW_TYPES)
    if (length(cds_types) > 0) {
      tracks <- c(tracks, list(list(
        name               = "CDS",
        separateFeaturesBy = "strand",
        position           = "both",
        dataType           = "feature",
        dataMethod         = "source",
        dataKeys           = "CDS",
        showLabels         = FALSE
      )))
    }

    # ncRNA track -- outside the CDS band
    ncrna_types <- intersect(types_present, .CGVIEW_NCRNA_TYPES)
    if (length(ncrna_types) > 0) {
      tracks <- c(tracks, list(list(
        name               = "ncRNA",
        separateFeaturesBy = "strand",
        position           = "outside",
        dataType           = "feature",
        dataMethod         = "source",
        dataKeys           = "ncRNA",
        showLabels         = FALSE
      )))
    }

    # Any remaining feature types get their own inside track
    covered <- c(.CGVIEW_ARROW_TYPES, .CGVIEW_NCRNA_TYPES)
    for (ot in setdiff(types_present, covered)) {
      tracks <- c(tracks, list(list(
        name               = ot,
        separateFeaturesBy = "none",
        position           = "inside",
        dataType           = "feature",
        dataMethod         = "source",
        dataKeys           = ot,
        showLabels         = FALSE
      )))
    }
  }

  # Variants track -- outermost ring (outside backbone)
  if (!is.null(mutations_df) && nrow(mutations_df) > 0) {
    tracks <- c(tracks, list(list(
      name               = "Variants",
      separateFeaturesBy = "none",
      position           = "outside",
      dataType           = "feature",
      dataMethod         = "source",
      dataKeys           = "mutations",
      showLabels         = FALSE
    )))
  }

  tracks
}

# ---- Utility helpers -------------------------------------------------------

#' @keywords internal
.cgview_feature_name <- function(row, type_fallback) {
  for (col in c("Name", "gene", "Alias", "ID", "locus_tag", "name")) {
    v <- if (col %in% names(row)) row[[col]] else NA
    if (!is.null(v) && length(v) > 0 && !is.na(v) && nzchar(as.character(v)))
      return(as.character(v))
  }
  paste0(type_fallback, "_", as.integer(row$start))
}

#' @keywords internal
.cgview_strand_int <- function(s) {
  if (is.null(s) || length(s) == 0 || is.na(s)) return(1L)
  if (is.numeric(s)) return(if (s < 0) -1L else 1L)
  if (identical(as.character(s), "-") || identical(as.character(s), "-1")) -1L
  else 1L
}

#' @keywords internal
.cgview_nice_interval <- function(len) {
  # Return a round number approximately 1/5 of the genome length
  raw <- len / 5
  mag <- 10 ^ floor(log10(raw))
  cands <- c(1, 2, 5, 10) * mag
  as.integer(cands[which.min(abs(cands - raw))])
}
