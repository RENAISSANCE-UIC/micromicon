# use_cases_visualization_plot_roi.R
# Internal implementation for plot_roi().
# Public API lives in controllers_plot_controller.R.

#' Internal worker for plot_roi
#'
#' All geometry, repulsion, and ggplot assembly resides here. `plot_roi()` is
#' the validated, documented user-facing entry point.
#'
#' @keywords internal
#' @noRd
plot_roi_impl <- function(roi,
                          title        = NULL,
                          arrow_height = 0.18,
                          head_prop    = 0.35,
                          neck_prop    = 0.60,
                          label_size   = 4.0,
                          colors       = NULL) {

  # Geometry helpers -----------------------------------------------------------

  .arrow_pts <- function(x0, x1, y, h, head_prop, neck_prop, dir) {
    L      <- x1 - x0
    head   <- head_prop * L
    x_neck <- x1 - head
    h_tip  <- h / 2
    h_body <- h_tip * neck_prop
    pts <- data.frame(
      x = c(x0,       x_neck,  x_neck,  x1, x_neck,  x_neck,  x0,       x0),
      y = c(y+h_body, y+h_body, y+h_tip, y,  y-h_tip, y-h_body, y-h_body, y+h_body)
    )
    if (dir == "left") {
      xm    <- (x0 + x1) / 2
      pts$x <- 2 * xm - pts$x
    }
    pts
  }

  # Build flat data frame from roi ---------------------------------------------

  .to_plotdf <- function(roi) {
    df    <- as.data.frame(roi, stringsAsFactors = FALSE)
    label <- rep(NA_character_, nrow(df))
    for (col in c("Name", "gene", "Alias", "type")) {
      if (col %in% names(df)) {
        v     <- as.character(df[[col]])
        label <- ifelse(is.na(label) | label == "", v, label)
      }
    }
    label[is.na(label)] <- "?"
    data.frame(
      seqnames   = as.character(df$seqnames),
      gene       = label,
      type       = if ("type" %in% names(df)) as.character(df$type) else "feature",
      xmin       = df$start,
      xmax       = df$end,
      mid        = (df$start + df$end) / 2,
      strand_num = ifelse(df$strand == "+", 1L, -1L),
      y          = 0,
      stringsAsFactors = FALSE
    )
  }

  # Polygon assembly -----------------------------------------------------------

  .make_polys <- function(plotdf, h, head_prop, neck_prop) {
    rows <- lapply(seq_len(nrow(plotdf)), function(i) {
      dir <- if (plotdf$strand_num[i] >= 0) "right" else "left"
      p   <- .arrow_pts(plotdf$xmin[i], plotdf$xmax[i], plotdf$y[i],
                        h, head_prop, neck_prop, dir)
      p$poly_id <- i
      p$gene    <- plotdf$gene[i]
      p$type    <- plotdf$type[i]
      p
    })
    do.call(rbind, rows)
  }

  # 2-D label repulsion --------------------------------------------------------
  #
  # Prefer x-movement (clamped to plot bounds). When x cannot resolve overlap
  # even after clamping, bump the label displaced further from its origin to
  # the next y-tier and reset its x to the natural midpoint so it resolves
  # cleanly within the new tier.

  .repel_labels <- function(x0, widths, x_lo, x_hi, y_base, y_step,
                             max_iter = 500) {
    n    <- length(x0)
    if (n <= 1L) return(data.frame(lx = x0, ly = rep(y_base, n)))
    half <- widths / 2
    x    <- x0
    y    <- rep(y_base, n)

    for (iter in seq_len(max_iter)) {
      moved <- FALSE
      for (cur_y in sort(unique(y))) {
        idx  <- which(y == cur_y)
        if (length(idx) < 2L) next
        sord <- idx[order(x[idx])]

        for (k in seq_len(length(sord) - 1L)) {
          a <- sord[k];  b <- sord[k + 1L]
          overlap <- (half[a] + half[b]) - (x[b] - x[a])
          if (overlap <= 1e-9) next

          push <- overlap / 2 + 1e-9
          xa   <- max(x[a] - push, x_lo + half[a])
          xb   <- min(x[b] + push, x_hi - half[b])

          if ((half[a] + half[b]) - (xb - xa) > 1e-9) {
            # x exhausted — bump label with greater displacement from origin
            if (abs(xb - x0[b]) >= abs(xa - x0[a])) {
              y[b] <- cur_y + y_step;  x[b] <- x0[b]
            } else {
              y[a] <- cur_y + y_step;  x[a] <- x0[a]
            }
          } else {
            x[a] <- xa;  x[b] <- xb
          }
          moved <- TRUE
        }
      }
      if (!moved) break
    }
    data.frame(lx = x, ly = y)
  }

  # Axis tick formatter  (bp / kbp / Mbp) -------------------------------------

  .bp_fmt <- function(x) {
    r <- diff(range(x, na.rm = TRUE))
    if (r >= 1e6) paste0(round(x / 1e6, 2), " Mbp")
    else if (r >= 1e3) paste0(round(x / 1e3, 1), " kbp")
    else paste0(x, " bp")
  }

  # ggplot2-style HCL palette (base R only) ------------------------------------

  .default_pal <- function(n) {
    if (n == 0L) return(character(0L))
    hues <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
    grDevices::hcl(h = hues, c = 100, l = 65)
  }

  # Main -----------------------------------------------------------------------

  plotdf     <- .to_plotdf(roi)
  arrow_poly <- .make_polys(plotdf, arrow_height, head_prop, neck_prop)

  # Colour map: one colour per gene row
  genes   <- plotdf$gene
  n_genes <- nrow(plotdf)
  pal     <- if (is.null(colors)) .default_pal(n_genes) else rep_len(colors, n_genes)
  if (!is.null(names(pal))) {
    color_map           <- setNames(pal[genes], genes)
    na_idx              <- is.na(color_map)
    color_map[na_idx]   <- .default_pal(sum(na_idx))
  } else {
    color_map <- setNames(pal, genes)
  }

  # Backbone bounds
  R     <- max(plotdf$xmax) - min(plotdf$xmin)
  x_pad <- R * 0.05
  x_lo  <- min(plotdf$xmin) - x_pad
  x_hi  <- max(plotdf$xmax) + x_pad
  bb_h  <- arrow_height * 0.10

  # Label repulsion
  bp_per_char <- R / 70 * (label_size / 3.2)
  label_w     <- nchar(plotdf$gene) * bp_per_char
  y_base      <- arrow_height / 2 + arrow_height * 0.55
  y_step      <- arrow_height * 0.65

  repelled <- .repel_labels(plotdf$mid, label_w, x_lo, x_hi,
                             y_base, y_step, max_iter = 500)

  label_df <- data.frame(
    lx    = repelled$lx,
    ly    = repelled$ly,
    mid_x = plotdf$mid,
    gene  = plotdf$gene
  )

  leader_threshold      <- bp_per_char * 0.10
  label_df$draw_leader  <- abs(label_df$lx - label_df$mid_x) > leader_threshold |
                           label_df$ly > y_base

  # y limits driven by actual label positions
  y_lo <- -(arrow_height / 2) - arrow_height * 0.25
  y_hi <- max(label_df$ly) + label_size * 0.012 * arrow_height * 8

  # ggplot assembly ------------------------------------------------------------
  ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = x_lo, xmax = x_hi,
                      ymin = -bb_h, ymax = bb_h,
                      fill = "grey78", color = NA) +
    ggplot2::geom_polygon(
      data = arrow_poly,
      ggplot2::aes(x = x, y = y, group = poly_id, fill = gene),
      color = "grey20", linewidth = 0.25
    ) +
    ggplot2::geom_segment(
      data = label_df[label_df$draw_leader, , drop = FALSE],
      ggplot2::aes(x    = mid_x, xend = lx,
                   y    = arrow_height / 2 + arrow_height * 0.08,
                   yend = ly - arrow_height * 0.08),
      color = "grey55", linewidth = 0.25, lineend = "round"
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = lx, y = ly, label = gene),
      size = label_size, vjust = 0, hjust = 0.5
    ) +
    ggplot2::scale_fill_manual(values = color_map, guide = "none") +
    ggplot2::scale_x_continuous(labels = .bp_fmt) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi),
                             expand = FALSE) +
    ggplot2::labs(
      x     = paste0("Position (", unique(plotdf$seqnames)[1L], ")"),
      y     = NULL,
      title = title
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.line.y        = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position    = "right"
    )
}
