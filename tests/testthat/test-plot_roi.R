# test-plot_roi.R — testthat tests for plot_roi() / plot_roi_impl()

skip_if_not_installed("ggplot2")

# =============================================================================
# Shared fixture
# =============================================================================

make_sample_roi <- function() {
  data.frame(
    seqnames = "U00096",
    start    = c(1971030, 1972836, 1973360, 1975329, 1976252,
                 1977266, 1977847, 1978518, 1979753, 1980188, 1981587),
    end      = c(1972691, 1973339, 1975324, 1976255, 1977139,
                 1977844, 1978197, 1978966, 1980181, 1981612, 1982387),
    strand   = c("-", "-", "-", "-", "-", "-", "-", "-", "+", "-", "-"),
    type     = "CDS",
    Name     = c("tar", "cheW", "cheA", "motB", "motA",
                 "flhC", "flhD", "pgaptmp_001971", "uspC", "otsA", "otsB"),
    Alias    = c("pgaptmp_001964", "pgaptmp_001965", "pgaptmp_001966",
                 "pgaptmp_001967", "pgaptmp_001968", "pgaptmp_001969",
                 "pgaptmp_001970", "pgaptmp_001971", "pgaptmp_001972",
                 "pgaptmp_001973", "pgaptmp_001974"),
    Note     = c(
      "methyl-accepting chemotaxis protein II",
      "chemotaxis protein CheW",
      "chemotaxis protein CheA",
      "flagellar motor protein MotB",
      "flagellar motor stator protein MotA",
      "flagellar transcriptional regulator FlhC",
      "flagellar transcriptional regulator FlhD",
      "IS1-like element IS1A family transposase",
      "universal stress protein UspC",
      "alpha,alpha-trehalose-phosphate synthase",
      "trehalose-phosphatase"
    ),
    stringsAsFactors = FALSE
  )
}

roi <- make_sample_roi()

# =============================================================================
# Return type
# =============================================================================

test_that("plot_roi() returns a ggplot object", {
  p <- plot_roi(roi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_roi() builds without error on the sample ROI", {
  expect_no_error(ggplot2::ggplot_build(plot_roi(roi)))
})

# =============================================================================
# Polygon geometry
# =============================================================================

test_that("arrow polygon data has 8 vertices x n_features rows", {
  p      <- plot_roi(roi)
  poly_d <- p$layers[[2]]$data          # layer 2 = geom_polygon
  n_feat <- nrow(roi)
  expect_equal(nrow(poly_d), n_feat * 8L)
})

test_that("each gene gets a unique poly_id", {
  p      <- plot_roi(roi)
  poly_d <- p$layers[[2]]$data
  expect_equal(length(unique(poly_d$poly_id)), nrow(roi))
})

# =============================================================================
# Colour handling
# =============================================================================

test_that("unnamed color vector is recycled to n_features", {
  two_colors <- c("#FF0000", "#0000FF")
  expect_no_error(plot_roi(roi, colors = two_colors))
})

test_that("named color vector maps correctly", {
  named <- c(tar = "#FF0000", cheW = "#00FF00")
  p <- plot_roi(roi, colors = named)
  scale_vals <- p$scales$scales[[1]]$palette(nrow(roi))
  expect_true(length(scale_vals) == nrow(roi))
})

test_that("NULL colors yields a ggplot without error", {
  expect_s3_class(plot_roi(roi, colors = NULL), "ggplot")
})

# =============================================================================
# Parameter pass-through
# =============================================================================

test_that("title is set in the plot", {
  p <- plot_roi(roi, title = "motA region")
  expect_equal(p$labels$title, "motA region")
})

test_that("arrow_height changes polygon y-extent", {
  p_small <- plot_roi(roi, arrow_height = 0.05)
  p_large <- plot_roi(roi, arrow_height = 0.40)
  y_range_small <- diff(range(p_small$layers[[2]]$data$y))
  y_range_large <- diff(range(p_large$layers[[2]]$data$y))
  expect_lt(y_range_small, y_range_large)
})

# =============================================================================
# Input validation errors
# =============================================================================

test_that("NULL roi raises an error", {
  expect_error(plot_roi(NULL), class = "rlang_error")
})

test_that("empty data frame raises an error", {
  expect_error(
    plot_roi(data.frame(seqnames = character(), start = integer(),
                        end = integer(), strand = character())),
    class = "rlang_error"
  )
})

test_that("missing required column raises an error", {
  bad <- roi[, setdiff(names(roi), "strand")]
  expect_error(plot_roi(bad), class = "rlang_error")
})

test_that("non-numeric start/end raises an error", {
  bad        <- roi
  bad$start  <- as.character(bad$start)
  expect_error(plot_roi(bad), class = "rlang_error")
})

test_that("non-positive arrow_height raises an error", {
  expect_error(plot_roi(roi, arrow_height = -0.1), class = "rlang_error")
  expect_error(plot_roi(roi, arrow_height =  0),   class = "rlang_error")
})

test_that("non-character colors raises an error", {
  expect_error(plot_roi(roi, colors = 1:3), class = "rlang_error")
})

# =============================================================================
# Warnings
# =============================================================================

test_that("head_prop > 1 raises a warning", {
  expect_warning(plot_roi(roi, head_prop = 1.5))
})

test_that("neck_prop > 1 raises a warning", {
  expect_warning(plot_roi(roi, neck_prop = 1.2))
})

# =============================================================================
# Edge cases
# =============================================================================

test_that("single-feature ROI plots without error", {
  expect_s3_class(plot_roi(roi[1L, ]), "ggplot")
})

test_that("ROI with missing Name falls back to Alias then type", {
  no_name        <- roi[1L, ]
  no_name$Name   <- NA_character_
  p              <- plot_roi(no_name)
  label_data     <- p$layers[[4]]$data   # geom_text layer (layers: rect, polygon, segment, text)
  expect_equal(label_data$gene, "pgaptmp_001964")
})

test_that("ROI with mixed strands plots without error", {
  mixed <- roi
  mixed$strand <- rep_len(c("+", "-"), nrow(mixed))
  expect_s3_class(plot_roi(mixed), "ggplot")
})
