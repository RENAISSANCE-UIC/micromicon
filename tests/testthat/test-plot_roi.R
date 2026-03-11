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

# =============================================================================
# genome_entity dispatch — fixtures
# =============================================================================

# Minimal feature table with two genes: "cheA" (CDS) and "cheA_dup" (CDS),
# plus a "gene" shadow for "cheA" to verify type=="gene" exclusion.
make_entity_features <- function() {
  data.frame(
    seqname = "U00096",
    start   = c(1973360L, 1973360L, 2000000L),
    end     = c(1975324L, 1975324L, 2001500L),
    strand  = c("-", "-", "+"),
    type    = c("CDS", "gene", "CDS"),
    Name    = c("cheA", "cheA", "cheA_dup"),
    ID      = c("cds-cheA", "gene-cheA", "cds-cheA_dup"),
    stringsAsFactors = FALSE
  )
}

make_test_entity <- function(feats = make_entity_features()) {
  new_genome_entity(
    features_df  = feats,
    indices_list = list(
      seqnames        = "U00096",
      locus_tag_index = integer(),
      gene_index      = list(),
      cds_hash        = new.env(hash = TRUE, parent = emptyenv())
    )
  )
}

# =============================================================================
# plot_roi.default
# =============================================================================

test_that("plot_roi.default errors with a helpful message", {
  expect_error(plot_roi(42L), class = "rlang_error")
  expect_error(plot_roi(list(a = 1)), class = "rlang_error")
})

# =============================================================================
# genome_entity dispatch — argument validation
# =============================================================================

test_that("providing neither gene nor coords raises an error", {
  entity <- make_test_entity()
  expect_error(plot_roi(entity), class = "rlang_error")
})

test_that("providing gene AND coords raises an error", {
  entity <- make_test_entity()
  expect_error(
    plot_roi(entity, gene = "cheA", start = 1000L, end = 2000L),
    class = "rlang_error"
  )
})

test_that("coordinate mode with only start (missing end) raises an error", {
  entity <- make_test_entity()
  expect_error(
    plot_roi(entity, start = 1000L),
    class = "rlang_error"
  )
})

test_that("coordinate mode with only end (missing start) raises an error", {
  entity <- make_test_entity()
  expect_error(
    plot_roi(entity, end = 2000L),
    class = "rlang_error"
  )
})

# =============================================================================
# genome_entity dispatch — gene search errors
# =============================================================================

test_that("gene with zero matches raises an error", {
  entity <- make_test_entity()
  expect_error(
    plot_roi(entity, gene = "nonexistent_gene_xyz"),
    class = "rlang_error"
  )
})

test_that("zero-match error message suggests search_features", {
  entity <- make_test_entity()
  expect_error(
    plot_roi(entity, gene = "nonexistent_gene_xyz"),
    regexp = "search_features"
  )
})

test_that("gene matching 'gene'-type rows only raises zero-match error", {
  # Feature table has one row with type=="gene" and Name=="cheA"
  # After filtering type=="gene", zero rows remain — should error.
  gene_only_feats <- data.frame(
    seqname = "U00096",
    start   = 1973360L,
    end     = 1975324L,
    strand  = "-",
    type    = "gene",
    Name    = "cheA",
    ID      = "gene-cheA",
    stringsAsFactors = FALSE
  )
  entity <- make_test_entity(gene_only_feats)
  expect_error(plot_roi(entity, gene = "cheA"), class = "rlang_error")
})

test_that("multiple CDS matches raises an error listing the hits", {
  entity <- make_test_entity()
  err <- tryCatch(
    plot_roi(entity, gene = "cheA"),
    error = function(e) e
  )
  expect_s3_class(err, "rlang_error")
  # Error message should mention the count of matches
  expect_match(conditionMessage(err), "2 features matched")
  # Error message should show hit coordinates
  expect_match(conditionMessage(err), "cheA")
})

test_that("multiple-match error shows feature type in output", {
  entity <- make_test_entity()
  err <- tryCatch(
    plot_roi(entity, gene = "cheA"),
    error = function(e) e
  )
  expect_match(conditionMessage(err), "CDS")
})

# =============================================================================
# genome_entity dispatch — gene mode (single match, mocked get_roi_features)
# =============================================================================

test_that("gene mode returns ggplot for a single-match gene", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  p <- plot_roi(entity, gene = "cheA_dup")
  expect_s3_class(p, "ggplot")
})

test_that("gene mode auto-populates title with gene name", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  p <- plot_roi(entity, gene = "cheA_dup")
  expect_match(p$labels$title, "cheA_dup")
})

test_that("gene mode respects explicit title override", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  p <- plot_roi(entity, gene = "cheA_dup", title = "My custom title")
  expect_equal(p$labels$title, "My custom title")
})

test_that("gene mode emits a flank message", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  expect_message(
    plot_roi(entity, gene = "cheA_dup"),
    regexp = "flank"
  )
})

test_that("gene mode default flank is 5000", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()
  captured_args <- list()

  local_mocked_bindings(
    get_roi_features = function(entity, contig, start, end, ...) {
      captured_args <<- list(start = start, end = end)
      make_sample_roi()
    },
    .package = "micromicon"
  )

  suppressMessages(plot_roi(entity, gene = "cheA_dup"))

  feat_start <- make_entity_features()$start[3L]  # cheA_dup
  feat_end   <- make_entity_features()$end[3L]
  expect_equal(captured_args$start, max(1L, feat_start - 5000L))
  expect_equal(captured_args$end,   feat_end + 5000L)
})

test_that("gene mode respects explicit flank override", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()
  captured_args <- list()

  local_mocked_bindings(
    get_roi_features = function(entity, contig, start, end, ...) {
      captured_args <<- list(start = start, end = end)
      make_sample_roi()
    },
    .package = "micromicon"
  )

  suppressMessages(plot_roi(entity, gene = "cheA_dup", flank = 1000L))

  feat_start <- make_entity_features()$start[3L]
  feat_end   <- make_entity_features()$end[3L]
  expect_equal(captured_args$start, max(1L, feat_start - 1000L))
  expect_equal(captured_args$end,   feat_end + 1000L)
})

# =============================================================================
# genome_entity dispatch — coordinate mode (mocked get_roi_features)
# =============================================================================

test_that("coordinate mode returns ggplot", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  p <- suppressMessages(
    plot_roi(entity, contig = "U00096", start = 1970000L, end = 1985000L)
  )
  expect_s3_class(p, "ggplot")
})

test_that("coordinate mode auto-populates title with contig:start-end", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  p <- suppressMessages(
    plot_roi(entity, contig = "U00096", start = 1970000L, end = 1985000L)
  )
  expect_match(p$labels$title, "U00096")
  expect_match(p$labels$title, "1970000")
  expect_match(p$labels$title, "1985000")
})

test_that("coordinate mode with no contig emits a notification", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()

  local_mocked_bindings(
    get_roi_features = function(...) make_sample_roi(),
    .package = "micromicon"
  )

  expect_message(
    plot_roi(entity, start = 1970000L, end = 1985000L),
    regexp = "contig"
  )
})

test_that("coordinate mode default flank is 0", {
  skip_if_not_installed("ggplot2")
  entity <- make_test_entity()
  captured_args <- list()

  local_mocked_bindings(
    get_roi_features = function(entity, contig, start, end, ...) {
      captured_args <<- list(start = start, end = end)
      make_sample_roi()
    },
    .package = "micromicon"
  )

  suppressMessages(
    plot_roi(entity, contig = "U00096", start = 1970000L, end = 1985000L)
  )
  expect_equal(captured_args$start, 1970000L)
  expect_equal(captured_args$end,   1985000L)
})
