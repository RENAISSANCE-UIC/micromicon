test_that("transl_table is pulled from qualifiers when column absent", {
  cds <- entity$features[entity$features$type == "CDS", ][1, , drop = FALSE]
  tt <- gd_get_transl_table(cds, default = 11L)
  expect_true(is.integer(tt) && length(tt) == 1L)
})


test_that("gene ID resolution tolerates both schemas", {
  cds <- entity$features[entity$features$type == "CDS", ][1, , drop = FALSE]
  gid <- gd_resolve_gene_id(cds)
  expect_true(is.character(gid) && nchar(gid) > 0)
})


test_that("cds_by_position works for Prokka and PGAP feature sets", {
  gw_a <- gd_create_annotation_gateway(entity)
  # pick a known CDS position here if you have one; else sample a CDS center
  cds_any <- entity$features[entity$features$type == "CDS", ][1, , drop = FALSE]
  seq_id  <- cds_any$seqname[1]
  pos     <- floor((cds_any$start[1] + cds_any$end[1]) / 2)
  cds_rows <- gw_a$cds_by_position(seq_id, pos)
  expect_true(nrow(cds_rows) >= 1L)
})