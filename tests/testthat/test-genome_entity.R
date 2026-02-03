test_that("genome_entity constructor works", {
  entity <- new_genome_entity()

  expect_s3_class(entity, "genome_entity")
  expect_true(all(c("sequences", "features", "metadata", "indices") %in% names(entity)))
})

test_that("genome_entity validation works", {
  entity <- new_genome_entity()

  expect_silent(validate_genome_entity(entity))

  # Test invalid entity
  bad_entity <- list()
  class(bad_entity) <- "genome_entity"
  expect_error(validate_genome_entity(bad_entity), "Missing required components")
})

test_that("genome_entity print method works", {
  entity <- new_genome_entity()

  expect_output(print(entity), "genome_entity")
  expect_output(print(entity), "Sequences:")
})

test_that("genome_entity summary method works", {
  entity <- new_genome_entity()

  expect_output(summary(entity), "genome_entity Summary")
  result <- summary(entity)
  expect_type(result, "list")
  expect_true("source" %in% names(result))
})
