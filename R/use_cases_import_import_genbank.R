#' Import GenBank Use Case
#'
#' @description
#' Business logic for importing GenBank format data into a genome_entity.
#' This use case orchestrates the transformation of raw GenBank data
#' (provided by a gateway) into a validated genome_entity object.
#'
#' Gateway Contract:
#' The gateway must provide a read(path) method that returns a list of records,
#' where each record is a list with:
#' - metadata: list with fields (locus, accession, length_bp, mol_type, topology, etc.)
#' - features: data.frame with columns (type, start, end, strand, gene, locus_tag, product, etc.)
#' - sequence: character string
#'
#' @param gateway Gateway object with read(path) method
#' @param file_path Path to GenBank file
#' @param options List of options (currently unused, for future extensibility)
#'
#' @return A validated genome_entity object
#' @export
#' @examples
#' \dontrun{
#' # With a real gateway
#' parser <- get_genbank_parser()
#' gateway <- create_genbank_gateway(parser)
#' entity <- execute_import_genbank(gateway, "data.gbk")
#'
#' # With a mock gateway (for testing)
#' mock_gateway <- list(
#'   read = function(path) {
#'     list(
#'       list(
#'         metadata = list(locus = "TEST", accession = "TEST001", length_bp = 100),
#'         features = data.frame(type = "gene", start = 1, end = 100),
#'         sequence = "ATCG"
#'       )
#'     )
#'   }
#' )
#' entity <- execute_import_genbank(mock_gateway, "test_ampC.gbk")
#' }
execute_import_genbank <- function(gateway, file_path, options = list()) {
  # Step 1: Read raw data via gateway
  gbk_list <- gateway$read(file_path)

  if (!is.list(gbk_list) || length(gbk_list) == 0) {
    cli::cli_abort("Gateway returned empty or invalid data")
  }

  # Step 2: Transform to entity structure
  # Initialize containers
  all_features <- list()
  records_list <- list()
  seqname_map <- character()

  # Process each record
  for (i in seq_along(gbk_list)) {
    rec <- gbk_list[[i]]
    meta <- rec$metadata

    # Determine seqname (use accession if available, else locus, else generic name)
    seqname <- if (!is.null(meta$accession) && !is.na(meta$accession) && nzchar(meta$accession)) {
      meta$accession
    } else if (!is.null(meta$locus) && !is.na(meta$locus) && nzchar(meta$locus)) {
      meta$locus
    } else {
      paste0("record_", i)
    }

    seqname_map[i] <- seqname

    # Build records data.frame row
    records_list[[i]] <- data.frame(
      record_id = i,
      seqname = seqname,
      locus = if (!is.null(meta$locus) && !is.na(meta$locus)) meta$locus else NA_character_,
      length_bp = if (!is.null(meta$length_bp) && !is.na(meta$length_bp)) meta$length_bp else nchar(rec$sequence),
      mol_type = if (!is.null(meta$mol_type) && !is.na(meta$mol_type)) meta$mol_type else NA_character_,
      topology = if (!is.null(meta$topology) && !is.na(meta$topology)) meta$topology else NA_character_,
      division = if (!is.null(meta$division) && !is.na(meta$division)) meta$division else NA_character_,
      date = if (!is.null(meta$date) && !is.na(meta$date)) meta$date else NA_character_,
      definition = if (!is.null(meta$definition) && !is.na(meta$definition)) meta$definition else NA_character_,
      accession = if (!is.null(meta$accession) && !is.na(meta$accession)) meta$accession else NA_character_,
      version = if (!is.null(meta$version) && !is.na(meta$version)) meta$version else NA_character_,
      organism = if (!is.null(meta$organism) && !is.na(meta$organism)) meta$organism else NA_character_,
      taxonomy = if (!is.null(meta$taxonomy) && !is.na(meta$taxonomy)) meta$taxonomy else NA_character_,
      stringsAsFactors = FALSE
    )

    # Process features
    if (!is.null(rec$features) && nrow(rec$features) > 0) {
      feat_df <- rec$features
      # Add record linkage
      feat_df$record_id <- i
      feat_df$seqname <- seqname

      # Ensure required columns exist
      if (!"start" %in% names(feat_df)) feat_df$start <- NA_integer_
      if (!"end" %in% names(feat_df)) feat_df$end <- NA_integer_
      if (!"strand" %in% names(feat_df)) feat_df$strand <- NA_character_
      if (!"type" %in% names(feat_df)) feat_df$type <- NA_character_
      if (!"gene" %in% names(feat_df)) feat_df$gene <- NA_character_
      if (!"locus_tag" %in% names(feat_df)) feat_df$locus_tag <- NA_character_
      if (!"product" %in% names(feat_df)) feat_df$product <- NA_character_

      all_features[[i]] <- feat_df
    }
  }

  # Step 3: Combine data
  records_df <- do.call(rbind, records_list)

  # Combine features
  if (length(all_features) > 0) {
    features_df <- do.call(rbind, all_features)
  } else {
    # Empty features
    features_df <- data.frame(
      record_id = integer(),
      seqname = character(),
      type = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      gene = character(),
      locus_tag = character(),
      product = character(),
      stringsAsFactors = FALSE
    )
  }

  # Build sequences (named character vector)
  dna_raw <- sapply(gbk_list, function(r) r$sequence)
  names(dna_raw) <- seqname_map

  # Step 4: Build indices
  indices_list <- execute_build_indices(features_df, records_df)

  # Step 5: Create entity
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = dna_raw,
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = features_df,
    metadata_df = records_df,
    indices_list = indices_list
  )

  # Step 5.5: Add source tracking for format conversion warnings
  attr(entity, "import_source") <- "genbank"

  # Step 6: Validate
  validate_genome_entity(entity)

  # Step 7: Return
  entity
}
