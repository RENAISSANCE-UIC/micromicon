#' Converter Functions (Transformation Layer)
#'
#' @description
#' Functions for converting between different genome data representations.
#' In Clean Architecture terms, these are transformation use cases.
#'
#' Note: In Phase 6, these functions remain in R/converters.R but should
#' eventually be moved to R/use_cases/transform/ for full Clean Architecture
#' compliance. They are used by controllers and other parts of the system.
#'
#' @name converters
#' @keywords internal
NULL

#' Convert GenBank list to genome_entity
#'
#' @description
#' Converts the output from read_gbk() (list of records) to a genome_entity object.
#'
#' @param gbk_list List of GenBank records from read_gbk()
#'
#' @return A genome_entity object
#' @export
#'
#' @examples
#' \dontrun{
#' gbk_list <- read_gbk("example.gbk")
#' entity <- gbk_to_entity(gbk_list)
#' }
gbk_to_entity <- function(gbk_list) {
  if (!is.list(gbk_list) || length(gbk_list) == 0) {
    cli::cli_abort("gbk_list must be a non-empty list")
  }

  # Initialize containers
  all_features <- list()
  records_list <- list()
  seqname_map <- character()

  # Process each record
  for (i in seq_along(gbk_list)) {
    rec <- gbk_list[[i]]
    meta <- rec$metadata

    # Determine seqname (use accession if available, else locus, else generic name)
    seqname <- if (!is.na(meta$accession) && nzchar(meta$accession)) {
      meta$accession
    } else if (!is.na(meta$locus) && nzchar(meta$locus)) {
      meta$locus
    } else {
      paste0("record_", i)
    }

    seqname_map[i] <- seqname

    # Build records data.frame row
    records_list[[i]] <- data.frame(
      record_id = i,
      seqname = seqname,
      locus = if (!is.na(meta$locus)) meta$locus else NA_character_,
      length_bp = if (!is.na(meta$length_bp)) meta$length_bp else nchar(rec$sequence),
      mol_type = if (!is.na(meta$mol_type)) meta$mol_type else NA_character_,
      topology = if (!is.na(meta$topology)) meta$topology else NA_character_,
      division = if (!is.na(meta$division)) meta$division else NA_character_,
      date = if (!is.na(meta$date)) meta$date else NA_character_,
      definition = if (!is.na(meta$definition)) meta$definition else NA_character_,
      accession = if (!is.na(meta$accession)) meta$accession else NA_character_,
      version = if (!is.na(meta$version)) meta$version else NA_character_,
      organism = if (!is.na(meta$organism)) meta$organism else NA_character_,
      taxonomy = if (!is.na(meta$taxonomy)) meta$taxonomy else NA_character_,
      replicon_type = if (!is.na(meta$replicon_type)) meta$replicon_type else NA_character_,
      replicon_name = if (!is.na(meta$replicon_name)) meta$replicon_name else NA_character_,
      stringsAsFactors = FALSE
    )

    # Process features
    if (nrow(rec$features) > 0) {
      feat_df <- rec$features
      # Add record linkage
      feat_df$record_id <- i
      feat_df$seqname <- seqname

      # Ensure all required columns exist
      if (!"start" %in% names(feat_df)) feat_df$start <- NA_integer_
      if (!"end" %in% names(feat_df)) feat_df$end <- NA_integer_
      if (!"strand" %in% names(feat_df)) feat_df$strand <- NA_integer_
      if (!"type" %in% names(feat_df)) feat_df$type <- NA_character_
      if (!"location_string" %in% names(feat_df)) feat_df$location_string <- NA_character_
      if (!"location_type" %in% names(feat_df)) feat_df$location_type <- NA_character_
      if (!"gene" %in% names(feat_df)) feat_df$gene <- NA_character_
      if (!"locus_tag" %in% names(feat_df)) feat_df$locus_tag <- NA_character_
      if (!"product" %in% names(feat_df)) feat_df$product <- NA_character_
      if (!"protein_id" %in% names(feat_df)) feat_df$protein_id <- NA_character_
      if (!"translation" %in% names(feat_df)) feat_df$translation <- NA_character_
      if (!"ranges" %in% names(feat_df)) feat_df$ranges <- I(vector("list", nrow(feat_df)))
      if (!"qualifiers" %in% names(feat_df)) feat_df$qualifiers <- I(vector("list", nrow(feat_df)))

      all_features[[i]] <- feat_df
    }
  }

  # Combine records
  records_df <- do.call(rbind, records_list)

  # Combine features
  if (length(all_features) > 0) {
    features_df <- do.call(rbind, all_features)
  } else {
    features_df <- data.frame(
      record_id = integer(),
      seqname = character(),
      type = character(),
      start = integer(),
      end = integer(),
      strand = integer(),
      location_string = character(),
      location_type = character(),
      gene = character(),
      locus_tag = character(),
      product = character(),
      protein_id = character(),
      translation = character(),
      stringsAsFactors = FALSE
    )
    features_df$ranges <- I(list())
    features_df$qualifiers <- I(list())
  }

  # Build sequences (named character vector)
  dna_raw <- sapply(gbk_list, function(r) r$sequence)
  names(dna_raw) <- seqname_map

  # Build indices
  seqnames_unique <- unique(seqname_map)
  seqname_to_record <- setNames(seq_along(seqname_map), seqname_map)

  # Locus tag index
  locus_tag_index <- integer()
  if (nrow(features_df) > 0 && "locus_tag" %in% names(features_df)) {
    has_lt <- !is.na(features_df$locus_tag) & nzchar(features_df$locus_tag)
    if (any(has_lt)) {
      locus_tag_index <- setNames(which(has_lt), features_df$locus_tag[has_lt])
    }
  }

  # Gene index
  gene_index <- list()
  if (nrow(features_df) > 0 && "gene" %in% names(features_df)) {
    has_gene <- !is.na(features_df$gene) & nzchar(features_df$gene)
    if (any(has_gene)) {
      gene_names <- unique(features_df$gene[has_gene])
      gene_index <- lapply(gene_names, function(g) {
        which(features_df$gene == g)
      })
      names(gene_index) <- gene_names
    }
  }

  # Create genome_entity
  entity <- new_genome_entity(
    sequences = list(
      dna_raw = dna_raw,
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features = list(
      df = features_df,
      granges = NULL
    ),
    metadata = list(
      records = records_df,
      source = "genbank",
      import_date = Sys.time(),
      gff_used = character(),
      fasta_used = character(),
      import_errors = list()
    ),
    indices = list(
      seqnames = seqnames_unique,
      seqname_to_record = seqname_to_record,
      locus_tag_index = locus_tag_index,
      gene_index = gene_index
    )
  )

  validate_genome_entity(entity)
  entity
}

#' Convert GFF3 + FASTA to genome_entity
#'
#' @description
#' Converts GFF3 and FASTA files to a genome_entity object.
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param auto_harmonize Logical; if TRUE, attempt to harmonize seqname mismatches
#' @param verbose Logical; if TRUE, print progress messages
#'
#' @return A genome_entity object
#' @export
#'
#' @examples
#' \dontrun{
#' entity <- gff_fasta_to_entity("annotation.gff3", "genome.fasta")
#' }
gff_fasta_to_entity <- function(gff_path, fasta_path, auto_harmonize = TRUE, verbose = TRUE) {

  if (!file.exists(gff_path)) {
    cli::cli_abort("GFF path not found: {gff_path}")
  }
  if (!file.exists(fasta_path)) {
    cli::cli_abort("FASTA path not found: {fasta_path}")
  }

  import_errors <- list()

  # Load FASTA
  if (verbose) cli::cli_inform("Loading FASTA sequence...")
  if (has_bioconductor()) {
    dna_bio <- Biostrings::readDNAStringSet(fasta_path)
    dna_raw <- as.character(dna_bio)
  } else {
    # Fallback: simple FASTA parser
    if (verbose) cli::cli_inform("Bioconductor not available; using simple FASTA parser")
    fasta_result <- .parse_fasta_simple(fasta_path)
    dna_raw <- fasta_result$sequences
    dna_bio <- NULL
  }

  # Load GFF3
  if (verbose) cli::cli_inform("Loading GFF3 annotation...")
  if (has_bioconductor()) {
    gff <- try(rtracklayer::import(gff_path, format = "gff3"), silent = TRUE)
    if (inherits(gff, "try-error")) {
      import_errors$gff_direct <- attr(gff, "condition")
      cli::cli_abort("Failed to import GFF3: {conditionMessage(import_errors$gff_direct)}")
    }
  } else {
    # Fallback: simple GFF3 parser
    if (verbose) cli::cli_inform("Bioconductor not available; using simple GFF3 parser")
    gff_result <- .parse_gff3_simple(gff_path)
    gff <- gff_result$features
  }

  # Convert GRanges to data.frame if needed
  if (has_bioconductor() && inherits(gff, "GRanges")) {
    features_df <- .granges_to_dataframe(gff)
  } else {
    features_df <- gff
  }

  # Build records metadata from FASTA headers and GFF
  seqnames_fasta <- names(dna_raw)
  records_list <- lapply(seq_along(seqnames_fasta), function(i) {
    seqname <- seqnames_fasta[i]
    data.frame(
      record_id = i,
      seqname = seqname,
      locus = seqname,
      length_bp = nchar(dna_raw[i]),
      mol_type = "DNA",
      topology = NA_character_,
      division = NA_character_,
      date = NA_character_,
      definition = NA_character_,
      accession = seqname,
      version = NA_character_,
      organism = NA_character_,
      taxonomy = NA_character_,
      replicon_type = NA_character_,
      replicon_name = seqname,
      stringsAsFactors = FALSE
    )
  })
  records_df <- do.call(rbind, records_list)

  # Build indices
  seqname_to_record <- setNames(seq_along(seqnames_fasta), seqnames_fasta)

  # Locus tag index
  locus_tag_index <- integer()
  if (nrow(features_df) > 0 && "locus_tag" %in% names(features_df)) {
    has_lt <- !is.na(features_df$locus_tag) & nzchar(features_df$locus_tag)
    if (any(has_lt)) {
      locus_tag_index <- setNames(which(has_lt), features_df$locus_tag[has_lt])
    }
  }

  # Gene index
  gene_index <- list()
  if (nrow(features_df) > 0 && "gene" %in% names(features_df)) {
    has_gene <- !is.na(features_df$gene) & nzchar(features_df$gene)
    if (any(has_gene)) {
      gene_names <- unique(features_df$gene[has_gene])
      gene_index <- lapply(gene_names, function(g) {
        which(features_df$gene == g)
      })
      names(gene_index) <- gene_names
    }
  }

  # Create genome_entity
  entity <- new_genome_entity(
    sequences = list(
      dna_raw = dna_raw,
      dna_bio = if (exists("dna_bio")) dna_bio else NULL,
      indexed_fa = NULL
    ),
    features = list(
      df = features_df,
      granges = if (has_bioconductor() && exists("gff") && inherits(gff, "GRanges")) gff else NULL
    ),
    metadata = list(
      records = records_df,
      source = "gff3_fasta",
      import_date = Sys.time(),
      gff_used = gff_path,
      fasta_used = fasta_path,
      import_errors = import_errors
    ),
    indices = list(
      seqnames = seqnames_fasta,
      seqname_to_record = seqname_to_record,
      locus_tag_index = locus_tag_index,
      gene_index = gene_index
    )
  )

  validate_genome_entity(entity)
  entity
}

#' Convert genome_entity to legacy GenBank list format
#'
#' @description
#' Converts a genome_entity back to the old read_gbk() list format
#' for backward compatibility.
#'
#' @param entity A genome_entity object
#'
#' @return List of GenBank records (old format)
#' @export
#'
#' @examples
#' \dontrun{
#' entity <- read_genome("example.gbk")
#' gbk_list <- entity_to_gbk_list(entity)
#' }
entity_to_gbk_list <- function(entity) {
  validate_genome_entity(entity)

  records_df <- entity$metadata$records
  features_df <- entity$features$df
  sequences <- entity$sequences$dna_raw

  # Split features by record_id
  gbk_list <- lapply(seq_len(nrow(records_df)), function(i) {
    rec <- records_df[i, ]

    # Extract features for this record
    if (nrow(features_df) > 0 && "record_id" %in% names(features_df)) {
      rec_features <- features_df[features_df$record_id == i, ]
      # Remove record_id and seqname columns
      rec_features$record_id <- NULL
      if ("seqname" %in% names(rec_features)) {
        rec_features$seqname <- NULL
      }
    } else {
      rec_features <- data.frame()
    }

    # Build metadata list
    metadata <- list(
      locus = rec$locus,
      length_bp = rec$length_bp,
      mol_type = rec$mol_type,
      topology = rec$topology,
      division = rec$division,
      date = rec$date,
      definition = rec$definition,
      accession = rec$accession,
      version = rec$version,
      organism = rec$organism,
      taxonomy = rec$taxonomy,
      replicon_type = rec$replicon_type,
      replicon_name = rec$replicon_name
    )

    # Get sequence
    seq <- sequences[rec$seqname]
    if (is.na(seq)) seq <- ""

    list(
      metadata = metadata,
      features = rec_features,
      sequence = seq
    )
  })

  gbk_list
}

#' Convert genome_entity to legacy genome_obj format
#'
#' @description
#' Converts a genome_entity to the old init_genome() object format
#' for backward compatibility.
#'
#' @param entity A genome_entity object
#'
#' @return Legacy genome_obj list
#' @export
#'
#' @examples
#' \dontrun{
#' entity <- read_genome(gff = "annotation.gff3", fasta = "genome.fasta")
#' genome_obj <- entity_to_legacy_genome_obj(entity)
#' }
entity_to_legacy_genome_obj <- function(entity) {
  validate_genome_entity(entity)

  if (!has_bioconductor()) {
    cli::cli_abort("Converting to legacy genome_obj requires Bioconductor packages")
  }

  # Get GRanges (create if needed)
  gff <- features(entity, format = "GRanges")

  # Get DNAStringSet (create if needed)
  fasta <- sequences(entity, format = "DNAStringSet")

  # Create FaFile if we have a fasta_used path (only present for GFF+FASTA workflows)
  fa <- NULL
  if ("fasta_used" %in% names(entity$metadata)) {
    fasta_used <- entity$metadata$fasta_used[1]
    if (!is.na(fasta_used) && nzchar(fasta_used)) {
      fasta_path <- fasta_used
      if (file.exists(fasta_path)) {
        # Index if needed
        if (!file.exists(paste0(fasta_path, ".fai"))) {
          Rsamtools::indexFa(fasta_path)
        }
        fa <- Rsamtools::FaFile(fasta_path)
      }
    }
  }

  # For multi-row metadata, take first row's values if columns exist
  # (these columns only exist for GFF+FASTA workflows, not GenBank)
  gff_used_val <- ""
  if ("gff_used" %in% names(entity$metadata) && nrow(entity$metadata) > 0) {
    gff_used_val <- entity$metadata$gff_used[1]
    if (is.na(gff_used_val)) gff_used_val <- ""
  }

  import_errs_val <- ""
  if ("import_errors" %in% names(entity$metadata) && nrow(entity$metadata) > 0) {
    import_errs_val <- entity$metadata$import_errors[1]
    if (is.na(import_errs_val)) import_errs_val <- ""
  }

  list(
    gff = gff,
    fasta = fasta,
    fa = fa,
    seqnames = entity$indices$seqnames,
    gff_used = gff_used_val,
    import_errs = import_errs_val
  )
}

# Helper: Simple FASTA parser (no dependencies)
.parse_fasta_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- grep("^>", lines)

  if (length(headers) == 0) {
    cli::cli_abort("No FASTA headers found in file")
  }

  sequences <- character(length(headers))
  names_vec <- character(length(headers))

  for (i in seq_along(headers)) {
    h_idx <- headers[i]
    # Extract name from header
    header_line <- lines[h_idx]
    name <- sub("^>\\s*", "", header_line)
    name <- strsplit(name, "\\s+")[[1]][1]  # Take first token
    names_vec[i] <- name

    # Extract sequence
    start <- h_idx + 1
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)

    if (start <= end) {
      seq_lines <- lines[start:end]
      seq <- paste(seq_lines, collapse = "")
      seq <- toupper(gsub("[^ACGTNacgtn]", "", seq))
      sequences[i] <- seq
    } else {
      sequences[i] <- ""
    }
  }

  names(sequences) <- names_vec
  list(sequences = sequences)
}

# Helper: Simple GFF3 parser (no dependencies)
.parse_gff3_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  # Remove comments
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(lines)]

  if (length(lines) == 0) {
    return(list(features = data.frame()))
  }

  # Parse each line
  features_list <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 9) return(NULL)

    # Parse attributes
    attrs_str <- fields[9]
    attrs <- .parse_gff3_attributes(attrs_str)

    data.frame(
      seqname = fields[1],
      source = fields[2],
      type = fields[3],
      start = as.integer(fields[4]),
      end = as.integer(fields[5]),
      score = fields[6],
      strand = fields[7],
      phase = fields[8],
      ID = if ("ID" %in% names(attrs)) attrs$ID else NA_character_,
      gene = if ("gene" %in% names(attrs)) attrs$gene else
             if ("Name" %in% names(attrs)) attrs$Name else NA_character_,
      locus_tag = if ("locus_tag" %in% names(attrs)) attrs$locus_tag else NA_character_,
      product = if ("product" %in% names(attrs)) attrs$product else NA_character_,
      protein_id = if ("protein_id" %in% names(attrs)) attrs$protein_id else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  features_df <- do.call(rbind, features_list[!sapply(features_list, is.null)])
  list(features = features_df)
}

# Helper: Parse GFF3 attributes string
.parse_gff3_attributes <- function(attrs_str) {
  if (is.na(attrs_str) || !nzchar(attrs_str)) return(list())

  pairs <- strsplit(attrs_str, ";")[[1]]
  attrs <- list()

  for (pair in pairs) {
    kv <- strsplit(pair, "=")[[1]]
    if (length(kv) == 2) {
      attrs[[kv[1]]] <- kv[2]
    }
  }

  attrs
}

# Helper: Convert GRanges to data.frame
.granges_to_dataframe <- function(gr) {
  df <- data.frame(
    seqname = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    strand = as.character(GenomicRanges::strand(gr)),
    stringsAsFactors = FALSE
  )

  # Add metadata columns
  mcols_df <- as.data.frame(S4Vectors::mcols(gr))
  if (ncol(mcols_df) > 0) {
    df <- cbind(df, mcols_df)
  }

  # Ensure required columns exist
  if (!"type" %in% names(df)) df$type <- NA_character_
  if (!"gene" %in% names(df)) df$gene <- NA_character_
  if (!"locus_tag" %in% names(df)) df$locus_tag <- NA_character_
  if (!"product" %in% names(df)) df$product <- NA_character_
  if (!"protein_id" %in% names(df)) df$protein_id <- NA_character_
  if (!"location_string" %in% names(df)) {
    # Synthesize location strings from coordinates
    df$location_string <- .synthesize_location_strings(df)
  }

  df
}

# Helper: Synthesize GenBank-style location strings from coordinates
.synthesize_location_strings <- function(df) {
  vapply(seq_len(nrow(df)), function(i) {
    loc <- paste0(df$start[i], "..", df$end[i])
    if (!is.na(df$strand[i]) && (df$strand[i] == "-" || df$strand[i] == -1)) {
      loc <- paste0("complement(", loc, ")")
    }
    loc
  }, character(1))
}
