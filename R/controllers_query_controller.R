# NOTE: extract_sequences_by_name and extract_sequences_by_coords are defined as
# S3 generics in generics_queries.R (as extract_by_name / extract_by_coords) with
# backward-compat aliases. The rich Bioconductor-based implementation lives in the
# extract_by_name.genome_entity method in that file.

# NOTE: search_features() is defined as an S3 generic in generics_queries.R.
# The genome_entity method works for genome_entity_gd as well via S3 dispatch
# inheritance (genome_entity_gd inherits from genome_entity).
# =============================================================================
# Query functions — moved from queries.R
# These are Bioconductor-dependent user-facing API functions.
# Each accepts genome_entity (converted to legacy format via entity_to_legacy_genome_obj)
# or the legacy list from init_genome() directly.
# =============================================================================

# Piping, tidyselect ops
#' @importFrom magrittr %>%

#' Extract Sequence from Genomic Coordinates
#'
#' @param genome_obj Genome object from init_genome() (expects $fa, $fasta (optional), $seqnames)
#' @param seqname Chromosome/contig name or numeric-like index (e.g., "1")
#' @param start Start coordinate (1-based, inclusive)
#' @param end End coordinate (1-based, inclusive)
#' @param strand Strand ("+", "-" or "*"); default "+"
#' @param auto_resolve If TRUE, resolve seqname to a valid FASTA header (default TRUE)
#' @param clamp If TRUE, clamp coordinates into contig bounds with a warning (default FALSE)
#' @return DNAStringSet of length 1 with the extracted sequence
#' @keywords internal
extract_sequence_by_coords <- function(genome_obj,
                                       seqname,
                                       start,
                                       end,
                                       strand = "+",
                                       auto_resolve = TRUE,
                                       clamp = FALSE) {
  # --- required namespaces ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("IRanges", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # --- basic checks ---
  stopifnot(is.list(genome_obj), "fa" %in% names(genome_obj), "seqnames" %in% names(genome_obj))
  if (!strand %in% c("+", "-", "*")) {
    cli::cli_abort("Strand must be one of '+', '-', or '*'.")
  }
  s <- as.integer(start); e <- as.integer(end)
  if (is.na(s) || is.na(e) || s < 1L || e < 1L || s > e) {
    cli::cli_abort("Start/end must be positive integers with start <= end.")
  }
  
  fa_headers <- genome_obj$seqnames
  
  # --- resolver: map seqname -> a valid FASTA header ---
  .resolve_seqname <- function(snm) {
    q <- as.character(snm)
    
    # 1) exact
    if (q %in% fa_headers) return(q)
    
    # 2) case-insensitive exact
    m <- match(tolower(q), tolower(fa_headers))
    if (!is.na(m)) return(fa_headers[m])
    
    # 3) numeric-like -> positional header (e.g., "1" -> first FASTA header)
    if (grepl("^[0-9]+$", q)) {
      i <- as.integer(q)
      if (!is.na(i) && i >= 1 && i <= length(fa_headers)) {
        cli::cli_alert_info("Interpreting seqname {.val {q}} as positional index -> {.val {fa_headers[i]}}")
        return(fa_headers[i])
      }
    }
    
    # 4) prefix match (case-insensitive)
    starts <- which(startsWith(tolower(fa_headers), tolower(q)))
    if (length(starts) == 1L) return(fa_headers[starts])
    
    # 5) fuzzy match (small distance)
    close <- try(agrep(q, fa_headers, ignore.case = TRUE, max.distance = 0.1), silent = TRUE)
    if (!inherits(close, "try-error") && length(close) == 1L) return(fa_headers[close])
    
    # Fail with a crisp diagnostic
    cli::cli_abort(c(
      "x" = "Unknown seqname.",
      "i" = paste0("Requested: ", q),
      "i" = paste0("Available FASTA headers: ", paste(fa_headers, collapse = ", "))
    ))
  }
  
  target <- if (isTRUE(auto_resolve)) .resolve_seqname(seqname) else as.character(seqname)
  
  # --- determine contig length for bound checks (prefer in-memory DNAStringSet if present) ---
  lengths <- NULL
  if (!is.null(genome_obj$fasta) && methods::is(genome_obj$fasta, "Biostrings::DNAStringSet")) {
    lengths <- Biostrings::width(genome_obj$fasta)
    names(lengths) <- names(genome_obj$fasta)
  } else if (!is.null(genome_obj$fasta) && methods::is(genome_obj$fasta, "DNAStringSet")) {
    lengths <- Biostrings::width(genome_obj$fasta)
    names(lengths) <- names(genome_obj$fasta)
  } else {
    # fallback to FA index
    fai <- try(Rsamtools::scanFaIndex(genome_obj$fa), silent = TRUE)
    if (!inherits(fai, "try-error")) {
      # scanFaIndex returns an index with fields; names(fai) are contig names; use $seqlengths
      lengths <- fai$seqlengths
      names(lengths) <- names(fai)
    }
  }
  
  # --- check and optionally clamp coordinates ---
  if (!is.null(lengths) && target %in% names(lengths)) {
    L <- as.integer(lengths[[target]])
    new_s <- s; new_e <- e
    if (!is.na(L)) {
      if (s < 1L || e > L) {
        if (!clamp) {
          cli::cli_abort(c(
            "x" = "Coordinates out of bounds.",
            "i" = sprintf("Requested %s:%d-%d (strand %s); contig length = %d", target, s, e, strand, L)
          ))
        } else {
          new_s <- max(1L, s); new_e <- min(L, e)
          cli::cli_warn(sprintf("Clamped coordinates to %s:%d-%d (contig length %d).", target, new_s, new_e, L))
        }
      }
    }
    s <- new_s; e <- new_e
  }
  
  # --- build ROI and extract ---
  roi <- GenomicRanges::GRanges(
    seqnames = target,
    ranges = IRanges::IRanges(start = s, end = e),
    strand = strand
  )
  
  # getSeq will automatically reverse-complement if strand == "-"
  seq <- Biostrings::getSeq(genome_obj$fa, roi)
  names(seq) <- sprintf("%s:%d-%d(%s)", target, s, e, strand)
  seq
}

# Genomic Context and Region Analysis Functions

#' Get Genomic Context Around Features
#'
#' @param genome_obj Genome object from init_genome()
#' @param features GRanges object or gene name pattern
#' @param flank_size Flanking region size (bp) on each side
#' @param feature_filter Optional regex pattern to filter context features
#' @return GRanges with features in the flanking regions
#' @keywords internal
get_genomic_context <- function(genome_obj,
                                features,
                                flank_size = 20000,
                                feature_filter = NULL) {

  # Check for required Bioconductor packages
  missing_pkgs <- character()
  required_pkgs <- c("GenomicRanges", "IRanges", "BiocGenerics")

  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    cli::cli_abort(c(
      "Missing required Bioconductor packages for get_genomic_context().",
      "x" = "Missing: {paste(missing_pkgs, collapse = ', ')}",
      "i" = "Install all required packages with:",
      " " = "if (!require('BiocManager')) install.packages('BiocManager')",
      " " = "BiocManager::install(c('GenomicRanges', 'Biostrings', 'IRanges', 'S4Vectors', 'BiocGenerics'))"
    ))
  }

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # Check that genome_obj has required components
  if (is.null(genome_obj$gff)) {
    cli::cli_abort(c(
      "genome_obj$gff is NULL.",
      "i" = "Did you create the genome object with init_genome()? This function requires GRanges data.",
      "i" = "Make sure Bioconductor packages were available when you created the genome object."
    ))
  }

  # If features is a string, find matching features using the same logic as extract_sequences_by_name
  if (is.character(features)) {
    gene_pattern <- features

    # Helper to match pattern in metadata fields
    .matches <- function(x, pat) {
      if (is.null(x)) return(rep(FALSE, length(genome_obj$gff)))
      y <- tryCatch(as.character(x), error = function(e) rep(NA_character_, length(x)))
      !is.na(y) & grepl(pat, y, perl = TRUE)
    }

    # Search Name first, then fallback to other fields
    name_idx <- .matches(genome_obj$gff$Name, gene_pattern)
    if (any(name_idx, na.rm = TRUE)) {
      features <- genome_obj$gff[name_idx]
    } else {
      # Fallback to gene, locus_tag, or product
      idx <- .matches(genome_obj$gff$gene, gene_pattern) |
        .matches(genome_obj$gff$locus_tag, gene_pattern) |
        .matches(genome_obj$gff$product, gene_pattern)
      idx[is.na(idx)] <- FALSE
      features <- genome_obj$gff[idx]
    }
  }

  if (length(features) == 0) {
    cli::cli_abort("No features provided or found")
  }

  # Expand features to include flanks
  ctx_region <- GenomicRanges::resize(features,
                       width = BiocGenerics::width(features) + 2 * flank_size,
                       fix = "center")

  # Find overlapping features
  ctx_feats <- IRanges::subsetByOverlaps(genome_obj$gff, ctx_region)
  
  # Apply filter if provided
  if (!is.null(feature_filter)) {
    ctx_feats <- ctx_feats[grepl(feature_filter, ctx_feats$Name, ignore.case = TRUE)]
  }
  
  cli::cli_inform("Found {length(ctx_feats)} feature{?s} in {flank_size} bp flanking region{?s}")
  
  return(ctx_feats)
}



#' Analyze Region of Interest
#'
#' @param genome_obj Genome object from init_genome() (expects $gff and $seqnames)
#' @param seqname Chromosome/contig name (can be numeric-like, e.g., "1")
#' @param start Start coordinate (1-based)
#' @param end End coordinate (1-based)
#' @param flank Flanking region size (default = 0)
#' @param feature_type Feature type to extract (default = "CDS")
#' @param tidy Return tibble (TRUE) or GRanges (FALSE)
#' @param auto_resolve If TRUE, resolve seqname against GFF seqlevels (default TRUE)
#' @param drop_extdb If TRUE, prefer human labels over extdb placeholders (default TRUE)
#' @return GRanges or tibble with features in ROI
#' @keywords internal
analyze_roi <- function(genome_obj,
                        seqname,
                        start,
                        end,
                        flank = 0,
                        feature_type = "CDS",
                        tidy = TRUE,
                        auto_resolve = TRUE,
                        drop_extdb = TRUE) {
  # --- namespaces ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("IRanges", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("BiocGenerics", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  gff <- genome_obj$gff
  
  # --------- resolver: map user seqname -> a level present in the GFF ---------
  .resolve_seqname_gff <- function(q, levels) {
    q <- as.character(q)
    # exact
    if (q %in% levels) return(q)
    # case-insensitive exact
    m <- match(tolower(q), tolower(levels))
    if (!is.na(m)) return(levels[m])
    # numeric-like -> positional level
    if (grepl("^[0-9]+$", q)) {
      i <- as.integer(q)
      if (!is.na(i) && i >= 1L && i <= length(levels)) {
        cli::cli_alert_info("Interpreting seqname {.val {q}} as positional index in GFF -> {.val {levels[i]}}")
        return(levels[i])
      }
    }
    # prefix
    starts <- which(startsWith(tolower(levels), tolower(q)))
    if (length(starts) == 1L) return(levels[starts])
    # fuzzy (light)
    close <- try(agrep(q, levels, ignore.case = TRUE, max.distance = 0.1), silent = TRUE)
    if (!inherits(close, "try-error") && length(close) == 1L) return(levels[close])
    cli::cli_abort(c(
      "x" = "Unknown seqname for GFF.",
      "i" = paste0("Requested: ", q),
      "i" = paste0("GFF seqlevels: ", paste(GenomeInfoDb::seqlevels(gff), collapse = ", "))
    ))
  }
  
  # --------- tidy normalization for PGAP fields ---------
  .tidy_features_pgappref <- function(gr, drop_extdb = TRUE) {
    n <- length(gr)
    # pull common columns safely
    mc <- S4Vectors::mcols(gr)
    getc <- function(x) if (x %in% names(mc)) as.character(mc[[x]]) else rep(NA_character_, n)
    
    gene       <- getc("gene")
    locus_tag  <- getc("locus_tag")
    name_raw   <- getc("Name")
    product    <- getc("product")
    type       <- getc("type")
    
    # Name priority: gene -> locus_tag -> Name (unless extdb) -> product -> ID
    name_out <- gene
    fill <- is.na(name_out) | name_out == ""
    name_out[fill] <- locus_tag[fill]
    fill <- is.na(name_out) | name_out == ""
    # use Name if not extdb:*, otherwise delay
    ok_name <- ifelse(!is.na(name_raw) & !grepl("^extdb:", name_raw), name_raw, NA_character_)
    name_out[fill] <- ok_name[fill]
    fill <- is.na(name_out) | name_out == ""
    name_out[fill] <- product[fill]
    fill <- is.na(name_out) | name_out == ""
    id <- getc("ID")
    name_out[fill] <- id[fill]
    
    # If still NA, synthesize
    name_out[is.na(name_out) | name_out == ""] <- paste0("feature_", which(is.na(name_out) | name_out == ""))
    
    # If drop_extdb, and we still have extdb in Name while locus_tag or gene exist, replace
    if (isTRUE(drop_extdb)) {
      is_ext <- grepl("^extdb:", name_out)
      replace_with <- ifelse(!is.na(gene) & gene != "", gene,
                             ifelse(!is.na(locus_tag) & locus_tag != "", locus_tag, name_out))
      name_out[is_ext] <- replace_with[is_ext]
    }
    
    tibble::tibble(
      seqnames = as.character(GenomicRanges::seqnames(gr)),
      start    = BiocGenerics::start(gr),
      end      = BiocGenerics::end(gr),
      width    = BiocGenerics::width(gr),
      strand   = as.character(GenomicRanges::strand(gr)),
      type     = as.character(type),
      Name     = name_out,
      Alias    = locus_tag,            # maps your desired Alias to locus_tag
      Note     = product               # maps your desired Note to product
    )
  }
  
  # ------------------- build ROI -------------------
  s <- as.integer(start); e <- as.integer(end); f <- as.integer(flank)
  if (is.na(s) || is.na(e) || s < 1L || e < 1L || s > e) {
    cli::cli_abort("Start/end must be positive integers with start <= end.")
  }
  
  # Resolve seqname against the GFF vocabulary
  target <- if (isTRUE(auto_resolve)) {
    .resolve_seqname_gff(seqname, GenomeInfoDb::seqlevels(gff))
  } else {
    as.character(seqname)
  }
  
  roi <- GenomicRanges::GRanges(
    seqnames = target,
    ranges   = IRanges::IRanges(start = max(1L, s - f), end = e + f)
  )
  
  cli::cli_inform(sprintf("Analyzing region: %s:%d-%d (\u00B1%d bp)",
                          as.character(seqname), s, e, f))
  
  # ------------------- overlap -------------------
  roi_feats <- IRanges::subsetByOverlaps(gff, roi)

  # Filter type if requested; drop region scaffolding
  if (!is.null(feature_type)) {
    mcols <- S4Vectors::mcols(roi_feats)
    if ("type" %in% names(mcols)) {
      roi_feats <- roi_feats[mcols$type == feature_type]
    }
  }
  
  cli::cli_inform(sprintf("Found %d %s feature(s) in ROI",
                          length(roi_feats), 
                          ifelse(is.null(feature_type), 
                                 "matching", 
                                 feature_type)))
  
  if (!isTRUE(tidy)) {
    return(roi_feats)
  }
  
  .tidy_features_pgappref(roi_feats, drop_extdb = drop_extdb)
}


#' Search Genome-wide for Features (Legacy Internal)
#'
#' @param genome_obj Genome object from init_genome()
#' @param pattern Regex pattern to search in Name or Note fields
#' @param field Field to search ("Name", "Note", "Alias", or "any")
#' @param feature_type Feature type filter (default = "CDS")
#' @param tidy Return tidy data frame (default = TRUE)
#' @return GRanges or tibble with matching features
#' @keywords internal
#'
search_features_legacy_internal <- function(genome_obj,
                            pattern,
                            field = "Note",
                            feature_type = "CDS",
                            tidy = TRUE) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  gff <- genome_obj$gff
  
  # Filter by feature type
  if (!is.null(feature_type)) {
    gff <- gff[gff$type == feature_type]
  }
  
  # Search in specified field(s)
  if (field == "any") {
    matches <- grepl(pattern, gff$Name) | 
      grepl(pattern, flatten_CharacterList(gff$Note)) |
      grepl(pattern, flatten_CharacterList(gff$Alias))
  } else if (field == "Note") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Note))
  } else if (field == "Alias") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Alias))
  } else if (field == "Name") {
    matches <- grepl(pattern, flatten_CharacterList(gff$Name))
  } else {
    matches <- grepl(pattern, gff[[field]])
  }
  
  result <- gff[matches]
  
  cli::cli_inform("Found {length(result)} feature{?s} matching '{pattern}' in {field} field{?s}")
  
  if (!tidy) {
    return(result)
  }
  
  return(tidy_features(result))
}

#' Extract FASTA Sequence from Region of Interest
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param flank Flanking region size in bp (default = 0)
#' @param strand Strand to extract ("+" or "-" or "*" for both); default = "+"
#' @param reverse_complement Apply reverse complement for minus strand (default = TRUE)
#' @param output_file Optional file path to write FASTA (default = NULL)
#' @param seq_name Custom sequence name for FASTA header (default = auto-generated)
#' @return DNAStringSet with extracted sequence(s)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Basic extraction
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000)
#' 
#' # With flanking sequence
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000, flank = 5000)
#' 
#' # Minus strand with reverse complement
#' roi_seq <- get_roi_fasta(genome, "1", 1828501, 1978000, strand = "-")
#' 
#' # Save to file
#' get_roi_fasta(genome, "1", 1828501, 1978000, 
#'               flank = 10000, 
#'               output_file = "roi_region.fasta")
#' }
get_roi_fasta <- function(genome_obj,
                          seqname,
                          start,
                          end,
                          flank = 0,
                          strand = "+",
                          reverse_complement = TRUE,
                          output_file = NULL,
                          seq_name = NULL) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # Input validation
  if (start > end) {
    cli::cli_abort("Start coordinate must be <= end coordinate")
  }

  if (flank < 0) {
    cli::cli_abort("Flank size must be >= 0")
  }

  if (!strand %in% c("+", "-", "*")) {
    cli::cli_abort("Strand must be '+', '-', or '*'")
  }
  
  # Calculate actual coordinates with flanking
  actual_start <- max(1, start - flank)
  actual_end <- end + flank

  # Get sequence length to validate end coordinate
  # Use fa (FaFile) if available, otherwise use fasta (DNAStringSet)
  if (!is.null(genome_obj$fa)) {
    si <- GenomeInfoDb::seqinfo(genome_obj$fa)
    seq_lengths <- GenomeInfoDb::seqlengths(si)
  } else if (!is.null(genome_obj$fasta)) {
    seq_lengths <- BiocGenerics::width(genome_obj$fasta)
    # Ensure names are preserved
    if (is.null(names(seq_lengths))) {
      names(seq_lengths) <- names(genome_obj$fasta)
    }
  } else {
    cli::cli_abort("No sequence data available in genome_obj (neither $fa nor $fasta)")
  }

  if (!seqname %in% names(seq_lengths)) {
    cli::cli_abort("Seqname '{seqname}' not found in FASTA. Available: {paste(names(seq_lengths), collapse = ', ')}")
  }

  # Adjust end if it exceeds chromosome length
  max_length <- seq_lengths[seqname]
  if (actual_end > max_length) {
    cli::cli_warn("End coordinate ({actual_end}) exceeds sequence length ({max_length}). Adjusting to {max_length}")
    actual_end <- max_length
  }
  
  # Generate sequence name if not provided
  if (is.null(seq_name)) {
    if (flank > 0) {
      seq_name <- sprintf("%s:%d-%d_flank%d(%s)", 
                          seqname, start, end, flank, strand)
    } else {
      seq_name <- sprintf("%s:%d-%d(%s)", 
                          seqname, start, end, strand)
    }
  }
  
  # Create GRanges object for extraction
  roi <- GenomicRanges::GRanges(
    seqnames = seqname,
    ranges = IRanges::IRanges(start = actual_start, end = actual_end),
    strand = strand
  )

  # Extract sequence
  cli::cli_inform("Extracting sequence from {seqname}:{start}-{end} (+/-{flank} bp, strand: {strand})")

  # Use fa (FaFile) if available, otherwise use fasta (DNAStringSet)
  if (!is.null(genome_obj$fa)) {
    seq <- Biostrings::getSeq(genome_obj$fa, roi)
  } else {
    # For DNAStringSet, extract the specific sequence then get subsequence
    chr_seq <- genome_obj$fasta[[seqname]]
    if (is.null(chr_seq)) {
      cli::cli_abort("Sequence '{seqname}' not found in DNAStringSet")
    }
    # Extract the region (accounting for strand)
    seq <- Biostrings::subseq(chr_seq, start = actual_start, end = actual_end)
    # Convert to DNAStringSet for consistency
    seq <- Biostrings::DNAStringSet(seq)
  }

  # Apply reverse complement if on minus strand and requested
  if (strand == "-" && reverse_complement) {
    cli::cli_inform("Applying reverse complement for minus strand")
    seq <- Biostrings::reverseComplement(seq)
  }

  # Set sequence name
  names(seq) <- seq_name

  # Report sequence info
  cli::cli_inform("Extracted {BiocGenerics::width(seq)} bp (requested region: {end - start + 1} bp, with flank: {actual_end - actual_start + 1} bp)")
  
  # Write to file if requested
  if (!is.null(output_file)) {
    cli::cli_inform("Writing sequence to: {output_file}")
    Biostrings::writeXStringSet(seq, filepath = output_file, format = "fasta")
  }
  
  return(seq)
}


#' Extract Multiple ROI Sequences as Batch
#'
#' @param genome_obj Genome object from init_genome()
#' @param roi_table Data frame with columns: seqname, start, end, and optionally strand, name
#' @param flank Flanking region size in bp (default = 0)
#' @param output_file Optional file path to write all sequences (default = NULL)
#' @return DNAStringSet with all extracted sequences
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Define multiple regions
#' regions <- tibble(
#'   seqname = c("1", "1", "1"),
#'   start = c(100000, 200000, 300000),
#'   end = c(150000, 250000, 350000),
#'   name = c("region1", "region2", "region3")
#' )
#' 
#' # Extract all at once
#' seqs <- get_roi_fasta_batch(genome, regions, flank = 1000)
#' 
#' # Save to file
#' get_roi_fasta_batch(genome, regions, 
#'                     flank = 1000, 
#'                     output_file = "all_regions.fasta")
#' }
get_roi_fasta_batch <- function(genome_obj,
                                roi_table,
                                flank = 0,
                                output_file = NULL) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # Validate input table
  required_cols <- c("seqname", "start", "end")
  if (!all(required_cols %in% names(roi_table))) {
    cli::cli_abort("roi_table must contain columns: {paste(required_cols, collapse = ', ')}")
  }
  
  cli::cli_inform("Extracting {nrow(roi_table)} ROI sequence{?s}...")
  
  # Extract sequences for each ROI
  seq_list <- lapply(seq_len(nrow(roi_table)), function(i) {
    row <- roi_table[i, ]
    
    # Get strand if available
    strand <- if ("strand" %in% names(row)) row$strand else "+"
    
    # Get custom name if available
    seq_name <- if ("name" %in% names(row)) {
      if (flank > 0) {
        sprintf("%s_%s:%d-%d_flank%d", row$name, row$seqname, row$start, row$end, flank)
      } else {
        sprintf("%s_%s:%d-%d", row$name, row$seqname, row$start, row$end)
      }
    } else {
      NULL  # Will auto-generate in get_roi_fasta
    }
    
    # Extract sequence
    seq <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = row$seqname,
      start = row$start,
      end = row$end,
      flank = flank,
      strand = strand,
      seq_name = seq_name
    )
    
    return(seq)
  })
  
  # Combine all sequences
  all_seqs <- do.call(c, seq_list)
  
  cli::cli_inform("Successfully extracted {length(all_seqs)} sequence{?s}")
  
  # Write to file if requested
  if (!is.null(output_file)) {
    cli::cli_inform("Writing all sequences to: {output_file}")
    writeXStringSet(all_seqs, filepath = output_file, format = "fasta")
  }
  
  return(all_seqs)
}


#' Extract FASTA for Features with Flanking Regions
#'
#' @param genome_obj Genome object from init_genome()
#' @param features GRanges object or gene name pattern
#' @param flank_upstream Upstream flanking size (default = 0)
#' @param flank_downstream Downstream flanking size (default = 0)
#' @param include_feature Include the feature itself (default = TRUE)
#' @param output_file Optional file path to write sequences (default = NULL)
#' @return DNAStringSet with extracted sequences
#' @export
#'
#' @examples
#' \dontrun{
#' # Get acrB genes with 5kb upstream and 2kb downstream
#' seqs <- get_feature_fasta(genome, "acrB", 
#'                          flank_upstream = 5000, 
#'                          flank_downstream = 2000)
#' 
#' # Get only upstream promoter region (exclude gene)
#' promoters <- get_feature_fasta(genome, "acrB",
#'                               flank_upstream = 1000,
#'                               flank_downstream = 0,
#'                               include_feature = FALSE)
#' }
get_feature_fasta <- function(genome_obj,
                              features,
                              flank_upstream = 0,
                              flank_downstream = 0,
                              include_feature = TRUE,
                              output_file = NULL) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # If features is a string, find matching features
  if (is.character(features)) {
    pattern <- features
    features <- genome_obj$gff[grepl(pattern, genome_obj$gff$Name)]
    
    if (length(features) == 0) {
      cli::cli_abort("No features found matching pattern: {pattern}")
    }
    
    cli::cli_inform("Found {length(features)} feature{?s} matching '{pattern}'")
  }
  
  # Fix seqnames if needed
  # Get actual sequence names from the data source
  if (!is.null(genome_obj$fa)) {
    actual_seqnames <- names(GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(genome_obj$fa)))
  } else if (!is.null(genome_obj$fasta)) {
    actual_seqnames <- names(genome_obj$fasta)
  } else {
    cli::cli_abort("No sequence data available")
  }

  # If there's only one sequence and features have different seqnames, fix them
  if (length(actual_seqnames) == 1) {
    cli::cli_inform("Fixing seqnames: '{as.character(unique(GenomicRanges::seqnames(features)))[1]}' -> '{actual_seqnames[1]}'")
    GenomeInfoDb::seqlevels(features) <- actual_seqnames
    GenomicRanges::seqnames(features) <- S4Vectors::Rle(rep(actual_seqnames[1], length(features)))
  } else if (length(actual_seqnames) == 0) {
    cli::cli_abort("FASTA/DNAStringSet has no sequence names set. This suggests a data structure issue.")
  }
  
  # Prepare extraction regions based on strand
  seq_list <- lapply(seq_along(features), function(i) {
    feat <- features[i]
    
    # Calculate coordinates accounting for strand
    if (as.character(GenomicRanges::strand(feat)) == "-") {
      # For minus strand, upstream is actually downstream in coordinates
      new_start <- BiocGenerics::start(feat) - flank_downstream
      new_end <- BiocGenerics::end(feat) + flank_upstream
    } else {
      # For plus strand
      new_start <- BiocGenerics::start(feat) - flank_upstream
      new_end <- BiocGenerics::end(feat) + flank_downstream
    }
    
    # If not including feature, adjust coordinates
    if (!include_feature) {
      if (as.character(GenomicRanges::strand(feat)) == "-") {
        new_start <- BiocGenerics::end(feat) + 1
        new_end <- BiocGenerics::end(feat) + flank_upstream
      } else {
        new_start <- BiocGenerics::start(feat) - flank_upstream
        new_end <- BiocGenerics::start(feat) - 1
      }
    }
    
    # Generate descriptive name
    feat_name <- if (!is.null(feat$Name)) {
      sprintf("%s_%s:%d-%d_up%d_down%d(%s)",
              feat$Name,
              as.character(GenomicRanges::seqnames(feat)),
              BiocGenerics::start(feat),
              BiocGenerics::end(feat),
              flank_upstream,
              flank_downstream,
              as.character(GenomicRanges::strand(feat)))
    } else {
      sprintf("%s:%d-%d_up%d_down%d(%s)",
              as.character(GenomicRanges::seqnames(feat)),
              BiocGenerics::start(feat),
              BiocGenerics::end(feat),
              flank_upstream,
              flank_downstream,
              as.character(GenomicRanges::strand(feat)))
    }
    
    # Extract sequence
    seq <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = as.character(GenomicRanges::seqnames(feat)),
      start = max(1, new_start),
      end = new_end,
      strand = as.character(GenomicRanges::strand(feat)),
      seq_name = feat_name
    )
    
    return(seq)
  })
  
  # Combine all sequences
  all_seqs <- do.call(c, seq_list)
  
  # Write to file if requested
  if (!is.null(output_file)) {
    cli::cli_inform("Writing all feature sequences to: {output_file}")
    writeXStringSet(all_seqs, filepath = output_file, format = "fasta")
  }
  
  return(all_seqs)
}


#' Quick wrapper: Extract ROI with Annotations
#'
#' @param genome_obj Genome object from init_genome()
#' @param seqname Chromosome/contig name
#' @param start Start coordinate
#' @param end End coordinate
#' @param flank Flanking region size (default = 0)
#' @param get_sequence Extract FASTA sequence (default = TRUE)
#' @param get_features Get overlapping features (default = TRUE)
#' @param output_fasta Optional FASTA output file
#' @return List with sequence (if requested) and features (if requested)
#' @keywords internal
extract_roi_complete <- function(genome_obj,
                                 seqname,
                                 start,
                                 end,
                                 flank = 0,
                                 get_sequence = TRUE,
                                 get_features = TRUE,
                                 output_fasta = NULL) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  result <- list()
  
  # Extract sequence if requested
  if (get_sequence) {
    cli::cli_alert_info("Extracting sequence")
    result$sequence <- get_roi_fasta(
      genome_obj = genome_obj,
      seqname = seqname,
      start = start,
      end = end,
      flank = flank,
      output_file = output_fasta
    )
  }
  
  # Get features if requested
  if (get_features) {
    cli::cli_alert_info("Finding overlapping features")
    result$features <- analyze_roi(
      genome_obj = genome_obj,
      seqname = seqname,
      start = start,
      end = end,
      flank = flank
    )
  }
  
  cli::cli_alert_success("Extraction complete")
  return(result)
}

#' Submit Protein BLAST Search
#'
#' Run Local BLASTP Against Protein Sequence
#'
#' @description
#' **IMPORTANT**: This function ONLY supports LOCAL BLAST databases.
#' Remote BLAST is not supported. You must have:
#'
#' 1. BLAST+ installed (blastp binary accessible)
#' 2. A local protein database (e.g., SwissProt, nr, RefSeq)
#' 3. BLASTDB environment variable set to the database directory
#'
#' @section Required External Dependencies:
#' - **BLAST+**: Command-line BLAST tools must be installed
#'   - Install: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#'   - Or via conda: `conda install -c bioconda blast`
#' - **Local Database**: Download from NCBI or build your own
#'   - SwissProt: ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
#'   - Extract to a directory and set BLASTDB environment variable
#'
#' @section Environment Variables:
#' Set `BLASTDB` to the directory containing your BLAST databases:
#' ```
#' Sys.setenv(BLASTDB = "/path/to/blast/databases")
#' ```
#' Or in your .Renviron file:
#' ```
#' BLASTDB=/path/to/blast/databases
#' ```
#'
#' @param sequence Protein sequence as character string or AAString
#' @param database Name of BLAST database (e.g., "swissprot", "nr")
#' @param dbdir Optional: explicit directory containing databases (overrides BLASTDB)
#' @param evalue E-value threshold (default 1e-10)
#' @param threads Number of threads to use (default: all available cores)
#' @param validate_db Logical; check if database exists before running (default TRUE)
#' @param max_hits Maximum number of hits to return (default 20)
#'
#' @return tibble with BLAST hits containing columns:
#'   - qseqid: query sequence ID
#'   - sacc: subject accession
#'   - stitle: subject title/description
#'   - pident: percent identity
#'   - length: alignment length
#'   - qcovs: query coverage per subject
#'   - bitscore: bit score
#'   - evalue: expect value
#'   - staxids: subject taxonomic IDs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Set BLASTDB environment variable
#' Sys.setenv(BLASTDB = "/path/to/blast/databases")
#'
#' # Run BLASTP
#' results <- blast_protein(
#'   sequence = "MPNFFIDRPIFAWVIAIIIMLAGGLAILK...",
#'   database = "swissprot",
#'   evalue = 1e-5,
#'   max_hits = 10
#' )
#'
#' # View results
#' print(results)
#' }
blast_protein <- function(sequence,
                          database = "swissprot",
                          dbdir = "",
                          evalue = 1e-10,
                          threads = parallel::detectCores(),
                          validate_db = TRUE,
                          max_hits = 20) {

  if (!requireNamespace("cli", quietly = TRUE)) cli::cli_abort("Package 'cli' is required")

  # Convert sequence to character if needed
  if (is(sequence, "AAString") || is(sequence, "DNAString")) {
    sequence <- as.character(sequence)
  }

  # Check for local database
  if (nzchar(dbdir)) {
    db_path <- file.path(dbdir, database)
  } else {
    db_path <- database
  }

  # Verify database exists if validation is requested
  if (validate_db) {
    # Try to find database files
    if (nzchar(dbdir)) {
      check_path <- file.path(dbdir, paste0(database, ".psq"))
    } else {
      # Check in current directory or BLASTDB path
      check_path <- paste0(database, ".psq")
      if (!file.exists(check_path) && nzchar(Sys.getenv("BLASTDB"))) {
        check_path <- file.path(Sys.getenv("BLASTDB"), paste0(database, ".psq"))
      }
    }

    if (!file.exists(check_path)) {
      cli::cli_abort(c(
        "Local BLAST database not found: {.val {database}}",
        i = "Searched for: {.path {check_path}}",
        i = "This function only supports LOCAL BLAST databases.",
        i = "Set {.env BLASTDB} environment variable or provide {.arg dbdir} parameter.",
        i = "Available databases should have files like {.file {database}.psq}"
      ))
    }
  }

  # Create temporary query file
  tmp_query <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp_query), add = TRUE)

  # Write sequence to file
  writeLines(c(">query", sequence), tmp_query)

  # Run BLASTP using BLAST gateway
  cli::cli_alert_info("Running local BLASTP against {.val {database}}")

  # If BLASTDB is set and no explicit dbdir provided, let BLAST use the env var
  use_dbdir <- if (nzchar(dbdir)) {
    dbdir
  } else if (nzchar(Sys.getenv("BLASTDB"))) {
    ""  # Let BLAST use BLASTDB env var to avoid path issues
  } else {
    ""
  }

  # Create gateway and run BLASTP
  gateway <- create_blast_gateway()
  results <- gateway$blastp_capture(
    query_faa = tmp_query,
    db = database,
    dbdir = use_dbdir,
    evalue = evalue,
    threads = threads,
    validate_db = FALSE,  # Already validated above
    more_args = c("-max_target_seqs", as.character(max_hits))
  )

  results
}


#' Parse BLAST XML Results
#'
#' @param xml_text XML text from BLAST results
#' @param top_n Number of top hits to return (default = 10)
#' @return Tibble with top BLAST hits
#' @keywords internal
parse_blast_xml <- function(xml_text, top_n = 10) {
  
  doc <- read_xml(xml_text)
  hits <- xml_find_all(doc, "//Hit")
  
  if (length(hits) == 0) {
    cli::cli_alert_warning("No BLAST hits found")
    return(tibble())
  }
  
  result <- tibble(
    accession = xml_text(xml_find_all(hits, ".//Hit_accession")),
    title = xml_text(xml_find_all(hits, ".//Hit_def")),
    length = as.integer(xml_text(xml_find_all(hits, ".//Hit_len")))
  )
  
  if (!is.null(top_n)) {
    result <- head(result, top_n)
  }
  
  return(result)
}



#' Flatten CharacterList to Character Vector
#'
#' @param x CharacterList object
#' @return Character vector with comma-separated values
flatten_CharacterList <- function(x) {
  sapply(x, function(v) {
    if (length(v)) paste(v, collapse = ",") else NA_character_
  })
}


#' Convert GRanges to Tidy Data Frame
#'
#' @param gr GRanges object
#' @param sort_by Sort by coordinate (default = TRUE)
#' @return Tibble with key annotation fields
#' @keywords internal
tidy_features <- function(gr, sort_by = "start") {
  
  if (length(gr) == 0) {
    return(tibble())
  }
  
  # Sort if requested
  if (!is.null(sort_by) && sort_by == "start") {
    gr <- gr[order(start(gr))]
  }
  
  # Convert to data frame and tidy up
  df <- as.data.frame(gr)
  
  # Flatten CharacterList columns if they exist
  if ("Alias" %in% names(df)) {
    df$Alias <- flatten_CharacterList(gr$Alias)
  }
  if ("Note" %in% names(df)) {
    df$Note <- flatten_CharacterList(gr$Note)
  }
  
  # Select key columns (adjust based on your GFF structure)
  key_cols <- c("seqnames", "start", "end", "width", "strand", "type", 
                "Name", "Alias", "Note")
  available_cols <- intersect(key_cols, names(df))
  
  result <- df %>%
    dplyr::select(dplyr::all_of(available_cols)) %>%
    tibble::as_tibble()
  
  return(result)
}


#' Compare Protein Sequences
#'
#' @param seq1 First sequence (AAString or character)
#' @param seq2 Second sequence (AAString or character)
#' @param type Alignment type (default = "global")
#' @return List with percent identity and alignment object
#' @keywords internal
compare_sequences <- function(seq1, seq2, type = "global") {
  
  if (!requireNamespace("pwalign", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package 'pwalign' is required",
      "i" = "Install with: BiocManager::install('pwalign')"
    ))
  }
  
  aln <- pwalign::pairwiseAlignment(seq1, seq2, type = type)
  pid <- pwalign::pid(aln)
  
  cli::cli_inform("Percent identity: {round(pid, 2)}%")
  
  list(
    percent_identity = pid,
    alignment = aln
  )
}

#' Complete Gene Analysis Workflow
#'
#' @param genome_obj Genome object from init_genome()
#' @param gene_pattern Gene name pattern
#' @param flank_size Flanking region size for context
#' @param blast Perform BLAST search (default = FALSE)
#' @return List with sequences, context, and optional BLAST results
#' @keywords internal
analyze_gene <- function(genome_obj,
                         gene_pattern,
                         flank_size = 20000,
                         blast = FALSE) {

  # Check for required Bioconductor packages
  missing_pkgs <- character()
  required_pkgs <- c("GenomicRanges", "Biostrings", "IRanges", "S4Vectors", "BiocGenerics")

  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    cli::cli_abort(c(
      "Missing required Bioconductor packages for analyze_gene().",
      "x" = "Missing: {paste(missing_pkgs, collapse = ', ')}",
      "i" = "Install all required packages with:",
      " " = "if (!require('BiocManager')) install.packages('BiocManager')",
      " " = "BiocManager::install(c('GenomicRanges', 'Biostrings', 'IRanges', 'S4Vectors', 'BiocGenerics'))"
    ))
  }

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    cli::cli_inform("Converting genome_entity to legacy format for analyze_gene()...")
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  cli::cli_h1("Gene Analysis Pipeline")

  # 1. Extract sequences
  dna_seqs <- extract_sequences_by_name(genome_obj, gene_pattern, translate = FALSE)
  aa_seqs <- extract_sequences_by_name(genome_obj, gene_pattern, translate = TRUE)
  
  cli::cli_alert_info("Sequence lengths:")
  print(BiocGenerics::width(aa_seqs))
  
  # 2. Get genomic context
  cli::cli_h2("Genomic Context")
  context <- get_genomic_context(genome_obj, gene_pattern, flank_size)
  
  # 3. Compare sequences if multiple found
  comparisons <- NULL
  if (length(aa_seqs) > 1) {
    cli::cli_h2("Sequence Comparisons")
    comparisons <- list()
    for (i in 1:(length(aa_seqs) - 1)) {
      for (j in (i + 1):length(aa_seqs)) {
        comp_name <- paste(names(aa_seqs)[i], "vs", names(aa_seqs)[j])
        comparisons[[comp_name]] <- compare_sequences(aa_seqs[i], aa_seqs[j])
      }
    }
  }
  
  # 4. Optional BLAST
  blast_results <- NULL
  if (blast && length(aa_seqs) > 0) {
    cli::cli_h2("BLAST Analysis")
    blast_results <- lapply(seq_along(aa_seqs), function(i) {
      cli::cli_inform("BLASTing {names(aa_seqs)[i]}...")
      blast_protein(aa_seqs[[i]])
    })
    names(blast_results) <- names(aa_seqs)
  }
  
  # Return comprehensive results
  list(
    dna_sequences = dna_seqs,
    protein_sequences = aa_seqs,
    genomic_context = tidy_features(context),
    sequence_comparisons = comparisons,
    blast_results = blast_results
  )
}


translate_cds <- function(dna, strand = "+", drop_terminal_stop = TRUE,
                          bact_gc = Biostrings::getGeneticCode("11")) {
  
  # specified bacterial/archaeal genetic code (11) 
  
  # Normalize orientation
  dna2 <- if (as.character(strand) == "-") Biostrings::reverseComplement(dna) else dna
  
  # Trim to nearest full codon to avoid fuzzy tail at ROI edges
  rem <- width(dna2) %% 3L
  if (rem != 0L) dna2 <- dna2[1:(width(dna2) - rem)]
  
  aa <- Biostrings::translate(dna2, genetic.code = bact_gc, if.fuzzy.codon = "X")
  aa <- as.character(aa)
  
  # Strip a terminal stop, keep internal stops as QC
  if (drop_terminal_stop) aa <- sub("\\*$", "", aa)
  aa
}

roi_cds_to_faa <- function(genome_obj, roi_tbl,
                           bact_gc = Biostrings::getGeneticCode("11"),
                           outfile = tempfile(fileext = ".faa")) {

  # Auto-convert genome_entity objects to legacy format
  if (inherits(genome_obj, "genome_entity")) {
    genome_obj <- entity_to_legacy_genome_obj(genome_obj)
  }

  # process data
  cds <- roi_tbl %>%
    filter(type == "CDS") %>%
    mutate(
      dna = pmap(list(seqnames, start, end, strand),
                 ~ extract_sequence_by_coords(genome_obj, ..1, ..2, ..3, ..4)),
      aa_raw = map2_chr(dna, strand, ~ as.character(Biostrings::translate(
        if (as.character(..2) == "-") Biostrings::reverseComplement(..1) else ..1,
        genetic.code = bact_gc, if.fuzzy.codon = "X"
      ))),
      aa = sub("\\*$", "", aa_raw),  # drop terminal stop only
      aa_len = nchar(aa),
      has_stop = str_detect(aa_raw, "\\*") & !str_detect(aa_raw, "\\*$")
    )
  
  # Header preference: Name, else Alias, else coord signature
  hdr <- cds$Name
  hdr[is.na(hdr) | hdr == ""] <- cds$Alias[is.na(hdr) | hdr == ""]
  fallback <- paste0(cds$seqnames, ":", cds$start, "-", cds$end, "(", cds$strand, ")")
  hdr[is.na(hdr) | hdr == ""] <- fallback[is.na(hdr) | hdr == ""]
  
  # Ensure uniqueness if the GFF repeats IDs
  if (any(duplicated(hdr))) {
    dup_idx <- which(duplicated(hdr) | duplicated(hdr, fromLast = TRUE))
    hdr[dup_idx] <- paste0(hdr[dup_idx], "_", seq_along(dup_idx))
  }
  
  lines <- unlist(map2(hdr, cds$aa, ~ c(paste0(">", .x), .y)))
  write_lines(lines, outfile)
  
  cds %>% mutate(qseqid = hdr, faa_path = outfile)
}




