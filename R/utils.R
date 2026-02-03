#' Legacy Utility Functions
#'
#' @description
#' Legacy utility functions from the original codebase. These are maintained for
#' backward compatibility but should not be used in new code.
#'
#' Note: The library() calls have been removed. Package dependencies are managed
#' through DESCRIPTION file (Imports/Suggests). Functions that require specific
#' packages should use requireNamespace() to check availability.
#'
#' @name legacy-utils
#' @keywords internal
NULL

# Setup and Data Loading Functions ====

# A tiny infix helper to mirror rlang's %||%, avoiding hard dependency.
`%||%` <- function(x, y) if (is.null(x)) y else x
# Helpers

# Scrubber for common prefixes. Extend as your corpus dictates.
.scrub_prefixes <- function(x) {
  x <- as.character(x)
  x <- sub("^lcl\\|", "", x, perl = TRUE)
  x <- sub("^chr",    "", x, ignore.case = TRUE)
  x
}

# Return FASTA names and lengths from either DNAStringSet or FaFile
.get_fasta_index <- function(fa_obj) {
  if (inherits(fa_obj, "DNAStringSet")) {
    nm <- names(fa_obj) %||% character(length(fa_obj))
    len <- Biostrings::width(fa_obj)
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }
  if (inherits(fa_obj, "FaFile")) {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop("Rsamtools is required to inspect FaFile indices.")
    }
    idx <- Rsamtools::scanFaIndex(fa_obj)
    nm <- names(idx)
    len <- vapply(idx, function(rec) rec$seqlength, integer(1))
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }
  # Fallback if you only have headers
  list(names = character(), lengths = integer())
}

# Pre-filter the GFF to drop rows with missing start/end
clean_gff_for_import <- function(gff_path, 
                                 drop_invalid = TRUE, 
                                 verbose = TRUE) {
  if (!file.exists(gff_path)) stop("File not found: ", gff_path)
  
  msg <- function(...) if (verbose) message(...)
  
  msg("Reading GFF lines...")
  lines <- readLines(gff_path, warn = FALSE)
  
  # Remove comment/directive lines and blank lines
  is_comment <- grepl("^\\s*#", lines)
  is_blank   <- !nzchar(trimws(lines))
  lines <- lines[!(is_comment | is_blank)]
  
  if (length(lines) == 0) stop("No feature lines found after removing comments: ", gff_path)
  
  # Split into fields by tab
  parts <- strsplit(lines, "\t", fixed = TRUE)
  lens  <- lengths(parts)
  
  # Keep rows with at least 9 fields; truncate extras
  keep_len <- lens >= 9
  dropped_len <- sum(!keep_len)
  
  if (dropped_len > 0 && verbose) {
    msg("Dropping ", dropped_len, " lines with fewer than 9 tab-separated fields.")
  }
  
  parts <- parts[keep_len]
  # Build a 9-column matrix
  mat <- do.call(rbind, lapply(parts, function(x) x[1:9]))
  df  <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  
  # Coerce start/end to integers; treat "." as NA
  to_int <- function(x) {
    x[x %in% c("", ".", "NA", "NaN")] <- NA_character_
    suppressWarnings(as.integer(x))
  }
  df$start <- to_int(df$start)
  df$end   <- to_int(df$end)
  
  # Drop invalid rows if requested
  keep <- rep(TRUE, nrow(df))
  invalid_coord <- is.na(df$start) | is.na(df$end)
  inverted      <- !invalid_coord & (df$start > df$end)
  
  if (drop_invalid) {
    dropped <- sum(invalid_coord | inverted)
    if (dropped > 0 && verbose) {
      msg("Dropping ", dropped, " rows with missing or inverted coordinates.")
    }
    df <- df[!(invalid_coord | inverted), , drop = FALSE]
  } else if (any(invalid_coord | inverted)) {
    warning("There are rows with missing or inverted coordinates; rtracklayer::import() may fail.")
  }
  
  if (nrow(df) == 0) stop("No valid rows remain after filtering.")
  
  # Write to a temporary GFF3 with the required directive
  tmp <- tempfile(fileext = ".gff3")
  con <- file(tmp, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("##gff-version 3", con)
  
  # Write rows with tab separation, no quotes, no headers
  write.table(
    df, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  if (verbose) msg("Clean GFF written to: ", tmp)
  tmp
}


#' Harmonize GFF seqlevels to FASTA headers
#'
#' Attempts, seriatim:
#' 1) Identity: if all GFF seqlevels already match FASTA headers.
#' 2) Positional numeric mapping: if GFF levels look like "1","3","4",... and counts match FASTA.
#' 3) Region-guided mapping: infer mapping from top-level 'region' features, preserving their order.
#' 4) (Optional) Length-guided mapping: match region widths to FASTA contig lengths greedily.
#'
#' Returns a modified genome_obj with gff seqlevels renamed to match genome_obj$seqnames.
#' Emits informative cli messages; uses cli::warn on partial successes; cli::abort only when hopeless.


harmonize_gff_seqlevels <- function(genome_obj,
                                    use_length_fallback = FALSE) {
  
  
  # --- Required namespaces (fail fast, clear) ---
  stopifnot(requireNamespace("cli", quietly = TRUE))
  stopifnot(requireNamespace("rtracklayer", quietly = TRUE))
  stopifnot(requireNamespace("Biostrings", quietly = TRUE))
  stopifnot(requireNamespace("Rsamtools", quietly = TRUE))
  stopifnot(requireNamespace("GenomeInfoDb", quietly = TRUE))
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("S4Vectors", quietly = TRUE))
  
  # --- Basic object checks ---
  stopifnot(
    is.list(genome_obj),
    all(c("gff","fa","seqnames") %in% names(genome_obj))
  )
  
  
  gff <- genome_obj$gff
  fa_headers <- genome_obj$seqnames
  
  # ---- Compatibility shims ---------------------------------------------------
  .set_seqlevels <- function(gr, keep, pruning = "coarse") {
    # Try modern form with pruning.mode, else fall back
    ok <- TRUE
    res <- try({
      GenomeInfoDb::seqlevels(gr, pruning.mode = pruning) <- keep
      gr
    }, silent = TRUE)
    if (inherits(res, "try-error")) {
      ok <- FALSE
    } else {
      return(res)
    }
    # Fallback (older GenomeInfoDb)
    GenomeInfoDb::seqlevels(gr) <- keep
    gr
  }
  
  .rename_seqlevels <- function(gr, map_named_vector, pruning = "coarse") {
    # Try modern form with pruning.mode, else fall back
    res <- try({
      GenomeInfoDb::renameSeqlevels(gr, map_named_vector, pruning.mode = pruning)
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    # Fallback (older GenomeInfoDb)
    GenomeInfoDb::renameSeqlevels(gr, map_named_vector)
  }
  # ---------------------------------------------------------------------------
  
  cli::cli_h2("Reconciling GFF seqlevels with FASTA headers")
  cli::cli_li("FASTA headers: {.val {fa_headers}}")
  cli::cli_li("GFF seqlevels: {.val {GenomeInfoDb::seqlevels(gff)}}")
  
  gff_levels <- GenomeInfoDb::seqlevels(gff)
  
  # 1) Identity
  if (all(gff_levels %in% fa_headers)) {
    keep <- intersect(fa_headers, gff_levels)
    gff <- .set_seqlevels(gff, keep)
    cli::cli_alert_success("GFF is already harmonized with FASTA. Using existing labels.")
    genome_obj$gff <- gff
    return(genome_obj)
  }
  
  # helper for renaming + validation
  .try_map <- function(map_named_vector, label = NULL) {
    gff2 <- .rename_seqlevels(gff, map_named_vector)
    if (all(GenomeInfoDb::seqlevels(gff2) %in% fa_headers)) {
      if (!is.null(label)) {
        cli::cli_li("Applied mapping ({.emph {label}}): {.field {paste(names(map_named_vector), '->', unname(map_named_vector), collapse=', ')}}")
      } else {
        cli::cli_li("Applied mapping: {.field {paste(names(map_named_vector), '->', unname(map_named_vector), collapse=', ')}}")
      }
      return(gff2)
    }
    NULL
  }
  
  # 2) Positional numeric mapping when counts match
  looks_numeric <- all(grepl("^[0-9]+$", gff_levels))
  if (looks_numeric && length(gff_levels) == length(fa_headers)) {
    cli::cli_alert_info(
      "Attempting positional numeric mapping: {.val {gff_levels}} -> {.val {fa_headers}}"
    )
    # Map each *existing* GFF level to the FASTA header at the same position
    # e.g., c("1","3","4","5","7") -> c("Chromosome","Plasmid_1_(...)","Plasmid_2_(...)","Plasmid_3_(...)","Plasmid_4_(...)")
    map <- stats::setNames(fa_headers, gff_levels)
    gff2 <- .try_map(map, label = "positional")
    if (!is.null(gff2)) {
      cli::cli_alert_success("Positional numeric mapping succeeded.")
      genome_obj$gff <- gff2
      return(genome_obj)
    } else {
      cli::cli_warn("Positional numeric mapping failed. Will try region-guided mapping next.")
    }
  } else {
    cli::cli_alert_info("Skipping positional numeric mapping (levels not purely numeric or counts differ).")
  }
  
  # 3) Region-guided mapping (order of region rows -> order of FASTA headers)
  has_type <- "type" %in% names(S4Vectors::mcols(gff))
  if (has_type) {
    is_region <- S4Vectors::mcols(gff)$type == "region"
    if (any(is_region)) {
      regions <- gff[is_region]
      reg_levels <- unique(as.character(GenomicRanges::seqnames(regions)))
      cli::cli_alert_info("Found {.strong region} rows. Region seqnames: {.val {reg_levels}}")
      
      if (length(reg_levels) == length(fa_headers)) {
        map <- stats::setNames(fa_headers, reg_levels)  # e.g., "1"->"Chromosome"
        cli::cli_alert_info("Attempting region-guided mapping.")
        gff2 <- .try_map(map, label = "region-guided")
        if (!is.null(gff2)) {
          cli::cli_alert_success("Region-guided mapping succeeded.")
          genome_obj$gff <- gff2
          return(genome_obj)
        } else {
          cli::cli_warn("Region-guided mapping did not converge.")
        }
      } else {
        cli::cli_warn("Region count does not match FASTA count. Expected {.val {length(fa_headers)}}, got {.val {length(reg_levels)}}.")
      }
    } else {
      cli::cli_alert_info("No region features available in GFF; skipping region-guided mapping.")
    }
  } else {
    cli::cli_alert_info("No 'type' column present in GFF metadata; cannot use region-guided mapping.")
  }
  
  # 4) Optional length-guided greedy matching
  if (isTRUE(use_length_fallback)) {
    cli::cli_h3("Length-guided fallback mapping")
    
    # FASTA lengths
    if (!is.null(genome_obj$fasta)) {
      fa_lengths <- Biostrings::width(genome_obj$fasta)
      names(fa_lengths) <- genome_obj$seqnames
    } else {
      fai <- Rsamtools::scanFaIndex(genome_obj$fa)
      fa_lengths <- fai$seqlengths
      names(fa_lengths) <- names(fai)
    }
    
    if (has_type && any(S4Vectors::mcols(gff)$type == "region")) {
      regions <- gff[S4Vectors::mcols(gff)$type == "region"]
      gff_lengths <- GenomicRanges::width(regions)
      names(gff_lengths) <- as.character(GenomicRanges::seqnames(regions))
      
      fa_left <- fa_lengths
      map <- character(length(gff_lengths))
      names(map) <- names(gff_lengths)
      for (k in names(gff_lengths)) {
        best <- names(which.min(abs(fa_left - gff_lengths[[k]])))
        map[[k]] <- best
        fa_left <- fa_left[names(fa_left) != best]
      }
      gff2 <- .try_map(map, label = "length-guided")
      if (!is.null(gff2)) {
        cli::cli_alert_success("Length-guided mapping succeeded.")
        genome_obj$gff <- gff2
        return(genome_obj)
      } else {
        cli::cli_warn("Length-guided mapping failed to produce a consistent vocabulary.")
      }
    } else {
      cli::cli_warn("Cannot compute region lengths without region features; skipping length fallback.")
    }
  }
  
  cli::cli_abort(c(
    "x" = "Unable to harmonize GFF seqlevels to FASTA headers.",
    "i" = paste0("GFF levels: ", paste(gff_levels, collapse = ", ")),
    "i" = paste0("FASTA: ", paste(fa_headers, collapse = ", ")),
    "i" = "Inspect 'region' rows and FASTA order, or enable use_length_fallback=TRUE."
  ))
}

#' Initialize Genome Resources
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param auto_index Logical; automatically index FASTA if needed
#' @return List with gff (GRanges), fasta (DNAStringSet), and fa (FaFile) objects
#' @export
#' 
#' 
