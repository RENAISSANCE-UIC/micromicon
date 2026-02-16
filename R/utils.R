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

# Setup and Data Loading Functions

# Note: %||% operator is defined in R/operators.R

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
      cli::cli_abort("Rsamtools is required to inspect FaFile indices")
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
  if (!file.exists(gff_path)) cli::cli_abort("File not found: {gff_path}")
  
  msg <- function(...) if (verbose) message(...)
  
  msg("Reading GFF lines...")
  lines <- readLines(gff_path, warn = FALSE)
  
  # Remove comment/directive lines and blank lines
  is_comment <- grepl("^\\s*#", lines)
  is_blank   <- !nzchar(trimws(lines))
  lines <- lines[!(is_comment | is_blank)]

  if (length(lines) == 0) cli::cli_abort("No feature lines found after removing comments: {gff_path}")
  
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
    cli::cli_warn("There are rows with missing or inverted coordinates; rtracklayer::import() may fail")
  }

  if (nrow(df) == 0) cli::cli_abort("No valid rows remain after filtering")
  
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
  
  # --- Basic object checks 
  stopifnot(
    is.list(genome_obj),
    all(c("gff","fa","seqnames") %in% names(genome_obj))
  )
  
  
  gff <- genome_obj$gff
  fa_headers <- genome_obj$seqnames
  
  # ---- Compatibility shims 
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
  # 
  
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



# Nullish fallback for scalars/length-1 values:
# Returns `b` iff a is NULL, length-0, or (length-1 AND NA); otherwise returns `a` unchanged.
scalar_coalesce <- function(a, b) {
  if (is.null(a))            return(b)
  if (length(a) == 0L)       return(b)
  if (length(a) == 1L && is.na(a)) return(b)
  a
}

# Fallback that ONLY treats NULL/length-0 as missing. Useful for sentinels that may be 0 or NA.
null_coalesce <- function(a, b) {
  if (is.null(a) || length(a) == 0L) b else a
}

# Element-wise coalesce (vector-safe). Keeps the type of the first non-NULL input.
# Recycles RHS vectors; ignores NULL args; treats NA as missing.
vec_coalesce <- function(...) {
  args <- list(...)
  # drop NULLs early
  args <- args[!vapply(args, is.null, logical(1))]
  if (length(args) == 0L) return(NULL)
  # target length
  lens <- vapply(args, length, integer(1))
  L <- max(lens)
  # recycle
  args <- lapply(args, function(x) if (length(x) == 0L) rep(NA, L) else rep_len(x, L))
  out <- args[[1]]
  if (length(args) > 1L) {
    for (k in 2:length(args)) {
      idx <- is.na(out)
      out[idx] <- args[[k]][idx]
    }
  }
  out
}

# Return the first non-NA element (or NA if none). Handy for "representative" choices.
first_non_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA)
  i <- which(!is.na(x))
  if (length(i)) x[i[1]] else x[1]
}

#' Convert genome_entity_gd events to a data frame
#'
#' Converts the nested event list from a genome_entity_gd object into a tidy
#' data frame with one row per event. Tags can either be expanded into individual
#' columns (tag_*) or kept as a list-column for programmatic access.
#'
#' @param gd A genome_entity_gd object
#' @param expand_tags Logical; if TRUE, flatten tags into separate tag_* columns;
#'   if FALSE, keep as a list-column
#' @param include_raw Logical; include the raw_line field
#' @param include_hash Logical; include the event hash field
#' @param types Character vector of event types to filter (e.g., c("SNP", "DEL"))
#' @param kinds Character vector of event kinds to filter: "mutation", "evidence", "validation"
#' @param n Maximum number of rows to return (default Inf for all)
#' @param cols Character vector of column names to prioritize in output order
#' @param stringsAsFactors Logical; passed to as.data.frame()
#' @return A data frame (or tibble if available) with one row per event
#' @export
#' @examples
#' \dontrun{
#'   # Get all events with tags as list-column
#'   ev_tbl <- gd_events_table(gd, expand_tags = FALSE)
#'
#'   # Get only mutations with expanded tags
#'   mut_tbl <- gd_events_table(gd, kinds = "mutation", expand_tags = TRUE)
#'
#'   # Filter to specific types
#'   snps <- gd_events_table(gd, types = c("SNP", "DEL"), n = 100)
#' }
gd_events_table <- function(
    gd,
    expand_tags = TRUE,
    include_raw = FALSE,
    include_hash = TRUE,
    types = NULL,
    kinds = NULL,
    n = Inf,
    cols = NULL,
    stringsAsFactors = FALSE
) {
  stopifnot(inherits(gd, "genome_entity_gd"))
  ev <- gd$events
  if (!length(ev)) return(.coerce_tbl(data.frame()))
  
  # Optional filters
  if (!is.null(types)) {
    keep <- vapply(ev, function(x) isTRUE(x$type %in% types), logical(1))
    ev <- ev[keep]
  }
  if (!is.null(kinds)) {
    keep <- vapply(ev, function(x) isTRUE(x$kind %in% kinds), logical(1))
    ev <- ev[keep]
  }
  if (!length(ev)) return(.coerce_tbl(data.frame()))
  
  # Helpers -
  .flatten_tags <- function(tags) {
    if (is.null(tags) || length(tags) == 0) return(list())
    out <- list()
    for (nm in names(tags)) {
      v <- tags[[nm]]
      if (length(v) >= 1) {
        v1  <- v[1]
        num <- suppressWarnings(as.numeric(v1))
        out[[paste0("tag_", nm)]] <-
          if (!is.na(num) && is.character(v1) && grepl("^[0-9.+-eE]+$", v1)) num else as.character(v1)
      }
    }
    out
  }
  
  # NEW: sanitize any length-0 to NA; collapse any length>1 to first element
  
  .scalarize_row <- function(r) {
    for (nm in names(r)) {
      v <- r[[nm]]
      if (is.list(v)) next                      # <-- keep list-cols intact
      if (is.null(v) || length(v) == 0L) {
        r[[nm]] <- NA
      } else if (length(v) > 1L) {
        r[[nm]] <- v[1]
      }
    }
    r
  }
  
  
  # Row assembly 
  rows <- lapply(ev, function(x) {
    row <- list(
      type     = x$type,
      kind     = x$kind,
      id       = x$id,
      rank     = x$rank,
      seq_id   = x$seq_id,
      position = x$position,
      contig   = x$contig,
      col6     = x$col6, col7 = x$col7, col8 = x$col8,
      
      snp_ref_base = x$snp_ref_base,
      snp_alt_base = x$snp_alt_base,
      del_size     = x$del_size,
      ins_seq      = x$ins_seq,
      ins_size     = x$ins_size,
      sub_size     = x$sub_size,
      sub_new_seq  = x$sub_new_seq,
      mob_repeat_name      = x$mob_repeat_name,
      mob_strand           = x$mob_strand,
      mob_duplication_size = x$mob_duplication_size,
      
      ev_frequency = x$ev_frequency,
      ev_quality   = x$ev_quality,
      ev_ref_cov_1 = x$ev_ref_cov_1, ev_ref_cov_2 = x$ev_ref_cov_2,
      ev_new_cov_1 = x$ev_new_cov_1, ev_new_cov_2 = x$ev_new_cov_2,
      ev_tot_cov_1 = x$ev_tot_cov_1, ev_tot_cov_2 = x$ev_tot_cov_2,
      ev_alignment_overlap = x$ev_alignment_overlap,
      ev_cov_minus = x$ev_cov_minus,
      ev_cov_plus  = x$ev_cov_plus,
      ev_pos_start = x$ev_pos_start,
      ev_pos_end   = x$ev_pos_end,
      
      ra_insert_position = x$ra_insert_position,
      ra_ref_base        = x$ra_ref_base,
      ra_new_base        = x$ra_new_base,
      
      start              = x$start,
      end                = x$end,
      mc_start_range     = x$mc_start_range,
      mc_end_range       = x$mc_end_range,
      
      side_1_seq_id      = x$side_1_seq_id,
      side_1_position    = x$side_1_position,
      side_1_strand      = x$side_1_strand,
      side_2_seq_id      = x$side_2_seq_id,
      side_2_position    = x$side_2_position,
      side_2_strand      = x$side_2_strand,
      jc_overlap         = x$jc_overlap
    )
    if (include_hash) row$hash <- x$hash
    if (include_raw)  row$raw_line <- x$raw_line
    
    if (isTRUE(expand_tags)) {
      row <- c(row, .flatten_tags(x$tags))
    } else {
      # Force a single list-column for tags, uniformly across rows
      row$tags <- I(list(if (is.null(x$tags)) NULL else x$tags))
    }
    
    .scalarize_row(row)  # ensure 1×1 scalars for all atomic fields
  })
  
  # Align names and bind 
  all_names <- Reduce(union, lapply(rows, names))
  rows_aligned <- lapply(rows, function(r) { r[setdiff(all_names, names(r))] <- NA; r[all_names] })
  
  # As data.frame, then optional tibble conversion
  df <- as.data.frame(
    do.call(rbind, lapply(rows_aligned, function(r) as.data.frame(r, stringsAsFactors = stringsAsFactors))),
    stringsAsFactors = stringsAsFactors, check.names = FALSE
  )
  
  # Light type nudges
  intish <- c("id","rank","position","del_size","ins_size","sub_size",
              "mob_strand","mob_duplication_size","ev_ref_cov_1","ev_ref_cov_2",
              "ev_new_cov_1","ev_new_cov_2","ev_tot_cov_1","ev_tot_cov_2",
              "ev_alignment_overlap","ev_cov_minus","ev_cov_plus","ev_pos_start","ev_pos_end",
              "mc_start_range","mc_end_range","side_1_position","side_2_position","jc_overlap")
  numish <- c("ev_frequency","ev_quality")
  for (nm in intersect(intish, names(df))) df[[nm]] <- suppressWarnings(as.integer(df[[nm]]))
  for (nm in intersect(numish, names(df))) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  
  # Preferred column order (optional)
  if (!is.null(cols)) {
    keep <- intersect(cols, names(df))
    df <- df[, c(keep, setdiff(names(df), keep)), drop = FALSE]
  }
  
  # Optional head
  if (is.finite(n)) df <- utils::head(df, n)
  
  .coerce_tbl(df)
}

.coerce_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Preview genome_entity_gd events in a compact table
#'
#' Displays a formatted preview of events from a genome_entity_gd object with
#' sensible column ordering and filtering. Returns the table invisibly for
#' further manipulation.
#'
#' @param gd A genome_entity_gd object
#' @param types Character vector of event types to filter (e.g., c("SNP", "RA"))
#' @param kinds Character vector of event kinds to filter: "mutation", "evidence", "validation"
#' @param n Number of rows to display (default 20)
#' @param expand_tags Logical; if TRUE, flatten tags into separate tag_* columns
#' @param include_raw Logical; include the raw GD line in output
#' @return The previewed data frame/tibble (invisibly)
#' @export
#' @examples
#' \dontrun{
#'   # Preview all events
#'   view_events(gd)
#'
#'   # View only SNPs
#'   view_events(gd, types = "SNP", n = 50)
#'
#'   # View evidence events
#'   view_events(gd, kinds = "evidence")
#' }
view_events <- function(
    gd,
    types = NULL,
    kinds = NULL,
    n = 20,
    expand_tags = TRUE,
    include_raw = FALSE
) {
  key_cols <- c(
    "type","kind","id","rank","seq_id","position","contig",
    "snp_ref_base","snp_alt_base","del_size","ins_seq","ins_size","sub_size",
    "mob_repeat_name","mob_strand","mob_duplication_size",
    "ra_insert_position","ra_ref_base","ra_new_base",
    "start","end","mc_start_range","mc_end_range",
    "side_1_seq_id","side_1_position","side_1_strand",
    "side_2_seq_id","side_2_position","side_2_strand","jc_overlap",
    "ev_frequency","ev_quality"
  )
  tag_pref <- c("tag_gene_name","tag_locus_tag","tag_genes_overlapping",
                "tag_mutation_category","tag_snp_type","tag_frequency")
  
  tbl <- gd_events_table(
    gd,
    expand_tags = expand_tags,
    include_raw = include_raw,
    include_hash = TRUE,
    types = types,
    kinds = kinds,
    n = n,
    cols = c(key_cols, tag_pref, "hash")
  )
  
  # Order: descending id, then type
  ord <- order(-as.integer(tbl$id), tbl$type, na.last = TRUE)
  tbl <- tbl[ord, , drop = FALSE]
  
  print(tbl, n = min(nrow(tbl), n))
  invisible(tbl)
}


#' Preview mutation events only
#'
#' Convenience wrapper for \code{view_events()} that filters to mutation events
#' only (SNP, SUB, DEL, INS, MOB, AMP, CON, INV).
#'
#' @param gd A genome_entity_gd object
#' @param n Number of rows to display (default 20)
#' @param expand_tags Logical; if TRUE, flatten tags into separate tag_* columns
#' @return The previewed data frame/tibble (invisibly)
#' @export
#' @seealso \code{\link{view_events}}, \code{\link{view_evidence}}
#' @examples
#' \dontrun{
#'   view_mutations(gd)
#'   view_mutations(gd, n = 50, expand_tags = FALSE)
#' }
view_mutations <- function(gd, n = 20, expand_tags = TRUE) {
  view_events(gd, kinds = "mutation", n = n, expand_tags = expand_tags)
}

#' Preview evidence events only
#'
#' Convenience wrapper for \code{view_events()} that filters to evidence events
#' only (RA, MC, JC, UN) and displays evidence-specific columns (coverage,
#' frequency, quality scores).
#'
#' @param gd A genome_entity_gd object
#' @param n Number of rows to display (default 25)
#' @param expand_tags Logical; if TRUE, flatten tags into separate tag_* columns
#' @return The previewed data frame/tibble (invisibly)
#' @export
#' @seealso \code{\link{view_events}}, \code{\link{view_mutations}}
#' @examples
#' \dontrun{
#'   view_evidence(gd)
#'   view_evidence(gd, n = 100)
#' }
view_evidence <- function(gd, n = 25, expand_tags = TRUE) {
  key_cols <- c(
    "type","id","seq_id","position","contig",
    # RA
    "ra_insert_position","ra_ref_base","ra_new_base",
    # MC
    "start","end","mc_start_range","mc_end_range",
    # JC
    "side_1_seq_id","side_1_position","side_1_strand",
    "side_2_seq_id","side_2_position","side_2_strand","jc_overlap",
    # UN
    # (uses start/end already)
    # Evidence numerics
    "ev_frequency","ev_quality","ev_ref_cov_1","ev_ref_cov_2",
    "ev_new_cov_1","ev_new_cov_2","ev_tot_cov_1","ev_tot_cov_2",
    "ev_alignment_overlap","ev_cov_minus","ev_cov_plus"
  )
  tag_pref <- c("tag_frequency","tag_mutation_category","tag_gene_name","tag_locus_tag")
  
  tbl <- gd_events_table(
    gd,
    kinds = "evidence",
    expand_tags = expand_tags,
    include_raw = FALSE,
    include_hash = TRUE,
    n = n,
    cols = c(key_cols, tag_pref, "hash")
  )
  
  # Sort: type then id
  ord <- order(tbl$type, as.integer(tbl$id), na.last = TRUE)
  tbl <- tbl[ord, , drop = FALSE]
  
  print(tbl, n = min(nrow(tbl), n))
  invisible(tbl)
}


#' Diagnostic probe for RA event parsing
#'
#' Shows the raw GD lines alongside parsed fields for RA (Read Alignment) events,
#' useful for debugging parser behavior or understanding breseq output.
#'
#' @param gd A genome_entity_gd object
#' @param k Maximum number of RA events to examine (default 10)
#' @return A list of lists with raw_line and parsed fields for each RA event
#' @export
#' @examples
#' \dontrun{
#'   probe <- probe_ra(gd, k = 5)
#'   probe[[1]]$raw      # Original GD line
#'   probe[[1]]$parsed   # Parsed fields
#' }
probe_ra <- function(gd, k = 10) {
  idx <- which(vapply(gd$events, function(e) identical(e$type, "RA"), logical(1)))
  idx <- head(idx, k)
  lapply(idx, function(i) {
    e <- gd$events[[i]]
    list(
      i = i,
      raw = e$raw_line,
      parsed = list(
        seq_id = e$seq_id,
        position = e$position,
        insert_position = e$ra_insert_position,
        ref_base = e$ra_ref_base,
        new_base = e$ra_new_base
      )
    )
  })
}

# Return a tidy two-column data.frame: key | value
tags_as_kv <- function(tags) {
  if (is.null(tags) || !length(tags)) {
    return(data.frame(key = character(0), value = character(0), stringsAsFactors = FALSE))
  }
  # Flatten multi-valued tags; keep first value index for traceability
  keys   <- rep.int(names(tags), vapply(tags, length, integer(1)))
  values <- unlist(tags, use.names = FALSE)
  idx    <- sequence(vapply(tags, length, integer(1)))  # 1,2,... per key
  data.frame(
    key   = as.character(keys),
    index = as.integer(idx),
    value = as.character(values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Peek at tags for a specific event
#'
#' Displays the tags for a single event from an events table in a tidy
#' key-value format, useful for exploring tag contents.
#'
#' @param ev_tbl A data frame from \code{gd_events_table()} with a tags list-column
#' @param i Row index of the event to inspect (default 1)
#' @param n Maximum number of tag key-value pairs to display (default 20)
#' @return A data frame with columns: key, index, value
#' @export
#' @seealso \code{\link{gd_events_table}}
#' @examples
#' \dontrun{
#'   ev_tbl <- gd_events_table(gd, expand_tags = FALSE)
#'   peek_tags(ev_tbl, i = 1, n = 25)
#' }
peek_tags <- function(ev_tbl, i = 1L, n = 20L) {
  kv <- tags_as_kv(ev_tbl$tags[[i]])
  if (nrow(kv) == 0) return(kv)
  utils::head(kv[order(kv$key, kv$index), , drop = FALSE], n)
}


#' Extract first value of a tag into a column
#'
#' Extracts the first value of a multi-valued tag from the tags list-column
#' and creates a new scalar column. Returns NA if the tag is absent.
#'
#' @param tbl A data frame from \code{gd_events_table()} with a tags list-column
#' @param name Name of the tag to extract (without "tag_" prefix)
#' @param into Name for the new column (default: "tag_{name}")
#' @return The input data frame with the new column added
#' @export
#' @seealso \code{\link{tag_get_concat}}, \code{\link{gd_events_table}}
#' @examples
#' \dontrun{
#'   ev_tbl <- gd_events_table(gd, expand_tags = FALSE)
#'   ev_tbl <- tag_get_first(ev_tbl, "gene_name")
#'   ev_tbl$tag_gene_name  # New column with first gene_name value
#' }
tag_get_first <- function(tbl, name, into = paste0("tag_", name)) {
  tbl[[into]] <- vapply(tbl$tags, function(t) {
    if (is.null(t) || is.null(t[[name]]) || !length(t[[name]])) return(NA_character_)
    as.character(t[[name]][1])
  }, character(1))
  tbl
}

#' Concatenate all values of a tag into a column
#'
#' Extracts all values of a multi-valued tag from the tags list-column,
#' concatenates them with a separator, and creates a new column.
#'
#' @param tbl A data frame from \code{gd_events_table()} with a tags list-column
#' @param name Name of the tag to extract (without "tag_" prefix)
#' @param into Name for the new column (default: "tag_{name}")
#' @param sep Separator for concatenating multiple values (default "|")
#' @return The input data frame with the new column added
#' @export
#' @seealso \code{\link{tag_get_first}}, \code{\link{gd_events_table}}
#' @examples
#' \dontrun{
#'   ev_tbl <- gd_events_table(gd, expand_tags = FALSE)
#'   ev_tbl <- tag_get_concat(ev_tbl, "locus_tag", sep = ";")
#'   ev_tbl$tag_locus_tag  # New column with concatenated values
#' }
tag_get_concat <- function(tbl, name, into = paste0("tag_", name), sep = "|") {
  tbl[[into]] <- vapply(tbl$tags, function(t) {
    if (is.null(t) || is.null(t[[name]]) || !length(t[[name]])) return(NA_character_)
    paste0(as.character(t[[name]]), collapse = sep)
  }, character(1))
  tbl
}

#' Pretty-print a single event
#'
#' Prints a human-readable summary of a single event, with type-specific
#' formatting for RA, JC, MC, and UN events.
#'
#' @param gd A genome_entity_gd object
#' @param i Index of the event to print (1-based)
#' @return The event object (invisibly)
#' @export
#' @examples
#' \dontrun{
#'   print_event(gd, 10)  # Print event #10
#'   print_event(gd, gd$provenance$by_kind$evidence_idx[1])  # First evidence event
#' }
print_event <- function(gd, i) {
  e <- gd$events[[i]]
  cat(sprintf("[%s] id=%s kind=%s\n", e$type, e$id, e$kind))
  if (e$type == "RA") {
    cat(sprintf("  %s:%s  ref=%s new=%s ins=%s  freq=%s\n",
                e$seq_id, e$position, e$ra_ref_base, e$ra_new_base,
                e$ra_insert_position, e$ev_frequency))
  } else if (e$type == "JC") {
    cat(sprintf("  %s:%s (%+d)  <->  %s:%s (%+d)  overlap=%s\n",
                e$side_1_seq_id, e$side_1_position, e$side_1_strand,
                e$side_2_seq_id, e$side_2_position, e$side_2_strand,
                e$jc_overlap))
  } else if (e$type == "MC") {
    cat(sprintf("  %s:[%s..%s]  ranges: +%s / -%s\n",
                e$seq_id, e$start, e$end, e$mc_start_range, e$mc_end_range))
  } else if (e$type == "UN") {
    cat(sprintf("  %s:[%s..%s]\n", e$seq_id, e$start, e$end))
  } else {
    cat("  (mutation row)\n")
  }
  # print a few tags
  kv <- utils::head(tags_as_kv(e$tags), 8)
  if (nrow(kv)) {
    cat("  tags:\n")
    apply(kv, 1, function(r) cat(sprintf("    - %s[%s]: %s\n", r["key"], r["index"], r["value"])))
  }
  invisible(e)
}

#' Format frequencies as percentages
#'
#' Converts numeric frequencies (0-1 range) to percentage strings for display.
#'
#' @param x Numeric vector of frequencies (values between 0 and 1)
#' @param digits Number of decimal places (default 2)
#' @return Character vector of formatted percentages (e.g., "85.50%")
#' @export
#' @examples
#' \dontrun{
#'   format_freq(c(0.855, 0.12, NA))  # "85.50%" "12.00%" NA
#'   format_freq(0.333, digits = 1)   # "33.3%"
#' }
format_freq <- function(x, digits = 2) {
  ifelse(is.na(x), NA, sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

