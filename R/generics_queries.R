#' Search Features
#'
#' @description
#' Search for features matching criteria (type, name, coordinates).
#' Works on \code{genome_entity} objects and, by S3 dispatch inheritance,
#' on \code{genome_entity_gd} objects as well.
#'
#' @param x A genome_entity (or genome_entity_gd) object
#' @param type Character; filter by feature type (e.g., "CDS", "gene").
#'   Case-insensitive. Note: Available types vary by annotation source.
#'   Bacterial annotations (Prokka, breseq) typically contain only "CDS"
#'   features without separate "gene" features. Use \code{unique(x$features$type)}
#'   to see what's available in your data.
#' @param pattern Character string to search for in feature annotations.
#'   Searches across ID, Name, Alias, gene, locus_tag, and product fields
#'   (case-insensitive). Compatible with various GFF3 annotation sources
#'   including breseq, Prokka, and NCBI.
#' @param seqname Character; filter by sequence name
#' @param start Integer; filter features starting at or after this position
#' @param end Integer; filter features ending at or before this position
#' @param strand Character; filter by strand ("+", "-")
#' @param ... Additional arguments (currently unused)
#'
#' @return data.frame of matching features
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Check what feature types are available
#' unique(genome$features$type)
#'
#' # Find all CDS features (most common in bacterial genomes)
#' cds <- search_features(genome, type = "CDS")
#'
#' # Find CDS on chr1
#' cds_chr1 <- search_features(genome, type = "CDS", seqname = "chr1")
#'
#' # Find features matching pattern (searches gene names, products, etc.)
#' dna_features <- search_features(genome, pattern = "dna")
#' }
search_features <- function(x, ...) {
  UseMethod("search_features")
}

#' @export
search_features.genome_entity <- function(x, type = NULL, pattern = NULL,
                                         seqname = NULL, start = NULL,
                                         end = NULL, strand = NULL, ...) {
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[tolower(feats$type) == tolower(type), ]
  }

  if (!is.null(pattern)) {
    matches <- rep(FALSE, nrow(feats))
    for (field in c("ID", "Name", "Alias", "gene", "locus_tag", "product")) {
      if (field %in% names(feats)) {
        field_matches <- grepl(pattern, feats[[field]], ignore.case = TRUE)
        matches <- matches | field_matches
      }
    }
    feats <- feats[matches, ]
  }

  if (!is.null(seqname) && "seqname" %in% names(feats)) {
    feats <- feats[feats$seqname == seqname, ]
  }

  if (!is.null(start) && "start" %in% names(feats)) {
    feats <- feats[feats$start >= start, ]
  }

  if (!is.null(end) && "end" %in% names(feats)) {
    feats <- feats[feats$end <= end, ]
  }

  if (!is.null(strand) && "strand" %in% names(feats)) {
    feats <- feats[feats$strand == strand, ]
  }

  feats
}

#' @export
search_features.default <- function(x, ...) {
  cli::cli_abort("search_features() not implemented for class {.cls {class(x)[1]}}")
}


#' Extract Sequences by Coordinates
#'
#' @param x Object to extract from
#' @param ... Additional arguments
#' @export
extract_by_coords <- function(x, ...) {
  UseMethod("extract_by_coords")
}

#' @export
extract_by_coords.genome_entity <- function(x, seqname, start, end,
                                           strand = "+",
                                           translate = FALSE,
                                           names = NULL, ...) {
  validate_genome_entity(x)
  options <- list(strand = strand, translate = translate, names = names)
  execute_extract_sequences_by_coords(x, seqname, start, end, options)
}

#' @export
extract_by_coords.default <- function(x, ...) {
  cli::cli_abort("extract_by_coords() not implemented for class {.cls {class(x)[1]}}")
}


#' Extract Sequences by Name
#'
#' @description
#' Extract DNA (or protein) sequences for features matching a name pattern.
#' Uses Bioconductor's GRanges and FaFile for robust extraction, including
#' automatic seqname harmonization between GFF and FASTA.
#'
#' @param x A genome_entity (or genome_entity_gd) object
#' @param pattern Character string or regex pattern (perl-compatible) to match
#'   feature names. Searches Name field first; falls back to gene, locus_tag,
#'   and product if Name yields no matches.
#' @param translate Logical; if TRUE, translate extracted CDS to amino acids
#'   (default FALSE)
#' @param genetic_code NCBI genetic code ID (default "11" for bacterial/archaeal)
#' @param auto_harmonize Logical; if TRUE, automatically harmonize seqname
#'   vocabulary between GFF and FASTA when mismatches are detected (default TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A DNAStringSet (or AAStringSet if translate = TRUE) of extracted sequences
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Extract DNA sequence for dnaA
#' seqs <- extract_by_name(genome, "dnaA")
#'
#' # Translate to protein
#' proteins <- extract_by_name(genome, "dnaA", translate = TRUE)
#' }
extract_by_name <- function(x, ...) {
  UseMethod("extract_by_name")
}

#' @export
extract_by_name.genome_entity <- function(x, pattern,
                                         translate = FALSE,
                                         genetic_code = "11",
                                         auto_harmonize = TRUE,
                                         ...) {
  # Check for required Bioconductor packages
  missing_pkgs <- character()
  required_pkgs <- c("Biostrings", "GenomeInfoDb", "GenomicRanges", "S4Vectors", "Rsamtools")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }
  if (length(missing_pkgs) > 0) {
    cli::cli_abort(c(
      "Missing required Bioconductor packages for extract_by_name().",
      "x" = "Missing: {paste(missing_pkgs, collapse = ', ')}",
      "i" = "Install all required packages with:",
      " " = "if (!require('BiocManager')) install.packages('BiocManager')",
      " " = "BiocManager::install(c('GenomicRanges', 'Biostrings', 'IRanges', 'S4Vectors', 'GenomeInfoDb', 'Rsamtools'))"
    ))
  }

  # Convert to legacy format for Bioconductor GRanges path
  genome_obj <- entity_to_legacy_genome_obj(x)

  if (is.null(genome_obj$gff)) {
    cli::cli_abort(c(
      "genome_obj$gff is NULL - this function requires GRanges annotation data.",
      "i" = "Did you create the genome object with read_genome()?",
      "i" = "Make sure Bioconductor packages (GenomicRanges, rtracklayer) are installed.",
      "i" = "Try reinstalling: BiocManager::install(c('GenomicRanges', 'rtracklayer', 'Biostrings'))"
    ))
  }

  gff <- genome_obj$gff
  fa_headers <- genome_obj$seqnames

  # ---------- helpers ----------
  .matches <- function(v, pat) {
    if (is.null(v)) return(rep(FALSE, length(GenomicRanges::seqnames(gff))))
    y <- tryCatch(as.character(v), error = function(e) rep(NA_character_, length(v)))
    !is.na(y) & grepl(pat, y, perl = TRUE)
  }
  .prefer_name_then_fallback <- function(gff, pat) {
    name_idx <- .matches(gff$Name, pat)
    if (any(name_idx, na.rm = TRUE)) return(gff[name_idx])
    idx <- .matches(gff$gene, pat) |
      .matches(gff$locus_tag, pat) |
      .matches(gff$product, pat)
    idx[is.na(idx)] <- FALSE
    gff[idx]
  }
  .drop_extdb_dupes <- function(gr) {
    if (length(gr) == 0L) return(gr)
    is_extdb <- !is.na(gr$Name) & grepl("^extdb:", gr$Name)
    if (!any(is_extdb)) return(gr)
    key <- if ("locus_tag" %in% names(S4Vectors::mcols(gr))) gr$locus_tag else gr$gene
    if (is.null(key)) key <- gr$Name
    keep <- rep(TRUE, length(gr))
    if (!all(is.na(key))) {
      keys_ext <- key[is_extdb]
      keys_good <- key[!is_extdb]
      drop_ext <- is_extdb & key %in% keys_good
      keep <- !drop_ext
    } else {
      keep <- !is_extdb
    }
    gr[keep]
  }
  .best_labels <- function(gr) {
    n <- length(gr)
    out <- rep(NA_character_, n)
    cand <- list(gr$gene, gr$locus_tag, gr$Name, gr$product)
    for (v in cand) {
      v <- tryCatch(as.character(v), error = function(e) rep(NA_character_, n))
      fill <- is.na(out) & !is.na(v)
      out[fill] <- v[fill]
    }
    all_extdb <- !is.na(out) & grepl("^extdb:", out)
    if (any(all_extdb)) {
      repl <- if ("locus_tag" %in% names(S4Vectors::mcols(gr))) gr$locus_tag else gr$gene
      repl <- tryCatch(as.character(repl), error = function(e) rep(NA_character_, n))
      swap <- all_extdb & !is.na(repl)
      out[swap] <- repl[swap]
    }
    out[is.na(out)] <- paste0("feature_", seq_len(n))[is.na(out)]
    make.unique(out, sep = "_")
  }
  # ---------- /helpers ----------

  features <- .prefer_name_then_fallback(gff, pattern)
  if (length(features) == 0L) {
    cli::cli_warn("No features found matching pattern: {.val {pattern}}")
    return(NULL)
  }

  # Filter by feature type based on translate flag
  if (isTRUE(translate)) {
    if ("type" %in% names(S4Vectors::mcols(features))) {
      cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
      if (any(cds_ok)) features <- features[cds_ok]
    }
  } else {
    if ("type" %in% names(S4Vectors::mcols(features))) {
      gene_ok <- S4Vectors::mcols(features)$type %in% c("gene")
      if (any(gene_ok)) {
        features <- features[gene_ok]
      } else {
        cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
        if (any(cds_ok)) features <- features[cds_ok]
      }
    }
  }

  features <- .drop_extdb_dupes(features)
  cli::cli_alert_success("Selected {length(features)} feature{?s} for extraction.")

  # Seqname vocabulary check and optional auto-harmonization
  fa_headers <- genome_obj$seqnames
  current_vocab_ok <- all(as.character(GenomicRanges::seqnames(features)) %in% fa_headers)
  if (!current_vocab_ok && isTRUE(auto_harmonize)) {
    cli::cli_alert_info("Seqname vocabulary mismatch detected. Invoking harmonizer...")
    genome_obj <- harmonize_gff_seqlevels(genome_obj)
    gff <- genome_obj$gff
    features <- .prefer_name_then_fallback(gff, pattern)
    if (isTRUE(translate)) {
      if ("type" %in% names(S4Vectors::mcols(features))) {
        cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
        if (any(cds_ok)) features <- features[cds_ok]
      }
    } else {
      if ("type" %in% names(S4Vectors::mcols(features))) {
        gene_ok <- S4Vectors::mcols(features)$type %in% c("gene")
        if (any(gene_ok)) {
          features <- features[gene_ok]
        } else {
          cds_ok <- S4Vectors::mcols(features)$type %in% c("CDS")
          if (any(cds_ok)) features <- features[cds_ok]
        }
      }
    }
    features <- .drop_extdb_dupes(features)
    fa_headers <- genome_obj$seqnames
    current_vocab_ok <- all(as.character(GenomicRanges::seqnames(features)) %in% fa_headers)
    if (!current_vocab_ok) {
      bad <- setdiff(unique(as.character(GenomicRanges::seqnames(features))), fa_headers)
      cli::cli_abort(c(
        "x" = "Seqnames mismatch persists after harmonization.",
        "i" = paste0("Offending seqnames: ", paste(bad, collapse = ", ")),
        "i" = "Inspect GFF 'region' rows and FASTA headers for order or naming issues."
      ))
    } else {
      cli::cli_alert_success("Seqname vocabulary reconciled; proceeding.")
    }
  } else if (!current_vocab_ok) {
    bad <- setdiff(unique(as.character(GenomicRanges::seqnames(features))), fa_headers)
    cli::cli_abort(c(
      "x" = "Seqnames mismatch between GFF and FASTA.",
      "i" = paste0("Offending seqnames: ", paste(bad, collapse = ", ")),
      "i" = "Set auto_harmonize = TRUE or run harmonize_gff_seqlevels() beforehand."
    ))
  }

  # Extract sequences — prefer FaFile (GFF+FASTA) over DNAStringSet (GenBank)
  if (!is.null(genome_obj$fa)) {
    dna_seqs <- Biostrings::getSeq(genome_obj$fa, features)
  } else if (!is.null(genome_obj$fasta)) {
    seqs_list <- lapply(seq_along(features), function(i) {
      feat <- features[i]
      seqname <- as.character(GenomicRanges::seqnames(feat))
      start_pos <- GenomicRanges::start(feat)
      end_pos <- GenomicRanges::end(feat)
      strand_val <- as.character(GenomicRanges::strand(feat))
      seq <- Biostrings::subseq(genome_obj$fasta[[seqname]], start = start_pos, end = end_pos)
      if (strand_val == "-") seq <- Biostrings::reverseComplement(seq)
      seq
    })
    dna_seqs <- Biostrings::DNAStringSet(seqs_list)
  } else {
    cli::cli_abort(c(
      "Cannot extract sequences: both genome_obj$fa and genome_obj$fasta are NULL.",
      "i" = "The genome object may be corrupted or incomplete."
    ))
  }

  names(dna_seqs) <- .best_labels(features)

  if (isTRUE(translate)) {
    cli::cli_alert_info("Translating to amino acids using genetic code {.val {genetic_code}}")
    aa_seqs <- Biostrings::translate(
      dna_seqs,
      genetic.code = Biostrings::getGeneticCode(genetic_code)
    )
    return(aa_seqs)
  }

  dna_seqs
}

#' @export
extract_by_name.default <- function(x, ...) {
  cli::cli_abort("extract_by_name() not implemented for class {.cls {class(x)[1]}}")
}


# Backward compatibility aliases
#' @rdname extract_by_coords
#' @export
extract_sequences_by_coords <- extract_by_coords

#' @rdname extract_by_name
#' @export
extract_sequences_by_name <- extract_by_name

#' Predict Mutations
#'
#' @description
#' Extracts and presents mutations from a genome data object.
#'
#' @param gd A genome data object
#' @param ... Additional arguments passed to methods
#' @param min_freq Minimum frequency threshold (default: 0)
#' @param include_structural Include structural variants (default: TRUE)
#' @param join Join method for multi-gene annotations (default: "slash")
#'
#' @return A data frame of mutations with friendly preview
#' @export
predict_mutations <- function(gd, ...) {
  UseMethod("predict_mutations")
}

#' Public wrapper: predict mutations (opinionated defaults; args pass-through)
#'
#' Calls predict_mutations_int(), announces row count, and prints a preview:
#' - If dplyr is installed: prints tibble (contract unaffected).
#' - If dplyr is not installed: prints base data.frame.
#' Any preview/logging errors are swallowed; the table is still returned.
#'
#' @param gd genome_entity_gd
#' @param ... reserved for future options passed to *_int (kept for forward-compat)
#' @param min_freq numeric scalar, keep rows with freq >= min_freq (default 0)
#' @param include_structural logical, include JC/MC/etc. (default TRUE)
#' @param join one of c("slash","pipe","newline") for multi-item fields (default "slash")
#' @return data.frame/tibble with predicted mutations (superset schema)
#' @export
#' @rdname predict_mutations
#' @export
predict_mutations.genome_entity_gd <- function(gd, ...,
                                               min_freq = 0,
                                               include_structural = TRUE,
                                               join = c("slash", "pipe", "newline")) {
  gd_assert(gd, "gd")
  join <- match.arg(join)
  
  # Build once via internal, regardless of downstream printing path
  tbl <- predict_mutations_int(
    gd,
    min_freq = min_freq,
    include_structural = include_structural,
    join = join,
    ...
  )
  
  # Announce + print a friendly preview, but never fail the call if printing goes sideways
  tryCatch(
    {
      n <- NROW(tbl)
      # row count announcement
      cli::cli_inform(c(
        "i" = sprintf("predict_mutations: %d row%s.", n, if (n == 1L) "" else "s")
      ))
      
      # choose printing strategy
      has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
      n_show <- min(n, 25L)
      
      if (has_dplyr) {
        # tibble path (do not modify tbl; only coerce for printing)
        cli::cli_inform(c(
          "i" = sprintf("Showing top %d as a tibble.", n_show)
        ))
        print(utils::head(dplyr::as_tibble(tbl), n_show))
        tbl <- dplyr::as_tibble(tbl)
      } else {
        cli::cli_inform(c(
          "!" = "{.pkg dplyr} not detected; showing a base data.frame preview.",
          "i" = "Install with: {.code install.packages('dplyr')} to prefer tibble display."
        ))
        print(utils::head(tbl, n_show))
      }
    },
    error = function(e) {
      # Stay serene; log and proceed to return the table
      cli::cli_warn(c(
        "!" = "Non-fatal issue while printing the preview.",
        ">" = conditionMessage(e)
      ))
    }
  )
  
  # Always return the fully built table (but quietly)
  invisible(tbl)
}

#' Compute Mutation Effects
#'
#' @description
#' Enriches mutation table with molecular consequences. Automatically selects
#' the best implementation based on operating system:
#' - Linux/macOS: Uses parallel processing for speed
#' - Windows: Uses fast serial processing (parallel not supported)
#'
#' @param gd A genome data object
#' @param pm_tbl A data frame from \code{predict_mutations()}
#' @param ... Additional arguments passed to methods
#'
#' @return Enriched data frame with consequence columns
#' @export
compute_effects <- function(gd, pm_tbl, ...) {
  UseMethod("compute_effects")
}

#' @export
compute_effects.genome_entity_gd <- function(gd, pm_tbl, ...) {
  # Dispatch to parallel on Linux/macOS, fast serial on Windows
  if (.Platform$OS.type == "windows") {
    pm_enrich_consequences_fast(gd, pm_tbl, ...)
  } else {
    pm_enrich_consequences_parallel(gd, pm_tbl, ...)
  }
}


#' Get Features Overlapping a Region of Interest
#'
#' @description
#' Returns a tidy tibble of genomic features that overlap a specified region,
#' with optional flanking extension. Works on \code{genome_entity} objects and,
#' by S3 dispatch inheritance, on \code{genome_entity_gd} objects as well.
#'
#' @param x A genome_entity (or genome_entity_gd) object
#' @param seqname Character; chromosome/contig name (numeric-like strings such
#'   as \code{"1"} are resolved automatically against the GFF seqlevels)
#' @param start Integer; start coordinate (1-based, inclusive)
#' @param end Integer; end coordinate (1-based, inclusive)
#' @param flank Integer; extend the query region by this many bp on each side
#'   (default 0)
#' @param feature_type Character; keep only features of this type, e.g.
#'   \code{"CDS"} (default). Pass \code{NULL} to return all feature types.
#' @param auto_resolve Logical; attempt to reconcile seqname to a GFF seqlevel
#'   when an exact match is not found (default TRUE)
#' @param drop_extdb Logical; replace \code{extdb:*} placeholder names with
#'   gene/locus_tag labels when available (default TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A tibble with columns: seqnames, start, end, width, strand, type,
#'   Name, Alias, Note
#' @export
#'
#' @examples
#' \dontrun{
#' gd <- parse_gd_annotated("sample.gd", ref_dir = "reference/")
#'
#' # Features in a 2.4 Mb region with 10 kb flanks
#' feats <- get_roi_features(gd, "1", 2369501, 2613000, flank = 10000)
#'
#' # All feature types, no flanks
#' feats <- get_roi_features(gd, "1", 2369501, 2613000, feature_type = NULL)
#' }
get_roi_features <- function(x, ...) {
  UseMethod("get_roi_features")
}

#' @export
get_roi_features.genome_entity <- function(x, seqname, start, end,
                                           flank = 0,
                                           feature_type = "CDS",
                                           auto_resolve = TRUE,
                                           drop_extdb = TRUE,
                                           ...) {
  analyze_roi(
    genome_obj   = x,
    seqname      = seqname,
    start        = start,
    end          = end,
    flank        = flank,
    feature_type = feature_type,
    tidy         = TRUE,
    auto_resolve = auto_resolve,
    drop_extdb   = drop_extdb
  )
}

#' @export
get_roi_features.default <- function(x, ...) {
  cli::cli_abort("get_roi_features() not implemented for class {.cls {class(x)[1]}}")
}
