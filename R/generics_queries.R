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
#' @param contig Character; filter by contig/sequence name
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
#' cds_chr1 <- search_features(genome, type = "CDS", contig = "chr1")
#'
#' # Find features matching pattern (searches gene names, products, etc.)
#' dna_features <- search_features(genome, pattern = "dna")
#' }
search_features <- function(x, ...) {
  UseMethod("search_features")
}

#' @rdname search_features
#' @export
search_features.genome_entity <- function(x, type = NULL, pattern = NULL,
                                         contig = NULL, start = NULL,
                                         end = NULL, strand = NULL, ...) {
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[tolower(feats$type) == tolower(type), ]
  }

  if (!is.null(pattern)) {
    name_fields <- intersect(c("ID", "Name", "Alias", "gene", "locus_tag"), names(feats))
    all_fields  <- intersect(c("ID", "Name", "Alias", "gene", "locus_tag", "product"), names(feats))
    mask <- rep(FALSE, nrow(feats))

    # Tier 1: exact match (case-insensitive) on primary name fields
    for (f in name_fields) {
      v    <- as.character(feats[[f]])
      mask <- mask | (!is.na(v) & tolower(v) == tolower(pattern))
    }

    if (!any(mask)) {
      # Tier 2: word-boundary match across all fields
      wb_pat <- paste0("\\b", pattern, "\\b")
      for (f in all_fields) {
        v    <- as.character(feats[[f]])
        mask <- mask | (!is.na(v) & grepl(wb_pat, v, ignore.case = TRUE, perl = TRUE))
      }
    }

    if (!any(mask)) {
      # Tier 3: substring match (original behaviour)
      for (f in all_fields) {
        mask <- mask | grepl(pattern, feats[[f]], ignore.case = TRUE)
      }
    }

    feats <- feats[mask, ]
  }

  if (!is.null(contig) && "seqname" %in% names(feats)) {
    feats <- feats[feats$seqname == contig, ]
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
extract_by_coords.genome_entity <- function(x, contig, start, end,
                                           strand = "+",
                                           translate = FALSE,
                                           names = NULL, ...) {
  validate_genome_entity(x)
  options <- list(strand = strand, translate = translate, names = names)
  execute_extract_sequences_by_coords(x, contig, start, end, options)
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

#' @rdname extract_by_name
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

  # Extract sequences -- prefer FaFile (GFF+FASTA) over DNAStringSet (GenBank)
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
#' @keywords internal
extract_sequences_by_coords <- extract_by_coords

#' @rdname extract_by_name
#' @keywords internal
extract_sequences_by_name <- extract_by_name

#' Predict Variants
#'
#' @description
#' Parses a \code{genome_entity_gd} object -- loaded from a breseq
#' \code{annotated.gd} file -- and returns a tidy data frame of mutations.
#' Annotation tags written by \code{gdtools ANNOTATE} (gene names, codon
#' changes, \code{snp_type}, inactivation status) are read directly from the
#' file; no genome-sequence computation is performed here.
#'
#' For the full molecular-consequence workflow (reference/alternate CDS and
#' protein sequences), pipe the result into \code{\link{annotate_variants}()}.
#'
#' @param gd A \code{genome_entity_gd} object produced by
#'   \code{\link{read_variants}()}.
#' @param ... Reserved for future arguments.
#' @param min_freq Numeric scalar in \eqn{[0, 1]}.  Rows with allele frequency
#'   below this threshold are dropped.  Applied to RA evidence only; structural
#'   variants (frequency fixed at 100\%) are unaffected.  Default: \code{0}
#'   (keep all).
#' @param include_structural Logical; include structural variants (MOB, AMP,
#'   INV, CON) detected by MC/JC evidence.  Default: \code{TRUE}.
#' @param join How to render multi-valued annotation fields (multiple genes,
#'   multiple products).  One of \code{"slash"} (default, \code{" / "}),
#'   \code{"pipe"} (\code{" | "}), or \code{"newline"} (\code{"\n"}).
#'
#' @return The input \code{gd} object with \code{gd$variants_predicted}
#'   populated (a data frame with one row per predicted variant).  Assign the
#'   result back: \code{gd <- predict_variants(gd)}.  Columns:
#' \describe{
#'   \item{\code{evidence}}{Evidence type supporting the call: \code{"RA"},
#'     \code{"MC"}, \code{"JC"}, or a combination such as \code{"MC JC"}.}
#'   \item{\code{type}}{Three-letter mutation type: \code{SNP}, \code{DEL},
#'     \code{INS}, \code{SUB}, \code{AMP}, \code{INV}, \code{MOB},
#'     \code{CON}.}
#'   \item{\code{seq_id}}{Reference contig / chromosome identifier.}
#'   \item{\code{position}}{Formatted genomic position (thousands-separated).
#'     Insertions append the within-position offset as \code{pos:k}.}
#'   \item{\code{mutation}}{Human-readable change label, e.g.\ \code{A->T},
#'     \code{delta 1,234 bp}, \code{+25 bp}, \code{AMP x3}.}
#'   \item{\code{freq}}{Allele frequency as a percentage string, e.g.\
#'     \code{"45.2\%"}.  Structural variants are reported as
#'     \code{"100.0\%"}.}
#'   \item{\code{annotation}}{Position within the gene for coding mutations
#'     (e.g.\ \code{"coding (45/1200 nt)"}); amino-acid change for SNPs
#'     (e.g.\ \code{"K123E (AAA->GAA)"}); intergenic distance otherwise.}
#'   \item{\code{gene}}{Gene symbol or locus tag, with a strand arrow appended
#'     (\code{->} / \code{<-}) for structural variants.}
#'   \item{\code{description}}{Gene product description.}
#'   \item{\code{var_type}}{Molecular consequence class read from breseq
#'     annotation tags.  Values:
#'     \code{synonymous}, \code{nonsynonymous}, \code{nonsense},
#'     \code{noncoding}, \code{pseudogene} (SNPs);
#'     \code{frameshift}, \code{inframe_deletion}, \code{inframe_insertion},
#'     \code{inframe_substitution}, \code{inactivating} (indels and structural
#'     variants);
#'     \code{intergenic} (any mutation outside annotated features).}
#' }
#'
#' @seealso
#' \code{\link{annotate_variants}()} to add reference/alternate sequences;\cr
#' \code{\link{filter_variants}()} to subset by gene, frequency, type, or
#' consequence;\cr
#' \code{\link{read_variants}()} to load the \code{annotated.gd} file.
#'
#' @export
predict_variants <- function(gd, ...) {
  UseMethod("predict_variants")
}

#' Public wrapper: predict variants (opinionated defaults; args pass-through)
#'
#' Calls predict_variants_int(), announces row count, and prints a preview:
#' - If dplyr is installed: prints tibble (contract unaffected).
#' - If dplyr is not installed: prints base data.frame.
#' Any preview/logging errors are swallowed; the table is still returned.
#'
#' @param gd genome_entity_gd
#' @param ... reserved for future options passed to *_int (kept for forward-compat)
#' @param min_freq numeric scalar, keep rows with freq >= min_freq (default 0)
#' @param include_structural logical, include JC/MC/etc. (default TRUE)
#' @param join one of c("slash","pipe","newline") for multi-item fields (default "slash")
#' @return The input \code{gd} object with \code{gd$variants_predicted} set to
#'   the mutation table.  Assign the result back:
#'   \code{gd <- predict_variants(gd)}.  Subsequent calls short-circuit when
#'   the cached table is present.
#' @export
#' @rdname predict_variants
#' @export
predict_variants.genome_entity_gd <- function(gd, ...,
                                               min_freq = 0,
                                               include_structural = TRUE,
                                               join = c("slash", "pipe", "newline")) {
  gd_assert(gd, "gd")
  join <- match.arg(join)

  .default_call <- (min_freq == 0 && isTRUE(include_structural) && join == "slash")

  if (.default_call && !is.null(gd$variants_predicted)) {
    cli::cli_inform("i" = "predict_variants: using cached table ({nrow(gd$variants_predicted)} rows).")
    return(invisible(gd))
  }

  # Build once via internal, regardless of downstream printing path
  tbl <- predict_variants_int(
    gd,
    min_freq = min_freq,
    include_structural = include_structural,
    join = join,
    ...
  )

  # Store result into gd so the caller gets an updated object
  gd$variants_predicted <- tbl

  # Announce + print a friendly preview, but never fail the call if printing goes sideways
  tryCatch(
    {
      n <- NROW(tbl)
      # row count announcement
      cli::cli_inform(c(
        "i" = sprintf("predict_variants: %d row%s.", n, if (n == 1L) "" else "s")
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
  
  # Return the updated gd object (with $variants_predicted populated)
  invisible(gd)
}

#' @export
predict_variants.default <- function(gd, ...) {
  cli::cli_abort(c(
    "{.fn predict_variants} requires a {.cls genome_entity_gd} as its first argument.",
    "i" = "Load your variant calls first with {.fn read_variants}."
  ))
}

#' Annotate Variants with Molecular Sequences and Consequences
#'
#' @description
#' Enriches a \code{\link{predict_variants}()} table with reference and
#' alternate DNA and protein sequences, and carries the \code{var_type}
#' consequence forward as the \code{consequence} column.
#'
#' For coding-region mutations the full reference CDS and its translated
#' protein are extracted from the genome, the mutation is applied in-silico,
#' and the alternate CDS / protein are returned.  For intergenic mutations a
#' flanking DNA window is extracted instead.  Structural variants (MOB, AMP,
#' INV, CON) receive reference sequences only.
#'
#' Consequences (\code{consequence}) are \emph{not} recomputed here -- they are
#' passed through verbatim from the \code{var_type} column produced by
#' \code{\link{predict_variants}()}, which in turn reads breseq's own
#' annotation tags.  This preserves breseq's terminology throughout (e.g.\
#' \code{nonsynonymous}, not \code{missense}).
#'
#' Processing uses parallel workers where available (fork-based on Linux/macOS;
#' PSOCK cluster on Windows), falling back to sequential execution if cluster
#' creation fails.
#'
#' @param gd A \code{genome_entity_gd} object produced by
#'   \code{\link{read_variants}()}.
#' @param pm_tbl A data frame produced by \code{\link{predict_variants}()}.
#'   Must contain columns \code{type}, \code{seq_id}, \code{position},
#'   \code{mutation}, and \code{gene}.  The \code{var_type} column, if present,
#'   is forwarded to \code{consequence} without modification.
#' @param ... Additional arguments passed to the underlying implementation.
#'   Recognised arguments:
#'   \describe{
#'     \item{\code{flank}}{Integer; number of bases to extract on each side of
#'       an intergenic mutation for the DNA window.  Default: \code{50}.}
#'     \item{\code{quiet}}{Logical; suppress progress messages.
#'       Default: \code{FALSE}.}
#'   }
#'
#' @return The input \code{pm_tbl} with the following columns added:
#' \describe{
#'   \item{\code{dna_ref}}{Reference DNA sequence.  Full CDS for coding
#'     mutations; a \code{flank}-bp window centred on the variant for
#'     intergenic mutations.}
#'   \item{\code{dna_alt}}{Alternate DNA sequence with the mutation applied.
#'     \code{NA} for structural variants.}
#'   \item{\code{aa_ref}}{Full reference protein sequence (coding only;
#'     \code{NA} otherwise).}
#'   \item{\code{aa_alt}}{Full alternate protein sequence with the mutation
#'     applied, truncated at the first stop codon (coding only;
#'     \code{NA} otherwise).}
#'   \item{\code{codon_ref}}{Reference codon at the mutation site (SNPs only).}
#'   \item{\code{codon_alt}}{Alternate codon at the mutation site (SNPs only).}
#'   \item{\code{consequence}}{Consequence class, copied from \code{var_type}.
#'     See \code{\link{predict_variants}()} for the full vocabulary.}
#'   \item{\code{region}}{\code{"coding"} or \code{"intergenic"}.}
#'   \item{\code{strand}}{Gene strand (\code{"+"} or \code{"-"}) for coding
#'     mutations; \code{NA} otherwise.}
#'   \item{\code{qc_note}}{Quality-control notes: position offset warnings,
#'     alternative start-codon normalisations, fallback notices.
#'     \code{NA} when no issues were detected.}
#' }
#'
#' @seealso
#' \code{\link{predict_variants}()} to build the input table;\cr
#' \code{\link{filter_variants}()} to subset by consequence or gene;\cr
#' \code{\link{get_reference_dna}()}, \code{\link{get_variant_dna}()},
#' \code{\link{get_reference_aa}()}, \code{\link{get_variant_aa}()} to
#' extract sequences from the result.
#'
#' @examples
#' \dontrun{
#' gd       <- read_variants(ref, "annotated.gd")
#' variants <- predict_variants(gd)
#'
#' # Full workflow
#' ct <- annotate_variants(gd, variants)
#'
#' # Quiet, wider intergenic window
#' ct <- annotate_variants(gd, variants, flank = 200, quiet = TRUE)
#'
#' # Extract sequences for a gene of interest
#' get_reference_aa(ct, gene = "rpoB")
#' get_variant_aa(ct,  gene = "rpoB")
#'
#' # Filter to nonsynonymous only
#' filter_variants(ct, consequence = "nonsynonymous")
#' }
#'
#' @export
annotate_variants <- function(gd, pm_tbl = NULL, ...) {
  UseMethod("annotate_variants")
}

#' @export
annotate_variants.genome_entity_gd <- function(gd, pm_tbl = NULL, ...) {
  if (is.null(pm_tbl)) {
    if (is.null(gd$variants_predicted)) {
      cli::cli_abort(c(
        "{.fn annotate_variants} needs a variant table.",
        "i" = "Run {.code gd <- predict_variants(gd)} first, then call {.code annotate_variants(gd)}."
      ))
    }
    pm_tbl <- gd$variants_predicted
  }
  pm_enrich_consequences_parallel(gd, pm_tbl, ...)
}

#' @export
annotate_variants.default <- function(gd, pm_tbl = NULL, ...) {
  cli::cli_abort(c(
    "{.fn annotate_variants} requires a {.cls genome_entity_gd} as its first argument.",
    "i" = "Load your variant calls first with {.fn read_variants}, then call {.fn predict_variants} to get a variant table."
  ))
}


#' Retrieve the Predicted Variants Table
#'
#' @description
#' Returns the predicted variants table stored in \code{gd$variants_predicted}
#' as a \code{\link[tibble]{tibble}}.  The table is created by
#' \code{\link{predict_variants}()}.
#'
#' @param gd A \code{genome_entity_gd} object.
#' @param ... Currently unused.
#'
#' @return A tibble with one row per predicted variant (see
#'   \code{\link{predict_variants}()} for column descriptions).
#'
#' @seealso
#' \code{\link{predict_variants}()} to build the table;\cr
#' \code{\link{annotate_variants}()} to enrich with molecular consequences.
#'
#' @examples
#' \dontrun{
#' gd <- read_variants("annotated.gd", ref_dir = "reference/")
#' gd <- predict_variants(gd)
#' variants_table <- get_predicted_variants_table(gd)
#' }
#'
#' @export
get_predicted_variants_table <- function(gd, ...) {
  UseMethod("get_predicted_variants_table")
}

#' @rdname get_predicted_variants_table
#' @export
get_predicted_variants_table.genome_entity_gd <- function(gd, ...) {
  if (is.null(gd$variants_predicted)) {
    cli::cli_abort(c(
      "{.fn get_predicted_variants_table} found no variant table in {.arg gd}.",
      "i" = "Run {.code gd <- predict_variants(gd)} first to generate the table."
    ))
  }
  tibble::as_tibble(gd$variants_predicted)
}

#' @export
get_predicted_variants_table.default <- function(gd, ...) {
  cli::cli_abort(c(
    "{.fn get_predicted_variants_table} requires a {.cls genome_entity_gd} as its first argument.",
    "i" = "Load your variant calls with {.fn read_variants}, then run {.code gd <- predict_variants(gd)}."
  ))
}


#' Get Features Overlapping a Region of Interest
#'
#' @description
#' Returns a tidy tibble of genomic features that overlap a specified region,
#' with optional flanking extension. Works on \code{genome_entity} objects and,
#' by S3 dispatch inheritance, on \code{genome_entity_gd} objects as well.
#'
#' @param x A genome_entity (or genome_entity_gd) object
#' @param contig Character; contig/chromosome name (numeric-like strings such
#'   as \code{"1"} are resolved automatically against the GFF seqlevels)
#' @param start Integer; start coordinate (1-based, inclusive)
#' @param end Integer; end coordinate (1-based, inclusive)
#' @param flank Integer; extend the query region by this many bp on each side
#'   (default 0)
#' @param feature_type Character; keep only features of this type, e.g.
#'   \code{"CDS"} (default). Pass \code{NULL} to return all feature types.
#' @param auto_resolve Logical; attempt to reconcile contig to a GFF seqlevel
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

#' @rdname get_roi_features
#' @export
get_roi_features.genome_entity <- function(x, contig, start, end,
                                           flank = 0,
                                           feature_type = "CDS",
                                           auto_resolve = TRUE,
                                           drop_extdb = TRUE,
                                           ...) {
  analyze_roi(
    genome_obj   = x,
    contig       = contig,
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

# -- read_variants --------------------------------------------------------------

#' Read Variant Calls Against a Reference Genome
#'
#' @description
#' Parses a variant file against an existing reference genome and returns a
#' \code{genome_entity_gd} object. The reference genome must already be loaded
#' with \code{read_genome()} before calling \code{read_variants()}.
#'
#' Currently supports breseq annotated genome diff (\code{format = "gd"}). The
#' function is designed to be extensible: future formats (e.g., VCF) will be
#' added via the \code{format} argument without changing the calling convention.
#'
#' @details
#' ## Format requirements
#' - **"gd"**: Must be an *annotated* genome diff (produced by \code{gdtools ANNOTATE}).
#'   Raw breseq \code{output.gd} files are not accepted because they lack the
#'   gene-level annotation needed for reference cross-checking.
#'
#' @param x A \code{genome_entity} object (the reference genome).
#' @param path Character; path to the variant file.
#' @param format Character; variant format. Currently \code{"gd"} (annotated breseq
#'   genome diff). Will be extended to support \code{"vcf"} in a future release.
#' @param strict Logical; if \code{TRUE} (default), stop on structural
#'   inconsistencies; if \code{FALSE}, warn and continue.
#' @param ... Additional arguments passed to the format-specific parser.
#'
#' @return
#' A \code{genome_entity_gd} object. In interactive sessions, also prints a
#' formatted summary to the console:
#' \itemize{
#'   \item \strong{File} -- basename of the variant file
#'   \item \strong{Mutations} -- total count and breakdown by type
#'     (SNP, DEL, INS, SUB, MOB, AMP, CON, INV)
#'   \item \strong{Next} -- three suggested follow-on function calls
#' }
#'
#' @seealso
#' \code{\link{predict_variants}()} to build the tidy mutation table;
#' \code{\link{annotate_variants}()} to add molecular consequences;
#' \code{\link{filter_variants}()} to subset by gene, type, or consequence.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the reference genome first
#' ref <- read_genome("reference.gbk")
#'
#' # Then load variants against it
#' gd <- read_variants(ref, "annotated.gd")
#'
#' # Canonical next steps
#' pm  <- predict_variants(gd)
#' eff <- annotate_variants(gd, pm)
#' filter_variants(gd, pm, gene = "mutL")
#' }
read_variants <- function(x, path, format = "gd", ...) {
  UseMethod("read_variants")
}

#' @rdname read_variants
#' @export
read_variants.genome_entity <- function(x, path, format = "gd",
                                        strict = TRUE, ...) {
  format <- match.arg(format, choices = c("gd"))
  switch(format,
    "gd" = .read_variants_gd(x, path, strict = strict, ...),
    cli::cli_abort("Unsupported variant format: {.val {format}}")
  )
}

#' @export
read_variants.default <- function(x, ...) {
  cli::cli_abort(c(
    "{.fn read_variants} requires a {.cls genome_entity} as its first argument.",
    "i" = "Load your reference genome first with {.fn read_genome}, then call {.fn read_variants}."
  ))
}

.read_variants_gd <- function(entity, path, strict = TRUE, ...) {
  gd <- parse_gd_annotated(gd_path = path, entity = entity, strict = strict, ...)
  .read_variants_report(gd, path)
  gd
}

#' @keywords internal
#' @noRd
.read_variants_report <- function(gd, path) {
  if (!interactive()) return(invisible(NULL))

  w   <- .micromicon_console_width()
  col <- getOption("micromicon.color.code", 39L)

  # -- Mutation counts -------------------------------------------------------

  mut_idx    <- gd$provenance$by_kind$mutation_idx
  mut_events <- gd$events[mut_idx]
  n_mut      <- length(mut_events)
  mut_types  <- vapply(mut_events, `[[`, character(1L), "type")
  type_tbl   <- table(mut_types)

  get_n <- function(key) {
    v <- type_tbl[key]
    if (length(v) == 0L || is.na(v)) 0L else as.integer(v)
  }

  type_order <- c("SNP", "DEL", "INS", "SUB", "MOB", "AMP", "CON", "INV")
  parts <- character()
  for (tp in type_order) {
    n <- get_n(tp)
    if (n > 0L) parts <- c(parts, paste0(formatC(n, format = "d", big.mark = ","), " ", tp))
  }
  n_other <- sum(vapply(setdiff(names(type_tbl), type_order), get_n, integer(1L)))
  if (n_other > 0L) parts <- c(parts, paste0(n_other, " other"))

  mut_str <- paste0(
    formatC(n_mut, format = "d", big.mark = ","),
    if (length(parts) > 0L) paste0("  \u2014  ", paste(parts, collapse = "  \u00b7  ")) else ""
  )

  # -- Print -----------------------------------------------------------------

  lbl <- function(x) sprintf("  %-11s", x)

  .micromicon_rule_title("Variants loaded", w)
  cat(lbl("File"),      basename(path), "\n", sep = "")
  cat(lbl("Mutations"), mut_str, "\n", sep = "")
  cat("\n")

  fn1 <- sprintf("%-34s", "predict_variants(gd)")
  fn2 <- sprintf("%-34s", "annotate_variants(gd, pm)")
  fn3 <- sprintf("%-34s", "filter_variants(gd, pm)")
  cat(lbl("Next"), .micromicon_col256(fn1, col), "build the mutation table\n",        sep = "")
  cat(lbl(""),     .micromicon_col256(fn2, col), "add molecular consequences\n",      sep = "")
  cat(lbl(""),     .micromicon_col256(fn3, col), "subset by gene, type, or consequence\n", sep = "")

  .micromicon_rule_full(w)
  invisible(NULL)
}
