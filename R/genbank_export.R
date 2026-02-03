#' Export Functions for Genome Data
#'
#' @description
#' Functions for exporting genome data to various formats (FASTA, GFF3).
#' In Clean Architecture terms, these handle output formatting and file writing.
#'
#' Note: In Phase 6, these functions remain in R/genbank_export.R and are used by
#' R/controllers/export_controller.R. Some of this logic could be moved to use cases
#' and presenters for full Clean Architecture compliance in future phases.
#'
#' @name export
#' @keywords internal
NULL

#' Write sequences to FASTA format
#'
#' @description
#' Writes genome sequences to FASTA format. Accepts both genome_entity objects
#' and legacy GenBank list format for backward compatibility.
#'
#' @param x A genome_entity object or list of GenBank records
#' @param file Output FASTA file path
#' @param wrap_width Integer; wrap sequences at this width (default 80)
#' @param auto_label Logical; if TRUE (default), infer chromosome vs plasmids
#' @param include_length Logical; append length annotation to header (default TRUE)
#'
#' @return Invisibly returns output path
#' @export
write_gbk_fasta <- function(x, file, wrap_width = 80L,
                            auto_label = TRUE, include_length = TRUE) {
  # Handle both genome_entity and legacy formats
  if (inherits(x, "genome_entity")) {
    gbk_list <- entity_to_gbk_list(x)
  } else {
    gbk_list <- x
  }

  .write_gbk_fasta_internal(gbk_list, file, wrap_width, auto_label, include_length)
}

# Internal FASTA writer
.write_gbk_fasta_internal <- function(gbk_list, file, wrap_width = 80L,
                            auto_label = TRUE, include_length = TRUE) {
  if (length(gbk_list) == 0L) {
    cli::cli_warn("No records to write.")
    return(invisible(NULL))
  }
  
  # Helper: wrap sequence at specified width
  wrap_seq <- function(seq_str, width) {
    if (is.na(seq_str) || !nzchar(seq_str)) return("")
    if (width <= 0L) return(seq_str)
    pos <- seq(1L, nchar(seq_str), by = width)
    paste(substring(seq_str, pos, pos + width - 1L), collapse = "\n")
  }
  
  # Compute actual sequence lengths
  seq_lengths <- vapply(gbk_list, function(r) {
    seq <- r$sequence
    if (is.na(seq) || !nzchar(seq)) return(0L)
    nchar(seq)
  }, integer(1L))
  
  # Sort by length (descending) if auto_label enabled
  if (auto_label) {
    ord <- order(seq_lengths, decreasing = TRUE)
    gbk_list <- gbk_list[ord]
    seq_lengths <- seq_lengths[ord]
  }
  
  # Filter out records without sequence
  has_seq <- seq_lengths > 0L
  if (!any(has_seq)) {
    cli::cli_warn("No valid sequences to write.")
    return(invisible(NULL))
  }
  gbk_list <- gbk_list[has_seq]
  seq_lengths <- seq_lengths[has_seq]
  
  # Build FASTA entries
  entries <- character(length(gbk_list))
  plasmid_counter <- 1L
  
  for (i in seq_along(gbk_list)) {
    r <- gbk_list[[i]]
    meta <- r$metadata
    seq <- r$sequence
    
    # Determine replicon label
    if (auto_label) {
      if (i == 1L) {
        # First (longest) = chromosome
        # Check if definition mentions "chromosome"
        is_chr <- !is.na(meta$definition) && 
          grepl("chromosome", meta$definition, ignore.case = TRUE)
        replicon_label <- if (is_chr) "Chromosome" else "Chromosome"
      } else {
        # Subsequent = plasmids
        replicon_label <- sprintf("Plasmid_%d_(presumptive)", plasmid_counter)
        plasmid_counter <- plasmid_counter + 1L
      }
    } else {
      # Use existing metadata if present
      if (!is.na(meta$replicon_type) && !is.na(meta$replicon_name)) {
        replicon_label <- paste0(meta$replicon_type, "_", meta$replicon_name)
      } else if (!is.na(meta$accession)) {
        replicon_label <- meta$accession
      } else if (!is.na(meta$locus)) {
        replicon_label <- meta$locus
      } else {
        replicon_label <- sprintf("contig_%d", i)
      }
    }
    
    # Build header parts
    hdr_parts <- replicon_label
    if (!is.na(meta$organism)) {
      hdr_parts <- paste(hdr_parts, meta$organism)
    }
    if (!is.na(meta$topology)) {
      hdr_parts <- paste(hdr_parts, sprintf("[topology=%s]", meta$topology))
    }
    if (include_length) {
      hdr_parts <- paste(hdr_parts, sprintf("[length=%d]", seq_lengths[i]))
    }
    
    entries[i] <- paste0(">", hdr_parts, "\n", wrap_seq(seq, wrap_width))
  }
  
  writeLines(entries, file)
  cli::cli_inform("Wrote {length(entries)} sequence{?s} to {.path {file}}")
  invisible(file)
}


# -------- GFF3 export ====

#' Write features to GFF3 format
#'
#' @description
#' Writes genomic features to GFF3 format. Accepts both genome_entity objects
#' and legacy GenBank list format for backward compatibility.
#'
#' @param x A genome_entity object or list of GenBank records
#' @param file Output GFF3 file path
#' @param include_sequence_region Logical; include ##sequence-region directives (default TRUE)
#' @param include_types Character vector; feature types to include (default: all)
#'
#' @details
#' GFF3 spec: coordinates are 1-based inclusive. Implements NCBI/INSDC conventions
#' for feature hierarchies (gene -> mRNA/CDS), ID/Parent relationships, and
#' standard attributes (gbkey, gene_biotype, etc.).
#'
#' @return Invisibly returns output path
#' @export
write_gbk_gff3 <- function(x, file, include_sequence_region = TRUE,
                           include_types = NULL) {
  # Handle both genome_entity and legacy formats
  if (inherits(x, "genome_entity")) {
    gbk_list <- entity_to_gbk_list(x)
  } else {
    gbk_list <- x
  }

  .write_gbk_gff3_internal(gbk_list, file, include_sequence_region, include_types)
}

# Internal GFF3 writer
.write_gbk_gff3_internal <- function(gbk_list, file, include_sequence_region = TRUE,
                           include_types = NULL) {
  if (length(gbk_list) == 0L) {
    cli::cli_warn("No records to convert.")
    return(invisible(NULL))
  }
  
  # Helper: URL-encode special chars in GFF3 attributes
  gff_escape <- function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    x <- gsub("%", "%25", x, fixed = TRUE)
    x <- gsub(";", "%3B", x, fixed = TRUE)
    x <- gsub("=", "%3D", x, fixed = TRUE)
    x <- gsub("&", "%26", x, fixed = TRUE)
    x <- gsub(",", "%2C", x, fixed = TRUE)
    x <- gsub("\t", "%09", x, fixed = TRUE)
    x <- gsub("\n", "%0A", x, fixed = TRUE)
    x
  }
  
  # Helper: format GFF3 attributes (order matters for readability)
  format_attributes <- function(attrs_list) {
    # Standard order: ID, Name, Parent, Dbxref, Ontology_term, Note, then alphabetical
    priority <- c("ID", "Name", "Parent", "Dbxref", "Ontology_term", "Note")
    keys <- names(attrs_list)
    # Remove empty values
    attrs_list <- attrs_list[!vapply(attrs_list, function(x) is.na(x) || !nzchar(x), logical(1L))]
    if (length(attrs_list) == 0L) return(".")
    
    keys <- names(attrs_list)
    ordered_keys <- c(
      intersect(priority, keys),
      setdiff(keys, priority)
    )
    
    pairs <- vapply(ordered_keys, function(k) {
      v <- attrs_list[[k]]
      paste0(k, "=", gff_escape(v))
    }, character(1L))
    
    paste(pairs, collapse = ";")
  }
  
  # Helper: infer source column from qualifiers
  infer_source <- function(qualifiers) {
    # Check inference qualifier for method
    inf <- qualifiers[["inference"]]
    if (!is.null(inf) && !is.na(inf[[1]])) {
      # e.g., "COORDINATES: similar to AA sequence:RefSeq:NP_416675.1"
      if (grepl("Protein Homology", inf[[1]], ignore.case = TRUE)) return("Protein Homology")
      if (grepl("ab initio", inf[[1]], ignore.case = TRUE)) return("Ab initio prediction")
      if (grepl("RefSeq", inf[[1]], ignore.case = TRUE)) return("RefSeq")
    }
    "."  # unknown
  }
  
  # Helper: compute CDS phase
  compute_phase <- function(ranges_df, strand) {
    if (nrow(ranges_df) == 0L) return(".")
    # Phase for first segment is 0; subsequent depend on cumulative length
    phases <- rep("0", nrow(ranges_df))
    if (nrow(ranges_df) > 1L) {
      cumlen <- 0L
      for (i in seq_len(nrow(ranges_df))) {
        seg_len <- ranges_df$end[i] - ranges_df$start[i] + 1L
        phases[i] <- as.character(cumlen %% 3L)
        cumlen <- cumlen + seg_len
      }
    }
    phases
  }
  
  # Collect all output lines
  gff_lines <- list()
  seq_region_lines <- list()
  
  # Process each record
  for (rec_idx in seq_along(gbk_list)) {
    r <- gbk_list[[rec_idx]]
    meta <- r$metadata
    feats <- r$features
    
    # Seqid: use locus (the contig/chromosome identifier in LOCUS line)
    seqid <- if (!is.na(meta$locus) && nzchar(meta$locus)) {
      meta$locus
    } else if (!is.na(meta$accession)) {
      meta$accession
    } else {
      paste0("unknown_", rec_idx)
    }
    
    # Sequence length (from actual sequence)
    seq_len <- if (!is.na(r$sequence) && nzchar(r$sequence)) {
      nchar(r$sequence)
    } else {
      meta$length_bp %||% NA_integer_
    }
    
    # Add ##sequence-region directive
    if (include_sequence_region && !is.na(seq_len)) {
      seq_region_lines[[length(seq_region_lines) + 1L]] <- 
        sprintf("##sequence-region %s 1 %d", seqid, seq_len)
    }
    
    # Add top-level "region" feature (from "source" feature in GenBank)
    if (nrow(feats) > 0L) {
      src_idx <- which(feats$type == "source")
      if (length(src_idx) > 0L) {
        src_feat <- feats[src_idx[1], ]
        src_quals <- src_feat$qualifiers[[1]]
        
        # Extract taxon from db_xref
        taxon <- NA_character_
        if (!is.null(src_quals$db_xref)) {
          dbx <- src_quals$db_xref
          taxon_match <- grep("^taxon:", dbx, value = TRUE)
          if (length(taxon_match)) {
            taxon <- sub("^taxon:", "", taxon_match[1])
          }
        }
        
        # Region attributes
        region_attrs <- list(
          ID = sprintf("%s:1..%d", seqid, seq_len),
          Dbxref = if (!is.na(taxon)) paste0("taxon:", taxon) else NA_character_,
          Name = meta$definition %||% "ANONYMOUS",
          gbkey = "Src",
          genome = if (rec_idx == 1L) "chromosome" else "plasmid",
          mol_type = meta$mol_type %||% "genomic DNA"
        )
        
        region_line <- paste(
          seqid, "Local", "region", 1, seq_len, ".", "+", ".",
          format_attributes(region_attrs),
          sep = "\t"
        )
        gff_lines[[length(gff_lines) + 1L]] <- region_line
      }
    }
    
    # Filter features by type if requested
    if (!is.null(include_types)) {
      feats <- feats[feats$type %in% include_types, , drop = FALSE]
    }
    
    if (!nrow(feats)) next
    
    # Track gene IDs for Parent relationships
    gene_map <- list()  # locus_tag -> gene ID
    
    # First pass: create genes
    for (i in seq_len(nrow(feats))) {
      feat <- feats[i, ]
      if (feat$type != "gene") next
      
      quals <- feat$qualifiers[[1]]
      locus_tag <- feat$locus_tag
      if (is.na(locus_tag)) next
      
      gene_id <- paste0("gene-", locus_tag)
      gene_map[[locus_tag]] <- gene_id
      
      strand_char <- switch(as.character(feat$strand), "1" = "+", "-1" = "-", ".")
      
      gene_attrs <- list(
        ID = gene_id,
        Name = feat$gene %||% locus_tag,
        gbkey = "Gene",
        gene_biotype = "protein_coding",
        locus_tag = locus_tag
      )
      
      if (!is.na(feat$gene)) gene_attrs$gene <- feat$gene
      
      # Check for partial annotations
      if (!is.null(quals$partial)) {
        gene_attrs$partial <- "true"
      }
      
      gene_line <- paste(
        seqid, ".", "gene", feat$start, feat$end, ".", strand_char, ".",
        format_attributes(gene_attrs),
        sep = "\t"
      )
      gff_lines[[length(gff_lines) + 1L]] <- gene_line
    }
    
    # Second pass: create CDS and other features
    for (i in seq_len(nrow(feats))) {
      feat <- feats[i, ]
      feat_type <- feat$type
      
      # Skip source (already handled) and gene (handled above)
      if (feat_type %in% c("source", "gene")) next
      
      quals <- feat$qualifiers[[1]]
      ranges_df <- feat$ranges[[1]]
      locus_tag <- feat$locus_tag
      
      strand_char <- switch(as.character(feat$strand), "1" = "+", "-1" = "-", ".")
      source_col <- infer_source(quals)
      
      # Determine ID and Parent
      feat_id <- if (feat_type == "CDS" && !is.na(feat$protein_id)) {
        paste0("cds-", feat$protein_id)
      } else if (!is.na(locus_tag)) {
        paste0(tolower(feat_type), "-", locus_tag)
      } else {
        paste0(seqid, "_", feat_type, "_", i)
      }
      
      parent_id <- if (!is.na(locus_tag) && !is.null(gene_map[[locus_tag]])) {
        gene_map[[locus_tag]]
      } else {
        NA_character_
      }
      
      # Build attributes
      attrs <- list(
        ID = feat_id,
        Name = feat$protein_id %||% feat$product %||% feat$gene %||% NA_character_,
        Parent = parent_id,
        gbkey = toupper(feat_type)
      )
      
      # Add common qualifiers
      if (!is.na(locus_tag)) attrs$locus_tag <- locus_tag
      if (!is.na(feat$gene)) attrs$gene <- feat$gene
      if (!is.na(feat$product)) attrs$product <- feat$product
      if (!is.na(feat$protein_id)) attrs$protein_id <- feat$protein_id
      
      # Add inference
      if (!is.null(quals$inference)) {
        attrs$inference <- quals$inference[[1]]
      }
      
      # Parse and add GO terms
      if (!is.null(quals$note)) {
        note <- quals$note[[1]]
        # Extract GO terms from note
        go_funcs <- gregexpr("GO:[0-9]{7}", note, perl = TRUE)
        if (go_funcs[[1]][1] != -1L) {
          go_terms <- regmatches(note, go_funcs)[[1]]
          if (length(go_terms)) {
            attrs$Ontology_term <- paste(go_terms, collapse = ",")
          }
        }
      }
      
      # Add partial flag
      if (!is.null(quals$partial)) {
        attrs$partial <- "true"
      }
      
      # Add transl_table for CDS
      if (feat_type == "CDS" && !is.null(quals$transl_table)) {
        attrs$transl_table <- quals$transl_table[[1]]
      }
      
      # Compute phase for CDS
      phases <- if (feat_type == "CDS") {
        compute_phase(ranges_df, feat$strand)
      } else {
        rep(".", nrow(ranges_df))
      }
      
      # Write one line per segment (for join/order)
      for (seg_idx in seq_len(nrow(ranges_df))) {
        seg <- ranges_df[seg_idx, ]
        if (is.na(seg$start) || is.na(seg$end)) next
        
        phase_char <- phases[seg_idx]
        
        # Attributes only on first segment
        attrs_str <- if (seg_idx == 1L) {
          format_attributes(attrs)
        } else {
          format_attributes(list(ID = paste0(feat_id, ".seg", seg_idx), Parent = feat_id))
        }
        
        line <- paste(
          seqid, source_col, feat_type, seg$start, seg$end, ".", strand_char, phase_char, attrs_str,
          sep = "\t"
        )
        gff_lines[[length(gff_lines) + 1L]] <- line
      }
    }
  }
  
  if (length(gff_lines) == 0L) {
    cli::cli_warn("No features to write to GFF3.")
    return(invisible(NULL))
  }
  
  # Write output
  header <- "##gff-version 3"
  all_lines <- c(header, unlist(seq_region_lines), unlist(gff_lines))
  writeLines(all_lines, file)
  cli::cli_inform("Wrote {length(gff_lines)} feature line{?s} to {.path {file}}")
  invisible(file)
}

# -------- LOGGING ====
logger <- function(level = "INFO") {
  function(msg, ...) {
    if (Sys.getenv("LOG_LEVEL", "INFO") == "DEBUG" || level != "DEBUG") {
      cli::cli_inform(paste0("[", level, "] ", msg), ...)
    }
  }
}
# log_info <- logger("INFO")
# log_warn <- logger("WARN")
# log_debug <- logger("DEBUG")


# -------- Unit Tests (DISABLED - Move to tests/testthat/) ====
# Note: These tests were embedded in the source file and executed at load time.
# They should be moved to proper test files in tests/testthat/.
# Commented out to avoid package load errors.
#
# ex <- '
# LOCUS       TEST0001       1200 bp    DNA     circular  BCT  01-JAN-2020
# DEFINITION  Example plasmid with two CDS features and a joined gene.
# ACCESSION   TEST0001
# VERSION     TEST0001.1  GI:999999999
# KEYWORDS    .
# SOURCE      Synthetic construct
#   ORGANISM  Artificial sequence
#             other sequences; synthetic sequences.
# REFERENCE   1  (bases 1 to 1200)
# FEATURES             Location/Qualifiers
#      source          1..1200
#                      /organism="Artificial sequence"
#                      /plasmid="pExample-1"
#      gene            complement(join(100..200,300..350))
#                      /gene="fooA"
#                      /locus_tag="LTG_0001"
#      CDS             complement(join(100..200,300..350))
#                      /gene="fooA"
#                      /product="fused protein alpha"
#                      /protein_id="PROT_0001"
#                      /translation="MSTNPKPQRKTK"
#      CDS             700..>950
#                      /gene="barB"
#                      /product="hypothetical protein"
#                      /protein_id="PROT_0002"
# ORIGIN
#         1 acgatcgatc gatcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga
#        61 tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga
# //
# '
# tf <- tempfile(fileext = ".gbk")
# writeLines(ex, tf)
# gb <- read_gbk(tf)
#
# testthat::test_that("FASTA export works", {
#   gbk <- read_gbk(tf)
#   tmp <- tempfile(fileext = ".fna")
#   write_gbk_fasta(gbk, tmp)
#   lines <- readLines(tmp)
#   testthat::expect_true(any(grepl("^>", lines)))  # has header
#   testthat::expect_true(any(grepl("^[ACGTN]+$", lines)))  # has sequence
# })
#
# testthat::test_that("GFF3 export works", {
#   gbk <- read_gbk(tf)
#   tmp <- tempfile(fileext = ".gff3")
#   write_gbk_gff3(gbk, tmp)
#   lines <- readLines(tmp)
#   testthat::expect_equal(lines[1], "##gff-version 3")
#   testthat::expect_true(any(grepl("CDS", lines)))
# })

