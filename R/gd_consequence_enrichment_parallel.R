#' Parallel Mutation Consequence Enrichment
#'
#' @description
#' Uses parallel processing to enrich mutations across multiple cores.
#' Each gene is processed independently in parallel.
#'
#' @param gd A GenomeData object
#' @param pm_tbl A data.frame from `predict_mutations()`
#' @param flank Integer, flanking bases for intergenic regions (default: 50)
#' @param quiet Logical, suppress messages (default: FALSE)
#' @param use_parallel Logical, enable parallel processing (default: TRUE)
#' @param mc.cores Number of cores to use (default: parallel::detectCores() - 1)
#'
#' @return Enriched data.frame with consequence columns
#' @keywords internal
pm_enrich_consequences_parallel <- function(gd, pm_tbl, flank = 50L, quiet = FALSE,
                                           use_parallel = TRUE,
                                           mc.cores = parallel::detectCores() - 1) {
  gd_assert(gd, "gd")
  stopifnot(is.data.frame(pm_tbl))

  # Check if parallel is available
  if (!requireNamespace("parallel", quietly = TRUE)) {
    cli::cli_warn("parallel package not available, falling back to sequential processing")
    use_parallel <- FALSE
  }

  # Check platform and adjust settings
  is_windows <- .Platform$OS.type == "windows"

  if (use_parallel && is_windows && mc.cores > 1) {
    if (!quiet) {
      cli::cli_warn(c(
        "!" = "Parallel processing not supported on Windows",
        "i" = "Falling back to sequential processing - this may take a while",
        "i" = "Consider running on Linux/macOS for {mc.cores}x speedup"
      ))
    }
    use_parallel <- FALSE
  }

  if (!use_parallel) {
    mc.cores <- 1
  }

  if (!quiet) {
    if (use_parallel) {
      cli::cli_inform(c(
        "i" = "pm_enrich_consequences_parallel(): enriching {nrow(pm_tbl)} mutation{?s} using {mc.cores} core{?s}"
      ))
    } else {
      cli::cli_inform(c(
        "i" = "pm_enrich_consequences_parallel(): enriching {nrow(pm_tbl)} mutation{?s} (sequential mode)"
      ))
    }
  }

  # Initialize output columns
  out <- pm_tbl
  n_rows <- nrow(out)
  out$dna_ref <- rep(NA_character_, n_rows)
  out$dna_alt <- rep(NA_character_, n_rows)
  out$aa_ref <- rep(NA_character_, n_rows)
  out$aa_alt <- rep(NA_character_, n_rows)
  out$codon_ref <- rep(NA_character_, n_rows)
  out$codon_alt <- rep(NA_character_, n_rows)
  out$codon_new <- rep(NA_character_, n_rows)
  out$consequence <- rep(NA_character_, n_rows)
  out$region <- rep(NA_character_, n_rows)
  out$strand <- rep(NA_character_, n_rows)
  out$qc_note <- rep(NA_character_, n_rows)

  out$..row_id <- seq_len(n_rows)

  # Separate intergenic mutations from coding ones (by annotation)
  out$..is_intergenic <- grepl("intergenic", out$annotation, ignore.case = TRUE)

  # Split into coding vs intergenic
  coding_rows <- out[!out$..is_intergenic, , drop = FALSE]
  intergenic_rows <- out[out$..is_intergenic, , drop = FALSE]

  # Clean gene names for coding mutations only
  if (nrow(coding_rows) > 0) {
    coding_rows$..gene_clean <- gsub("\\s*[→←]\\s*$", "", coding_rows$gene)
    coding_rows$..gene_clean <- trimws(coding_rows$..gene_clean)

    # Handle multi-gene annotations (e.g., "geneA|geneB" or "geneA / geneB")
    # Take the first gene in the list
    coding_rows$..gene_clean <- gsub("\\|.*$", "", coding_rows$..gene_clean)  # Remove everything after |
    coding_rows$..gene_clean <- gsub("\\s*/\\s*.*$", "", coding_rows$..gene_clean)  # Remove everything after /
    coding_rows$..gene_clean <- trimws(coding_rows$..gene_clean)
  }

  # Split coding mutations by gene
  gene_groups <- if (nrow(coding_rows) > 0) {
    split(coding_rows, coding_rows$..gene_clean, drop = FALSE)
  } else {
    list()
  }

  # Create shared caches (will be copied to each worker)
  gene_dna_cache <- new.env(parent = emptyenv())
  gene_aa_cache <- new.env(parent = emptyenv())
  gene_meta_cache <- new.env(parent = emptyenv())

  # Process genes (parallel or sequential with progress bar)
  gene_names <- names(gene_groups)

  if (use_parallel) {
    # PARALLEL PROCESSING: Process each gene in parallel
    if (!quiet) {
      cli::cli_inform(c(
        "i" = "Processing {length(gene_names)} gene{?s} in parallel..."
      ))
      start_time <- Sys.time()
    }

    enriched_groups <- parallel::mclapply(
      gene_names,
      function(gene) {
        if (is.na(gene) || !nzchar(gene)) {
          return(gene_groups[[gene]])
        }

        .pf_enrich_gene_batch(
          gd,
          gene_groups[[gene]],
          flank,
          quiet = TRUE,  # Suppress messages in parallel
          gene_dna_cache,
          gene_aa_cache,
          gene_meta_cache
        )
      },
      mc.cores = mc.cores,
      mc.preschedule = FALSE  # Dynamic scheduling for load balancing
    )

    if (!quiet) {
      elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
      cli::cli_inform(c(
        "v" = "Parallel processing completed in {elapsed}s"
      ))
    }
  } else {
    # SEQUENTIAL PROCESSING: With progress bar
    enriched_groups <- vector("list", length(gene_names))

    if (!quiet) {
      pb <- txtProgressBar(min = 0, max = length(gene_names), style = 3)
    }

    for (i in seq_along(gene_names)) {
      gene <- gene_names[i]

      if (is.na(gene) || !nzchar(gene)) {
        enriched_groups[[i]] <- gene_groups[[gene]]
      } else {
        enriched_groups[[i]] <- .pf_enrich_gene_batch(
          gd,
          gene_groups[[gene]],
          flank,
          quiet = TRUE,
          gene_dna_cache,
          gene_aa_cache,
          gene_meta_cache
        )
      }

      if (!quiet) {
        setTxtProgressBar(pb, i)
      }
    }

    if (!quiet) {
      close(pb)
      cat("\n")  # New line after progress bar
    }
  }

  # Combine results
  for (i in seq_along(enriched_groups)) {
    enriched_group <- enriched_groups[[i]]

    # Skip errors from parallel workers (will be caught by intergenic fallback)
    if (inherits(enriched_group, "try-error") || !is.data.frame(enriched_group)) {
      next
    }

    for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                  "codon_ref", "codon_alt", "codon_new", "consequence",
                  "region", "strand", "qc_note")) {
      out[enriched_group$..row_id, col] <- enriched_group[[col]]
    }
  }

  # Handle mutations that weren't enriched (intergenic or failed gene enrichment)
  # These will have region = NA still
  unenriched_idx <- which(is.na(out$region))
  if (length(unenriched_idx) > 0) {
    for (i in unenriched_idx) {
      out[i, ] <- .pf_enrich_intergenic(gd, out[i, , drop = FALSE], flank, quiet)
    }
  }

  # Clean up
  out$..row_id <- NULL
  out$..gene_clean <- NULL
  out$..is_intergenic <- NULL

  if (!quiet) {
    n_coding <- sum(out$region == "coding", na.rm = TRUE)
    n_intergenic <- sum(out$region == "intergenic", na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Enriched {n_coding} coding and {n_intergenic} intergenic mutation{?s}"
    ))
  }

  out
}
