#' Build Indices for Fast Lookup
#'
#' @description
#' Use case that builds indices for fast feature lookup. This is pure business logic
#' that operates on data frames without any I/O or framework dependencies.
#'
#' Indices built:
#' - seqnames: Unique sequence names
#' - locus_tag_index: Map locus_tag to feature row number
#' - gene_index: Map gene name to feature row number(s)
#'
#' @param features_df data.frame with features (must have seqname column)
#' @param metadata_df data.frame with sequence metadata (must have seqname column)
#'
#' @return List with indices components
#' @keywords internal
execute_build_indices <- function(features_df, metadata_df) {
  # Extract unique seqnames from metadata
  seqnames_unique <- unique(metadata_df$seqname)

  # Build locus_tag index
  locus_tag_index <- integer()
  if (nrow(features_df) > 0 && "locus_tag" %in% names(features_df)) {
    has_lt <- !is.na(features_df$locus_tag) & nzchar(features_df$locus_tag)
    if (any(has_lt)) {
      locus_tag_index <- setNames(which(has_lt), features_df$locus_tag[has_lt])
    }
  }

  # Build gene index (multiple features can have same gene name)
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

  # Build cds_hash: O(1) environment mapping any gene identifier to its feature
  # row index. CDS features take priority over other types; first CDS match wins
  # when multiple CDS share the same identifier.
  cds_hash <- new.env(hash = TRUE, parent = emptyenv())
  if (nrow(features_df) > 0) {
    id_cols <- intersect(
      c("locus_tag", "gene", "Name", "Alias", "ID"),
      names(features_df)
    )
    is_cds <- if ("type" %in% names(features_df)) {
      !is.na(features_df$type) & features_df$type == "CDS"
    } else {
      rep(FALSE, nrow(features_df))
    }

    for (col in id_cols) {
      vals <- features_df[[col]]
      for (i in seq_len(nrow(features_df))) {
        val <- vals[[i]]
        # GFF3 list-columns (e.g. Alias) may be length-0 or multi-value vectors
        if (is.null(val) || length(val) == 0) next
        for (v in val) {
          if (is.na(v) || !nzchar(v)) next
          if (!exists(v, envir = cds_hash, inherits = FALSE)) {
            assign(v, i, envir = cds_hash)
          } else if (is_cds[[i]]) {
            existing_row <- get(v, envir = cds_hash, inherits = FALSE)
            if (!is_cds[[existing_row]]) {
              assign(v, i, envir = cds_hash)  # upgrade non-CDS to CDS
            }
            # first CDS match wins; do not overwrite CDS with later CDS
          }
        }
      }
    }
  }

  list(
    seqnames = seqnames_unique,
    locus_tag_index = locus_tag_index,
    gene_index = gene_index,
    cds_hash = cds_hash
  )
}
