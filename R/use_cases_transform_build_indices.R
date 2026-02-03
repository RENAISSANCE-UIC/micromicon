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
#' @export
#' @examples
#' features <- data.frame(
#'   seqname = c("chr1", "chr1", "chr2"),
#'   locus_tag = c("GENE001", "GENE002", "GENE003"),
#'   gene = c("geneA", "geneB", "geneC"),
#'   stringsAsFactors = FALSE
#' )
#'
#' metadata <- data.frame(
#'   seqname = c("chr1", "chr2"),
#'   length = c(1000, 500),
#'   stringsAsFactors = FALSE
#' )
#'
#' indices <- execute_build_indices(features, metadata)
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

  list(
    seqnames = seqnames_unique,
    locus_tag_index = locus_tag_index,
    gene_index = gene_index
  )
}
