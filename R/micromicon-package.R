#' micromicon: Genomic Analysis for Microbial Genomes
#'
#' @description
#' Tools for loading, querying, and analyzing microbial reference genomes and
#' variant calls (breseq genome diffs). Built around two core S3 objects:
#' \code{genome_entity} (reference genome) and \code{genome_entity_gd}
#' (genome + variant calls).
#'
#' ## Typical workflow
#'
#' ```r
#' # 1. Load a reference genome
#' ref <- read_genome("reference.gbk")
#'
#' # 2. Explore features and sequences
#' get_features(ref, type = "CDS")
#' get_gene_dna(ref, "dnaA")
#'
#' # 3. Load variant calls
#' gd <- read_variants(ref, "annotated.gd")
#'
#' # 4. Analyse mutations
#' pm  <- predict_mutations(gd)
#' eff <- compute_effects(gd, pm)
#' ```
#'
#' Call \code{\link{micromicon_functions}()} at the console for a full
#' categorised function listing.
#'
#' @keywords internal
"_PACKAGE"

# ── micromicon_functions ───────────────────────────────────────────────────────

#' List All User-Facing Functions
#'
#' Prints a categorised, formatted overview of every exported function in the
#' micromicon package. Useful for discovering what is available without leaving
#' the console.
#'
#' @return Called for its side effect (printing). Returns \code{invisible(NULL)}.
#' @export
#'
#' @examples
#' micromicon_functions()
micromicon_functions <- function() {

  cli::cli_h1("micromicon")

  cli::cli_h2("Reading data")
  cli::cli_ul(c(
    "{.fn read_genome}    \u2014 load a reference genome (.gbk, .gff3 + .fasta)",
    "{.fn read_variants}  \u2014 load variant calls against a reference (GD; VCF coming)"
  ))

  cli::cli_h2("Genomic features")
  cli::cli_ul(c(
    "{.fn get_features}      \u2014 all annotated features as a data.frame",
    "{.fn search_features}   \u2014 filter features by type, name pattern, or coordinates",
    "{.fn get_roi_features}  \u2014 features overlapping a genomic region",
    "{.fn gene_info}         \u2014 full details for a single gene"
  ))

  cli::cli_h2("Sequences")
  cli::cli_ul(c(
    "{.fn get_gene_dna}  \u2014 CDS nucleotide sequence for a gene",
    "{.fn get_gene_aa}   \u2014 amino acid sequence for a gene",
    "{.fn get_roi_dna}   \u2014 nucleotide sequence for any genomic region",
    "{.fn sequences}     \u2014 all contig/chromosome sequences"
  ))

  cli::cli_h2("Coordinate mapping")
  cli::cli_ul(c(
    "{.fn map_cds_to_genomic}   \u2014 CDS position \u2192 genomic position",
    "{.fn map_genomic_to_cds}   \u2014 genomic position \u2192 CDS position",
    "{.fn extract_by_name}      \u2014 extract sequence and features by gene name",
    "{.fn extract_by_coords}    \u2014 extract sequence and features by coordinates"
  ))

  cli::cli_h2("Genome metadata")
  cli::cli_ul(c(
    "{.fn get_genome_metadata}  \u2014 organism name, length, topology, etc.",
    "{.fn get_genome_seqnames}  \u2014 contig / chromosome names"
  ))

  cli::cli_h2("Variants & mutations")
  cli::cli_ul(c(
    "{.fn predict_mutations}        \u2014 tidy mutation table from a genome diff",
    "{.fn compute_effects}          \u2014 molecular consequences (aa change, consequence class)",
    "{.fn validate_variant_in_gene} \u2014 check whether a variant falls within a gene"
  ))

  cli::cli_h2("Export")
  cli::cli_ul(c(
    "{.fn export_genome}  \u2014 write genome as GFF3, FASTA, or GenBank"
  ))

  cli::cli_inform(c(
    "i" = "Use {.code ?<function>} for full documentation, e.g. {.code ?read_genome}.",
    "i" = "Use {.code ?micromicon} for the package overview."
  ))

  invisible(NULL)
}
