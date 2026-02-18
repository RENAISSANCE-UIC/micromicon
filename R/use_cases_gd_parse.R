#' Parse a breseq annotated.gd file with type-aware dispatch
#'
#' Parses a breseq annotated.gd file (produced by gdtools ANNOTATE) and creates
#' a genome_entity_gd object. The parser uses type-aware dispatch to handle
#' different mutation types (SNP, DEL, INS, SUB, MOB, AMP, CON, INV), evidence
#' types (RA, MC, JC, UN), and validation types (TSEQ, PFLP, etc.). Events are
#' automatically binned by kind (mutation/evidence/validation) for efficient filtering.
#'
#' @param gd_path Path to the annotated.gd file (not raw output.gd)
#' @param entity A genome_entity object (reference genome; fields will be hoisted)
#' @param strict Logical; if TRUE, stop on inconsistencies; if FALSE, warn and continue
#' @param fasta_path Optional path to FASTA file for provenance checksums
#' @param gff3_path Optional path to GFF3 file for provenance checksums
#' @param gbk_path Optional path to GenBank file for provenance checksums
#' @return A genome_entity_gd object with parsed events binned by kind in
#'   \code{provenance$by_kind} (mutation_idx, evidence_idx, validation_idx)
#' @export
#' @examples
#' \dontrun{
#'   ref <- read_genome("reference.gbk")
#'   gd <- parse_gd_annotated("annotated.gd", ref, strict = TRUE)
#'   length(gd$events)
#'   gd$provenance$by_kind$mutation_idx  # indices of mutation events
#' }
parse_gd_annotated <- function(gd_path, entity, strict = TRUE,
                               fasta_path = NULL, gff3_path = NULL,
                               gbk_path = NULL) {
  if (!file.exists(gd_path)) cli::cli_abort("File does not exist: {gd_path}")
  stopifnot(inherits(entity, "genome_entity"))

  contig_lengths <- setNames(as.numeric(entity$metadata$length_bp), entity$metadata$seqname)
  ref_manifest   <- reference_manifest_from_genome_entity(entity, fasta_path, gff3_path, gbk_path)

  raw <- parse_gd_raw(gd_path, contig_lengths, strict)

  gd_checksum <- gd_digest(gd_path, file = TRUE)

  obj <- new_genome_entity_gd(
    header    = raw$header,
    events    = raw$events,
    file      = list(path = normalizePath(gd_path), checksum = gd_checksum),
    entity    = entity,
    reference = ref_manifest,
    strict    = strict
  )

  obj$provenance$by_kind <- list(
    mutation_idx   = which(raw$kind_vec == "mutation"),
    evidence_idx   = which(raw$kind_vec == "evidence"),
    validation_idx = which(raw$kind_vec == "validation")
  )

  validate_genome_entity_gd(obj, strict = strict)
}
