# --- object constructor + validator 

#' Internal constructor: genome_entity_gd (with hoisted fields)
#'
#' Builds a genome_entity_gd that is shape-compatible with genome_entity by
#' hoisting core fields to the top level for backward compatibility.
#'
#' @keywords internal
new_genome_entity_gd <- function(header, events, file,
                                 entity, reference, strict = TRUE) {
  stopifnot(is.list(header), is.list(events), is.list(file),
            inherits(entity, "genome_entity"))

  # ---- Hoist core genome_entity fields for backward compatibility
  sequences <- entity$sequences %||% list()
  features  <- entity$features  %||% NULL
  metadata  <- entity$metadata  %||% NULL
  indices   <- entity$indices   %||% list()

  # ---- Rectify event schema (fail-fast if missing required fields) ----
  req <- c(
    "type","id","rank","seq_id","position","col6","col7","col8",
    "contig","tags","raw_line","hash",
    "is_evidence","evidence_class",
    "snp_alt_base","snp_ref_base","del_size","ins_seq","ins_size",
    "sub_size","sub_new_seq","mob_repeat_name","mob_strand",
    "mob_duplication_size","ev_frequency","ev_quality",
    "ev_ref_cov_1","ev_ref_cov_2","ev_new_cov_1","ev_new_cov_2",
    "ev_tot_cov_1","ev_tot_cov_2","ev_alignment_overlap",
    "ev_cov_minus","ev_cov_plus","ev_pos_start","ev_pos_end"
  )
  events <- lapply(events, function(ev) {
    missing <- setdiff(req, names(ev))
    if (length(missing)) {
      stop("Event missing fields: ", paste(missing, collapse = ", "))
    }
    ev
  })

  # ---- Provenance scaffold
  provenance <- list(
    gd_file   = file,
    reference_checksums = reference$checksums %||% list(),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    host       = Sys.info()[["nodename"]],
    validation = list(status = "unvalidated", messages = character())
  )

  # ---- Assemble object
  x <- list(
    # Hoisted parent fields (back-compat surface)
    sequences  = sequences,
    features   = features,
    metadata   = metadata,
    indices    = indices,
    # GD extensions
    header     = header,
    events     = events,
    reference  = reference,
    provenance = provenance,
    strict     = isTRUE(strict)
  )

  class(x) <- c("genome_entity_gd", "genome_entity")

  x
}


validate_genome_entity_gd <- function(x, strict = x$strict) {
  stopifnot(inherits(x, "genome_entity_gd"))
  msgs <- character(); status <- "ok"
  
  # Identity invariant: hashes must be unique
  hashes <- vapply(x$events, `[[`, character(1), "hash")
  if (anyDuplicated(hashes)) {
    msgs   <- c(msgs, "Duplicate event hashes detected.")
    status <- "error"
  }
  
  contig_lengths <- setNames(
    as.numeric(x$metadata$length_bp %||% numeric(0)),
    x$metadata$seqname %||% character(0)
  )
  
  for (ev in x$events) {
    contig <- ev$contig
    pos    <- ev$position
    
    # --- NEW: skip validation when contig is missing/unknown ---
    if (is.null(contig) || length(contig) == 0L || is.na(contig) || !nzchar(contig)) next
    
    if (!(contig %in% names(contig_lengths))) {
      msgs   <- c(msgs, sprintf("Contig not present in reference: %s", contig))
      status <- if (strict) "error" else "warn"
      next
    }
    
    if (!is.null(pos) && !is.na(pos)) {
      ln <- contig_lengths[[contig]]
      if (pos < 1L || pos > ln) {
        msgs <- c(msgs, sprintf("Position out of bounds for %s: %d not in [1..%d]", contig, pos, ln))
        status <- if (strict) "error" else "warn"
      }
    }
  }
  
  if (identical(status, "error") && strict) {
    stop(paste(c("Validation failed (strict mode).", msgs), collapse = "\n- "))
  } else if (identical(status, "warn") ||
             (identical(status, "error") && !strict)) {
    cli::cli_warn(paste(c("Validation warnings:", msgs), collapse = "\n- "))
  }
  
  x$provenance$validation$status   <- if (strict && identical(status, "ok")) "ok_strict" else status
  x$provenance$validation$messages <- unique(c(x$provenance$validation$messages, msgs))
  x
}

