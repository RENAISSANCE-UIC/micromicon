# --- object constructor + validator -------------------------------------------
#' @keywords internal
new_genome_entity_gd <- function(header, events, file, 
                                 entity, reference, strict = TRUE) {
  stopifnot(is.list(header), is.list(events), is.list(file), inherits(entity, "genome_entity"))
  
  # Required fields in every event (rectangular core + minimal flags)
  req <- c(
    "type","id","rank","seq_id","position","col6","col7","col8",
    "contig","tags","raw_line","hash",
    # at least these evidence flags exist (values may be NA)
    "is_evidence","evidence_class",
    # and descriptive mutation fields exist (values may be NA)
    "snp_alt_base","snp_ref_base","del_size","ins_seq","ins_size",
    "sub_size","sub_new_seq","mob_repeat_name","mob_strand","mob_duplication_size",
    # evidence descriptive columns (values may be NA)
    "ev_frequency","ev_quality",
    "ev_ref_cov_1","ev_ref_cov_2","ev_new_cov_1","ev_new_cov_2","ev_tot_cov_1","ev_tot_cov_2",
    "ev_alignment_overlap","ev_cov_minus","ev_cov_plus","ev_pos_start","ev_pos_end"
  )
  
  events <- lapply(events, function(ev) {
    missing <- setdiff(req, names(ev))
    if (length(missing)) stop("Event missing fields: ", paste(missing, collapse = ", "))
    ev
  })
  
  provenance <- list(
    gd_file   = file,
    reference_checksums = reference$checksums %||% list(),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    host       = Sys.info()[["nodename"]],
    validation = list(status = "unvalidated", messages = character())
  )
  
  structure(
    list(
      header     = header,
      events     = events,
      entity     = entity,
      reference  = reference,
      provenance = provenance,
      strict     = isTRUE(strict)
    ),
    class = c("genome_entity_gd","micromicon_gd")
  )
}

validate_genome_entity_gd <- function(x, strict = x$strict) {
  stopifnot(inherits(x, "genome_entity_gd"))
  msgs <- character(); status <- "ok"
  
  # Identity invariant: hashes must be unique
  hashes <- vapply(x$events, `[[`, character(1), "hash")
  if (anyDuplicated(hashes)) { msgs <- c(msgs, "Duplicate event hashes detected."); status <- "error" }
  
  # Contig & coordinate discipline vs Mode A
  contig_lengths <- setNames(as.numeric(x$entity$metadata$length_bp), x$entity$metadata$seqname)
  for (ev in x$events) {
    if (!ev$contig %in% names(contig_lengths)) {
      msgs <- c(msgs, sprintf("Contig in GD not present in reference: %s", ev$contig))
      status <- if (strict) "error" else "warn"
    } else if (!is.na(ev$position)) {
      ln <- contig_lengths[[ev$contig]]
      if (ev$position < 1L || ev$position > ln) {
        msgs <- c(msgs, sprintf("Position out of bounds for %s: %d not in [1..%d]", ev$contig, ev$position, ln))
        status <- if (strict) "error" else "warn"
      }
    }
  }
  
  if (identical(status, "error") && strict) {
    stop(paste(c("Validation failed (strict mode).", msgs), collapse = "\n- "))
  } else if (identical(status, "warn") || (identical(status, "error") && !strict)) {
    cli::cli_warn(paste(c("Validation warnings:", msgs), collapse = "\n- "))
  }
  
  x$provenance$validation$status   <- if (strict && identical(status, "ok")) "ok_strict" else status
  x$provenance$validation$messages <- unique(c(x$provenance$validation$messages, msgs))
  x
}

