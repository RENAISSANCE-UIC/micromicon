
#' Constructor for genome_entity_gd
#' @keywords internal
new_genome_entity_gd <- function(header, events, file, reference = NULL, strict = TRUE) {
  stopifnot(is.list(header), is.list(events), is.list(file))
  stopifnot(is.null(reference) || is.list(reference))
  
  # light normalization of events: ensure standard names exist
  events <- lapply(events, function(ev) {
    req <- c("type","id","rank","genome","position","ref_seq_fixed","alt_or_len","tags","raw_line","hash")
    missing <- setdiff(req, names(ev))
    if (length(missing)) stop("Event missing fields: ", paste(missing, collapse = ", "))
    ev
  })
  
  provenance <- list(
    gd_file = file,
    reference_checksums = if (!is.null(reference) && !is.null(reference$checksums)) reference$checksums else NULL,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    host       = Sys.info()[["nodename"]],
    validation = list(status = "unvalidated", messages = character())
  )
  
  structure(
    list(
      header     = header,
      events     = events,
      reference  = reference,  # may be NULL until caller binds it
      provenance = provenance,
      strict     = isTRUE(strict)
    ),
    class = c("genome_entity_gd", "micromicon_gd")
  )
}

#' Validate a genome_entity_gd against a Mode A reference (if provided)
#'
#' @param x genome_entity_gd
#' @param strict logical
#' @export
validate_genome_entity_gd <- function(x, strict = x$strict) {
  stopifnot(inherits(x, "genome_entity_gd"))
  
  msgs <- character()
  status <- "ok"
  
  # Check for unique hashes (identity invariant)
  hashes <- vapply(x$events, `[[`, character(1), "hash")
  if (anyDuplicated(hashes)) {
    msgs <- c(msgs, "Duplicate event hashes detected.")
    status <- "error"
  }
  
  # If reference provided, verify contig name/length discipline
  if (!is.null(x$reference) && !is.null(x$reference$contigs)) {
    contigs_df <- x$reference$contigs
    # Gather contig mentions from tags if present (e.g., seq_id) or infer from header if you have a map
    # You can refine this as your GD schema stabilizes.
    # Here we only check lengths if header enumerates contigs (optional).
    if (!all(c("name","length") %in% names(contigs_df))) {
      msgs <- c(msgs, "Reference contigs must have columns: name, length.")
      status <- "error"
    }
    
    # Example: if header includes CONTIG and LENGTH lists, compare.
    # This is a placeholder; wire to your real header keys when available.
    hdr_names  <- unname(unlist(x$header[["CONTIG"]] %||% list()))
    hdr_lengths <- suppressWarnings(as.numeric(unname(unlist(x$header[["LENGTH"]] %||% list()))))
    
    if (length(hdr_names) && length(hdr_lengths) && length(hdr_names) == length(hdr_lengths)) {
      ref_map <- setNames(contigs_df$length, contigs_df$name)
      for (j in seq_along(hdr_names)) {
        nm <- hdr_names[j]
        ln <- hdr_lengths[j]
        if (!nm %in% names(ref_map)) {
          msgs <- c(msgs, sprintf("Contig in GD not present in reference: %s", nm))
          status <- if (strict) "error" else "warn"
        } else if (!isTRUE(as.numeric(ref_map[[nm]]) == as.numeric(ln))) {
          msgs <- c(msgs, sprintf("Length mismatch for contig %s: GD=%s, REF=%s", nm, ln, ref_map[[nm]]))
          status <- if (strict) "error" else "warn"
        }
      }
    } else {
      # If GD header doesn't enumerate contigs/lengths, we decline to liftover or remap.
      # We only assert that we won't silently transform coordinates later.
      msgs <- c(msgs, "GD header lacks contig/length enumeration; no remapping will be performed.")
      if (strict) {
        status <- "warn"
      }
    }
  }
  
  if (identical(status, "error") && strict) {
    stop(paste(c("Validation failed (strict mode).", msgs), collapse = "\n- "))
  } else if (identical(status, "warn") || (identical(status, "error") && !strict)) {
    warning(paste(c("Validation warnings:", msgs), collapse = "\n- "))
  }
  
  x$provenance$validation$status    <- if (strict && identical(status, "ok")) "ok_strict" else status
  x$provenance$validation$messages  <- unique(c(x$provenance$validation$messages, msgs))
  x
}
