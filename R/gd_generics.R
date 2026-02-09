
#' @export
effects <- function(x, ...) UseMethod("effects")

#' @export
effects.genome_entity_gd <- function(x, reference = x$reference, ...) {
  # Late-bound: refuse to infer effects until validation has passed
  if (is.null(reference)) stop("No reference bound; effects() requires a validated reference.")
  if (!identical(x$provenance$validation$status, "ok_strict") &&
      !identical(x$provenance$validation$status, "ok")) {
    stop("Reference validation has not passed; refusing to infer effects.")
  }
  
  # Placeholder: return a data.frame with one row per event and a list-column
  # "per_gene" that can have multiple effect rows for overlapping CDS, etc.
  # This preserves multiplicity; a separate helper can flatten on request.
  data.frame(
    hash   = vapply(x$events, `[[`, character(1), "hash"),
    type   = vapply(x$events, `[[`, character(1), "type"),
    rank   = vapply(x$events, `[[`, integer(1), "rank"),
    genome = vapply(x$events, `[[`, integer(1), "genome"),
    position = vapply(x$events, `[[`, integer(1), "position"),
    per_gene = I(vector("list", length(x$events))), # fill later
    stringsAsFactors = FALSE
  )
}

#' @export
summary.genome_entity_gd <- function(object, ...) {
  n <- length(object$events)
  kinds <- table(vapply(object$events, `[[`, character(1), "type"))
  out <- list(
    n_events = n,
    by_type = sort(kinds, decreasing = TRUE),
    validation = object$provenance$validation
  )
  class(out) <- "summary_genome_entity_gd"
  out
}

#' @export
print.summary_genome_entity_gd <- function(x, ...) {
  cat(sprintf("Genome events: %d\n", x$n_events))
  cat("By type:\n")
  for (k in names(x$by_type)) {
    cat(sprintf("  - %s: %d\n", k, x$by_type[[k]]))
  }
  cat("Validation status:", x$validation$status, "\n")
  invisible(x)
}

#' @export
view_event <- function(x, ...) UseMethod("view_event")

#' @export
view_event.genome_entity_gd <- function(x, id = NULL, hash = NULL, ...) {
  if (!is.null(hash)) {
    idx <- which(vapply(x$events, `[[`, character(1), "hash") == hash)
  } else if (!is.null(id)) {
    idx <- which(vapply(x$events, `[[`, integer(1), "id") == as.integer(id))
  } else stop("Supply either id or hash.")
  
  if (!length(idx)) return(NULL)
  x$events[[idx[1]]]
}

#' @export
view_gene <- function(x, ...) UseMethod("view_gene")

#' @export
view_gene.genome_entity_gd <- function(x, gene, ...) {
  # Find any event with tags$gene containing the symbol
  idx <- which(vapply(x$events, function(ev) {
    tg <- ev$tags[["gene"]] %||% ev$tags[["gene_name"]] %||% character(0)
    any(tg == gene)
  }, logical(1)))
  x$events[idx]
}

#' @export
as_tibble <- function(x, ...) UseMethod("as_tibble")

#' @export
as_tibble.genome_entity_gd <- function(x, flatten = FALSE, ...) {
  df <- data.frame(
    type     = vapply(x$events, `[[`, character(1), "type"),
    id       = vapply(x$events, `[[`, integer(1), "id"),
    rank     = vapply(x$events, `[[`, integer(1), "rank"),
    genome   = vapply(x$events, `[[`, integer(1), "genome"),
    position = vapply(x$events, `[[`, integer(1), "position"),
    ref_seq  = vapply(x$events, `[[`, character(1), "ref_seq_fixed"),
    alt_or_len = vapply(x$events, `[[`, character(1), "alt_or_len"),
    hash     = vapply(x$events, `[[`, character(1), "hash"),
    raw_line = vapply(x$events, `[[`, character(1), "raw_line"),
    stringsAsFactors = FALSE
  )
  # Attach tags as a list-column by default (multiplicity preserved)
  df$tags <- I(lapply(x$events, `[[`, "tags"))
  
  if (!flatten) return(df)
  
  # Flatten tags into repeated rows: one row per (event x tag-value)
  out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    tg <- df$tags[[i]]
    if (!length(tg)) {
      cbind(df[i, setdiff(names(df), "tags"), drop = FALSE],
            tag_key = NA_character_, tag_value = NA_character_)
    } else {
      rows <- do.call(rbind, lapply(names(tg), function(k) {
        vals <- tg[[k]]
        if (!length(vals)) {
          data.frame(tag_key = k, tag_value = NA_character_, stringsAsFactors = FALSE)
        } else {
          data.frame(tag_key = k, tag_value = vals, stringsAsFactors = FALSE)
        }
      }))
      cbind(df[i, setdiff(names(df), "tags"), drop = FALSE], rows)
    }
  }))
  rownames(out) <- NULL
  out
}

#' @export
write_ndjson <- function(x, ...) UseMethod("write_ndjson")

#' @export
write_ndjson.genome_entity_gd <- function(x, path, include_header = TRUE, ...) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  
  # optional header record as first line
  if (isTRUE(include_header)) {
    hdr <- list(
      type = "header",
      header = x$header,
      provenance = x$provenance
    )
    writeLines(jsonlite::toJSON(hdr, auto_unbox = TRUE), con = con, sep = "\n")
  }
  
  # one event per line
  for (ev in x$events) {
    rec <- list(type = "event", event = ev)
    writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE), con = con, sep = "\n")
  }
  invisible(normalizePath(path))
}
