#' Pretty-printer for genome_entity_gd
#'
#' Prints the genome_entity view (via hoisted fields), then a compact genome_diff summary.
#' @param x A genome_entity_gd object.
#' @param ... Passed to parent print method.
#' @export
print.genome_entity_gd <- function(x, ...) {
  # 1) Safety checks
  if (!inherits(x, "genome_entity_gd")) {
    stop("print.genome_entity_gd: object is not a 'genome_entity_gd'")
  }

  # 2) Print the parent genome_entity view using hoisted fields
  #    NextMethod() will dispatch to print.genome_entity()
  NextMethod("print")

  cat("\n")  # visual spacer

  # 3) Summarize the GD and print under a genome_diff header
  s <- tryCatch(summary(x), error = function(e) e)

  cat("genome_diff\n")
  cat("===============\n\n")

  if (inherits(s, "error")) {
    # If summary failed for any reason, don't explode the printer.
    cat("(summary unavailable: ", conditionMessage(s), ")\n", sep = "")
    return(invisible(x))
  }

  # Defensive extraction
  n_events <- if (!is.null(s$n_events)) s$n_events else length(x$events)
  by_type  <- if (!is.null(s$by_type))  s$by_type  else table(character(0))
  status   <- tryCatch(s$validation$status, error = function(e) "unknown")

  cat(sprintf("Genome events: %d\n", as.integer(n_events)))
  cat("By type:\n")
  if (length(by_type)) {
    for (k in names(by_type)) {
      cat(sprintf("  - %s: %d\n", k, as.integer(by_type[[k]])))
    }
  } else {
    cat("  (none)\n")
  }
  cat("Validation status:", status, "\n")

  invisible(x)
}

summary.genome_entity_gd <- function(object, ...) {
  n <- length(object$events)
  kinds <- sort(table(vapply(object$events, `[[`, character(1), "type")), decreasing = TRUE)
  out <- list(n_events = n, by_type = kinds, validation = object$provenance$validation)
  class(out) <- "summary_genome_entity_gd"
  out
}






