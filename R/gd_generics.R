
# JUST STUBS RIGHT NOW!!!
#' Pretty-printer for genome_entity_gd
#'
#' Prints the inner genome_entity first, then a compact genome_diff summary.
#' @param x A genome_entity_gd object.
#' @param ... Passed to inner print/summary if needed.
print.genome_entity_gd <- function(x, ...) {
  # 1) Safety checks
  if (!inherits(x, "genome_entity_gd")) {
    stop("print.genome_entity_gd: object is not a 'genome_entity_gd'")
  }
  
  # 2) Print the embedded genome_entity using its own printer
  if (!is.null(x$entity) && inherits(x$entity, "genome_entity")) {
    print(x$entity)  # uses print.genome_entity()
  } else {
    cat("genome_entity\n")
    cat("===============\n")
    cat("(no embedded genome_entity available)\n")
  }
  
  cat("\n")  # visual spacer
  
  # 3) Summarize the GD and print under a genome_diff header
  #    We reuse your existing summary(.), then a small custom formatter.
  s <- tryCatch(summary(x), error = function(e) e)
  
  cat("genome_diff\n")
  cat("===============\n\n")
  
  if (inherits(s, "error")) {
    # If summary failed for any reason, don't explode the printer.
    cat("(summary unavailable: ", conditionMessage(s), ")\n", sep = "")
    return(invisible(x))
  }
  
  # s should be a summary_genome_entity_gd (per your code). We'll print it concisely here.
  # If you prefer to reuse your existing print.summary_genome_entity_gd(), replace this
  # section with: print(s); return(invisible(x))
  
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






