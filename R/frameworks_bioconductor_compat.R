#' Bioconductor Compatibility Layer
#'
#' @description
#' Functions for checking Bioconductor package availability and version compatibility.
#' This layer enables graceful degradation when Bioconductor packages are not available.
#'
#' @keywords internal


#' Check if Bioconductor packages are available
#'
#' @description
#' Checks whether the required Bioconductor packages (GenomicRanges, Biostrings,
#' rtracklayer, Rsamtools) are installed and can be loaded.
#'
#' @return Logical indicating if Bioconductor packages are available
#' @export
#'
#' @examples
#' if (has_bioconductor()) {
#'   message("Bioconductor is available")
#' }
has_bioconductor <- function() {
  requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("Biostrings", quietly = TRUE) &&
    requireNamespace("rtracklayer", quietly = TRUE) &&
    requireNamespace("Rsamtools", quietly = TRUE)
}


#' Check if a specific Bioconductor package is available
#'
#' @param package Character string naming the package to check
#' @return Logical indicating if the package is available
#' @keywords internal
has_bioc_package <- function(package) {
  requireNamespace(package, quietly = TRUE)
}


#' Get Bioconductor package versions
#'
#' @description
#' Returns version information for installed Bioconductor packages.
#'
#' @return Named list of package versions, or NULL if packages not available
#' @export
bioconductor_versions <- function() {
  if (!has_bioconductor()) {
    return(NULL)
  }

  packages <- c("GenomicRanges", "Biostrings", "rtracklayer", "Rsamtools")
  versions <- lapply(packages, function(pkg) {
    tryCatch(
      as.character(utils::packageVersion(pkg)),
      error = function(e) NA_character_
    )
  })

  setNames(versions, packages)
}


#' Check minimum Bioconductor version requirements
#'
#' @param min_versions Named list of minimum versions (e.g., list(GenomicRanges = "1.40.0"))
#' @return Logical indicating if all version requirements are met
#' @keywords internal
check_bioc_versions <- function(min_versions = list()) {
  if (!has_bioconductor()) {
    return(FALSE)
  }

  current <- bioconductor_versions()

  for (pkg in names(min_versions)) {
    if (!pkg %in% names(current)) {
      return(FALSE)
    }

    current_ver <- package_version(current[[pkg]])
    min_ver <- package_version(min_versions[[pkg]])

    if (current_ver < min_ver) {
      return(FALSE)
    }
  }

  TRUE
}


#' Suggest Bioconductor installation
#'
#' @description
#' Prints a helpful message about installing Bioconductor packages.
#'
#' @param feature Character string describing the feature that requires Bioconductor
#' @export
suggest_bioconductor <- function(feature = "this feature") {
  message(
    "Bioconductor packages are required for ", feature, ".\n",
    "To install Bioconductor packages, run:\n\n",
    '  if (!require("BiocManager", quietly = TRUE))\n',
    '      install.packages("BiocManager")\n',
    '  BiocManager::install(c("GenomicRanges", "Biostrings", "rtracklayer", "Rsamtools"))\n'
  )
}
