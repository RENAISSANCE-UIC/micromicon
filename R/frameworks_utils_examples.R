#' Get Path to Example Files
#'
#' @description
#' Helper function to locate example data files included with the package.
#' Useful for testing and documentation examples.
#'
#' @param filename Character string of the filename in inst/extdata
#'
#' @return Character string with full path to the example file
#' @export
#'
#' @examples
#' \dontrun{
#' # Get path to example GenBank file
#' gbk_file <- get_example_file("test_ampC.gbk")
#'
#' # Read the example file
#' genome <- read_genome(gbk_file)
#' }
get_example_file <- function(filename) {
  # Get file from installed package
  pkg_file <- system.file("extdata", filename, package = "micromicon")

  if (!file.exists(pkg_file) || pkg_file == "") {
    cli::cli_abort(c(
      "Example file not found: {filename}",
      "i" = "Available example files can be listed with: list.files(system.file('extdata', package = 'micromicon'))"
    ))
  }

  pkg_file
}
