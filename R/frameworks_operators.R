#' Null-coalescing Operator
#'
#' @description
#' The null-coalescing operator `%||%` returns the right-hand side if the
#' left-hand side is NULL, otherwise it returns the left-hand side.
#'
#' This is a common pattern for providing default values when a variable
#' might be NULL.
#'
#' @param x An object (can be NULL)
#' @param y A fallback value to use if x is NULL
#'
#' @return Returns `x` if it is not NULL, otherwise returns `y`
#'
#' @details
#' **Pronunciation**: This operator is commonly called "or-else" or "null-or".
#' The unwieldy "grapes-or-or-grapes" refers to the visual appearance of the
#' percent signs (%) looking like grapes, but "or-else" is much clearer.
#'
#' This operator is similar to rlang's `%||%` but implemented here to avoid
#' a hard dependency on rlang.
#'
#' @examples
#' # Basic usage
#' NULL %||% "default"        # Returns "default"
#' "value" %||% "default"     # Returns "value"
#'
#' # Common pattern for function arguments
#' my_function <- function(x = NULL) {
#'   x <- x %||% 10  # Use 10 if x is NULL
#'   x * 2
#' }
#'
#' @name null-coalesce
#' @rdname null-coalesce
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
