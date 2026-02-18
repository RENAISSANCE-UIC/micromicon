
len_bp <- function(x) {
  # DNAString / RNAString / AAString
  if (methods::is(x, "XString")) return(Biostrings::width(x))
  # XStringSet (take first by default; caller can choose otherwise)
  if (methods::is(x, "XStringSet")) return(Biostrings::width(x)[1L])
  # plain character
  if (is.character(x)) return(nchar(x[1L]))
  # raw or integer vectors (occasionally used to store sequences)
  if (is.raw(x) || is.integer(x)) return(length(x))
  # last resort: try character coercion
  return(nchar(as.character(x)[1L]))
}
