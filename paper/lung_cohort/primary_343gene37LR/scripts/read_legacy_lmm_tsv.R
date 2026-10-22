# The legacy writer left multiline lme4 messages unquoted in a TSV.
# A valid record always begins with the length-prefixed communication key.
read_legacy_lmm_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  stopifnot(length(lines) >= 2L, startsWith(lines[1L], "feature_id\t"))
  body <- lines[-1L]
  starts <- grepl("^[0-9]+:[^\t]*\t", body)
  stopifnot(starts[1L])
  records <- vapply(split(body, cumsum(starts)), paste, character(1L),
                    collapse = "\\n")
  read.delim(text = paste(c(lines[1L], records), collapse = "\n"),
             quote = "", comment.char = "", check.names = FALSE,
             stringsAsFactors = FALSE)
}
