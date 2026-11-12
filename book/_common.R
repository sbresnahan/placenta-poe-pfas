# Helpers shared by the book chapters: include a repository markdown file,
# shifting its headings so it sits under the chapter title, and turning
# ```mermaid blocks into raw <pre class="mermaid"> for the mermaid script.
html_escape <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

include_md <- function(path, shift = 0, drop_title = TRUE, drop_contents = TRUE) {
  x <- readLines(path, warn = FALSE)
  if (drop_title) { h1 <- grep("^# ", x)[1]; if (!is.na(h1)) x <- x[-h1] }
  if (drop_contents) {
    s <- grep("^Contents$", x)
    if (length(s)) { e <- which(x == "---" & seq_along(x) > s[1])[1]; x <- x[-(s[1]:e)] }
  }
  out <- character(0); in_mm <- FALSE; in_code <- FALSE
  for (l in x) {
    if (!in_code && grepl("^```mermaid", l)) { in_mm <- TRUE; out <- c(out, "", "```{=html}", '<pre class="mermaid">'); next }
    if (in_mm && grepl("^```", l)) { in_mm <- FALSE; out <- c(out, "</pre>", "```", ""); next }
    if (in_mm) { out <- c(out, html_escape(l)); next }
    if (grepl("^```", l)) in_code <- !in_code
    if (!in_code && shift > 0 && grepl("^#+ ", l)) l <- paste0(strrep("#", shift), l)
    out <- c(out, l)
  }
  cat(out, sep = "\n")
}
show_file <- function(path, lang = "bash") {
  cat(paste0("```", lang), readLines(path, warn = FALSE), "```", sep = "\n")
}
