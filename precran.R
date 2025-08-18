# precran.R — Run typical pre-CRAN checks
# Usage from the package root:
#   Rscript precran.R
# or source("precran.R") from an interactive session.

message("==== Pre-CRAN preflight starting ====")

# 1) Ensure required helper packages are available
required <- c(
  "devtools",     # document, build, check, win-builder helpers
  "rcmdcheck",    # programmatic R CMD check
  "spelling",     # spell-check Rd, vignettes, R
  "urlchecker",   # validate and fix URLs in Rd, R, vignettes
  "lintr",        # static code analysis
  "styler",       # code style checks
  "desc",         # DESCRIPTION helpers
  "rhub"          # R-hub v2 remote checks (Windows/Linux)
)

missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  lapply(required, requireNamespace, quietly = TRUE, character.only = TRUE)
})

# 2) Basic metadata sanity checks
message("\n==> Checking DESCRIPTION sanity")
d <- desc::desc(file = "DESCRIPTION")
fields_needed <- c("Package", "Version", "Title", "Description", "License")
missing_fields <- fields_needed[!fields_needed %in% d$fields()]
if (length(missing_fields)) {
  stop("Missing DESCRIPTION fields: ", paste(missing_fields, collapse = ", "))
}
if (grepl("TODO|TBD|FILL ME|FIXME", d$get("Description"))) {
  stop("DESCRIPTION 'Description' contains placeholders. Please fix.")
}
message("DESCRIPTION looks OK")

# 3) Regenerate docs and namespace
message("\n==> Generating Rd files and NAMESPACE via roxygen2")
devtools::document(quiet = TRUE)

# # 4) Optional: build README if a README.Rmd exists
# if (file.exists("README.Rmd")) {
#   message("\n==> Building README")
#   try(devtools::build_readme(), silent = TRUE)
# }

# 5) NEWS and CRAN-RELEASE presence check (non-fatal)
if (!file.exists("NEWS.md")) message("Note: Consider maintaining a NEWS.md")
if (!file.exists("CRAN-SUBMISSION") && !file.exists("CRAN-RELEASE")) {
  message("Note: Consider using usethis::use_cran_comments() and keeping CRAN comments/notes")
}

# 6) URL checks
message("\n==> Checking URLs")
bad <- urlchecker::url_check()
if (nrow(bad)) {
  print(bad)
  message("You can try: urlchecker::url_update() to auto-fix redirects.")
} else {
  message("All URLs OK")
}

# 7) Spell check
message("\n==> Spell checking package")
sp <- spelling::spell_check_package()
if (nrow(sp)) {
  print(sp[, c("file", "line", "word")])
  message("Add accepted words to inst/WORDLIST or fix typos.")
} else {
  message("No spelling issues found")
}

# 8) Lint R code
message("\n==> Linting R code")
lint_res <- lintr::lint_package()
if (length(lint_res)) {
  print(lint_res)
  message("Consider running styler::style_pkg() or address lintr suggestions.")
} else {
  message("No lints found")
}

# 9) Local CRAN-style check
message("\n==> Local R CMD check --as-cran")
# rcmdcheck returns an object with summaries; devtools::check() opens viewer; prefer rcmdcheck here
chk <- rcmdcheck::rcmdcheck(args = c("--as-cran"), error_on = "never", build_args = "--no-manual")
summary <- rcmdcheck::parse_check(chk$results)
print(chk)
if (length(chk$errors))  message("Errors: ", length(chk$errors))
if (length(chk$warnings)) message("Warnings: ", length(chk$warnings))
if (length(chk$notes))    message("Notes: ", length(chk$notes))

# 10) Build source tarball
message("\n==> Building source tarball")
pkgfile <- devtools::build(path = "dist", quiet = TRUE)
message("Built: ", pkgfile)

# 11) Remote checks with R-hub v2 (Windows + Linux)
message("\n==> R-hub v2 checks (Windows + Linux)")
# Requires a GitHub repo or local setup; see ?rhubv2
# If your package is on GitHub, set gh_url to your repo. Otherwise, rhub will infer from current git remote.
# Common platforms include: 'windows-x86_64-release', 'ubuntu-gcc-release'.
# Use rhub::rhub_check() with defaults or specify platforms explicitly.
try({
  #p <- rhub::rhub_platforms()
  rhub::rhub_check(platforms = "ubuntu-latest (r-release)")

  # pick <- function(rx) {
  #   x <- p$name[grepl(rx, p$name, ignore.case = TRUE)]
  #   if (length(x)) x[1] else character(0)
  # }
  #
  # plats <- unique(c(
  #   pick("^ubuntu.*r-release"),
  #   pick("^macos.*r-release"),
  #   pick("^windows.*r-release")
  # ))
  #
  # print(plats)  # sanity check; should be non-empty exact names
  # rhub::rhub_check(platforms = plats)
}, silent = TRUE)

# 12) Optional: Win-builder checks via devtools
message("\n==> Optional: win-builder checks")
message("Skip with options(precran.skip_winbuilder = TRUE) to bypass")
if (!isTRUE(getOption("precran.skip_winbuilder", FALSE))) {
  try(devtools::check_win_release(), silent = TRUE)
  try(devtools::check_win_devel(),   silent = TRUE)
}

message("\n==== Pre-CRAN preflight complete ====")
