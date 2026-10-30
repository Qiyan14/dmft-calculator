# Rebuilds the website (docs/) from the app (app/app.R).
#
# Run in R from this folder (the repository root) after changing app/app.R,
# then commit and push the updated docs/ folder:
#   install.packages("shinylive")   # once
#   source("build.R")

shinylive::export(
  appdir = "app",
  destdir = "docs",
  # The app needs only shiny, which Shinylive already includes.
  # Set this to TRUE if you add other R packages to the app.
  wasm_packages = FALSE,
  template_params = list(
    title = "Salivary biomarker DMFT calculator",
    include_in_head = paste0(
      '<meta name="description" content="Estimate the risk of high DMFT (DMFT > 13) ',
      'from six salivary biomarkers, age and BMI. Runs entirely in your browser.">'
    ),
    # Loading message shown while R starts in the browser
    include_after_body = paste(readLines("loading.html", encoding = "UTF-8"), collapse = "\n")
  )
)

# Not needed: the code-editor copy of the app that Shinylive also writes
unlink("docs/edit", recursive = TRUE)

# GitHub Pages: serve the files exactly as they are (no Jekyll processing)
file.create("docs/.nojekyll")
