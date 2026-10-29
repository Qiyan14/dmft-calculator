# Updates the website after a change to app/app.R only, without a full Shinylive rebuild.
# It rewrites the copy of app.R inside docs/app.json and changes the version tag that
# docs/index.html adds to its app.json request, so returning visitors load the new app.
#
# Run in R from this folder (the repository root), then commit and push app/ and docs/:
#   install.packages(c("jsonlite", "digest"))   # once
#   source("update_app_source.R")
# Use build.R instead when the Shinylive version or the app's packages change.

src_file <- "app/app.R"
src <- readChar(src_file, file.size(src_file), useBytes = TRUE)
Encoding(src) <- "UTF-8"
writeLines(jsonlite::toJSON(list(list(name = "app.R", content = src)), auto_unbox = TRUE),
           "docs/app.json", useBytes = TRUE, sep = "")

tag <- substr(digest::digest(file = src_file, algo = "sha1"), 1, 7)
idx <- readLines("docs/index.html", encoding = "UTF-8")
hit <- grepl("app\\.json\\?v=[0-9a-f]+", idx)
if (!any(hit)) stop("docs/index.html has no versioned app.json request; restore it with git checkout docs/index.html")
idx <- sub("app\\.json\\?v=[0-9a-f]+", paste0("app.json?v=", tag), idx)
writeLines(idx, "docs/index.html", useBytes = TRUE)
message("docs/app.json rewritten from ", src_file, "; app.json version tag set to ", tag)
