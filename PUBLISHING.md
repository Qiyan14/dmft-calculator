# Publishing and updating the website

The website at <https://qiyan14.github.io/dmft-calculator/> is the `docs/` folder of this repository, served by GitHub Pages (**Settings → Pages**: deploy from the branch **main**, folder **/docs**).

## Updating the calculator

1. Get a copy of the repository: in [GitHub Desktop](https://desktop.github.com/), choose **File → Clone repository** and pick `Qiyan14/dmft-calculator` (or run `git clone https://github.com/Qiyan14/dmft-calculator`).
2. Edit `app/app.R`, for example the coefficients in `PUBLISHED_MODEL`, and check it with `shiny::runApp("app")`.
3. In R, from the repository folder, run `source("update_app_source.R")`. This copies `app/app.R` into `docs/app.json` and updates the version tag in `docs/index.html` (see below). A full rebuild with `install.packages("shinylive")` (once) and `source("build.R")` is needed only for a new Shinylive version or added R packages; the first build downloads the Shinylive files, a few hundred MB, once.
4. Commit and push: in GitHub Desktop, **Commit to main** and then **Push origin**. The site updates within a few minutes; visitors may need to reload the page.

The site has about 190 files and the GitHub website accepts at most 100 files per upload, so update `docs/` with GitHub Desktop or git rather than the upload page.

## Limits

A GitHub Pages site can be up to 1 GB and has a soft bandwidth limit of 100 GB per month. A first visit downloads about 42 MB, mostly R itself, so the limit covers roughly 2,000 first-time visitors a month; repeat visits mostly use the browser's cache. GitHub Pages may not be used to run an online business or a paid software service, so if the calculator becomes a paid product, move it to commercial hosting.

## Selected age-and-sex model (October 2026)

The selected model (six biomarkers plus age and sex) is bundled in `app/app.R` as `PUBLISHED_MODEL`, with the full-precision coefficients from `dmft_selected_model.csv` (written by `analysis/Figures_DMFT.V3.R`), female (code 1) as the reference sex category, and the labels Female and Male for the sex input. Every visitor can use it without uploading a file. If the coefficients change, copy the new values into `PUBLISHED_MODEL` and run `Rscript tests/test_selected_model.R app/app.R` (after updating the values typed in the test) before publishing. Do not replace the model with the earlier age/teeth/BMI coefficients or a demographics-only fit.

The desktop copy of the app (`Risk_Calculator_Shiny_DMFT.V6.R`) is the same file except for one line: it also reads `dmft_selected_model.csv` when that file is placed beside it.

## Updating only the app source

The October 2026 update changed only `app/app.R` and the matching copy of it inside `docs/app.json`; the Shinylive runtime files in `docs/shinylive/` are unchanged. `docs/index.html` requests `app.json?v=<tag>` with caching turned off, so returning visitors receive the new app; the tag is the first seven characters of the SHA-1 checksum of `app/app.R` (in R, `substr(digest::digest(file = "app/app.R", algo = "sha1"), 1, 7)`). For a source-only update, run `source("update_app_source.R")` in R from the repository folder; it rewrites `docs/app.json` from `app/app.R` and sets the tag.

`build.R` rebuilds all of `docs/`, including `docs/index.html` from the Shinylive template, which removes the versioned request. After a full rebuild, restore the file with `git checkout docs/index.html` and run `source("update_app_source.R")`.
