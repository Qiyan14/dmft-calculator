# Publishing and updating the website

The website at <https://qiyan14.github.io/dmft-calculator/> is the `docs/` folder of this repository, served by GitHub Pages (**Settings → Pages**: deploy from the branch **main**, folder **/docs**).

## Updating the calculator

1. Get a copy of the repository: in [GitHub Desktop](https://desktop.github.com/), choose **File → Clone repository** and pick `Qiyan14/dmft-calculator` (or run `git clone https://github.com/Qiyan14/dmft-calculator`).
2. Edit `app/app.R`, for example the coefficients in `PUBLISHED_MODEL`, and check it with `shiny::runApp("app")`.
3. In R, from the repository folder, run `install.packages("shinylive")` (once) and then `source("build.R")`. This rebuilds `docs/`. The first build downloads the Shinylive files, a few hundred MB, once.
4. Commit and push: in GitHub Desktop, **Commit to main** and then **Push origin**. The site updates within a few minutes; visitors may need to reload the page.

The site has about 190 files and the GitHub website accepts at most 100 files per upload, so update `docs/` with GitHub Desktop or git rather than the upload page.

## Limits

A GitHub Pages site can be up to 1 GB and has a soft bandwidth limit of 100 GB per month. A first visit downloads about 42 MB, mostly R itself, so the limit covers roughly 2,000 first-time visitors a month; repeat visits mostly use the browser's cache. GitHub Pages may not be used to run an online business or a paid software service, so if the calculator becomes a paid product, move it to commercial hosting.


## Selected age-and-sex model update (October 2026)

The new results select six biomarkers plus age and sex. Run `analysis/Figures_DMFT.V3.R` in the matching preprocessing session to generate `dmft_selected_model.csv`. Put this verified file in `app/` before rebuilding if the selected model should be available to every visitor without uploading a file. The export includes the sex reference category; do not replace it with the previous age/teeth/BMI coefficients or a demographics-only fit.

Without this CSV, the site deliberately shows the new panel and development metrics while withholding its probability calculation. Visitors can still load a complete model CSV or train their own model. The current source-only update changes `app/app.R` and its matching embedded entry in `docs/app.json`; existing Shinylive runtime assets are unchanged. Future full rebuilds should continue to use `build.R`.
