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
