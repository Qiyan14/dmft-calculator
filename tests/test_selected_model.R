# Checks the calculator's selected DMFT model (Risk Calculator, both versions).
# Run with the path of the app file:
#   Rscript tests/test_selected_model.R app/app.R                        (repository root; browser version)
#   Rscript tests/test_selected_model.R Risk_Calculator_Shiny_DMFT.V6.R  (desktop version, with tests/ beside it)
# 1. The app's prediction equals the published equation (Table 3) computed
#    independently, for 10,000 random inputs covering zeros, both sexes and ages 18-100.
# 2. Saving the model to CSV and loading it back gives identical predictions and keeps
#    the sex labels.
# 3. A full-precision export in the format written by Equation_Coefficients_DMFT.V2.R /
#    Figures_DMFT.V3.R is accepted by the desktop loader and predicts identically.
# 4. The example inputs shown in the manuscript figures give the reported risk.
args <- commandArgs(trailingOnly = TRUE)
app_file <- if (length(args)) args[1] else "app.R"
options(dmft.app.test = TRUE)
suppressMessages(source(app_file, local = TRUE, encoding = "UTF-8"))

m <- PUBLISHED_MODEL
# --- Table 3, typed independently of the app's model object -----------------
b0    <- -59.61927668169686
b_bm  <- c(B511 = -0.039074625958201936, P475 = 0.062822476580071382, LAC = 0.079517386755282984,
           BG5 = -0.95827948205066982, BH = -0.89043902031579725, P950 = -0.88287872589894956)
b_age <- 0.1178751841582006
b_male <- 0.11282243181054701
table3_eta <- function(x, age, male) b0 + sum(b_bm * log2(x + 5e-7)) + b_age * age + b_male * male

set.seed(20261029)
n <- 10000
diffs <- numeric(n); dp <- numeric(n)
for (i in seq_len(n)) {
  x <- ifelse(runif(6) < 0.35, 0, exp(rnorm(6, 3, 3)))   # zeros and widely spread positives
  age <- sample(18:100, 1); sex <- sample(c("1", "2"), 1)
  pr <- predict_model(m, x, list(Age_DOS = age, Sex = sex))
  ref <- table3_eta(x, age, as.numeric(sex == "2"))
  diffs[i] <- abs(pr$eta - ref); dp[i] <- abs(pr$p - plogis(ref))
}
cat(sprintf("1. App vs Table 3 equation, %d random inputs: max |difference| in log-odds %.3g, in probability %.3g\n",
            n, max(diffs), max(dp)))
stopifnot(max(diffs) < 1e-10)

# --- Save and reload ---------------------------------------------------------
f <- tempfile(fileext = ".csv")
utils::write.csv(model_to_table(m), f, row.names = FALSE, na = "")
m2 <- model_from_table(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE), "test.csv")
x <- c(1, 1, 1, 0, 0, 0)
for (sex in c("1", "2")) stopifnot(abs(predict_model(m2, x, list(Age_DOS = 40, Sex = sex))$eta -
                                       predict_model(m, x, list(Age_DOS = 40, Sex = sex))$eta) < 1e-12)
stopifnot(identical(unname(m2$level_labels$Sex[c("1", "2")]), c("Female", "Male")))
cat("2. Model saved to CSV and reloaded: identical predictions; sex labels kept (1 = Female, 2 = Male)\n")

# --- Export in the analysis scripts' format ------------------------------------
exp_tab <- data.frame(
  Term = c("(Intercept)", m$markers$key, "Age_DOS", "Sex2", "(Threshold)", "(Delta)", "(Case group)",
           "(Reference group)", "(Reference level: Sex)"),
  Label = c("(Intercept)", m$markers$key, "Age_DOS", "Sex2", "Classification threshold", "log2 offset",
            "high DMFT (DMFT > 13)", "low DMFT (DMFT <= 13)", "1"),
  Coefficient_raw_input = c(b0, b_bm, b_age, b_male, m$threshold, 5e-7, NA, NA, NA),
  Input = c("", rep("log2(raw value + 5e-07)", 6), "value in original units", "indicator: Sex = 2", rep("", 5)),
  stringsAsFactors = FALSE)
d <- tempfile(); dir.create(d)
utils::write.csv(exp_tab, file.path(d, "dmft_selected_model.csv"), row.names = FALSE, na = "")
ld <- load_published_coefficients(m, d)
stopifnot(identical(ld$msg, "Complete selected-model coefficients and sex coding loaded."))
for (sex in c("1", "2")) stopifnot(abs(predict_model(ld$model, x, list(Age_DOS = 40, Sex = sex))$eta -
                                       table3_eta(x, 40, sex == "2")) < 1e-12)
stopifnot(identical(level_label(ld$model, "Sex", c("1", "2")), c("Female", "Male")))
cat("3. Analysis-script export accepted by the desktop loader; identical predictions and labels\n")

# --- Example shown in the manuscript figures -----------------------------------
ex <- predict_model(m, EXAMPLE$markers, list(Age_DOS = 40, Sex = "1"))
cat(sprintf("4. Example inputs (1, 1, 1, 0, 0, 0; age 40; female): eta = %.4f, risk = %.2f%% (threshold %.2f%%)\n",
            ex$eta, 100 * ex$p, 100 * m$threshold))
stopifnot(abs(ex$p - 0.9066) < 5e-5)
cat("All checks passed.\n")
