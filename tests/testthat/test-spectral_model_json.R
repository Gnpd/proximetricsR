data("NIRcannabis", package = "proximetricsR")

dat <- NIRcannabis[41:80, ]
X <- dat$spc
rownames(X) <- NULL
Y <- matrix(dat$THC, dimnames = list(41:80, "THC"))

method <- fit_plsr(6, "standard")
control <- calibration_control(validation_type = "kfold", number = 3, folds = "random", tuning_parameter = "rsq", seed = 42)
pretreats <- preprocess_recipe(
  prep_snv(),
  prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp"),
  device = "unspecified"
)

model <- calibrate(
  X, Y,
  data = dat, preprocess = pretreats, method = method, control = control,
  metadata = add_model_metadata(unit = "%"), verbose = FALSE
)

test_that("save_spectral_model requires a spectral_model object", {
  expect_error(save_spectral_model(list(), tempfile()), "'object' must be of class 'spectral_model'")
})

test_that("save_spectral_model/load_spectral_model round-trips predictions", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(model, file)

  reloaded <- load_spectral_model(file)
  expect_s3_class(reloaded, "spectral_model")
  expect_null(reloaded$formula)

  original_pred <- predict(model, newdata = X, ncomp = 1:6, verbose = FALSE)$predictions
  reloaded_pred <- predict(reloaded, newdata = X, ncomp = 1:6, verbose = FALSE)$predictions

  expect_equal(unname(reloaded_pred), unname(original_pred), tolerance = 1e-08)
})

test_that("load_spectral_model reconstructs identifying fields and metadata", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(model, file)
  reloaded <- load_spectral_model(file)

  expect_equal(reloaded$target_variable, model$target_variable)
  expect_equal(reloaded$predictor_variables, model$predictor_variables)
  expect_equal(reloaded$final_ncomp, model$final_ncomp)
  expect_equal(reloaded$preprocess$device, model$preprocess$device)
  expect_length(reloaded$preprocess$steps, length(model$preprocess$steps))
  expect_equal(reloaded$metadata$Unit, model$metadata$Unit)
})

test_that("load_spectral_model reconstructs the fitted spectral_fit numerically", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(model, file)
  reloaded <- load_spectral_model(file)

  expect_s3_class(reloaded$final_model$model, "spectral_fit")
  expect_equal(
    reloaded$final_model$model$coefficients,
    model$final_model$model$coefficients,
    tolerance = 1e-08
  )
  expect_equal(
    unname(reloaded$final_model$model$x_means),
    unname(model$final_model$model$x_means),
    tolerance = 1e-08
  )
  expect_equal(
    unname(reloaded$final_model$model$intercept),
    unname(model$final_model$model$intercept),
    tolerance = 1e-08
  )
  expect_equal(reloaded$final_model$model$method$ncomp, model$final_model$model$method$ncomp)
  expect_equal(reloaded$final_model$model$method$type, model$final_model$model$method$type)
})

test_that("load_spectral_model errors on a file that is not a spectral_model export", {
  file <- tempfile(fileext = ".json")
  writeLines('{"estimator_class": "SomethingElse"}', file)
  expect_error(load_spectral_model(file), "does not look like a proximetricsR spectral_model")
})

# A recipe with a *vector-valued* step parameter (band) and a step that narrows
# the grid, so the round trip covers both param simplification and processed_wavs.
trim_recipe <- preprocess_recipe(
  prep_wav_trim(band = c(1100, 1600)),
  prep_snv(),
  prep_derivative(m = 1, w = 11, p = 2, algorithm = "savitzky-golay"),
  device = "unspecified"
)
trim_model <- calibrate(
  X, Y,
  data = dat, preprocess = trim_recipe, method = fit_plsr(3, "standard"),
  control = calibration_control("none"), verbose = FALSE
)

test_that("load_spectral_model restores vector-valued step parameters as atomic", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(trim_model, file)
  reloaded <- load_spectral_model(file)

  band <- reloaded$preprocess$steps[[1]]$band
  # read back with simplifyVector = FALSE, a JSON array arrives as a list and
  # reaches the step executor as one, where min()/max() fail
  expect_true(is.numeric(band))
  expect_equal(band, c(1100, 1600))
  expect_equal(
    predict(reloaded, newdata = X, verbose = FALSE)$predictions,
    predict(trim_model, newdata = X, verbose = FALSE)$predictions
  )
})

test_that("save_spectral_model/load_spectral_model round-trips processed_wavs", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(trim_model, file)
  reloaded <- load_spectral_model(file)

  # one entry per step, plus the incoming grid
  expect_named(reloaded$processed_wavs, paste0("step_", 0:length(trim_recipe$steps)))
  expect_equal(
    lapply(reloaded$processed_wavs, as.numeric),
    lapply(trim_model$processed_wavs, as.numeric)
  )
})

test_that("the JSON round trip is bit-exact, not merely close", {
  # jsonlite's digits = NA writes 15 significant digits, which loses up to ~22
  # ULP on a float64; the writers use digits = I(17) so doubles read back identical.
  file <- tempfile(fileext = ".json")
  save_spectral_model(trim_model, file)
  reloaded <- load_spectral_model(file)

  expect_identical(
    reloaded$final_model$model$coefficients,
    trim_model$final_model$model$coefficients
  )
  expect_identical(
    predict(reloaded, newdata = X, verbose = FALSE)$predictions,
    predict(trim_model, newdata = X, verbose = FALSE)$predictions
  )
})

test_that("a reloaded model re-exports to an identical sklearn pipeline", {
  file <- tempfile(fileext = ".json")
  save_spectral_model(trim_model, file)
  reloaded <- load_spectral_model(file)

  strip_created <- function(j) sub('"created_at":"[^"]*"', "", j)
  expect_identical(
    strip_created(export_sklearn_model(reloaded)),
    strip_created(export_sklearn_model(trim_model))
  )
})

test_that("export_sklearn_model refuses a model with no processed_wavs", {
  # predict() recomputes the grid, so a model from a format_version 1 file
  # predicts fine but cannot be re-exported; without the guard that wrote a
  # structurally valid pipeline whose coefficients were all null
  stripped <- trim_model
  stripped$processed_wavs <- NULL

  expect_error(export_sklearn_model(stripped), "processed wavelength grid")
})
