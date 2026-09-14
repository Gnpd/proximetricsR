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
