data("NIRcannabis", package = "proximetricsR")

dat <- NIRcannabis[41:80, ]
X <- dat$spc
rownames(X) <- NULL
Y <- matrix(dat$THC, dimnames = list(41:80, "THC"))

base_recipe <- preprocess_recipe(
  prep_wav_trim(band = c(1100, 1600)),
  prep_snv(),
  prep_derivative(m = 1, w = 11, p = 5, algorithm = "savitzky-golay"),
  device = "unspecified"
)

model <- calibrate(
  X, Y,
  data = dat, preprocess = base_recipe, method = fit_plsr(5, "standard"),
  control = calibration_control("none"), verbose = FALSE
)

test_that("export_sklearn_model requires a spectral_model object", {
  expect_error(export_sklearn_model(list()), "'object' must be of class 'spectral_model'")
})

test_that("export_sklearn_model rejects nwp-based models", {
  nwp_model <- calibrate(
    X, Y,
    data = dat, preprocess = preprocess_recipe(prep_snv()),
    method = fit_plsr(5, "nwp"), control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(nwp_model), 'type = "nwp"')
})

test_that("export_sklearn_model rejects nwp preprocessing steps", {
  nwp_recipe <- preprocess_recipe(
    prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp"),
    device = "unspecified"
  )
  nwp_model <- calibrate(
    X, Y,
    data = dat, preprocess = nwp_recipe, method = fit_plsr(5, "standard"),
    control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(nwp_model), 'algorithm = "nwp"')
})

test_that("export_sklearn_model rejects prep_derivative(algorithm = 'gap-segment')", {
  recipe <- preprocess_recipe(
    prep_derivative(m = 1, w = 3, p = 5, algorithm = "gap-segment"),
    device = "unspecified"
  )
  gap_model <- calibrate(
    X, Y,
    data = dat, preprocess = recipe, method = fit_plsr(5, "standard"),
    control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(gap_model), "gap-segment")
})

test_that("export_sklearn_model rejects prep_transform(to = 'reflectance')", {
  recipe <- preprocess_recipe(
    prep_transform(to = "reflectance"),
    device = "proxiscout"
  )
  reflectance_model <- calibrate(
    X, Y,
    data = dat, preprocess = recipe, method = fit_plsr(5, "standard"),
    control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(reflectance_model), 'to = "reflectance"')
})

test_that("export_sklearn_model rejects prep_wav_trim(trim_constant_edges = TRUE)", {
  recipe <- preprocess_recipe(
    prep_wav_trim(band = c(1100, 1600), trim_constant_edges = TRUE),
    device = "unspecified"
  )
  trim_model <- calibrate(
    X, Y,
    data = dat, preprocess = recipe, method = fit_plsr(5, "standard"),
    control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(trim_model), "trim_constant_edges = TRUE")
})

test_that("export_sklearn_model rejects prep_resample", {
  recipe <- preprocess_recipe(
    prep_resample(grid = c(1100, 1600, 5)),
    device = "proximate"
  )
  resample_model <- calibrate(
    X, Y,
    data = dat, preprocess = recipe, method = fit_plsr(5, "standard"),
    control = calibration_control("none"), verbose = FALSE
  )
  expect_error(export_sklearn_model(resample_model), "prep_resample")
})

test_that("export_sklearn_model produces the expected Pipeline JSON shape", {
  json <- export_sklearn_model(model)
  doc <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_equal(doc$estimator_class, "Pipeline")
  steps <- doc$params$steps
  expect_length(steps, length(base_recipe$steps) + 1)

  step_classes <- sapply(steps, function(s) s[[2]]$estimator_class)
  expect_equal(
    unlist(step_classes),
    c("RangeCut", "StandardNormalVariate", "SavitzkyGolay", "NIRWiseLinearModel")
  )

  model_step <- steps[[length(steps)]][[2]]
  expect_equal(model_step$params$fit_method, "plsr")
  expect_equal(model_step$params$type, "standard")
  expect_length(model_step$attributes$coef_, length(model$final_model$model$x_means))

  expect_equal(doc$metadata$domain, "sklearn")
  expect_equal(doc$metadata$source, "proximetricsR")
})

test_that("export_sklearn_model can write to a file", {
  file <- tempfile(fileext = ".json")
  result <- export_sklearn_model(model, file = file)
  expect_true(file.exists(file))
  expect_equal(as.character(result), paste(readLines(file), collapse = "\n"))
})
