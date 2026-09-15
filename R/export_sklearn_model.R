#' @title Export a spectral_model as an sklearn/chemotools Pipeline (JSON)
#' @name export_sklearn_model
#'
#' @description
#'
#' \loadmathjax
#'
#' Serializes a \code{\link{spectral_model}} object into openmodels-shaped JSON
#' whose \code{estimator_class} values are real scikit-learn/chemotools class
#' names, so it can be loaded directly in Python as a working
#' \code{sklearn.pipeline.Pipeline} via
#' \href{https://github.com/Gnpd/openmodels}{openmodels}, with
#' \href{https://github.com/paucablop/chemotools}{chemotools} and
#' \href{https://github.com/Gnpd/proximetricsr-estimators}{proximetricsr-estimators}
#' registered as custom estimator providers:
#'
#' \preformatted{
#' from openmodels import SerializationManager, SklearnSerializer
#' from chemotools.utils.discovery import all_estimators as chemotools_estimators
#' from proximetricsr_estimators.utils.discovery import all_estimators as pr_estimators
#'
#' manager = SerializationManager(
#'     SklearnSerializer(custom_estimators=[chemotools_estimators, pr_estimators])
#' )
#' model = manager.load("exported_from_R.json")
#' predictions = model.predict(X)
#' }
#'
#' @usage
#' export_sklearn_model(object, file = NULL)
#'
#' @param object an object of class \code{spectral_model}, as returned by
#' \code{\link{calibrate}}.
#' @param file an optional character string with the path (including file name)
#' where the JSON output should be written. If \code{NULL} (default), no file is
#' written and the JSON string is returned.
#'
#' @return If \code{file = NULL} (default), the JSON string is returned visibly
#' so it can be inspected or assigned to a variable. If \code{file} is specified,
#' the JSON is written to that file and returned invisibly.
#'
#' @details
#' Unlike \code{\link{save_spectral_model}}, this export is restricted to
#' preprocessing/method combinations that have a working chemotools/scikit-learn
#' equivalent reconstructible via a plain \code{Pipeline.predict(X)} call. A
#' clear, named error is raised for:
#' \itemize{
#'   \item \code{prep_derivative()}/\code{prep_smooth()} steps using
#'     \code{algorithm = "nwp"} (BUCHI NIRWise PLUS exact-match preprocessing
#'     math has no chemotools equivalent) -- note this restriction does *not*
#'     apply to \code{fit_plsr()}/\code{fit_xlsr()}'s own \code{type = "nwp"}:
#'     predictions from an \code{"nwp"}-type model are numerically identical to
#'     a \code{"modified"}-type model fitted on the same data (the two only
#'     differ in the internal score-space representation, not in
#'     \code{coefficients}/\code{intercept}), so \code{type = "nwp"} models
#'     export the same as \code{"modified"} ones;
#'   \item \code{prep_derivative(algorithm = "gap-segment")} (its chemotools
#'     equivalent, \code{NorrisWilliams}, precomputes an internal kernel that
#'     has not yet been verified to reproduce correctly from R);
#'   \item \code{prep_transform(to = "reflectance")} (mirrors
#'     \code{\link{proxiscout_write_model}}'s existing behaviour of only
#'     supporting the reflectance-to-absorbance direction);
#'   \item \code{prep_wav_trim(trim_constant_edges = TRUE)} (data-dependent, no
#'     chemotools equivalent);
#'   \item \code{prep_resample()} (its chemotools equivalent,
#'     \code{XAxisInterpolator}, requires scikit-learn metadata routing --
#'     the incoming spectra's x-axis must be passed explicitly to every
#'     \code{predict()} call -- unlike every other step here, so it cannot be
#'     used in a plain \code{Pipeline.predict(X)} workflow).
#' }
#' Use \code{\link{save_spectral_model}}/\code{\link{load_spectral_model}} instead
#' for a full-fidelity R-native round trip of any proximetricsR model, without
#' these restrictions.
#'
#' The PLS/XLS regression step is exported as a
#' \code{proximetricsr_estimators.regression.NIRWiseLinearModel}: prediction from
#' an already-fitted proximetricsR model is always affine
#' (\code{(X - x_means) \%*\% t(coefficients) + intercept}), regardless of fitting
#' algorithm, so this single class covers every supported
#' \code{fit_plsr()}/\code{fit_xlsr()} combination -- \code{fit_method}/\code{type}/
#' \code{min_w}/\code{max_w} are kept as constructor parameters purely for
#' provenance, matching \code{\link{fit_plsr}}/\code{\link{fit_xlsr}}.
#'
#' \code{prep_transform(to = "absorbance")} is exported as
#' \code{chemotools.physics.IntensityConversion(input_unit = "reflectance",
#' output_unit = "pseudoabsorbance")}, not \code{output_unit = "absorbance"}:
#' proximetricsR's \code{to = "absorbance"} computes \mjeqn{A = -\log_{10}(R)}{A =
#' -log10(R)} directly from reflectance, which is what chemotools calls
#' "pseudoabsorbance" (its "absorbance"/"transmittance" pair instead follows the
#' Beer-Lambert transmittance convention, a different physical quantity).
#'
#' @seealso \code{\link{calibrate}}, \code{\link{save_spectral_model}},
#' \code{\link{proxiscout_write_model}}
#'
#' @examples
#' \donttest{
#' data("NIRcannabis")
#' recipe <- preprocess_recipe(
#'   prep_wav_trim(band = c(1100, 1600)),
#'   prep_snv(),
#'   device = "unspecified"
#' )
#' model <- calibrate(CBDA ~ spc,
#'   data = NIRcannabis, preprocess = recipe,
#'   method = fit_plsr(5, type = "standard"),
#'   control = calibration_control("none"), verbose = FALSE
#' )
#' json <- export_sklearn_model(model)
#' }
#' @author Leonardo Ramirez-Lopez
#' @export
export_sklearn_model <- function(object, file = NULL) {
  if (!inherits(object, "spectral_model")) {
    stop("'object' must be of class 'spectral_model'.")
  }
  model <- object$final_model$model
  if (is.null(model)) {
    stop("'object' does not contain a fitted model (object$final_model$model is NULL).")
  }
  if (!is.null(file) && (!is.character(file) || length(file) != 1)) {
    stop("'file' must be a single character string, if provided.")
  }

  steps <- object$preprocess$steps
  step_docs <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    x_axis_in <- as.numeric(object$processed_wavs[[paste0("step_", i - 1)]])
    translated <- .translate_prep_step(steps[[i]], x_axis_in)
    step_name <- paste0("step", i, "_", gsub("^prep_", "", steps[[i]]$method))
    step_docs[[i]] <- list(step_name, translated)
  }

  model_step <- list(
    "model",
    list(
      estimator_class = "NIRWiseLinearModel",
      params = list(
        fit_method = model$method$fit_method,
        type = model$method$type,
        ncomp = object$final_ncomp,
        min_w = model$method$min_w,
        max_w = model$method$max_w
      ),
      attributes = list(
        x_means_ = unname(model$x_means),
        coef_ = unname(model$coefficients[object$final_ncomp, ]),
        intercept_ = unname(model$intercept)[1],
        n_features_in_ = length(model$x_means),
        feature_names_in_ = object$predictor_variables
      ),
      attribute_types = list(
        x_means_ = "ndarray",
        coef_ = "ndarray",
        intercept_ = "float",
        n_features_in_ = "int",
        feature_names_in_ = "ndarray"
      ),
      attribute_dtypes = list(
        x_means_ = "float64",
        coef_ = "float64",
        feature_names_in_ = "object"
      )
    )
  )

  all_steps <- c(step_docs, list(model_step))
  # openmodels' generic deserializer only recurses into a nested (name, estimator)
  # pair when its type is tagged as a *list* of ("str", <EstimatorClass>) pairs in
  # param_types -- this must be supplied explicitly here since these dicts are
  # hand-authored in R rather than produced by openmodels' own Python-side
  # introspection (which derives it automatically from Pipeline.get_params()).
  step_types <- lapply(all_steps, function(s) list("str", s[[2]]$estimator_class))

  doc <- list(
    estimator_class = "Pipeline",
    params = list(
      steps = all_steps,
      memory = NULL,
      verbose = FALSE
    ),
    param_types = list(
      steps = step_types,
      memory = "NoneType",
      verbose = "bool"
    ),
    metadata = list(
      producer_name = "sklearn",
      producers = list(
        sklearn = "unknown",
        chemotools = "unknown",
        proximetricsr_estimators = "unknown"
      ),
      domain = "sklearn",
      openmodels_format_version = 2L,
      openmodels_version = "unknown",
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      dependency_versions = list(numpy = "unknown", scipy = "unknown"),
      source = "proximetricsR",
      proximetricsR_version = as.character(utils::packageVersion("proximetricsR")),
      r_version = R.version.string
    )
  )

  wired <- .to_wire(doc)
  json <- toJSON(wired, auto_unbox = TRUE, null = "null", na = "null", digits = NA)

  if (!is.null(file)) {
    writeLines(json, con = file)
    return(invisible(json))
  }
  json
}

#' @title Translate one proximetricsR preprocessing step to a chemotools estimator dict
#' @description internal function used by \code{\link{export_sklearn_model}}
#' @param step an object of class \code{preprocessing}, one entry of
#' \code{object$preprocess$steps}.
#' @param x_axis_in numeric vector of the wavelength/wavenumber grid entering this
#' step (\code{object$processed_wavs[[paste0("step_", i - 1)]]}), used only by
#' steps whose chemotools equivalent needs an explicit x-axis (currently
#' \code{prep_wav_trim}).
#' @return A list with \code{estimator_class} and \code{params} (and, where
#' relevant, \code{param_types}/\code{param_dtypes}) describing the equivalent
#' chemotools transformer.
#' @keywords internal
.translate_prep_step <- function(step, x_axis_in) {
  n_features_in <- length(x_axis_in)

  switch(step$method,
    prep_snv = list(
      estimator_class = "StandardNormalVariate",
      params = list(),
      # StandardNormalVariate.fit() sets no fitted state beyond n_features_in_
      # (mean/std are computed per-row at transform time, not stored).
      attributes = list(n_features_in_ = n_features_in),
      attribute_types = list(n_features_in_ = "int")
    ),
    prep_derivative = {
      if (identical(step$algorithm, "nwp")) {
        stop(
          "export_sklearn_model() does not support prep_derivative(algorithm = \"nwp\") ",
          "(no chemotools equivalent). Use save_spectral_model()/load_spectral_model() ",
          "for a full-fidelity R-native export instead."
        )
      } else if (identical(step$algorithm, "savitzky-golay")) {
        # SavitzkyGolay.fit() only mirrors window_length/polyorder/deriv into
        # window_length_/polyorder_/deriv_ (no other computed state) + n_features_in_.
        list(
          estimator_class = "SavitzkyGolay",
          params = list(window_length = step$w, polyorder = step$p, deriv = step$m),
          attributes = list(
            window_length_ = step$w, polyorder_ = step$p, deriv_ = step$m,
            n_features_in_ = n_features_in
          ),
          attribute_types = list(
            window_length_ = "int", polyorder_ = "int", deriv_ = "int",
            n_features_in_ = "int"
          )
        )
      } else if (identical(step$algorithm, "gap-segment")) {
        # Not yet supported: chemotools.derivative.NorrisWilliams precomputes a
        # kernel_ (a convolution of an internal smoothing kernel and a derivative
        # kernel, see NorrisWilliams.fit()) that isn't just its constructor params
        # mirrored with a trailing underscore -- unlike every other step mapped
        # here, replicating it correctly needs verification against chemotools'
        # exact internal kernel formulas that hasn't been done yet.
        stop(
          "export_sklearn_model() does not yet support ",
          "prep_derivative(algorithm = \"gap-segment\") (its chemotools equivalent, ",
          "NorrisWilliams, precomputes an internal kernel not yet replicated here). ",
          "Use save_spectral_model()/load_spectral_model() for a full-fidelity ",
          "R-native export instead."
        )
      } else {
        stop("Unsupported prep_derivative algorithm '", step$algorithm, "'.")
      }
    },
    prep_smooth = {
      if (identical(step$algorithm, "savitzky-golay")) {
        # chemotools.smooth.SavitzkyGolayFilter precomputes a convolution kernel via
        # scipy.signal.savgol_coeffs(window_length, polyorder, deriv=0, use="conv"),
        # rather than calling scipy's savgol_filter at transform time. We reuse the
        # existing sgf() helper (see proxiscout_write_model.R, which already emits
        # Savitzky-Golay coefficients for the ProxiScout device JSON format), which
        # implements the same standard pseudo-inverse-of-Vandermonde-matrix formula.
        # Cross-checked numerically (transliterating sgf() to Python/numpy and
        # comparing against scipy.signal.savgol_coeffs(..., use="conv") for several
        # (window, polyorder) pairs): the m=0 (smoothing) coefficient vector is
        # symmetric, so it is identical to its own reversal and matches scipy's
        # "conv" kernel exactly -- unlike derivative orders (m>0), there is no
        # conv-vs-dot ordering ambiguity to get wrong for this particular mapping.
        # The rev() call below is therefore a no-op for m=0, kept only so this
        # matches the general pattern (and would be needed if this were ever
        # extended to export a nonzero-derivative smoothing kernel).
        kernel <- rev(as.vector(sgf(p = step$p, n = step$w, m = 0)))
        list(
          estimator_class = "SavitzkyGolayFilter",
          params = list(window_length = step$w, polyorder = step$p),
          attributes = list(
            window_length_ = step$w,
            polyorder_ = step$p,
            kernel_ = kernel,
            # _half_ = (window_length_-1)%/%2 is a *private* attribute
            # chemotools._BaseFIRFilter.transform() reads directly at transform
            # time (chemotools/smooth/_base.py::_apply_filter_1d). openmodels'
            # own attribute extraction skips leading-underscore attributes, so a
            # plain Python-native round trip of a real SavitzkyGolayFilter would
            # actually drop this too -- it only works here because we set it
            # directly, bypassing that extraction step.
            `_half_` = (step$w - 1L) %/% 2L,
            n_features_in_ = n_features_in
          ),
          attribute_types = list(
            window_length_ = "int", polyorder_ = "int",
            kernel_ = "ndarray", `_half_` = "int", n_features_in_ = "int"
          ),
          attribute_dtypes = list(kernel_ = "float64")
        )
      } else if (identical(step$algorithm, "moving-average")) {
        # MeanFilter.fit() only mirrors window_length into window_length_ (the mean
        # is computed at transform time via scipy's uniform_filter1d) + n_features_in_.
        list(
          estimator_class = "MeanFilter",
          params = list(window_length = step$w),
          attributes = list(window_length_ = step$w, n_features_in_ = n_features_in),
          attribute_types = list(window_length_ = "int", n_features_in_ = "int")
        )
      } else {
        stop("Unsupported prep_smooth algorithm '", step$algorithm, "'.")
      }
    },
    prep_detrend = list(
      estimator_class = "PolynomialCorrection",
      params = list(order = step$p, indices = NULL),
      # indices = NULL fits the polynomial to every point (matching prospectr::
      # detrend's whole-spectrum behaviour), i.e. indices_ = 0:(n_features_in_-1)
      # (0-based, matching PolynomialCorrection.fit()'s own `list(range(0, len(X[0])))`).
      attributes = list(
        indices_ = seq(0L, n_features_in - 1L),
        n_features_in_ = n_features_in
      ),
      attribute_types = list(indices_ = "ndarray", n_features_in_ = "int"),
      attribute_dtypes = list(indices_ = "int32")
    ),
    prep_transform = {
      if (!identical(step$to, "absorbance")) {
        stop(
          "export_sklearn_model() only supports prep_transform(to = \"absorbance\"); ",
          "'to = \"reflectance\"' has no supported chemotools export path (mirrors ",
          "proxiscout_write_model()'s existing behaviour)."
        )
      }
      # proximetricsR's "absorbance" here is A = -log10(R) computed directly from
      # reflectance: chemotools calls this conversion pair "pseudoabsorbance"
      # (reflectance-based), distinct from its "absorbance"/"transmittance" pair
      # (Beer-Lambert, transmittance-based) -- see export_sklearn_model()'s details.
      # IntensityConversion.fit() sets no fitted state beyond n_features_in_.
      list(
        estimator_class = "IntensityConversion",
        params = list(input_unit = "reflectance", output_unit = "pseudoabsorbance"),
        attributes = list(n_features_in_ = n_features_in),
        attribute_types = list(n_features_in_ = "int")
      )
    },
    prep_wav_trim = {
      if (isTRUE(step$trim_constant_edges)) {
        stop(
          "export_sklearn_model() does not support ",
          "prep_wav_trim(trim_constant_edges = TRUE) (data-dependent, no chemotools ",
          "equivalent)."
        )
      }
      if (length(step$band) == 0) {
        stop("export_sklearn_model() requires a non-empty 'band' in prep_wav_trim().")
      }
      # chemotools.feature_selection.RangeCut resolves start/end to indices via
      # nearest-value lookup (chemotools._axis_mixin.XAxisMixin._find_index:
      # argmin(abs(axis - target))) and stores x_axis_/wavenumbers_ as
      # x_axis[start_index_:end_index_] (Python half-open slice) -- replicated
      # here rather than deferring index resolution to a (nonexistent) R-side fit().
      start_index <- which.min(abs(x_axis_in - min(step$band))) - 1L
      end_index <- which.min(abs(x_axis_in - max(step$band))) - 1L
      selected_axis <- x_axis_in[(start_index + 1L):end_index]
      list(
        estimator_class = "RangeCut",
        params = list(
          start = min(step$band),
          end = max(step$band),
          x_axis = x_axis_in
        ),
        param_types = list(x_axis = "ndarray"),
        param_dtypes = list(x_axis = "float64"),
        attributes = list(
          start_index_ = start_index,
          end_index_ = end_index,
          x_axis_ = selected_axis,
          wavenumbers_ = selected_axis,
          n_features_in_ = n_features_in
        ),
        attribute_types = list(
          start_index_ = "int", end_index_ = "int",
          x_axis_ = "ndarray", wavenumbers_ = "ndarray", n_features_in_ = "int"
        ),
        attribute_dtypes = list(x_axis_ = "float64", wavenumbers_ = "float64")
      )
    },
    prep_resample = stop(
      "export_sklearn_model() does not support prep_resample(): its chemotools ",
      "equivalent (XAxisInterpolator) requires scikit-learn metadata routing at ",
      "predict time (the incoming spectra's x-axis must be passed explicitly to ",
      "every predict() call), unlike every other step here, so it is not usable in ",
      "a plain Pipeline.predict(X) workflow. Use save_spectral_model()/",
      "load_spectral_model() for a full-fidelity R-native export instead."
    ),
    stop("Unsupported preprocessing step '", step$method, "' for export_sklearn_model().")
  )
}
