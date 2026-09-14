#' @title Save and load a spectral_model as a portable JSON file
#' @name spectral_model_json
#' @aliases save_spectral_model
#' @aliases load_spectral_model
#'
#' @description
#'
#' \loadmathjax
#'
#' \code{save_spectral_model} writes a fitted \code{\link{spectral_model}} object to
#' a portable, versioned, human-readable JSON file. \code{load_spectral_model} reads
#' it back into an object usable with \code{\link{predict.spectral_model}}.
#'
#' @usage
#' save_spectral_model(object, file, metadata = NULL)
#'
#' load_spectral_model(file)
#'
#' @param object an object of class \code{spectral_model}, as returned by
#' \code{\link{calibrate}}.
#' @param file a string with the path to the JSON file to write to (for
#' \code{save_spectral_model}) or read from (for \code{load_spectral_model}).
#' @param metadata an optional named list of additional free-form metadata to
#' store alongside the model (e.g. \code{list(title = "...", description = "...")}).
#' Merged into the file's \code{metadata} block.
#'
#' @return For \code{save_spectral_model}, invisibly returns \code{NULL}; called for
#' its side effect of writing \code{file}. For \code{load_spectral_model}, an object
#' of class \code{spectral_model} sufficient for \code{\link{predict.spectral_model}}.
#'
#' @details
#' This stores the fitted \code{\link{spectral_fit}} model (\code{object$final_model$model}),
#' its \code{\link{preprocess_recipe}}, the cross-validation tuning grid used to pick
#' the optimal number of components (if any), and identifying fields
#' (\code{target_variable}, \code{predictor_variables}, \code{final_ncomp}, and
#' \code{object$metadata} if present). It intentionally does **not** store the full
#' calibration/cross-validation audit trail (per-fold predictions in \code{model_cv},
#' \code{calibration_statistics_all}, \code{detected_outliers}, \code{initial_fit},
#' \code{input_data}) -- only what is needed to describe and predict from the model.
#'
#' The reloaded object always predicts via \code{\link{predict.spectral_model}}'s
#' non-formula path: pass \code{newdata} to \code{predict()} as a matrix of spectra,
#' or a data.frame with a \code{$spc} matrix column -- the original calibration
#' \code{formula} (if \code{object} was fitted via \code{calibrate.formula}) is not
#' preserved.
#'
#' The JSON shape (\code{estimator_class}/\code{params}/\code{attributes}/
#' \code{metadata}) mirrors the wire format used by the Python
#' \href{https://github.com/Gnpd/openmodels}{openmodels} library, adapted to R (R has
#' no equivalent of numpy's dtype ambiguity, so this format tracks matrix
#' dimensions/dimnames explicitly instead of a parallel type/dtype map).
#'
#' For a JSON export that a Python installation of openmodels/
#' \href{https://github.com/paucablop/chemotools}{chemotools} can load directly as a
#' real \code{sklearn.pipeline.Pipeline}, see \code{\link{export_sklearn_model}}
#' instead. That format is restricted to preprocessing/method combinations that have
#' a Python equivalent (it excludes \code{algorithm = "nwp"}), whereas
#' \code{save_spectral_model}/\code{load_spectral_model} support every proximetricsR
#' model without restriction.
#'
#' @seealso \code{\link{calibrate}}, \code{\link{predict.spectral_model}},
#' \code{\link{export_sklearn_model}}
#'
#' @examples
#' \donttest{
#' data("NIRcannabis")
#' model <- calibrate(CBDA ~ spc,
#'   data = NIRcannabis, preprocess = preprocess_recipe(prep_snv()),
#'   method = fit_plsr(5), control = calibration_control("none"), verbose = FALSE
#' )
#' file <- tempfile(fileext = ".json")
#' save_spectral_model(model, file)
#' reloaded <- load_spectral_model(file)
#' predict(reloaded, newdata = NIRcannabis[1:5, ], verbose = FALSE)
#' }
#' @author Leonardo Ramirez-Lopez
#' @export
save_spectral_model <- function(object, file, metadata = NULL) {
  if (!inherits(object, "spectral_model")) {
    stop("'object' must be of class 'spectral_model'.")
  }
  if (missing(file) || !is.character(file) || length(file) != 1) {
    stop("'file' must be a single character string with the path to write to.")
  }
  if (!is.null(metadata) && !is.list(metadata)) {
    stop("'metadata' must be a list, if provided.")
  }

  model <- object$final_model$model
  if (is.null(model)) {
    stop("'object' does not contain a fitted model (object$final_model$model is NULL).")
  }

  model_class <- switch(model$method$fit_method,
    plsr = "ProximetricsPLS",
    xlsr = "ProximetricsXLS",
    stop("Unknown fit method '", model$method$fit_method, "'.")
  )

  preprocess_params <- lapply(object$preprocess$steps, function(step) {
    list(
      estimator_class = step$method,
      params = step[setdiff(names(step), c("method", "compatible_devices"))]
    )
  })

  model_cv_grid <- object$final_model$model_cv$grid

  doc <- list(
    estimator_class = "ProximetricsSpectralModel",
    params = list(
      target_variable = object$target_variable,
      predictor_variables = object$predictor_variables,
      final_ncomp = object$final_ncomp,
      device = object$preprocess$device,
      preprocess = preprocess_params
    ),
    attributes = list(
      model = list(
        estimator_class = model_class,
        params = model$method[setdiff(names(model$method), "fit_method")],
        attributes = model[setdiff(names(model), "method")]
      ),
      model_cv_grid = model_cv_grid
    ),
    metadata = c(
      list(
        producer_name = "proximetricsR",
        producer_version = as.character(utils::packageVersion("proximetricsR")),
        domain = "proximetricsR",
        format_version = 1L,
        created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        dependency_versions = list(
          R = R.version.string,
          jsonlite = as.character(utils::packageVersion("jsonlite"))
        )
      ),
      if (!is.null(object$metadata)) list(model_metadata = .to_wire(object$metadata)),
      metadata
    )
  )

  wired <- .to_wire(doc)
  json <- toJSON(wired, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  writeLines(json, con = file)
  invisible(NULL)
}

#' @rdname spectral_model_json
#' @export
load_spectral_model <- function(file) {
  if (missing(file) || !is.character(file) || length(file) != 1) {
    stop("'file' must be a single character string with the path to read from.")
  }

  raw <- fromJSON(file, simplifyVector = FALSE)
  doc <- .from_wire(raw)

  if (!identical(doc$estimator_class, "ProximetricsSpectralModel")) {
    stop(
      "'file' does not look like a proximetricsR spectral_model JSON file ",
      "(expected estimator_class 'ProximetricsSpectralModel', got '",
      doc$estimator_class, "')."
    )
  }

  model_doc <- doc$attributes$model
  fit_method <- switch(model_doc$estimator_class,
    ProximetricsPLS = "plsr",
    ProximetricsXLS = "xlsr",
    stop("Unknown model estimator_class '", model_doc$estimator_class, "'.")
  )

  method <- c(list(fit_method = fit_method), model_doc$params)
  class(method) <- c(paste0("fit_", fit_method), "fit_constructor")

  model <- c(list(method = method), model_doc$attributes)
  class(model) <- c("spectral_fit", "list")

  preprocess_steps <- lapply(doc$params$preprocess, function(step) {
    out <- c(list(method = step$estimator_class), step$params)
    class(out) <- c("preprocessing", "list")
    out
  })
  preprocess <- structure(
    list(
      steps = preprocess_steps,
      device = doc$params$device,
      preprocessing_order = paste(
        gsub("^prep_", "", sapply(preprocess_steps, `[[`, "method")),
        collapse = " > "
      )
    ),
    class = c("preprocess_recipe", "list")
  )

  results <- list(
    target_variable = doc$params$target_variable,
    predictor_variables = unlist(doc$params$predictor_variables),
    final_model = list(
      model = model,
      model_cv = if (!is.null(doc$attributes$model_cv_grid)) {
        list(grid = doc$attributes$model_cv_grid)
      } else {
        NULL
      }
    ),
    final_ncomp = doc$params$final_ncomp,
    preprocess = preprocess
  )

  model_metadata <- NULL
  for (item in doc$metadata) {
    if (is.list(item) && !is.null(item$model_metadata)) {
      model_metadata <- item$model_metadata
    }
  }
  if (!is.null(doc$metadata$model_metadata)) {
    model_metadata <- doc$metadata$model_metadata
  }
  results$metadata <- model_metadata

  class(results) <- c("spectral_model", "list")
  results
}

# --- Internal wire-format helpers -------------------------------------------
# .to_wire()/.from_wire() give R matrices and named vectors an explicit,
# self-describing JSON representation (data + dim + dimnames, or names + values)
# rather than relying on jsonlite's default (lossy for dimnames/names, and
# version-dependent) auto-simplification of arrays back into vectors/matrices.
# Plain unnamed scalars/vectors and ordinary named lists round-trip through
# jsonlite without any special handling and are passed through unchanged.

#' @noRd
.to_wire <- function(x) {
  if (is.matrix(x)) {
    list(
      `__r_kind__` = "matrix",
      data = as.vector(x),
      dim = as.integer(dim(x)),
      rownames = rownames(x),
      colnames = colnames(x)
    )
  } else if (is.list(x)) {
    stats::setNames(lapply(x, .to_wire), names(x))
  } else if (!is.null(x) && !is.null(names(x)) && length(x) >= 1) {
    list(
      `__r_kind__` = "named_vector",
      names = names(x),
      values = unname(x)
    )
  } else {
    x
  }
}

#' @noRd
.from_wire <- function(x) {
  if (is.list(x) && !is.null(x[["__r_kind__"]])) {
    kind <- x[["__r_kind__"]]
    if (identical(kind, "matrix")) {
      dims <- unlist(x$dim)
      m <- matrix(unlist(x$data), nrow = dims[1], ncol = dims[2])
      rn <- x$rownames
      cn <- x$colnames
      if (!is.null(rn) || !is.null(cn)) {
        dimnames(m) <- list(
          if (is.null(rn)) NULL else unlist(rn),
          if (is.null(cn)) NULL else unlist(cn)
        )
      }
      return(m)
    } else if (identical(kind, "named_vector")) {
      return(stats::setNames(unlist(x$values), unlist(x$names)))
    }
  }
  if (is.list(x)) {
    return(stats::setNames(lapply(x, .from_wire), names(x)))
  }
  x
}
