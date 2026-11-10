# Average assays that map to the same HGNC gene.
collapse_assays_to_genes <- function(assays) {
  if (ncol(assays) == 0L || !all(vapply(assays, is.numeric, logical(1)))) {
    stop("Assay columns must be numeric, with at least one assay.")
  }

  gene_names <- toupper(trimws(names(assays)))
  if (anyNA(gene_names) || any(!nzchar(gene_names))) {
    stop("Assay columns must have nonempty HGNC gene names.")
  }

  assay_matrix <- as.matrix(assays)
  genes <- unique(gene_names)
  values <- vapply(genes, function(gene) {
    means <- rowMeans(assay_matrix[, gene_names == gene, drop = FALSE],
                      na.rm = TRUE)
    means[is.nan(means)] <- NA_real_
    means
  }, numeric(nrow(assay_matrix)))

  matrix(values, nrow = nrow(assay_matrix),
         dimnames = list(NULL, genes))
}


# Match training features and fill unavailable values with saved medians.
prepare_clock_input <- function(gene_matrix, clock) {
  features <- clock$features
  if (!is.character(features) || length(features) == 0L ||
      anyNA(features) || any(!nzchar(features)) || anyDuplicated(features)) {
    stop("clock$features must contain unique, nonempty gene names.")
  }

  medians <- clock$training_medians[features]
  if (!is.numeric(medians) || length(medians) != length(features) ||
      any(!is.finite(medians))) {
    stop("Saved training medians must be named and finite for all features.")
  }

  limit <- clock$max_absent_feature_fraction
  if (!is.numeric(limit) || length(limit) != 1L ||
      !is.finite(limit) || limit < 0 || limit > 1) {
    stop("clock$max_absent_feature_fraction must be between 0 and 1.")
  }

  measured <- intersect(features, colnames(gene_matrix))
  absent_fraction <- 1 - length(measured) / length(features)
  if (absent_fraction > limit) {
    stop(sprintf("%.1f%% of clock features are absent (limit: %.1f%%).",
                 100 * absent_fraction, 100 * limit))
  }

  n_samples <- nrow(gene_matrix)
  defaults <- matrix(rep(medians, each = n_samples), nrow = n_samples,
                     dimnames = list(NULL, features))
  model_matrix <- defaults
  imputed <- matrix(TRUE, nrow = n_samples, ncol = length(features),
                    dimnames = list(NULL, features))

  if (length(measured) > 0L) {
    model_matrix[, measured] <- gene_matrix[, measured, drop = FALSE]
    imputed[, measured] <- !is.finite(gene_matrix[, measured, drop = FALSE])
  }
  model_matrix[imputed] <- defaults[imputed]

  list(matrix = model_matrix,
       n_imputed = rowSums(imputed),
       absent_fraction = absent_fraction)
}


#' Apply the inflammaging clock
#'
#' Estimates inflammaging-associated biological age
#' Protein measurements are matched by HGNC gene symbols.
#' Multiple assays corresponding to the same gene are averaged within each sample.
#'Missing values are imputed using medians derived from the training dataset.
#'
#' @param data A data frame with samples in rows and measurements in columns.
#' @param clock A pretrained clock object. By default,
#'   the model included in the package is loaded.
#' @param id_column Optional column containing sample IDs.
#' @param age_column Optional column containing chronological age.
#'
#' @return A data frame containing predicted age, the number
#'   of imputed features, the fraction of absent features,
#'   and age acceleration when chronological age is provided.
#'
#' @export
#' @examples
# clock <- list(
#   features = c("IL6", "TNF"),
#   training_medians = c(IL6 = 1, TNF = 2),
#   max_absent_feature_fraction = 0.5,
#   model_type = "elastic_net",
#   intercept = 30,
#   coefficients = c(IL6 = 2, TNF = 1),
#   age_calibration = c(2, 0.9)
# )
#' \dontrun{
#' clock <- load_inflammaging_clock()
#'
#' results <- inflammaging_clock(
#'   samples,
#'   clock = clock,
#'   id_column = "sample_id",
#'   age_column = "age"
#' )
#' }
inflammaging_clock <- function(data, clock = load_inflammaging_clock(),
                               id_column = NULL, age_column = NULL) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a nonempty data frame.")
  }

  metadata_columns <- c(id_column, age_column)
  if (!all(metadata_columns %in% names(data))) {
    stop("The specified ID or age column is missing from data.")
  }
  if (!is.null(age_column) && !is.numeric(data[[age_column]])) {
    stop("The age column must be numeric.")
  }

  assays <- data[, !(names(data) %in% metadata_columns), drop = FALSE]
  gene_matrix <- collapse_assays_to_genes(assays)
  input <- prepare_clock_input(gene_matrix, clock)

  model_type <- clock$model_type
  if (!is.character(model_type) || length(model_type) != 1L ||
      is.na(model_type)) {
    stop("clock$model_type must be a single model-type string.")
  }

  if (model_type %in% c("ridge", "lasso", "elastic_net")) {
    beta <- clock$coefficients[clock$features]
    if (!is.numeric(beta) || length(beta) != length(clock$features) ||
        any(!is.finite(beta)) || !is.numeric(clock$intercept) ||
        length(clock$intercept) != 1L || !is.finite(clock$intercept)) {
      stop("Linear models require a finite intercept and named coefficients for every feature.")
    }
    predicted_age <- as.numeric(clock$intercept + input$matrix %*% beta)
  } else {
    if (is.null(clock$fit)) {
      stop("The saved model must provide clock$fit for this model type.")
    }
    predicted_age <- as.numeric(stats::predict(clock$fit, newdata = input$matrix))
  }

  if (length(predicted_age) != nrow(data)) {
    stop("The model returned an unexpected number of predictions.")
  }

  ids <- if (is.null(id_column)) seq_len(nrow(data)) else data[[id_column]]
  scores <- data.frame(
    sample_id = ids,
    predicted_age = predicted_age,
    imputed_features = input$n_imputed,
    absent_feature_fraction = input$absent_fraction
  )

  if (!is.null(age_column)) {
    calibration <- clock$age_calibration
    if (!is.numeric(calibration) || length(calibration) != 2L ||
        any(!is.finite(calibration))) {
      stop("clock$age_calibration must contain a finite intercept and slope.")
    }
    expected_age <- calibration[1] + calibration[2] * data[[age_column]]
    scores$age_acceleration <- predicted_age - expected_age
  }

  scores
}
