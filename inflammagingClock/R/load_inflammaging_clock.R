#' Load the bundled inflammaging clock
#'
#' @return The clock object stored in `inflammaging_clock.rds`.
#' @export
#' @examples clock <- load_inflammaging_clock()
#' \dontrun{
#' clock <- load_inflammaging_clock()
#' }
load_inflammaging_clock <- function() {
  model_path <- system.file(
    "extdata", "inflammaging_clock.rds",
    package = "inflammagingClock"
  )

  if (!nzchar(model_path)) {
    stop(
      "Bundled model not found. Add inst/extdata/inflammaging_clock.rds",
      "to the package source and reinstall the package.",
      call. = FALSE
    )
  }

  readRDS(model_path)
}
