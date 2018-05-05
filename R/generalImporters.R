#' Import narrowPeak file
#'
#' #https://charlesjb.github.io/How_to_import_narrowPeak/
#'
#' @param path Character corresponding to path to narrowPeak formatted file.
#' @return A GRanges object.
#' @export
import_narrowpeak <- function(path) {
  cols <- c(signalValue = "numeric",
            pValue = "numeric",
            qValue = "numeric",
            peak = "integer")
  gr <- rtracklayer::import(path,
                            format = "bed",
                            extraCols = cols)
  gr
}
