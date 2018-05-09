#' Get fragment sizes for a whole bam or a region
#'
#' Will only work on PE data.
#' @import foreach
#' @import Rsamtools
#' @import GenomicAlignments
#' @export
get_frag_sizes <- function(bamfiles, regions = NULL) {
  if (is.null(names(bamfiles))) {
    bamfiles <- as.list(bamfiles)
    names(bamfiles) <- unlist(bamfiles)
  }
  pm <- ScanBamParam(what= c("isize", "qname", "mate_status"),
                     flag = scanBamFlag(isProperPair = T))
  if (!is.null(regions)) bamWhich(pm) <- regions
  res <- foreach(i = bamfiles) %dopar% {
    scanBam(BamFile(i, asMates = T), param = pm) %>%
      lapply(tibble::as_tibble) %>%
      dplyr::bind_rows(.id = "Locus")
  }
  res %>% magrittr::set_names(names(bamfiles)) %>%
    dplyr::bind_rows(.id = "Bam")
}

## VPLOT!!!!!!

## tornado plots


