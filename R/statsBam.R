#' Count aligned reads in a bam file.
#'
#' Convenience wrapper for a few
#'
#' When primary_align_only is set to FALSE and no extra arguments
#' are given this will use an index if one exists.
#'
#' @return A list of integer values corresponding to counts.
#' @import Rsamtools
#' @export
reads_aligned <- function(bamfiles, idx, primary_align_only = F,
                          return_tibble = T, millions = F, ...) {
  if (!is.null(names(bamfiles))) names(bamfiles) <- names(bamfiles)
  extra_args <- list(...)
  use_extra <- (length(extra_args) > 0)
  bai_files_exist <- lapply(bamfiles, paste0, ".bai") %>%
    lapply(file.exists) %>% unlist %>% all
  if (primary_align_only | use_extra) {
    if (!primary_align_only) primary_align_only <- NA # return all aln
    bv <- BamViews(bamfiles)
    bf <- scanBamFlag(isUnmappedQuery = F,
                      isSecondaryAlignment = !primary_align_only,
                      ...)
    pm <- ScanBamParam(bf)
    res <- lapply(countBam(bv, param = pm), `[[`, "records")
  } else {
    if (!bai_files_exist) stop("Please index your bams.")
    res <- lapply(bamfiles,
           FUN = function(x) sum(idxstatsBam(x)$mapped))
  }
  if (!return_tibble) {
    return(ifelse(millions, lapply(res, FUN = function(x) x/1e6), res) %>%
             magrittr::set_names(names(res)))
  }

  tibble::tibble(Name = names(res), mapped = unlist(res)) %>%
    dplyr::mutate(mapped = ifelse(millions, mapped/1e6, mapped))
}
