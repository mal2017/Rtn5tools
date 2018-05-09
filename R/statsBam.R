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
  if (is.null(names(bamfiles))) {
    bamfiles <- as.list(bamfiles)
    names(bamfiles) <- unlist(bamfiles)
  }
  extra_args <- list(...)
  use_extra <- (length(extra_args) > 0)
  bai_files_exist <- lapply(bamfiles, paste0, ".bai") %>%
    lapply(file.exists) %>% unlist %>% all
  # option 1 : use extra flags
  if (primary_align_only | use_extra) {
    if (!primary_align_only) primary_align_only <- NA # return all aln
    bv <- BamViews(bamfiles)
    bf <- scanBamFlag(isUnmappedQuery = F,
                      isSecondaryAlignment = !primary_align_only,
                      ...)
    pm <- ScanBamParam(bf)
    res <- lapply(countBam(bv, param = pm), `[[`, "records")
  } else { # option 2 : take all alignments, no flags
    if (!bai_files_exist) stop("Please index your bams.")
    res <- lapply(bamfiles,
           FUN = function(x) sum(idxstatsBam(x)$mapped))
  }
  if (millions) res <- lapply(res,magrittr::divide_by,1e6)
  if (!return_tibble) return(res)
  tibble(name = names(res), mapped = unlist(res))
}
