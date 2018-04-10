#' find_motif_instances
#'
#' Convenience/wrapper function for finding motifs with motifmatchr.
#'
#' @importFrom GenomicRanges mcols
#' @import motifmatchr
#' @export
find_motif_instances <- function(object, pwms,
                                 genome = 'BSgenome.Hsapiens.UCSC.hg19') {
  tfbs <- matchMotifs(pwms, object,
                      genome = genome,
                      out = 'positions')
  tfbs <- names(tfbs) %>%
    lapply(function(x) {a <- tfbs[[x]];
           mcols(a) <- data.frame(mcols(a), pwm = rep(x, length(a)));
           a}) %>%
    do.call("c",.) %>%
    return
}
