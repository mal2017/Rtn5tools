#' tn5_cuts
#'
#' Find tn5 insertions in a single region in a single sample.
#'
#' @export
#' @importFrom Rsamtools ScanBamParam
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import dplyr
tn5_cuts <- function(obj, bamfile, bai = bamfile) {
  stopifnot(is(obj, 'GRanges') & length(obj) == 1)
  seqlevels(obj) <- seqlevelsInUse(obj)
  reads <- readGAlignments(bamfile,
                           index = bai,
                           param = ScanBamParam(which = obj))
  seqlevels(reads) <- seqlevelsInUse(reads)
  stopifnot(reads@strand@values %in% c("+","-"))
  readgr <- GRanges(reads)
  plus <- subset(readgr, readgr@strand@values == "+") %>%
    shift(4) %>% resize(1, fix = "start")
  minus <- subset(readgr, readgr@strand@values == "-") %>%
    shift(-5) %>% resize(1, fix ="end")
  cov <- function(x) {
    coverage(x) %>%
      unlist %>%
      .[start(obj):end(obj)] %>%
      as.vector
  }
  plus %<>% subsetByOverlaps(obj, ignore.strand = T)
  minus %<>% subsetByOverlaps(obj, ignore.strand = T)

  pos <- (1:width(obj)) - ceiling(width(obj)/2)

  if (length(plus) == 0) {
    res_plus <- tibble(pos = pos, cuts= rep(0,width(obj)), cut_strand="forward")
  } else {
    res_plus <- cov(plus) %>%
      tibble(pos = pos, cuts = ., cut_strand="forward")
  }
  if (length(minus) == 0) {
    res_minus <- tibble(pos = pos, cuts= rep(0,width(obj)), cut_strand="rev")
  } else {
    res_minus <- cov(minus) %>%
      tibble(pos = pos, cuts = ., cut_strand="rev")
  }
  bind_rows(res_plus, res_minus) %>%
    mutate(region = as.character(obj)) -> ins
  class(ins) %<>% c(.,'tn5_cuts')
  return(ins)
}


#' tn5_cuts_in_regions
#'
#' Find tn5 insertions across multiple regions in a single sample.
#'
#' @export
#' @import GenomicRanges
#' @import foreach
#' @import doParallel
#' @import dplyr
tn5_cuts_in_regions <- function(regions, bamfile, cores = 1, bai = bamfile) {
  stopifnot(is(regions, 'GRanges'))
  if (cores > 1) {
    registerDoParallel(cores)
  }
  stopifnot(file.exists(bamfile))
  stopifnot(file.exists(paste0(bai,'.bai')))
  ins <- foreach(i = 1:length(regions), .combine = bind_rows) %dopar% {tn5_cuts(regions[i], bamfile, bai)}
  class(ins) %<>% c(.,'tn5_cuts')
  return(ins)
}


#' tn5cuts_to_freq
#'
#' From a tibble of cuts, derive a tibble of region-wise
#' cutting frequencies for each strand. Frequency for each strand
#' is with respect to all cuts within each region, regardless of strand.
#'
#' @import dplyr
#' @export
tn5cuts_to_freq <- function(insertions) {
  stopifnot(is(insertions, 'tn5_cuts'))
  insertions %>%
    group_by(region) %>%
    mutate(cut_freq_in_region = cuts/sum(cuts)) %>%
    ungroup %>%
    group_by(pos,cut_strand) %>%
    summarize(mean_cut_freq = mean(cut_freq_in_region)) -> res
  class(res) %<>% c(.,'cuts_to_frequency')
  return(res)
}

#' make a plot of the output from the  tn5cuts_to_freq function
#' @import ggplot2
#' @export
plot_insertion_freq <- function(freq_tbl) {
  stopifnot(is(freq_tbl, 'cuts_to_frequency'))
  ggplot(freq_tbl,aes(x = pos, y = mean_cut_freq, group = cut_strand)) +
    geom_line(aes(color=cut_strand)) +
    scale_color_manual(values = c(forward = "darkorange", rev = "darkgreen"),
                       name = "Strand") +
    theme_classic() +
    ylab("Mean Insertion Frequency") +
    xlab("Position Relative to Center") +
    theme(axis.title.y=element_text(face="italic", color='darkgray'),
          axis.title.x=element_text(face="italic", color='darkgray'),
          legend.text = element_text(colour="darkgray", face='italic'),
          legend.title = element_text(colour="darkgray", face = 'bold'),
          legend.justification = c(1, 1), legend.position = c(1, 1))

}

