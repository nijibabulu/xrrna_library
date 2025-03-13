seqs <- purrr::map(fs::dir_ls("work/utrs_igs_v3/", glob="**.fa"),
                   Biostrings::readDNAStringSet)

annotation_str <- function(lens) {
  n_oligos <- purrr::map_int(lens$length, \(x) IRanges::slidingWindows(IRanges::IRanges(1, x), 200, 50) |> _[[1]] |> length()) |> sum()
  glue::glue("n = {nrow(lens)}\nmedian = {median(lens$length)}\nmean = {signif(mean(lens$length), 2)}\nrange = ({min(lens$length)}, {max(lens$length)})\ntotal length = {sum(lens$length)}\ncalculated num oligos = {n_oligos}")
}

lens <- purrr::map_int(seqs, \(seq) length(seq[[1]])) |> tibble::enframe( value = "length")
label <- annotation_str(lens)
lens |>
  ggplot2::ggplot(ggplot2::aes(length)) +
  ggplot2::geom_histogram() + 
  ggplot2::labs(y = "count") +
  ggplot2::theme_bw() + 
  ggplot2::annotate("text", x = 1500, y=200, 
                    label = label)


lens_filt <- lens |> dplyr::filter(length >= 50)
label <- annotation_str(lens_filt)
lens_filt |>
  ggplot2::ggplot(ggplot2::aes(length)) +
  ggplot2::geom_histogram() + 
  ggplot2::labs(y = "count") +
  ggplot2::theme_bw() + 
  ggplot2::annotate("text", x = 3000, y=200, 
                    label = label)
