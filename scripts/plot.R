files <- fs::dir_ls("work/utrs_igs_fulls_v4/", glob = "**.fa")
seqs <- purrr::map(files, Biostrings::readDNAStringSet)
file_data <- fs::path_file(files) |>
  stringr::str_split("[.]") |>
  purrr::map(rlang::set_names,
             c("id", "version_segment", "feature_type", "fa")) |>
  purrr::map(as.list) |>
  purrr::map(tibble::as_tibble) |>
  purrr::list_rbind() |>
  tidyr::separate(version_segment, c("version", "segment"), sep = "_") |>
  tidyr::separate(segment, c("start", "end"), convert = TRUE) |>
  dplyr::mutate(length = end - start)

summaries <- file_data |>
  dplyr::filter(feature_type != "junction") |>
  dplyr::summarize(
    n = dplyr::n(),
    total_length = sum(length),
    median_length = median(length),
    mean_length = mean(length),
    label = glue::glue(
      "n = {n}, total length = {total_length}\nmedian length = {median_length}\nmean_length = {prettyNum(mean_length, digits = 2)}"
    ),
    .by = feature_type
  )

file_data |>
  dplyr::filter(feature_type != "junction") |>
  dplyr::mutate(type = forcats::fct_rev(
    ifelse(feature_type != "complete" &
             length > 2000, "removed", "kept")
  )) |>
  ggplot2::ggplot() +
  ggplot2::geom_histogram(ggplot2::aes(length, fill = type)) +
  ggplot2::facet_wrap( ~ feature_type, scales = "free") +
  ggplot2::scale_y_continuous(trans = "pseudo_log") +
  ggplot2::geom_text(
    ggplot2::aes(label = label, x = Inf, y = Inf),
    data = summaries,
    hjust = 1,
    vjust = 1
  ) + 
  ggplot2::theme_bw()

annotation_str <- function(lens) {
  n_oligos <- purrr::map_int(
    lens$length,
    \(x) IRanges::slidingWindows(IRanges::IRanges(1, x), 200, 50) |> _[[1]] |> length()
  ) |> sum()
  glue::glue(
    "n = {nrow(lens)}\nmedian = {median(lens$length)}\nmean = {signif(mean(lens$length), 2)}\nrange = ({min(lens$length)}, {max(lens$length)})\ntotal length = {sum(lens$length)}\ncalculated num oligos = {n_oligos}"
  )
}

lens <- purrr::map_int(seqs, \(seq) length(seq[[1]])) |> tibble::enframe(value = "length")
label <- annotation_str(lens)
lens |>
  ggplot2::ggplot(ggplot2::aes(length)) +
  ggplot2::geom_histogram() +
  ggplot2::labs(y = "count") +
  ggplot2::theme_bw() +
  ggplot2::annotate("text",
                    x = 1500,
                    y = 200,
                    label = label)


lens_filt <- lens |> dplyr::filter(length >= 50)
label <- annotation_str(lens_filt)
lens_filt |>
  ggplot2::ggplot(ggplot2::aes(length)) +
  ggplot2::geom_histogram() +
  ggplot2::labs(y = "count") +
  ggplot2::theme_bw() +
  ggplot2::annotate("text",
                    x = 3000,
                    y = 200,
                    label = label)
