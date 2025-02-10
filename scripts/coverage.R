df <- readr::read_csv("work/utrs_igs_v2/virus_stats.csv")

df |>
  dplyr::count(family) |>
  dplyr::mutate(family = forcats::fct_reorder(family, n)) |>
  ggplot2::ggplot(ggplot2::aes(n, family)) + 
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = n), nudge_x = 3) +
  ggplot2::theme_minimal()

df |>
  dplyr::count(order) |>
  dplyr::mutate(order = forcats::fct_reorder(order, n)) |>
  ggplot2::ggplot(ggplot2::aes(n, order)) + 
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = n), nudge_x = 3) +
  ggplot2::theme_minimal()


df |>
  dplyr::count(class) |>
  dplyr::mutate(class = forcats::fct_reorder(class, n)) |>
  ggplot2::ggplot(ggplot2::aes(n, class)) + 
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = n), nudge_x = 10) +
  ggplot2::theme_minimal()

df |>
  dplyr::count(phylum) |>
  dplyr::mutate(phylum = forcats::fct_reorder(phylum, n)) |>
  ggplot2::ggplot(ggplot2::aes(n, phylum)) + 
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = n), nudge_x = 10) +
  ggplot2::theme_minimal()
