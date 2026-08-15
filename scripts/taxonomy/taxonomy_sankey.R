# Taxonomy Sankey/alluvial plotting helpers for MAGOR_16S TaxAnalysis notebooks.
# Uses only packages already required by the notebooks: phyloseq, ggplot2 and scales.

.magor_clean_taxon <- function(x, rank) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | !nzchar(x) |
    grepl("^(na|nan|none|null)$", x, ignore.case = TRUE)
  x[bad] <- paste("Unclassified", tools::toTitleCase(rank))
  x
}

.magor_is_unclassified <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x))) |
    grepl("unclass|unassign|unknown|uncultured", x, ignore.case = TRUE)
}

.magor_bezier_band <- function(x0, x1, source_ymin, source_ymax,
                               target_ymin, target_ymax, n = 36L) {
  t <- seq(0, 1, length.out = n)
  smooth <- 3 * t^2 - 2 * t^3
  x <- x0 + (x1 - x0) * t
  upper <- source_ymax + (target_ymax - source_ymax) * smooth
  lower <- source_ymin + (target_ymin - source_ymin) * smooth

  data.frame(
    x = c(x, rev(x)),
    y = c(upper, rev(lower))
  )
}

magor_taxonomy_sankey <- function(
    ps,
    ranks = c("phylum", "class", "order", "family", "genus", "species"),
    top_n = 20L,
    include_unclassified = FALSE,
    min_label_fraction = 0.002,
    title = NULL
) {
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("The scales package is required.")
  }

  tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps))) {
    otu <- t(otu)
  }

  weights <- rowSums(otu, na.rm = TRUE)
  shared_taxa <- intersect(rownames(tax), names(weights))
  if (!length(shared_taxa)) {
    stop("No shared taxa were found between the taxonomy and abundance tables.")
  }

  tax <- tax[shared_taxa, , drop = FALSE]
  weights <- weights[shared_taxa]
  ranks <- ranks[ranks %in% colnames(tax)]
  if (length(ranks) < 2L) {
    stop("At least two requested taxonomic ranks must be present.")
  }

  for (rank in ranks) {
    tax[[rank]] <- .magor_clean_taxon(tax[[rank]], rank)
  }

  keep <- is.finite(weights) & weights > 0
  if (!include_unclassified) {
    keep <- keep & !.magor_is_unclassified(tax[[ranks[1L]]])
  }
  tax <- tax[keep, , drop = FALSE]
  weights <- weights[keep]
  if (!length(weights) || sum(weights) <= 0) {
    stop("No positive classified abundance remained for the Sankey plot.")
  }

  top_taxa <- vector("list", length(ranks))
  names(top_taxa) <- ranks

  for (rank in ranks) {
    totals <- tapply(weights, tax[[rank]], sum, na.rm = TRUE)
    totals <- sort(totals, decreasing = TRUE)
    selected <- names(totals)[seq_len(min(as.integer(top_n), length(totals)))]
    top_taxa[[rank]] <- selected

    other_label <- paste("Other", tools::toTitleCase(rank))
    tax[[rank]][!tax[[rank]] %in% selected] <- other_label
  }

  path_data <- tax[, ranks, drop = FALSE]
  path_data$weight <- as.numeric(weights)
  paths <- stats::aggregate(weight ~ ., data = path_data, FUN = sum)
  paths <- paths[is.finite(paths$weight) & paths$weight > 0, , drop = FALSE]
  paths$fraction <- paths$weight / sum(paths$weight)
  paths$path_id <- seq_len(nrow(paths))

  first_rank <- ranks[1L]
  first_rank_totals <- stats::aggregate(
    fraction ~ fill_group,
    data = data.frame(
      fill_group = paths[[first_rank]],
      fraction = paths$fraction,
      stringsAsFactors = FALSE
    ),
    FUN = sum
  )
  first_rank_totals <- first_rank_totals[
    order(-first_rank_totals$fraction, first_rank_totals$fill_group),
    , drop = FALSE
  ]
  fill_groups <- first_rank_totals$fill_group
  fill_values <- grDevices::hcl.colors(length(fill_groups), palette = "Dark 3")
  names(fill_values) <- fill_groups
  fill_values[grepl("^Other ", names(fill_values))] <- "#B8B8B8"
  fill_values[grepl("^Unclassified ", names(fill_values))] <- "#E0E0E0"

  node_list <- lapply(seq_along(ranks), function(rank_index) {
    rank <- ranks[rank_index]
    node_data <- stats::aggregate(
      fraction ~ taxon,
      data = data.frame(
        taxon = paths[[rank]],
        fraction = paths$fraction,
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
    node_data <- node_data[
      order(-node_data$fraction, node_data$taxon),
      , drop = FALSE
    ]
    node_data$ymax <- cumsum(node_data$fraction)
    node_data$ymin <- node_data$ymax - node_data$fraction

    dominant_data <- stats::aggregate(
      fraction ~ taxon + fill_group,
      data = data.frame(
        taxon = paths[[rank]],
        fill_group = paths[[first_rank]],
        fraction = paths$fraction,
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
    dominant_data <- dominant_data[
      order(dominant_data$taxon, -dominant_data$fraction, dominant_data$fill_group),
      , drop = FALSE
    ]
    dominant_data <- dominant_data[!duplicated(dominant_data$taxon), , drop = FALSE]

    node_data$fill_group <- dominant_data$fill_group[
      match(node_data$taxon, dominant_data$taxon)
    ]
    node_data$rank <- rank
    node_data$rank_index <- rank_index
    node_data$xmin <- rank_index - 0.11
    node_data$xmax <- rank_index + 0.11
    node_data$label <- ifelse(
      node_data$fraction >= min_label_fraction,
      node_data$taxon,
      ""
    )
    node_data
  })
  nodes <- do.call(rbind, node_list)
  rownames(nodes) <- NULL

  lookup_ymin <- setNames(nodes$ymin, paste(nodes$rank, nodes$taxon, sep = "\r"))

  flow_polygon_list <- list()
  flow_summary_list <- list()
  polygon_index <- 1L

  for (rank_index in seq_len(length(ranks) - 1L)) {
    source_rank <- ranks[rank_index]
    target_rank <- ranks[rank_index + 1L]

    links <- stats::aggregate(
      fraction ~ source + target + fill_group,
      data = data.frame(
        source = paths[[source_rank]],
        target = paths[[target_rank]],
        fill_group = paths[[first_rank]],
        fraction = paths$fraction,
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
    links$link_id <- seq_len(nrow(links))

    source_order <- links[
      order(links$source, links$target, links$fill_group, links$link_id),
      , drop = FALSE
    ]
    source_order$source_offset <- ave(
      source_order$fraction,
      source_order$source,
      FUN = function(z) cumsum(z) - z
    )
    source_base <- lookup_ymin[
      paste(source_rank, source_order$source, sep = "\r")
    ]
    source_order$source_ymin <- as.numeric(source_base) + source_order$source_offset
    source_order$source_ymax <- source_order$source_ymin + source_order$fraction

    target_order <- links[
      order(links$target, links$source, links$fill_group, links$link_id),
      , drop = FALSE
    ]
    target_order$target_offset <- ave(
      target_order$fraction,
      target_order$target,
      FUN = function(z) cumsum(z) - z
    )
    target_base <- lookup_ymin[
      paste(target_rank, target_order$target, sep = "\r")
    ]
    target_order$target_ymin <- as.numeric(target_base) + target_order$target_offset
    target_order$target_ymax <- target_order$target_ymin + target_order$fraction

    links$source_ymin <- source_order$source_ymin[
      match(links$link_id, source_order$link_id)
    ]
    links$source_ymax <- source_order$source_ymax[
      match(links$link_id, source_order$link_id)
    ]
    links$target_ymin <- target_order$target_ymin[
      match(links$link_id, target_order$link_id)
    ]
    links$target_ymax <- target_order$target_ymax[
      match(links$link_id, target_order$link_id)
    ]
    links$source_rank <- source_rank
    links$target_rank <- target_rank
    links$rank_index <- rank_index
    flow_summary_list[[rank_index]] <- links

    for (link_row in seq_len(nrow(links))) {
      polygon <- .magor_bezier_band(
        x0 = rank_index + 0.11,
        x1 = rank_index + 1L - 0.11,
        source_ymin = links$source_ymin[link_row],
        source_ymax = links$source_ymax[link_row],
        target_ymin = links$target_ymin[link_row],
        target_ymax = links$target_ymax[link_row]
      )
      polygon$polygon_id <- polygon_index
      polygon$fill_group <- links$fill_group[link_row]
      flow_polygon_list[[polygon_index]] <- polygon
      polygon_index <- polygon_index + 1L
    }
  }

  flow_polygons <- do.call(rbind, flow_polygon_list)
  flow_summary <- do.call(rbind, flow_summary_list)
  rownames(flow_polygons) <- NULL
  rownames(flow_summary) <- NULL

  if (is.null(title)) {
    title <- sprintf("Taxonomic Sankey: top %d taxa at each rank", as.integer(top_n))
  }

  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = flow_polygons,
      ggplot2::aes(x = x, y = y, group = polygon_id, fill = fill_group),
      alpha = 0.48,
      colour = NA
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill_group
      ),
      colour = "white",
      linewidth = 0.22,
      alpha = 0.96
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(
        x = rank_index,
        y = (ymin + ymax) / 2,
        label = label
      ),
      size = 2.55,
      lineheight = 0.9,
      check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(ranks),
      labels = tools::toTitleCase(ranks),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Top ", as.integer(top_n),
        " taxa are retained independently at every rank; remaining lineages are collapsed into rank-specific Other nodes."
      ),
      x = NULL,
      y = "Total relative abundance",
      fill = tools::toTitleCase(first_rank)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.45, "lines"),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(10, 25, 10, 10)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )

  list(
    plot = plot,
    paths = paths,
    nodes = nodes,
    links = flow_summary,
    top_taxa = top_taxa,
    palette = fill_values
  )
}

save_magor_taxonomy_sankey <- function(
    ps,
    output_dir,
    file_prefix = "taxonomy_sankey_top20",
    ranks = c("phylum", "class", "order", "family", "genus", "species"),
    top_n = 20L,
    include_unclassified = FALSE,
    min_label_fraction = 0.002,
    title = NULL,
    width = 24,
    height = 14,
    dpi = 300
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  result <- magor_taxonomy_sankey(
    ps = ps,
    ranks = ranks,
    top_n = top_n,
    include_unclassified = include_unclassified,
    min_label_fraction = min_label_fraction,
    title = title
  )

  png_file <- file.path(output_dir, paste0(file_prefix, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_prefix, ".pdf"))
  path_file <- file.path(output_dir, paste0(file_prefix, "_paths.tsv"))
  node_file <- file.path(output_dir, paste0(file_prefix, "_nodes.tsv"))
  link_file <- file.path(output_dir, paste0(file_prefix, "_links.tsv"))

  ggplot2::ggsave(
    filename = png_file,
    plot = result$plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = pdf_file,
    plot = result$plot,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )

  utils::write.table(
    result$paths,
    path_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    result$nodes,
    node_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    result$links,
    link_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  result$files <- list(
    png = png_file,
    pdf = pdf_file,
    paths = path_file,
    nodes = node_file,
    links = link_file
  )
  result
}
