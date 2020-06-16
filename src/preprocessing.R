
unify_row_data <- function(row_data_list) {
  do.call(rbind, row_data_list) %>%
    as_tibble() %>%
    unique() %>%
    arrange(channel_name, marker_name)
}

unify_col_data <- function(col_data_list) {
  col_union <- lapply(col_data_list, function(x) as_tibble(colData(x)))
  bind_rows(col_union, .id="cell_type") %>%
    select(sample_id, starts_with("patient"), everything()) %>%
    select(-file_name) %>%
    mutate_at(vars(matches("Age|percent|Score")), as.numeric) %>%
    mutate_at(vars(-matches("Age|percent|Score")), as_factor)
}

subsample_experiments <- function(x_list, p_keep = 0.05) {
  for (i in seq_along(x_list)) {
    D <- ncol(x_list[[i]])
    sample_ix <- sample(D, D * p_keep, replace=FALSE)
    x_list[[i]] <- x_list[[i]][, sample_ix]
  }
  x_list
}

data_list <- function(pattern) {
  global <- ls(envir=.GlobalEnv)
  cell_types <- global[grep(pattern, global)]
  x <- lapply(cell_types, get)
  names(x) <- cell_types
  x
}

quantile_transform <- function(exper, max_rank = 2000) {
  x <- assay(exper)
  for (j in seq_len(nrow(exper))) {
    x[j, ] <- pmin(rank(x[j, ], ties = "random"), max_rank) / max_rank
  }

  assay(exper) <- x
  exper
}

polygonize <- function(im) {
  polys <- stars::st_as_stars(im) %>%
    sf::st_as_sf(merge = TRUE) %>%
    sf::st_cast("POLYGON")

  colnames(polys)[1] <- "cellLabelInImage"
  polys %>%
    mutate(geometry = sf::st_buffer(geometry, dist = 0)) %>%
    group_by(cellLabelInImage) %>%
    summarise(n_polys = n(), .groups = "drop") %>%
    dplyr::select(-n_polys)
}

backgroundProp <- function(x, ...) {
  if (nrow(x) == 0) { # case of no neighbors
    return (tibble(immuneGroup = NA, props = NA))
  }

  props <- table(x$cellLabelInImage %in% c(0, 1))
  tibble(background = names(props), props = props / sum(props))
}

typeProps <- function(x, ...) {
  if (nrow(x) == 0) { # case of no neighbors
    return (tibble(cell_type = NA, props = NA))
  }

  props <- table(x$cell_type, useNA = "ifany")
  tibble(cellType = names(props), props = as.numeric(props) / sum(props)) %>%
    filter(props != 0)
}

cell_type <- function(exper) {
  colData(exper) %>%
    as.data.frame() %>%
    dplyr::select(tumor_group, immune_group) %>%
    mutate(
      cell_type = paste0(tumor_group, immune_group),
      cell_type = gsub("not immune", "", cell_type),
      cell_type = gsub("Immune", "", cell_type),
      ) %>%
    .[["cell_type"]] %>%
    as_factor()
}

#' Apply fun to Graph Neighborhoods
#'
#' @param cell_id The ID of the cell to extract a local neighborhood around.
#' @param G The graph object giving the connections between cell_ids.
#' @param polys A spatial data.frame with a column (geometry) giving the spatial
#'   geometry of each cell.
#' @param fun A function that can be applied to a data.frame whose rows are
#'   pixels and whose columns give features of those pixels (e.g., immune
#'   group).
#' @return result A tibble mapping the cell to statistics calculated by fun.
graph_stats_cell <- function(cell_id, G, polys, fun, ...) {
  ball <- igraph::neighbors(G, as.character(cell_id))
  cell_stats <- polys %>%
    filter(cellLabelInImage %in% names(ball)) %>%
    group_map(fun)

  cell_stats[[1]] %>%
    dplyr::mutate(cellLabelInImage = cell_id) %>%
    dplyr::select(cellLabelInImage, everything())
}

extract_graph <- function(geometries, K = 5) {
  nb <- spdep::knn2nb(
    spdep::knearneigh(sf::st_centroid(geometries[1:nrow(geometries), ]), K)
  )
  labels <- unique(geometries$cellLabelInImage)
  dists <- sapply(geometries$geometry, sf::st_centroid) %>%
    t()

  relations_data <- list()
  for (i in seq_along(nb)) {
    relations_data[[i]] <- tibble(
      from = labels[i],
      to = labels[nb[[i]]]
    )

    relations_data[[i]]$dist <- pracma::distmat(
      dists[i, ], dists[nb[[i]], ]
    ) %>%
      as.numeric()
  }

  relations_data <- bind_rows(relations_data)
  igraph::graph_from_data_frame(relations_data, labels)
}

#' Apply fun to Local Neighborhoods
#'
#' @param cell_id The ID of the cell to extract a local neighborhood around.
#' @param im The raster object giving the pixel-level information about the
#'   sample.
#' @param polys A spatial data.frame with a column (geometry) giving the spatial
#'   geometry of each cell.
#' @param fun A function that can be applied to a data.frame whose rows are
#'   pixels and whose columns give features of those pixels (e.g., immune
#'   group).
#' @param buffer_radius The size of the window around cell_id, to use to subset
#'   the raster on which to apply fun.
#' @param plot_masks If you want to see what the subsets of cells looks like,
#'   you can use this.
#' @return result A tibble mapping the cell to statistics calculated by fun.
raster_stats_cell <- function(cell_id, im, polys, fun, buffer_radius=90,
                              plot_masks=TRUE) {
  sub_poly <- polys %>%
    filter(cellLabelInImage == cell_id) %>%
    .[["geometry"]] %>%
    sf::st_centroid() %>%
    sf::st_buffer(dist=buffer_radius)

  im_ <- mask(im, as_Spatial(sub_poly))
  if (plot_masks) {
    plot(im_)
  }

  melted_im <- as.matrix(im_) %>%
    melt(na.rm=TRUE, value.name = "cellLabelInImage") %>%
    left_join(polys, by = "cellLabelInImage") %>%
    group_map(fun)

  melted_im[[1]] %>%
    mutate(cellLabelInImage = cell_id) %>%
    dplyr::select(cellLabelInImage, everything())
}

#' Wrapper for Local Statistics
#'
#' @param cell_ids A vector of cell IDs on which to apply a function to
#' @param type Either "raster" or "graph". Specifies the types of neighborhoods
#'   (image or graph) on which to compute statistics.
loop_stats <- function(cell_ids, type="raster", ...) {
  cell_fun <- ifelse(type == "raster", raster_stats_cell, graph_stats_cell)

  result <- list()
  for (i in seq_along(cell_ids)) {
    result[[i]] <- cell_fun(cell_ids[i], ...)
  }

  bind_rows(result)
}

generate_model <- function(n_ft) {
  keras_model_sequential() %>%
    layer_dense(units = 32, input_shape = n_ft) %>%
    layer_activation('relu') %>%
    layer_dense(units = 32, input_shape = 32) %>%
    layer_activation('relu') %>%
    layer_dense(units = 32, input_shape = 32) %>%
    layer_activation('relu') %>%
    layer_dense(units = 32, input_shape = 32) %>%
    layer_activation('relu') %>%
    layer_dense(units = 32, input_shape = 32) %>%
    layer_activation('relu') %>%
    layer_dense(units = 2) %>%
    compile(optimizer = optimizer_adam(lr=1e-2), loss = "mae")
}

load_mibi <- function(data_dir, n_paths = NULL) {
  load(file.path(data_dir, "mibiSCE.rda"))
  tiff_paths <- list.files(file.path(data_dir, "TNBC_shareCellData"), "*.tif*", full=T)

  if (is.null(n_paths)) {
    n_paths = length(tiff_paths)
  }

  tiff_paths <- tiff_paths[1:n_paths]
  sample_names <- stringr::str_extract(tiff_paths, "[0-9]+")
  summary(mibi.sce)
  colData(mibi.sce)$cell_type <- cell_type(mibi.sce)
  list(
    tiffs = tiff_paths,
    mibi = mibi.sce[, colData(mibi.sce)$SampleID %in% sample_names]
  )
}

spatial_subsample <- function(tiff_paths, exper, qsize=500) {
  ims <- list()
  for (i in seq_along(tiff_paths)) {
    print(paste0("cropping ", i, "/", length(tiff_paths)))
    r <- raster::raster(tiff_paths[[i]])
    ims[[i]] <- raster::crop(r, raster::extent(1, qsize, 1, qsize))
  }

  names(ims) <- stringr::str_extract(tiff_paths, "[0-9]+")
  cur_cells <- sapply(ims, raster::unique) %>%
    melt() %>%
    dplyr::rename(cellLabelInImage = "value", SampleID = "L1") %>%
    tidyr::unite(sample_by_cell, SampleID, cellLabelInImage, remove=F)

  scell <- colData(exper) %>%
    as.data.frame() %>%
    dplyr::select(SampleID, cellLabelInImage) %>%
    tidyr::unite(sample_by_cell, SampleID, cellLabelInImage) %>%
    .[["sample_by_cell"]]

  list(
    ims = ims,
    exper = exper[, scell %in% cur_cells$sample_by_cell]
  )
}

sample_proportions <- function(SampleID, cluster) {
  tab <- table(SampleID, cluster)
  props <- tab / rowSums(tab)
  props[hclust(dist(props))$order, ]
}

subgraphs <- function(G, order=3) {
  ids <- V(G)
  SG <- list()
  for (i in seq_along(ids)) {
    ball <- do.call(c, igraph::ego(G, ids[[i]], order=order))
    SG[[i]] <- igraph::induced_subgraph(G, ball)

  }
  names(SG) <- ids

  SG
}

entropies <- function(G, clusters_) {
  ents <- list()
  for (g in seq_along(G)) {
    counts <- table(clusters_[names(V(G[[g]]))])
    fq <- counts / sum(counts)
    ents[[g]] <- -sum(fq * log(fq))
  }

  do.call(c, ents)
}

cluster_props <- function(G, clusters_) {
  props <- list()
  for (g in seq_along(G)) {

    if (length(intersect(names(V(G[[g]])), names(clusters_))) == 0) { # hack for tiny data, not needed in real application
      for (j in 1:5) {
        G[[g]] <- set.vertex.attribute(G[[g]], "name", j, names(clusters_)[j])
      }
    }

    counts <- table(clusters_[names(V(G[[g]]))])
    props[[g]] <- counts / sum(counts)
  }

  props_df <- bind_rows(props) %>%
    mutate_all(function(z) { z[is.na(z)] <- 0; z })
  colnames(props_df) <- paste0("prop_", colnames(props_df))
  props_df
}

avg_dists <- function(G) {
  dists <- list()

  for (g in seq_along(G)) {
    dists[[g]] <- mean(edge_attr(G[[g]], "dist"))
  }

  do.call(c, dists)
}

plot_fits <- function(x, y, glmnet_fit, rf_fit) {
  y_hat <- predict(glmnet_fit, x)
  plot(y, y_hat, ylim = range(y), xlim = range(y))
  abline(0, 1)
  y_hat <- predict(rf_fit, x)
  points(y, y_hat, col = "red")
  data.frame(y = y, y_hat = y_hat)
}

fit_wrapper <- function(x, y) {
  glmnet_fit <- glmnet::cv.glmnet(x, y)
  plot(glmnet_fit)
  rf_fit <- caret::train(x, y)
  print(rf_fit)
  list(rf = rf_fit, glmnet = glmnet_fit)
}

extract_graphs <- function(ims, sample_names) {
  graphs <- list()
  for (i in seq_along(sample_names)) {
    print(sprintf("graph %s/%s", i, length(sample_names)))
    poly <- polygonize(ims[[sample_names[i]]]) %>%
      filter(cellLabelInImage > 1)
    graphs[[i]] <- extract_graph(poly)
  }

  graphs
}

spatial_wrapper <- function(sample_names, graphs, cluster_ids) {
  spatial <- list()
  for (i in seq_along(sample_names)) {
    print(sprintf("spatial stats %s/%s", i, length(sample_names)))
    SG <- subgraphs(graphs[[i]])
    ptrn <- paste0("^", sample_names[i], "_")
    clusters_ <- cluster_ids[grepl(ptrn, names(cluster_ids))]
    names(clusters_) <- gsub(ptrn, "", names(clusters_))

    spatial[[sample_names[i]]] <- tibble(
      scell = paste0(sample_names[i], "_", names(V(graphs[[i]]))),
      entropy = entropies(SG, clusters_),
      avg_dists = avg_dists(SG),
      cluster_props(SG, clusters_)
    )
  }

  bind_rows(spatial) %>%
    mutate_all(function(z) { z[is.na(z)] <- 0; z})
}
