# Original code from:
# https://github.com/nicolasdugue/DirectedLouvain

normalize_to_edges <- function(x) {
  # Returns:
  # - edges_num: numeric matrix (u, v, w) with u,v in 1..N, w numeric
  # - nodes: character vector of original node names (order used for mapping)

  # --- igraph ---
  if (inherits(x, "igraph")) {
    vnames <- igraph::V(x)$name
    # Use vertex names if available; otherwise numeric vertex IDs coerced to character
    el <- igraph::as_edgelist(x, names = !is.null(vnames))
    u_char <- as.character(el[,1])
    v_char <- as.character(el[,2])

    w <- igraph::E(x)$weight
    if (is.null(w)) w <- rep(1, nrow(el))

    # Nodes are all unique endpoints, in first-occurrence order
    nodes <- unique(c(u_char, v_char))
    u <- match(u_char, nodes)
    v <- match(v_char, nodes)

    edges_num <- cbind(u, v, as.numeric(w))
    return(list(edges_num = edges_num, nodes = nodes))
  }

  # --- data.frame edge list ---
  if (is.data.frame(x)) {
    if (ncol(x) < 2) stop("Edge list data.frame must have at least two columns (from, to)")
    u_char <- as.character(x[[1]])
    v_char <- as.character(x[[2]])
    w <- if (ncol(x) >= 3) as.numeric(x[[3]]) else rep(1, nrow(x))

    nodes <- unique(c(u_char, v_char))
    u <- match(u_char, nodes)
    v <- match(v_char, nodes)

    edges_num <- cbind(u, v, w)
    return(list(edges_num = edges_num, nodes = nodes))
  }

  # --- adjacency matrix ---
  if (is.matrix(x)) {
    nr <- nrow(x); nc <- ncol(x)
    if (nr != nc) stop("Adjacency matrix must be square")

    rn <- rownames(x); cn <- colnames(x)
    if (!is.null(rn) && !is.null(cn)) {
      if (!setequal(rn, cn)) stop("Row and column names must refer to the same node set")
      # Align columns to row order for consistent indexing
      x <- x[, rn, drop = FALSE]
      nodes <- rn
    } else if (!is.null(rn)) {
      x <- x[, rn, drop = FALSE]; nodes <- rn
    } else if (!is.null(cn)) {
      x <- x[cn, , drop = FALSE]; nodes <- cn
    } else {
      nodes <- as.character(seq_len(nr))
    }

    nz <- which(x != 0, arr.ind = TRUE)
    if (nrow(nz) == 0) {
      edges_num <- matrix(numeric(0), ncol = 3)
      return(list(edges_num = edges_num, nodes = nodes))
    }

    u <- nz[,1]
    v <- nz[,2]
    w <- as.numeric(x[nz])

    edges_num <- cbind(u, v, w)
    return(list(edges_num = edges_num, nodes = nodes))
  }

  stop("Unsupported input type: must be igraph, data.frame edge list, or adjacency matrix (matrix)")
}

directed_louvain <- function(x,
                             precision = 0.0001,
                             gamma = 1.0,
                             # reproducibility = FALSE,
                             renumbering = FALSE,
                             randomized = TRUE,
                             level = "final",
                             verbose = FALSE) {

  norm <- normalize_to_edges(x)

  if (is.null(lvl)) stop("level must have a value.")
  if (level == "final") {
    lvl <- -2
  } else if (identical(level, "hierarchical")) {
    stop("hierarchical output not yet supported.")
  } else {
    lvl <- as.integer(level - 1)
  }

  part <- directed_louvain_fast(norm$edges_num, precision, gamma,
                           #reproducibility,
                           renumbering, randomized, lvl, verbose)
  names(part) <- norm$nodes
  part
}

directed_modularity <- function(x, partition) {
  norm <- normalize_to_edges(x)
  nodes <- norm$nodes

  if (!is.null(names(partition))) {
    idx <- match(nodes, names(partition))
    if (any(is.na(idx)))
      stop("Partition names must cover all node names")
    part_vec <- partition[idx]
  } else {
    if (length(partition) != length(nodes))
      stop("Partition length must equal number of nodes")
    part_vec <- partition
  }

  directed_modularity_fast(norm$edges_num, part_vec)
}
