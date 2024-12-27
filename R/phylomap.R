#' Turns a phylogeny into a graph projected onto geographical space
#'
#' This function takes a phylogeny object (class phylo) and the geographical
#' coordinates of the leaves, infers the internal node locations, and then
#' converts the phylogeny into an sfnetworks spatial graph.
#'
#' @import ape
#' @import phytools
#' @import dplyr
#' @import magrittr
#' @import tidyr
#' @import tidygraph
#' @import sfnetworks
#' @import sf
#'
#' @param tree Phylogenetic tree of class `phylo`.
#' @param data An object of class `data.frame` or `matrix` with coordinates of the tips in the tree in columns 'lon' and 'lat' and rownames = tip labels
#' @param crs The coordinate reference system of the lon and lat data. Default is 4326, the EPSG code for World Geodetic System 1984. See https://epsg.io/4326.
#' @return An `sfnetwork` object.
#'
#' @examples
#' library(phylomapR)
#' library(ape)
#' library(phytools)
#' library(ggraph)
#' library(ggplot2)
#' library(magrittr)
#' library(dplyr)
#' library(sfnetworks)
#' library(tidygraph)
#' random_tree_rooted<-rtree(26,rooted = T,tip.label = letters)
#' random_xy<-data.frame(
#'   lon=fastBM(random_tree_rooted,bounds=c(15,25),sig2=2),
#'   lat=fastBM(random_tree_rooted,bounds=c(0,10),sig2=2)
#' )
#' # Note that fastBM produces named vectors so we don't have to add rownames ourselves
#' test_rooted<-phylomap(random_tree_rooted,random_xy)
#' ### basic plot ----
#' test_rooted %>%
#'   ggraph('sf')+
#'   geom_node_sf()+
#'   geom_edge_sf()+
#'   theme_bw()
#'
#' ### add graphy info to nodes and edges and plot ----
#' test_rooted %>%
#'   activate(nodes) %>%
#'   mutate(centrality=centrality_betweenness(normalized = F)) %>%
#'   activate(edges) %>%
#'   mutate(betweenness=centrality_edge_betweenness()) %>%
#'   ggraph('sf')+
#'   geom_edge_sf(aes(alpha=betweenness),linewidth=2)+
#'   geom_node_sf(aes(colour=centrality,size=centrality))+
#'   theme_bw()+
#'   scale_colour_viridis_c()
#'
#' @export
phylomap <- function(tree,data,crs=4326) {
  # Checks ----
  nt <- length(tree$tip.label)
  if (nt != nrow(data) ||
      is.null(rownames(data)) ||
      any(!rownames(data) %in% tree$tip.label) ||
      any(!tree$tip.label %in% rownames(data)))
    stop("Tree and data don't match. tree$tip.label needs to correspond to rownames(data)")
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data needs to be either a matrix or a data.frame.")
  if (any(!c("lon", "lat") %in% colnames(data)))
    stop("Data needs to have coordinate columns named 'lon' and 'lat'.")
  if( any(data[,"lon"] < -180 | data[,"lon"] > 180) || any(data[,"lat"] < -90 | data[,"lat"] > 90))
    warning("Coordinates are not in decimal degrees. Make sure you specify the correct crs!")

  # Processing ----
  data <- data[tree$tip.label, ]
  tree$node.label <- paste0("I", c(1:(Nnode(tree))))
  xy <- data[, c("lon", "lat")]

  # Find internal node coordinates with fastAnc ----
  xy_ancestral <- apply(xy, 2, fastAnc, tree = tree)

  xy_interpolated <-
    data.frame(
      node = rownames(xy_ancestral),
      name = tree$node.label,
      xy_ancestral) %>%
    bind_rows(
      data.frame(
        node = as.character(1:Ntip(tree)),
        name = tree$tip.label,
        lon = xy[, "lon"],
        lat = xy[, "lat"]
      )
    )

  # Making a tidy graph ----
  graph_interpolated<-
    as_tbl_graph(tree) %>%
    activate(nodes) %>%
    left_join(xy_interpolated)

  # Convert the graph into a spatial network ----
  graph_interpolated_sfnet<-
    graph_interpolated %>%
    as_sfnetwork(coords=c("lon","lat"),
                 crs=crs) %>%
    # we can add the line geometry automatically!
    to_spatial_explicit()
  # Below is the old way of adding the line geometry
  # make_edges_explicit() %>%
  # make_edges_valid()

  return(graph_interpolated_sfnet$explicit)
}
