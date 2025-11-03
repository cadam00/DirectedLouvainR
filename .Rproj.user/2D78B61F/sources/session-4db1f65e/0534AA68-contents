\name{directed_modularity}
\alias{directed_modularity}
\title{
Directed Lovain clustering
}
\description{
Get partitions for a directed Lovain clustering for a graph.
}
\usage{
directed_modularity(x, partition)
}
\arguments{
\item{x}{
Either an \code{igraph} object from \pkg{igraph} package, an edge list
(\code{data.frame} object) or an adjacency matrix (\code{matrix} object). If an
edge list, then it can have either the first two columns, \code{"from"} and
\code{"to"}, or three columns, containing \code{"weight"}.
}
\item{partition}{
\code{integer} vector object, with length equal to the number of the graph
nodes.
}
}
\details{
Directed modularity, as described in Dugué and Perez (2015, 2022).
This implementation has taken parts of code from the supplementary code for the
above paper, located in \cr
\href{https://github.com/nicolasdugue/DirectedLouvain}{
https://github.com/nicolasdugue/DirectedLouvain}.
}

\value{
Single modularity numeric value.
}

\references{
Dugué, N., & Perez, A. (2022). Direction matters in complex networks: A
theoretical and applied study for greedy modularity optimization. \emph{Physica
A: Statistical Mechanics and its Applications}, 603, 127798.
\doi{10.1016/j.physa.2022.127798}

Dugué, N., & Perez, A. (2015). Directed Louvain: maximizing modularity in
directed networks (Doctoral dissertation, Université d'Orléans).
}

\seealso{
\code{
\link[DirectedLouvainR]{directed_modularity}
}
}

\examples{
# Edge list as data.frame
el <- data.frame(
  from   = c(1, 2, 3),
  to     = c(2, 3, 1),
  weight = c(1, 2, 3)
)

part <- directed_louvain(el)
directed_modularity(el, part)

# Edge list as data.frame
elc <- data.frame(
  from   = c("A", "B", "C"),
  to     = c("B", "C", "A"),
  weight = c(1, 2, 3),
  stringsAsFactors = FALSE
)

part <- directed_louvain(elc)
directed_modularity(elc, part)

# Adjacency matrix
adj <- matrix(c(0,1,0,
                0,0,2,
                3,0,0), nrow=3, byrow=TRUE)

rownames(adj) <- c("x","y","z")
colnames(adj) <- c("z","x","y")  # deliberately permuted

part <- directed_louvain(adj)   # columns realigned to row order
directed_modularity(adj, part)

# igraph
library(igraph)

g <- make_ring(4)
E(g)$weight <- c(1,2,3,4)
V(g)$name <- c("n1","n2","n3","n4")

part <- directed_louvain(g)
directed_modularity(g, part)
}

