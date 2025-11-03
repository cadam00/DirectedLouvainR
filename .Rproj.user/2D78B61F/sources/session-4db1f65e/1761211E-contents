\name{directed_louvain}
\alias{directed_louvain}
\title{
Directed Lovain clustering
}
\description{
Get partitions for a directed Lovain clustering for a graph.
}
\usage{
directed_louvain(x, precision = 0.0001, gamma = 1.0, renumbering = FALSE,
                 randomized = TRUE, level = "final", verbose = FALSE)
}
\arguments{
\item{x}{
Either an \code{igraph} object from \pkg{igraph} package, an edge list
(\code{data.frame} object) or an adjacency matrix (\code{matrix} object). If an
edge list, then it can have either the first two columns, \code{"from"} and
\code{"to"}, or three columns, containing \code{"weight"}.
}
\item{precision}{
A numeric value representing the minimal improvement in modularity required to
continue computing.
}
\item{gamma}{
Size of the generated clusters and it is the resolution parameter. The relative
significance of intra-community edges is scaled. Smaller clusters are associated
with a higher gamma and larger clusters with a lower one. The default is for the
original modularity definition and equal to \code{1}.
}
\item{renumbering}{
\code{logical} object restarting partition numbering from \code{1}, until the
maximum number of communities detected. The default is \code{FALSE},
corresponding to the order that they are detected.
}
\item{randomized}{
\code{logical} object for randomized node updates, escaping locally optimal
solutions. It means allowing to shuffle the order in which nodes are considered
for moves. The default is \code{TRUE}, meaning randomized node updates.
}
\item{level}{
Either \code{character} or \code{integer} object for controlling the number of
levels (steps) of the algorithm. \code{level = "final"} means that the final
level is returned, according to the precision above. If \code{level = 1}, then
the first level (step) of the algorithm is returned, if \code{level = 2}, then
the second level (step) etc. The default is \code{"final"}.
}
\item{verbose}{
\code{logical} object for printing additional information on console.
}
}
\details{
Algorithm of directed Louvain algorithm, utilizing the directed modularity
definition, as described in Dugué and Perez (2015, 2022). This implementation
has taken parts of code from the supplementary code for the above paper, located
in \cr \href{https://github.com/nicolasdugue/DirectedLouvain}{
https://github.com/nicolasdugue/DirectedLouvain}.

\code{directed_louvain} refers to a more \code{R}-friendly version for the
inputs. The internal \cr \code{DirectedLouvainR:::directed_louvain_fast}
accepts for \code{x} a \code{matrix} object with two first columns as the node
IDs and the third optional the weights.
\code{DirectedLouvainR:::directed_louvain_fast} is not listed here, because its
inputs are not normalized and so its use is more "dangerous", but it is a bit
faster because of the lack of pre-processing and IDs are advised to start from
1 until the number of nodes.
}

\value{
Partition label number per node.
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

directed_louvain(el)

# Edge list as data.frame
elc <- data.frame(
  from   = c("A", "B", "C"),
  to     = c("B", "C", "A"),
  weight = c(1, 2, 3),
  stringsAsFactors = FALSE
)

directed_louvain(elc)

# Adjacency matrix
adj <- matrix(c(0,1,0,
                0,0,2,
                3,0,0), nrow=3, byrow=TRUE)

rownames(adj) <- c("x","y","z")
colnames(adj) <- c("z","x","y")  # deliberately permuted

directed_louvain(adj)   # columns realigned to row order

# igraph
library(igraph)

g <- make_ring(4)
E(g)$weight <- c(1,2,3,4)
V(g)$name <- c("n1","n2","n3","n4")

directed_louvain(g)
}

