#include <Rcpp.h>
#include "community.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector directed_louvain_fast(NumericMatrix edges,
                                    double precision = 0.0001,
                                    double gamma = 1.0, // bool reproducibility = false,
                                    bool renumbering = false,
                                    bool randomized = true,
                                    int level = -2,
                                    bool verbose = false) {
  // Build Community object
  int nrow = edges.nrow();
  int ncol = edges.ncol();
  if (ncol != 2 && ncol != 3)
    stop("edges must have 2 or 3 columns: (src, dst) or (src, dst, weight)");

  unsigned int n = 0;
  for (int i = 0; i < nrow; i++) {
    n = std::max(n, (unsigned int) edges(i,0));
    n = std::max(n, (unsigned int) edges(i,1));
  }
  // n is already the number of nodes, since R IDs are 1-based
  std::vector<std::vector<std::pair<unsigned int,double>>> adj(n);

  for (int i = 0; i < nrow; i++) {
    unsigned int u = (unsigned int) edges(i,0) - 1;
    unsigned int v = (unsigned int) edges(i,1) - 1;
    if (u >= n || v >= n) Rcpp::stop("Edge index out of range");
    double w = (ncol == 3) ? edges(i,2) : 1.0;
    adj[u].push_back(std::make_pair(v, w));
  }

  Graph g(adj, // reproducibility,
          renumbering, false, true);
  Community c(g, precision, gamma, // reproducibility,
              renumbering, randomized);

  c.run(verbose, level, "");

  std::map<unsigned int, unsigned int> lvl;
  if (level == -2) {
    lvl = c.get_last_level();
  } else {
    lvl = c.get_level(level < 0 ? 0 : level);
  }

  IntegerVector res(lvl.size());
  CharacterVector names(lvl.size());
  unsigned int i = 0;
  for (auto &kv : lvl) {
    res[i] = kv.second + 1;
    names[i] = std::to_string(kv.first + 1);
    i++;
  }
  res.attr("names") = names;
  return res;
}


// [[Rcpp::export]]
double directed_modularity_fast(NumericMatrix edges, IntegerVector partition) {
  int nrow = edges.nrow();
  int ncol = edges.ncol();
  if (ncol != 2 && ncol != 3)
    stop("edges must have 2 or 3 columns: (src, dst) or (src, dst, weight)");

  // number of nodes = max node ID
  int n_nodes = 0;
  for (int i = 0; i < nrow; i++) {
    n_nodes = std::max(n_nodes, (int) edges(i,0));
    n_nodes = std::max(n_nodes, (int) edges(i,1));
  }

  if (partition.size() != n_nodes)
    stop("Partition length must equal number of nodes");

  // compute in/out degrees
  std::vector<double> kout(n_nodes, 0.0), kin(n_nodes, 0.0);
  double m = 0.0;
  for (int i = 0; i < nrow; i++) {
    int u = edges(i,0) - 1;
    int v = edges(i,1) - 1;
    double w = (ncol == 3) ? edges(i,2) : 1.0;
    kout[u] += w;
    kin[v]  += w;
    m += w;
  }

  // aggregate by community
  int n_comm = *std::max_element(partition.begin(), partition.end());
  std::vector<double> tot_out(n_comm, 0.0), tot_in(n_comm, 0.0), in(n_comm, 0.0);

  for (int i = 0; i < nrow; i++) {
    int u = edges(i,0) - 1;
    int v = edges(i,1) - 1;
    double w = (ncol == 3) ? edges(i,2) : 1.0;
    if (partition[u] == partition[v]) {
      in[partition[u]-1] += w;
    }
  }

  for (int u = 0; u < n_nodes; u++) {
    int c = partition[u] - 1;
    tot_out[c] += kout[u];
    tot_in[c]  += kin[u];
  }

  // compute modularity
  double Q = 0.0;
  for (int c = 0; c < n_comm; c++) {
    Q += in[c]/m - (tot_out[c] * tot_in[c]) / (m*m);
  }

  return Q;
}

