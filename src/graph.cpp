#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cassert>
#include <Rcpp.h>

#include <graph.hpp>

void init_attributes(
    Graph &g,
    std::vector<std::vector<std::pair<unsigned int,double>>> &LOUT,
    std::vector<std::vector<std::pair<unsigned int,double>>> &LIN,
    bool verbose
) {
  // Count arcs
  g.arcs = 0;
  for (const auto &nbrs : LOUT)
    g.arcs += static_cast<unsigned int>(nbrs.size());

  // Prepare containers
  g.outdegrees.resize(g.nodes);
  g.indegrees.resize(g.nodes);
  g.outcoming_arcs.resize(g.arcs);
  g.incoming_arcs.resize(g.arcs);

  // Only allocate weights arrays if g.weighted was set by the constructor
  if (g.weighted) {
    g.outcoming_weights.resize(g.arcs);
    g.incoming_weights.resize(g.arcs);
  } else {
    g.outcoming_weights.clear();
    g.incoming_weights.clear();
  }

  // Fill CSR: out
  unsigned long outpos = 0;
  for (unsigned int i = 0; i < g.nodes; ++i) {
    const auto &nbrs = LOUT[i];
    for (size_t k = 0; k < nbrs.size(); ++k) {
      g.outcoming_arcs[outpos + k] = nbrs[k].first;
      if (g.weighted) g.outcoming_weights[outpos + k] = nbrs[k].second;
    }
    outpos += static_cast<unsigned long>(nbrs.size());
    g.outdegrees[i] = outpos;
  }

  // Fill CSR: in
  unsigned long inpos = 0;
  for (unsigned int i = 0; i < g.nodes; ++i) {
    const auto &nbrs = LIN[i];
    for (size_t k = 0; k < nbrs.size(); ++k) {
      g.incoming_arcs[inpos + k] = nbrs[k].first;
      if (g.weighted) g.incoming_weights[inpos + k] = nbrs[k].second;
    }
    inpos += static_cast<unsigned long>(nbrs.size());
    g.indegrees[i] = inpos;
  }

  // Total weight = sum of weighted out-degrees
  g.total_weight = 0.0;
  for (unsigned int i = 0; i < g.nodes; ++i)
    g.total_weight += g.weighted_out_degree(i);

  if (verbose) {
    Rcpp::Rcerr << "init_attributes: nodes=" << g.nodes
              << " arcs=" << g.arcs
              << " weighted=" << (g.weighted ? "true" : "false")
              << " total_weight=" << g.total_weight << std::endl;
  }
}


// Default constructor
Graph::Graph() {
  this->nodes         = 0;
  this->arcs          = 0;
  this->total_weight  = 0;
  this->weighted      = false;

  this->outcoming_arcs.clear();
  this->incoming_arcs.clear();
  this->outcoming_weights.clear();
  this->incoming_weights.clear();
  this->outdegrees.clear();
  this->indegrees.clear();
}

// Constructor from weighted adjacency list
Graph::Graph(std::vector<std::vector<std::pair<unsigned int, double>>> adjacency_list,
             bool renumbering, bool verbose, bool weighted) {
  std::vector<std::vector<std::pair<unsigned int,double>>> LIN(adjacency_list.size());

  this->correspondance.clear();
  for (unsigned int i = 0; i < adjacency_list.size(); ++i)
    correspondance.push_back(i);

  this->nodes = adjacency_list.size();

  // Build incoming links from outgoing links
  for (unsigned int i = 0 ; i < adjacency_list.size(); ++i) {
    for (auto v : adjacency_list[i]) {
      LIN[v.first].push_back(std::make_pair(i,v.second));
    }
  }

  init_attributes(*this, adjacency_list, LIN, verbose);
}

// Constructor from unweighted adjacency list
Graph::Graph(std::vector<std::vector<unsigned int>> adjacency_list,
             bool renumbering, bool verbose, bool weighted) {
  std::vector<std::vector<std::pair<unsigned int,double>>> LOUT(adjacency_list.size());
  std::vector<std::vector<std::pair<unsigned int,double>>> LIN(adjacency_list.size());

  this->correspondance.clear();
  for (unsigned int i = 0; i < adjacency_list.size(); ++i)
    correspondance.push_back(i);

  this->nodes = adjacency_list.size();

  // Build incoming links from outgoing links
  for (unsigned int i = 0 ; i < adjacency_list.size(); ++i) {
    for (auto v : adjacency_list[i]) {
      LOUT[i].push_back(std::make_pair(v,1));
      LIN[v].push_back(std::make_pair(i,1));
    }
  }

  init_attributes(*this, LOUT, LIN, verbose);
}

// Copy constructor
Graph::Graph(const Graph &g) {
  this->weighted          = g.weighted;
  this->nodes             = g.nodes;
  this->arcs              = g.arcs;
  this->total_weight      = g.total_weight;
  this->outcoming_arcs    = g.outcoming_arcs;
  this->incoming_arcs     = g.incoming_arcs;
  this->outdegrees        = g.outdegrees;
  this->indegrees         = g.indegrees;
  this->outcoming_weights = g.outcoming_weights;
  this->incoming_weights  = g.incoming_weights;
  this->correspondance    = g.correspondance;
}

// Write graph to binary file
void Graph::write(std::string outfile) {
  std::ofstream foutput;
  foutput.open(outfile, std::fstream::out | std::fstream::binary);
  assert(foutput.rdstate() == std::ios::goodbit);

  foutput.write((char*)(&this->nodes), sizeof(unsigned int));
  foutput.write((char*)(&this->outdegrees[0]), sizeof(unsigned long) * this->nodes);
  foutput.write((char*)(&this->outcoming_arcs[0]), sizeof(unsigned int) * this->arcs);
  if (this->weighted)
    foutput.write((char*)(&this->outcoming_weights[0]), sizeof(double) * this->arcs);
  foutput.write((char*)(&this->indegrees[0]), sizeof(unsigned long) * this->nodes);
  foutput.write((char*)(&this->incoming_arcs[0]), sizeof(unsigned int) * this->arcs);
  if (this->weighted)
    foutput.write((char*)(&this->incoming_weights[0]), sizeof(double) * this->arcs);
  foutput.write((char*)(&this->correspondance[0]), sizeof(unsigned long) * this->nodes);
  foutput.close();
}

// Load graph from binary file
void Graph::load(std::string filename, bool verbose) {
  std::ifstream finput;
  finput.open(filename, std::fstream::in | std::fstream::binary);
  assert(finput.rdstate() == std::ios::goodbit);
  this->outcoming_weights.clear();
  this->incoming_weights.clear();

  finput.read((char*)&this->nodes, sizeof(unsigned int));
  if (verbose)
    Rcpp::Rcerr << "number of nodes: " << this->nodes << std::endl;

  this->outdegrees.resize(this->nodes);
  finput.read((char*)&this->outdegrees[0], this->nodes * sizeof(unsigned long));

  this->arcs = outdegrees[this->nodes - 1];
  this->outcoming_arcs.resize(this->arcs);
  finput.read((char*)(&this->outcoming_arcs[0]), this->arcs * sizeof(unsigned int));

  if (this->weighted) {
    this->outcoming_weights.resize(arcs);
    finput.read((char*)&this->outcoming_weights[0], this->arcs * sizeof(double));
  }

  this->indegrees.resize(this->nodes);
  finput.read((char*)&this->indegrees[0], this->nodes * sizeof(unsigned long));

  if (verbose)
    Rcpp::Rcerr << "number of arcs:" << this->indegrees[this->nodes - 1] << std::endl;

  this->incoming_arcs.resize(this->arcs);
  finput.read((char*)(&this->incoming_arcs[0]), this->arcs * sizeof(unsigned int));

  if (this->weighted) {
    this->incoming_weights.resize(this->arcs);
    finput.read((char*)&this->incoming_weights[0], this->arcs * sizeof(double));
  }

  this->correspondance.resize(this->nodes);
  finput.read((char*)(&this->correspondance[0]), this->nodes * sizeof(unsigned long));

  this->total_weight = 0.f;
  for (unsigned int i = 0; i < this->nodes; ++i) {
    this->total_weight += this->weighted_out_degree(i);
  }

  if (verbose)
    Rcpp::Rcerr << "total weight: " << this->total_weight << std::endl;
  finput.close();
}

// Display graph in edgelist format
void Graph::display() const {
  for (unsigned int node = 0; node < this->nodes; ++node) {
    size_t p = this->out_neighbors(node);
    Rcpp::Rcout << this->correspondance[node] << ":";
    for (unsigned int i = 0; i < out_degree(node); ++i) {
      if (this->weighted)
        Rcpp::Rcout << " (" << this->outcoming_arcs[p + i] << " " << this->outcoming_weights[p + i] << ")";
      else
        Rcpp::Rcout << " " << this->outcoming_arcs[p+i];
    }
    Rcpp::Rcout << std::endl;
  }
}

double Graph::count_selfloops(unsigned int node) {
  assert(node < this->nodes);
  size_t p = this->out_neighbors(node);
  for (unsigned int i=0 ; i < this->out_degree(node) ; ++i) {
    if (this->outcoming_arcs[p+i]==node) {
      if (this->weighted)
        return this->outcoming_weights[p+i];
      else
        return 1.;
    }
  }
  return 0.;
}

double Graph::weighted_out_degree(unsigned int node) {
  assert(node < this->nodes);
  if (!this->weighted)
    return this->out_degree(node);
  else {
    size_t p = this->out_neighbors(node);
    double res = 0;
    for (unsigned int i=0 ; i < this->out_degree(node) ; ++i)
      res += this->outcoming_weights[p+i];
    return res;
  }
}

double Graph::weighted_in_degree(unsigned int node) {
  assert(node < this->nodes);
  if (!this->weighted)
    return this->in_degree(node);
  else {
    size_t p = this->in_neighbors(node);
    double res = 0;
    for (unsigned int i=0 ; i < this->in_degree(node) ; ++i)
      res += this->incoming_weights[p+i];
    return res;
  }
}

// Friend and static functions are in a separate include
// #include "impl/graph_friend_static.inc"
