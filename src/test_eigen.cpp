#include <vector>
#include <array>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector Ax_tree(const Rcpp::IntegerMatrix &edge,
                            const Rcpp::NumericVector &lengths,
                            const Rcpp::NumericVector &x,
                            int nNodes) {
  Rcpp::NumericVector result(nNodes);

  int nEdges = edge.nrow();
  for (int i = 0; i < nEdges; i++) {
    int u = edge(i, 0) - 1; // convert to zero-based
    int v = edge(i, 1) - 1;

    double w = lengths[i];

    // adjacency contribution
    result[u] += w * x[v];
    result[v] += w * x[u];
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector Lx_tree_weighted(const Rcpp::IntegerMatrix &edge,
                               const Rcpp::NumericVector &w,
                               const Rcpp::NumericVector &x,
                               int nNodes) {
  Rcpp::NumericVector result(nNodes);

  int nEdges = edge.nrow();
  for (int i = 0; i < nEdges; i++) {

    int u = edge(i, 0) - 1;  // zero-based
    int v = edge(i, 1) - 1;
    double weight = w[i];

    // contribution to Laplacian:
    // L(u) += w * (x[u] - x[v])
    // L(v) += w * (x[v] - x[u])

    double diff_uv = x[u] - x[v];

    result[u] += weight * diff_uv;
    result[v] -= weight * diff_uv;  // = weight * (x[v] - x[u])
  }

  return result;
}



// some chatGPT magic:



// =============================================================
// Memoized (cached) precalc object stored in a static variable
// =============================================================

static Rcpp::List precalc_cache;
static bool precalc_initialized = false;


// =============================================================
// Precompute edge arrays (zero-based) + degree vector
// =============================================================
//' precompute laplacian stuff
//' @param edge edge
//' @param w weights
//' @param nNodes nnodes
//' @export
// [[Rcpp::export]]
List precalc_Laplacian_mem(const IntegerMatrix &edge,
                           const NumericVector &w,
                           int nNodes)
{
  int m = edge.nrow();
  IntegerVector u(m), v(m);
  NumericVector weights(m);
  NumericVector deg(nNodes);


  for (int i = 0; i < m; ++i) {
    int uu = edge(i, 0) - 1;
    int vv = edge(i, 1) - 1;
    double wi = w[i];
    u[i] = uu;
    v[i] = vv;
    weights[i] = wi;
    deg[uu] += wi;
    deg[vv] += wi;
  }


  precalc_cache = List::create(
    Named("u") = u,
    Named("v") = v,
    Named("w") = weights,
    Named("deg") = deg
  );


  precalc_initialized = true;
  return precalc_cache;
}

// =============================================================
// Compute Lx using the memorized precalc object
// =============================================================
//' lx_mem
//' @param x vector
//' @export
// [[Rcpp::export]]
NumericVector Lx_mem(const NumericVector &x)
{
  if (!precalc_initialized)
    stop("Precalc object not initialized. Call precalc_Laplacian_mem() first.");


  IntegerVector u = precalc_cache["u"];
  IntegerVector v = precalc_cache["v"];
  NumericVector w = precalc_cache["w"];
  NumericVector deg = precalc_cache["deg"];


  int nNodes = deg.size();
  int m = u.size();


  if ((int)x.size() != nNodes)
    stop("Length of x does not match precalc node count.");


  NumericVector result(nNodes);


  // Start with D * x
  for (int i = 0; i < nNodes; ++i)
    result[i] = deg[i] * x[i];


  // Subtract A * x term
  for (int i = 0; i < m; ++i) {
    int uu = u[i];
    int vv = v[i];
    double wi = w[i];


    result[uu] -= wi * x[vv];
    result[vv] -= wi * x[uu];
  }


  return result;
}


//' lx_mem_fast
//' @param x vector
//' @export
// [[Rcpp::export]]
NumericVector Lx_mem_fast(const NumericVector &x) {
  if (!precalc_initialized)
    stop("Call precalc_Laplacian_mem() first.");

  IntegerVector u_ = precalc_cache["u"];
  IntegerVector v_ = precalc_cache["v"];
  NumericVector w_ = precalc_cache["w"];
  NumericVector deg_ = precalc_cache["deg"];

  int n = deg_.size();
  int m = u_.size();

  const int* u = u_.begin();
  const int* v = v_.begin();
  const double* w = w_.begin();
  const double* xv = x.begin();
  const double* deg = deg_.begin();

  NumericVector out(n);
  double* y = out.begin();

  // D * x
  for (int i = 0; i < n; i++)
    y[i] = deg[i] * xv[i];

  // -A * x
  for (int i = 0; i < m; i++) {
    int uu = u[i];
    int vv = v[i];
    double wi = w[i];

    y[uu] -= wi * xv[vv];
    y[vv] -= wi * xv[uu];
  }

  return out;
}


// =============================================================
// Raw-pointer CSR Laplacian multiplier (fastest single-threaded)
// =============================================================


static bool csr_initialized = false;


// CSR storage
static std::vector<int> csr_indptr; // size n+1
static std::vector<int> csr_indices; // neighbor indices
static std::vector<double> csr_weights; // corresponding weights
static std::vector<double> csr_deg; // degree vector
static int csr_n = 0;


// =============================================================
// Build CSR adjacency structure
// =============================================================
// [[Rcpp::export]]
void precalc_Laplacian_csr(const IntegerMatrix &edge,
                           const NumericVector &w,
                           int nNodes)
{
  csr_n = nNodes;
  int m = edge.nrow();


  // prepare temporary neighbor lists
  std::vector<std::vector<std::pair<int,double>>> adj(nNodes);
  csr_deg.assign(nNodes, 0.0);


  for (int i = 0; i < m; i++) {
    int u = edge(i,0) - 1;
    int v = edge(i,1) - 1;
    double wi = w[i];


    adj[u].push_back({v, wi});
    adj[v].push_back({u, wi});


    csr_deg[u] += wi;
    csr_deg[v] += wi;
  }


  // convert to CSR
  csr_indptr.resize(nNodes + 1);
  csr_indptr[0] = 0;
  for (int i = 0; i < nNodes; i++)
    csr_indptr[i+1] = csr_indptr[i] + adj[i].size();


  int total = csr_indptr[nNodes];
  csr_indices.resize(total);
  csr_weights.resize(total);


  for (int i = 0; i < nNodes; i++) {
    int start = csr_indptr[i];
    for (size_t k = 0; k < adj[i].size(); k++) {
      csr_indices[start + k] = adj[i][k].first;
      csr_weights[start + k] = adj[i][k].second;
    }
  }


  csr_initialized = true;
}


// =============================================================
// Compute Lx = Deg(x) - A x using CSR (FAST)
// =============================================================
// [[Rcpp::export]]
NumericVector Lx_csr(const NumericVector &x)
{
  if (!csr_initialized)
    stop("Call precalc_Laplacian_csr() first.");


  if ((int)x.size() != csr_n)
    stop("x has incorrect length.");


  const double* xv = x.begin();
  NumericVector out(csr_n);
  double* y = out.begin();


  const int* indptr = csr_indptr.data();
  const int* indices = csr_indices.data();
  const double* w = csr_weights.data();
  const double* deg = csr_deg.data();


  for (int i = 0; i < csr_n; i++) {
    double yi = deg[i] * xv[i];
    int start = indptr[i];
    int end = indptr[i+1];
    for (int k = start; k < end; k++) {
      yi -= w[k] * xv[ indices[k] ];
    }
    y[i] = yi;
  }


  return out;
}
