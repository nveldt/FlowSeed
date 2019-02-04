/**
 * @file pprpush_weighted_mex.cc
 * Implement a PPR clustering scheme.
 * 
 * mex pprclus_mex.cc CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims
 */


#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>

#include <unordered_set>
#include <unordered_map>
#define tr1ns std
#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif

#include <mex.h>

struct sparsevec {
  typedef tr1ns::unordered_map<mwIndex,double> map_type;
  map_type map;
  /** Get an element and provide a default value when it doesn't exist
   * This command does not insert the element into the vector
   */
  double get(mwIndex index, double default_value=0.0) {
    map_type::iterator it = map.find(index);
    if (it == map.end()) {
      return default_value;
    } else {
      return it->second;
    }
  }

  /** Compute the sum of all the elements
   * Implements compensated summation
   */
  double sum() {
    double s=0.;
    for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
      s += it->second;
    }
    return s;
  }
  
  /** Compute the max of the element values 
   * This operation returns the first element if the vector is empty.
   */
  mwIndex max_index() {
    mwIndex index=0;
    double maxval=std::numeric_limits<double>::min();
    for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
      if (it->second>maxval) { maxval = it->second; index = it->first; }
    }
    return index;
  }
};

struct sparserow {
    mwSize n, m; 
    mwIndex *ai;
    mwIndex *aj;
    double *a;
    double *d;
};

double machine_eps() {
    double eps = 1.;
    for (size_t i=0; i<52; ++i) {
        eps = eps/2.;
    }
    return eps;
}
double sr_degree(sparserow *s, mwIndex u) {
    if (s->d == NULL) {
        return (s->ai[u+1] - s->ai[u]);
    } else {
        return (s->d[u]);
    }
}

template <class Queue>
int compute_local_pagerank(sparserow *s, sparsevec& r, sparsevec& p, 
    double alpha, double epsilon, int max_push_count, Queue& q) 
{
  sparsevec r1;
  for (sparsevec::map_type::iterator it=r.map.begin(),itend=r.map.end();
        it!=itend;++it){
    //mexPrintf("r[%i] = %f, %f, eps*deg\n", (int)(it->first), it->second, epsilon*sr_degree(s,it->first));
    if (it->second > epsilon*sr_degree(s,it->first)) {
      q.push(it->first);
      r1.map[it->first] = it->second - epsilon*sr_degree(s,it->first);
    }
  }

  double epsilonf = epsilon*(1+machine_eps());
  int push_count = 0;
  while (q.size()>0 && push_count < max_push_count) {
    push_count += 1;
    mwIndex u = q.front();
    q.pop();
    double du = sr_degree(s, u);
    double moving_probability = r1.map[u]; 
    if (moving_probability <= 0) { continue; } // nothing to update here!
    r.map[u] = epsilon*(double)du;
    r1.map[u] = 0;
    p.map[u] += moving_probability;

    // mexPrintf("step %5i node %llu degree %lf move prob %f\n",
    //  push_count, u, du, moving_probability);

    for (mwIndex nzi=s->ai[u]; nzi<s->ai[u+1]; nzi++) {
      double neighbor_update = alpha*moving_probability*(s->a[nzi])/(du);

      mwIndex x = s->aj[nzi];
      double dx = sr_degree(s, x);
      double rxold = r.get(x);
      double rxnew = rxold + neighbor_update;
      r.map[x] = rxnew;
      if (rxnew > epsilonf*dx && rxold <= epsilonf*dx) {
        r1.map[x] = std::max(neighbor_update - (dx*epsilon-rxold), 0.);
        q.push(x);
      } else if (rxnew > epsilonf*dx) {
        r1.map[x] += neighbor_update;
      }
    }
  }
  
  return (push_count);
}



struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

void cluster_from_sweep(sparserow* G, sparsevec& p, 
      std::vector<mwIndex>& cluster, double *outcond, double* outvolume,
      double *outcut)
{
  // now we have to do the sweep over p in sorted order by value
  typedef std::vector< std::pair<int, double> > vertex_prob_type;
  vertex_prob_type prpairs(p.map.begin(), p.map.end());
  std::sort(prpairs.begin(), prpairs.end(), greater2nd());

  // compute cutsize, volume, and conductance
  std::vector<double> conductance(prpairs.size());
  std::vector<mwIndex> volume(prpairs.size());
  std::vector<mwIndex> cutsize(prpairs.size());

  size_t i=0;
  tr1ns::unordered_map<int,size_t> rank;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    rank[it->first] = i;
  }
  //printf("support=%i\n",prpairs.size());
  mwIndex total_degree = G->ai[G->m];
  mwIndex curcutsize = 0;
  mwIndex curvolume = 0;
  i=0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    mwIndex v = it->first;
    mwIndex deg = G->ai[v+1]-G->ai[v];
    mwIndex change = deg;
    for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
      mwIndex nbr = G->aj[nzi];
      if (rank.count(nbr) > 0) {
        if (rank[nbr] < rank[v]) {
          change -= 2;
        }
      }
    }
    curcutsize += change;
    //if (curvolume + deg > target_vol) {
      //break;
    //}
    curvolume += deg;
    volume[i] = curvolume;
    cutsize[i] = curcutsize;
    if (curvolume == 0 || total_degree-curvolume==0) {
      conductance[i] = 1;
    } else {
      conductance[i] = (double)curcutsize/
                        (double)std::min(curvolume,total_degree-curvolume);
    }
    //printf("%5i : cut=%6i vol=%6i prval=%8g cond=%f\n", i, curcutsize, curvolume, it->second, conductance[i]);
  }
  // we stopped the iteration when it finished, or when it hit target_vol
  size_t lastind = i;
  double mincond = std::numeric_limits<double>::max();
  size_t mincondind = 0; // set to zero so that we only add one vertex 
  for (i=0; i<lastind; i++) {
    if (conductance[i] < mincond) {
      mincond = conductance[i];
      mincondind = i;
    }
  }
  //printf("mincond=%f mincondind=%i\n", mincond, mincondind);
  if (lastind == 0) {
    // add a case 
    mincond = 0.0;
  }
  i = 0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend && i<mincondind+1; ++it, ++i) {
    cluster.push_back(it->first);
  }
  if (outcond) { *outcond = mincond; }
  if (outvolume) { *outvolume = volume[mincondind]; }
  if (outcut) { *outcut = cutsize[mincondind]; }
}

struct local_pagerank_stats {
    double conductance;
    double volume;
    double support;
    double steps;
    double nedges;
    double eps;
    double cut;
};

/** Cluster will contain a list of all the vertices in the cluster
 * @param set the set of starting vertices to use
 * @param alpha the value of alpha in the PageRank computation
 * @param target_vol the approximate number of edges in the cluster
 * @param p the pagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template <class Queue>
int hypercluster_pagerank_multiple(sparserow* G, 
    sparsevec& v, double alpha, double eps, 
    sparsevec& p, sparsevec &r, Queue& q,
    std::vector<mwIndex>& cluster, local_pagerank_stats *stats)
{
  // reset data
  p.map.clear();
  r.map.clear();
  q.empty();
  
  assert(eps > 0);
  assert(alpha < 1.0); assert(alpha > 0.0);

  for (sparsevec::map_type::iterator it=v.map.begin(),itend=v.map.end();
        it!=itend;++it) {
    r.map[it->first] = (1.-alpha)*it->second;
  }
  
  double pr_eps = eps;
  if (stats) { stats->eps = pr_eps; }
  
  //printf("find_cluster: target_vol=%7lli alpha=%5.3ld pr_eps=%ld\n", target_vol, alpha, pr_eps);
  
  // calculate an integer number of maxsteps
  double maxsteps = 100./(pr_eps*(1.-alpha));
  maxsteps = std::min(maxsteps, 0.5*(double)std::numeric_limits<int>::max());
      
  int nsteps = compute_local_pagerank(G, r, p, alpha, pr_eps, (int)maxsteps, q);
  int support = r.map.size(); 
  //if (stats) { stats->nedges = nedges; } 
  if (stats) { stats->steps = nsteps; }
  if (stats) { stats->support = support; }
  
  //mexPrintf("setsize=%zu, nsteps=%i, support=%i\n", set.size(), nsteps, support);

  // reset the residual to the pagerank vector and scale the residual
  // this preserves the PageRank solution to return
  if (nsteps > 0) {
    r = p;
  } else {
    // in this case p = 0, and we should sweep over the residual anyway
  }

  

  // scale the probablities by their degree
  for (sparsevec::map_type::iterator it=r.map.begin(),itend=r.map.end();
    it!=itend;++it) {
    it->second *= 1.0/(double)std::max((double)sr_degree(G,it->first),pr_eps);  
  }
  double *outcond = NULL;
  double *outvolume = NULL;
  double *outcut = NULL;
  if (stats) { outcond = &stats->conductance; }
  if (stats) { outvolume = &stats->volume; }
  if (stats) { outcut = &stats->cut; }
  cluster_from_sweep(G, r, cluster, outcond, outvolume, outcut);
  return (0);
}

void pprgrow(sparserow* G, sparsevec& v, std::vector< mwIndex >& set, double alpha,
    double pr_eps, double* fcond, double* fcut,
    double* fvol, sparsevec& p)
{
    sparsevec r;
    std::queue<mwIndex> q;
    local_pagerank_stats stats;
    std::vector<mwIndex> bestclus;
    hypercluster_pagerank_multiple(G, v, alpha, pr_eps, 
        p, r, q, bestclus, &stats);
    set = bestclus;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
}

void copy_array_to_sparse_vector(const mxArray* v, sparsevec& vec, double* degs)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    
    if (mxIsSparse(v)) {

        mxAssert(mxIsDouble(v), "array type is not double");

        double *p = mxGetPr(v);
        mwIndex *inds = mxGetIr(v);
        mwIndex *indend = mxGetJc(v);
        double vsum = 0.;

        //mexPrintf("nnz = %i\n", (int)indend[1]);

        for (mwIndex nzi = indend[0]; nzi < indend[1]; ++nzi) {
            //mexPrintf("%i %f\n", (int)inds[nzi], p[nzi]);
            vec.map[inds[nzi]] = p[nzi];
            vsum += p[nzi];
        }
        //mexPrintf("vsum = %.16f\n", vsum);
        mxAssert(vsum >= 1-10*indend[1]*(2.2204460492503131e-16), "v not stochastic");
        mxAssert(vsum <= 1+10*indend[1]*(2.2204460492503131e-16), "v not stochastic");

    } else {
        size_t n = mxGetNumberOfElements(v);
        double *p = mxGetPr(v);

        double svol = 0.;

        for (size_t i=0; i<n; ++i) {
            double elem = p[i];
            mxAssert(elem >= 1, "Only positive integer elements allowed");
            mwIndex ind = (mwIndex)elem - 1;
            svol += degs[ind];
        }

        for (size_t i=0; i<n; ++i) {
            double elem = p[i];
            mwIndex ind = (mwIndex)elem - 1;
            vec.map[ind] = degs[ind]/svol;
        }
    }
}


// USAGE
// [bestset,cond,cut,vol,y] = pprpush_weighted_mex(A,degs,set,eps,alpha)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mxAssert(nrhs > 2 && nrhs < 6, "3-5 inputs required.");

    const mxArray* mat = prhs[0];
    const mxArray* deg = prhs[1];
    const mxArray* set = prhs[2];

    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");

    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    
    if (nlhs > 1) { plhs[1] = cond; }
    if (nlhs > 2) { plhs[2] = cut; }
    if (nlhs > 3) { plhs[3] = vol; }

    mxAssert(nlhs <= 5, "Too many output arguments");

    double alpha = 0.99;
    if (nrhs >= 5) {
        alpha = mxGetScalar(prhs[4]);
    }
    mxAssert(alpha >= 0. && alpha < 1, "alpha must be 0 <= alpha < 1");

    // use a strange sentinal
    double eps = 1.e-5;
    if (nrhs >= 4) {
        eps = mxGetScalar(prhs[3]);
    }

    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);
    r.aj = mxGetIr(mat);
    r.a = mxGetPr(mat);
    r.d = mxGetPr(deg);

    sparsevec v;
    copy_array_to_sparse_vector( set, v, r.d ); // this code needs r.d
                                                // in case the input 
                                                // is a list of ids
    std::vector< mwIndex > cluster;
    sparsevec ppr;
    pprgrow(&r, v, cluster, alpha, eps,
        mxGetPr(cond), mxGetPr(cut), mxGetPr(vol), ppr);

    if (nlhs > 0) { 
        mxArray* cassign = mxCreateDoubleMatrix(cluster.size(),1,mxREAL);
        plhs[0] = cassign;

        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<cluster.size(); ++i) {
            ci[i] = (double)(cluster[i] + 1);
        }
    }

    if (nlhs > 4) { // sets output "y" to the vector computed
        mxArray* pprvec = mxCreateDoubleMatrix(r.n,1,mxREAL);
        plhs[4] = pprvec;
        double *ci = mxGetPr(pprvec);
        for (sparsevec::map_type::iterator it=ppr.map.begin(),itend=ppr.map.end();
             it!=itend;++it) {
            ci[it->first] = it->second;
        }
    }
}