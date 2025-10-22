// pareto.c
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "pareto.h"
#include "eval.h"

/* ------------- Objective adapter -------------
   Convert all objectives to a "minimize" view:
   - f1: conflicts          (minimize as-is)
   - f2: utilization        (maximize) -> negate to minimize
   - f3: workload variance  (minimize as-is)
*/
static inline double obj_as_min(const Timetable* S, int k){
    switch(k){
        case 0: return S->f1_conflicts;           // minimize
        case 1: return -S->f2_utilization;        // maximize -> minimize by negation
        default: return S->f3_workload_var;       // minimize
    }
}

/* ------------- Non-dominated sort -------------
   O(N^2) rank assignment is fine for moderate N.
   Assumes 'dominates(a,b)' is declared in pareto.h and
   compares Timetable a vs b across the three objectives
   with correct min/max semantics (we adapt here already).
*/
void non_dominated_sort(ParetoView* V){
    if(!V || V->size <= 0) return;
    for(int i=0;i<V->size;++i){
        int r = 0;
        for(int j=0;j<V->size;++j){
            if(i==j) continue;
            if(dominates(&V->pop[j], &V->pop[i])) ++r;
        }
        V->rank[i] = r; // 0 means non-dominated
    }
}

/* -------- Portable qsort() comparator context --------
   We need to sort indices of a *single* Pareto front by a given
   objective k (in "minimize" space). Use static globals so this
   works on all platforms (no qsort_r).
*/
static ParetoView* gV = NULL;
static int gK = 0;

static int cmp_obj_asc_wrap(const void* a, const void* b){
    int ia = *(const int*)a;
    int ib = *(const int*)b;
    double va = obj_as_min(&gV->pop[ia], gK);
    double vb = obj_as_min(&gV->pop[ib], gK);
    if(va < vb) return -1;
    if(va > vb) return  1;
    return 0;
}

/* ------------- Crowding distance (NSGA-II style) -------------
   - Computed per front (same rank).
   - Boundary points on ANY objective get INFINITY crowding.
   - Interior points sum normalized (v_next - v_prev)/range for k=0..2.
*/
void compute_crowding(ParetoView* V){
    if(!V || V->size <= 0) return;

    // reset crowding
    for(int i=0;i<V->size;++i) V->crowding[i] = 0.0;

    // find max rank to iterate fronts
    int rmax = 0;
    for(int i=0;i<V->size;++i) if(V->rank[i] > rmax) rmax = V->rank[i];

    // scratch index buffer (at most V->size elements in a front)
    int* idx = (int*)malloc(sizeof(int) * (V->size > 0 ? V->size : 1));
    if(!idx) return;

    for(int r=0; r<=rmax; ++r){
        // collect indices for this front
        int n = 0;
        for(int i=0;i<V->size;++i) if(V->rank[i] == r) idx[n++] = i;
        if(n == 0) continue;

        // trivial fronts: keep extremes
        if(n <= 2){
            for(int t=0;t<n;++t) V->crowding[idx[t]] = INFINITY;
            continue;
        }

        // process each objective k
        for(int k=0;k<3;++k){
            // compute range in "minimize" space
            double vmin = DBL_MAX, vmax = -DBL_MAX;
            for(int t=0;t<n;++t){
                double v = obj_as_min(&V->pop[idx[t]], k);
                if(v < vmin) vmin = v;
                if(v > vmax) vmax = v;
            }
            double range = vmax - vmin;
            if(range < 1e-12) range = 1.0; // avoid division by ~0 (all equal -> zero contribution)

            // sort this front by objective k ascending
            gV = V; gK = k;
            qsort(idx, n, sizeof(int), cmp_obj_asc_wrap);

            // boundary points get infinity to preserve extremes
            V->crowding[idx[0]]   = INFINITY;
            V->crowding[idx[n-1]] = INFINITY;

            // interior points: accumulate normalized perimeter distance
            for(int t=1; t<n-1; ++t){
                int im = idx[t-1], ic = idx[t], ip = idx[t+1];
                double vprev = obj_as_min(&V->pop[im], k);
                double vnext = obj_as_min(&V->pop[ip], k);
                double inc   = (vnext - vprev) / range;  // >= 0
                if(!isinf(V->crowding[ic]))              // keep INF if boundary on another objective
                    V->crowding[ic] += inc;
            }
        }
    }

    free(idx);
}

/* ------------- Binary tournament -------------
   Rank first (lower is better), then crowding (higher is better).
*/
int binary_tournament(ParetoView* V, int a, int b){
    if(V->rank[a] < V->rank[b]) return a;
    if(V->rank[b] < V->rank[a]) return b;
    return (V->crowding[a] > V->crowding[b]) ? a : b;
}
