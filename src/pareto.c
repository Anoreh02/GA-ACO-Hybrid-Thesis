#include <float.h>
#include <stdlib.h>
#include "pareto.h"
#include "eval.h"

// Simple O(N^2) non-dominated rank
void non_dominated_sort(ParetoView* V){
    for(int i=0;i<V->size;++i){
        int r=0;
        for(int j=0;j<V->size;++j) if(i!=j && dominates(&V->pop[j], &V->pop[i])) ++r;
        V->rank[i]=r;
    }
}

void compute_crowding(ParetoView* V){
    // Minimal placeholder (you can implement full 3-objective crowding later)
    for(int i=0;i<V->size;++i) V->crowding[i]=0.0;
}

int binary_tournament(ParetoView* V, int a, int b){
    if(V->rank[a] < V->rank[b]) return a;
    if(V->rank[b] < V->rank[a]) return b;
    return (V->crowding[a] > V->crowding[b]) ? a : b;
}
