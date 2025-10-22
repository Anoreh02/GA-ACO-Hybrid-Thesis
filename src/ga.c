#include <stdlib.h>
#include <string.h>
#include "ga.h"
#include "eval.h"
#include "pareto.h"
#include "util.h"

static int L2C(const Instance* I, int L){
    return I->lecture_to_course[L];
}

static void random_feasible_assign(Instance* I, Timetable* S, int L){
    int c = L2C(I,L);
    int tries=0;
    while(tries++<1000){
        int t = rand() % I->num_timeslots;
        if(I->feasible[c][t]){
            S->lec_to_timeslot[L] = t;
            S->lec_to_room[L]     = rand() % I->num_rooms; // improve later
            return;
        }
    }
    S->lec_to_timeslot[L] = -1;
    S->lec_to_room[L]     = 0;
}

static void init_random(Instance* I, Timetable* S){
    for(int L=0; L<I->total_lectures; ++L){
        random_feasible_assign(I,S,L);
    }
    evaluate(I,S);
}

static void crossover(const Timetable* A, const Timetable* B, Timetable* C, Instance* I){
    int cut = rand() % I->total_lectures;
    for(int L=0; L<I->total_lectures; ++L){
        const Timetable* src = (L<cut)?A:B;
        C->lec_to_timeslot[L] = src->lec_to_timeslot[L];
        C->lec_to_room[L]     = src->lec_to_room[L];
    }
}

static void mutate(Instance* I, Timetable* S, double pm){
    for(int L=0; L<I->total_lectures; ++L){
        if(unif01() < pm) random_feasible_assign(I,S,L);
    }
}

void ga_run(Instance* I, GAParams P, Timetable* out_best, int k_export, Timetable* out_front){
    Timetable* pop  = calloc(P.pop_size, sizeof(Timetable));
    Timetable* next = calloc(P.pop_size, sizeof(Timetable));

    for(int i=0;i<P.pop_size;++i) init_random(I, &pop[i]);
    Timetable best = pop[0];

    for(int it=0; it<P.iters; ++it){
        int*    rank = calloc(P.pop_size, sizeof(int));
        double* crowd= calloc(P.pop_size, sizeof(double));
        ParetoView V = {.pop=pop, .size=P.pop_size, .rank=rank, .crowding=crowd};
        non_dominated_sort(&V); compute_crowding(&V);

        for(int i=0;i<P.pop_size;++i) if(dominates(&pop[i], &best)) best = pop[i];

        for(int i=0;i<P.pop_size;++i){
            int a = rand()%P.pop_size, b = rand()%P.pop_size;
            int p1 = binary_tournament(&V,a,b);
            a = rand()%P.pop_size; b = rand()%P.pop_size;
            int p2 = binary_tournament(&V,a,b);

            crossover(&pop[p1], &pop[p2], &next[i], I);
            mutate(I, &next[i], P.pm);
            evaluate(I, &next[i]);
        }
        free(rank); free(crowd);

        Timetable* tmp=pop; pop=next; next=tmp;
    }

    *out_best = best;

    // Export some rank-0 individuals as seeds
    int*    rank = calloc(P.pop_size, sizeof(int));
    double* crowd= calloc(P.pop_size, sizeof(double));
    ParetoView V = {.pop=pop, .size=P.pop_size, .rank=rank, .crowding=crowd};
    non_dominated_sort(&V);
    int count=0;
    for(int i=0;i<P.pop_size && count<k_export;++i){
        if(V.rank[i]==0) out_front[count++] = pop[i];
    }
    for(;count<k_export;++count) out_front[count]=best;

    free(rank); free(crowd);
    free(pop); free(next);
}
