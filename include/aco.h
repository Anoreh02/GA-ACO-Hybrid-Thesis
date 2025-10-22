#pragma once
#include "model.h"

typedef struct {
    int ants;
    int iters;
    double alpha, beta;   // tau^alpha * eta^beta
    double rho;           // evaporation (local)
    double init_tau;      // initial pheromone
    double selective_p;   // top p% candidate timeslots per lecture
} ACOParams;

void aco_refine(Instance* I, ACOParams P, int seeds, const Timetable* seed_solutions,
                Timetable* out_best);
