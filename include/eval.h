#pragma once
#include "model.h"

typedef struct {
    double w1, w2, w3; // optional if you still use a scalarized debug score
} EvalWeights;

void evaluate(Instance* I, Timetable* S);

// Pareto dominance for (f1 minimize, f2 maximize, f3 minimize)
// Keep if you still want multi-objective mechanics alongside CB-CTT score
int  dominates(const Timetable* A, const Timetable* B);

// Convenience: compute single CB-CTT soft score (sum of soft penalties)
int  ctt_soft_score(const Timetable* S);
