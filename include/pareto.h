#pragma once
#include "model.h"

typedef struct {
    Timetable* pop;
    int size;
    int*    rank;      // lower better
    double* crowding;  // higher better
} ParetoView;

void non_dominated_sort(ParetoView* V);
void compute_crowding(ParetoView* V);
int  binary_tournament(ParetoView* V, int a, int b);
