#pragma once
#include "model.h"

typedef struct {
    int pop_size;
    int iters;
    double pm; // mutation rate per lecture
} GAParams;

void ga_run(Instance* I, GAParams P, Timetable* out_best, int k_export, Timetable* out_front);
