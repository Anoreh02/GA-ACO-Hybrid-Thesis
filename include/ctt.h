#pragma once
#include "model.h"

// Loads a CB-CTT .ctt file into Instance.
// Returns 0 on success, <0 on failure.
int load_ctt(const char* path, Instance* I);
