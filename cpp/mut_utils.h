/*
    Utils: Auxiliary functions for imperf.cpp
    @file utils.h
    @author Akshay Avvaru
    @version 0.1 06/08/2020
*/
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <bitset>
#include <fstream>
#include <chrono>
#include <string>
#include "utils.h"

using namespace std;
using namespace chrono;
using namespace utils;


// global variables
unordered_map<string, string> rClassMap;
unordered_map<string, utils::RepeatTracker> globalRepeatTracker;
// running individually for each motif is faster than
// running than in parllel
const uint motif_size = 4;
const float fraction_mutations = 0.1;
const bool debug = 1;

namespace mut_utils {

}