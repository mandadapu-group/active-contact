// created by Amaresh Sahu
// amaresh.sahu@berkeley.edu
// 28 Sept. 2018


// required external libraries
#include <cmath>
#include <math.h>
#include <random>
#include <vector>
#include <cassert>
#include <iostream>
#include <stdlib.h>

// physical values
#define KBT             4.19     // pN.nm (30 celsius)
#define K_BENDING       120.     // pN.nm (~26 k_B T), note k_bending = 2 * kappa
#define LAMBDA          4.0e-3   // pN/nm
#define LENGTH          25133.   // nm (2 * pi * R_0)
#define VISCOSITY       7.972e-4 // pg/nm/us
#define FINAL_TIME      7.0e6    // us

// simulation constants
#define NUM_MODES       50
#define NUM_SMALL_MODES 20
#define EPSILON         1.0e-5

// active particle parameters
#define NUM_PARTICLES   7
#define TAU_R           5.e5     // us
#define TAU_P           5.e4     // us
#define TRAVERSAL_TIME  5.e5     // us, time for particle to traverse vesicle
#define PARTICLE_RADIUS 250.     // nm
#define VESICLE_RADIUS  4.e3     // nm
#define PARTICLE_SPEED  1.6e-2   // nm/us, particle speed

// option flags
#define SHOW_OUTPUT     true
#define AVERAGE_MODES   true
#define INCLUDE_ACTIVE  true


