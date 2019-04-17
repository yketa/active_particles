#ifndef CMKDE_H
#define CMKDE_H

// ---- STRUCTURES ----

typedef struct SampleData{
  // Structure containing all data relevant to computation in this library.
  // Reflected by active_particles.mkde._SampleData.
  //
  // data sample
  int n;              // number of data points
  int d;              // number of dimensions
  double **data;      // array of data points
  // computed asymptotic mean integrated squared error (AMISE) and its derivatives with respect to the bandwidths
  double *AMISE;      // AMISE
  double *gradAMISE;  // gradient of AMISE
  double **hessAMISE; // Hessian matrix of AMISE
  // optimised bandwidths
  double *h;          // bandwidths
} SampleData;

// ---- PROTOTYPES ----

// AMISE AND ITS GRADIENT AND HESSIAN MATRIX

void AMISE(SampleData *sd);

void gradAMISE(SampleData *sd);

void hessAMISE(SampleData *sd);

// SIMPLIFYING FUNCTIONS

double phi(double x);

#endif
