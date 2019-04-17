// Extension cmkde.c (together with cmkde.h) provides functions to efficiently
// compute the asymptotic mean squared integrated error (AMISE), its gradient
// and Hessian matrix. It is meant to be used by the active_particles.mkde
// Python module. (see
// https://yketa.github.io/UBC_2018_Wiki/#Multivariate%20kernel%20density%20estimation)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "cmkde.h"

#define M_PI 3.14159265358979323846264338327

// AMISE AND ITS GRADIENT AND HESSIAN MATRIX

void AMISE(SampleData *sd){
  // Compute asymptotic mean integrated squared error (AMISE) at sd->AMISE.
  //
  // sd [pointer] : Data structure.

  double sumBC = 0;

  double deltaijk;

  double sumsq;
  double Bij;
  double Bijcoeff[2];
  Bijcoeff[0] = -(2*(double)sd->d + 4);
  Bijcoeff[1] = pow((double)sd->d, 2) + 2*(double)sd->d;
  double Cij;

  for (int i = 0; i < sd->n; i++){
    for (int j = i + 1; j < sd->n; j++){

      sumsq = 0;
      Cij = 1;

      for (int k = 0; k < sd->d; k++){

        deltaijk = (sd->data[i][k] - sd->data[j][k])/sd->h[k];
        sumsq += pow(deltaijk, 2);

        Cij *= phi(deltaijk);

      }

      Bij = pow(sumsq, 2) + Bijcoeff[0]*sumsq + Bijcoeff[1];
      sumBC += Bij*Cij;

    }
  }

  double A = 1;
  for (int k = 0; k < sd->d; k++){
    A /= sd->h[k];
  }

  sd->AMISE[0] = A*(1/(pow(2*sqrt(M_PI), (double)sd->d)*(double)sd->n)
    + sumBC/(2*(double)sd->n*((double)sd->n - 1)));

}

void gradAMISE(SampleData *sd){
  // Compute gradient of AMISE at sd->gradAMISE.
  //
  // sd [pointer] : Data structure.

  double sumBC = 0;
  double *sumpBCBpC;
  sumpBCBpC = (double *)malloc(sd->d*sizeof(double));
  for (int p = 0; p < sd->d; p++){
    sumpBCBpC[p] = 0;
  }

  double deltaijk;

  double sumsq;
  double Bij;
  double *pBij;
  pBij = (double *)malloc(sd->d*sizeof(double));
  double Bijcoeff[2];
  Bijcoeff[0] = -(2*(double)sd->d + 4);
  Bijcoeff[1] = pow((double)sd->d, 2) + 2*(double)sd->d;
  double Cij;
  double *pCij;
  pCij = (double *)malloc(sd->d*sizeof(double));

  for (int i = 0; i < sd->n; i++){
    for (int j = i + 1; j < sd->n; j++){

      sumsq = 0;
      Cij = 1;

      for (int k = 0; k < sd->d; k++){

        deltaijk = (sd->data[i][k] - sd->data[j][k])/sd->h[k];
        sumsq += pow(deltaijk, 2);

        pBij[k] = pow(deltaijk, 2)/sd->h[k];
        Cij *= phi(deltaijk);
        pCij[k] = pow(deltaijk, 2)/sd->h[k];

      }

      Bij = pow(sumsq, 2) + Bijcoeff[0]*sumsq + Bijcoeff[1];
      sumBC += Bij*Cij;

      for (int p = 0; p < sd->d; p++){

        pBij[p] *= -4*sumsq + 4*(double)sd->d + 8;
        pCij[p] *= Cij;

        sumpBCBpC[p] += pBij[p]*Cij + Bij*pCij[p];

      }

    }
  }

  double A = 1;
  for (int k = 0; k < sd->d; k++){
    A /= sd->h[k];
  }
  double pA;

  for (int p = 0; p < sd->d; p++){

    pA = -A/sd->h[p];

    sd->gradAMISE[p] = pA/(pow(2*sqrt(M_PI), (double)sd->d)*(double)sd->n)
      + (pA*sumBC + A*sumpBCBpC[p])/(2*(double)sd->n*((double)sd->n - 1));

  }

  free(sumpBCBpC);
  free(pBij);
  free(pCij);

}

void hessAMISE(SampleData *sd){
  // Compute Hessian matrix of AMISE at sd->hessAMISE.
  //
  // sd [pointer] : Data structure.

  double sumBC = 0;
  double *sumpBCBpC;
  sumpBCBpC = (double *)malloc(sd->d*sizeof(double));
  double **sumpqBCpBqCqBpCBpqC;
  sumpqBCpBqCqBpCBpqC = (double **)malloc(sd->d*sizeof(double *));
  for (int p = 0; p < sd->d; p++){
    sumpBCBpC[p] = 0;
    sumpqBCpBqCqBpCBpqC[p] = (double *)malloc(sd->d*sizeof(double));
    for (int q = 0; q < sd->d; q++){
      sumpqBCpBqCqBpCBpqC[p][q] = 0;
    }
  }

  double deltaijk;
  double deltaijp;
  double deltaijq;

  double sumsq;
  double Bij;
  double *pBij;
  pBij = (double *)malloc(sd->d*sizeof(double));
  double **pqBij;
  pqBij = (double **)malloc(sd->d*sizeof(double *));
  for (int p = 0; p < sd->d; p++){
    pqBij[p] = (double *)malloc(sd->d*sizeof(double));
  }
  double Bijcoeff[2];
  Bijcoeff[0] = -(2*(double)sd->d + 4);
  Bijcoeff[1] = pow((double)sd->d, 2) + 2*(double)sd->d;
  double Cij;
  double *pCij;
  pCij = (double *)malloc(sd->d*sizeof(double));
  double **pqCij;
  pqCij = (double **)malloc(sd->d*sizeof(double *));
  for (int p = 0; p < sd->d; p++){
    pqCij[p] = (double *)malloc(sd->d*sizeof(double));
  }

  for (int i = 0; i < sd->n; i++){
    for (int j = 0; j < sd->n; j++){

      sumsq = 0;
      Cij = 1;

      for (int k = 0; k < sd->d; k++){

        deltaijk = (sd->data[i][k] - sd->data[j][k])/sd->h[k];
        sumsq += pow(deltaijk, 2);

        pBij[k] = pow(deltaijk, 2)/sd->h[k];
        Cij *= phi(deltaijk);
        pCij[k] = pow(deltaijk, 2)/sd->h[k];

      }

      Bij = pow(sumsq, 2) + Bijcoeff[0]*sumsq + Bijcoeff[1];
      sumBC += Bij*Cij;

      for (int p = 0; p < sd->d; p++){

        pBij[p] *= -4*sumsq + 4*(double)sd->d + 8;
        pCij[p] *= Cij;

        sumpBCBpC[p] += pBij[p]*Cij + Bij*pCij[p];

      }

      for (int p = 0; p < sd->d; p++){

        deltaijp = (sd->data[i][p] - sd->data[j][p])/sd->h[p];

        for (int q = p; q < sd->d; q++){

          deltaijq = (sd->data[i][q] - sd->data[j][q])/sd->h[q];

          pqBij[p][q] = 8*pow(deltaijp, 2)*pow(deltaijq, 2)
            /(sd->h[p]*sd->h[q]);
          pqCij[p][q] = pow(deltaijp, 2)*pow(deltaijq, 2)*Cij
            /(sd->h[p]*sd->h[q]);
          if ( p == q ){
            pqBij[p][q] += pow(deltaijp/sd->h[p], 2)
              *(12*sumsq - 6*(2*(double)sd->d + 4));
            pqCij[p][q] -= 3*pow(deltaijp/sd->h[p], 2)*Cij;
          }

          sumpqBCpBqCqBpCBpqC[p][q] += pqBij[p][q]*Cij + pBij[p]*pCij[q]
            + pBij[q]*pCij[p] + Bij*pqCij[p][q];

        }

      }

    }
  }

  double A = 1;
  for (int k = 0; k < sd->d; k++){
    A /= sd->h[k];
  }
  double pA;
  double qA;
  double pqA;

  for (int p = 0; p < sd->d; p++){

    pA = -A/sd->h[p];

    for (int q = p; q < sd->d; q++){

      qA = -A/sd->h[q];
      pqA = -pA/sd->h[q];
      if ( p == q ){
        pqA *= 2;
      }

      sd->hessAMISE[p][q] =
        pqA/(pow(2*sqrt(M_PI), (double)sd->d)*(double)sd->n)
        + (pqA*sumBC + pA*sumpBCBpC[q] + qA*sumpBCBpC[p]
          + A*sumpqBCpBqCqBpCBpqC[p][q])/(2*(double)sd->n*((double)sd->n - 1));
      sd->hessAMISE[q][p] = sd->hessAMISE[p][q];

    }

  }

  free(sumpBCBpC);
  free(pBij);
  free(pCij);
  for (int p = 0; p < sd->d; p++){
    free(sumpqBCpBqCqBpCBpqC[p]);
    free(pqBij[p]);
    free(pqCij[p]);
  }
  free(sumpqBCpBqCqBpCBpqC);
  free(pqBij);
  free(pqCij);

}

// SIMPLIFYING FUNCTIONS

double phi(double x){
  // Standard normal in distribution evaluated in x.

  return exp(-0.5*pow(x, 2))/sqrt(2*M_PI);
}
