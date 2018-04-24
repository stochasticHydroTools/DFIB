/*
  TetrahedronSphereVolume.c
 

  Created by yuanxun bao on 7/16/14.

*/

#include <stdio.h>
#include <math.h>
#include <mex.h>

double pivot[3] = {0, 0, 0};

double dotProduct(double u[], double v[]){
  double d;
  d = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
  return d;
}

void crossProduct(double u[], double v[], double w[]){
  
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
  
}

void TetrahedronSphereVolume(double *V, double *T, int nV, int nT, double *vol){
  
  int t[3];
  double X1[3], X2[3], X3[3];
  double v1[3], v2[3], v3[3];
  double w[3];
  int j,k;

  
  for(k = 0; k < nT; k++){
    
    /* indices of the 3 vertices of the l-th trianle */
    t[0] = (int) T[k] - 1;
    t[1] = (int) T[k + nT] - 1;
    t[2] = (int) T[k + 2*nT] - 1;
    
    /* fill the 3 vertices from V */
    for (j = 0; j<3; j++){
      X1[j] = V[ j*nV + t[0]];
      X2[j] = V[ j*nV + t[1]];
      X3[j] = V[ j*nV + t[2]];
    }
    
    /* get 3 vector of a tetrahedron */
    for (j = 0; j<3; j++){
      v1[j] = X1[j] - pivot[j];
      v2[j] = X2[j] - pivot[j];
      v3[j] = X3[j] - pivot[j];
    }
    
    crossProduct(v2,v3,w);
    (*vol) += 1./6 * dotProduct(v1,w);
    
    /* printf("%f \n", *vol); */
  
  }
  
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  
  int numTriangle, numVertex;
  double *Vertex, *Triangle;
  double *Volume;
  
  
  Vertex      = mxGetPr(prhs[0]);
  Triangle    = mxGetPr(prhs[1]);
  numVertex   = mxGetScalar(prhs[2]);
  numTriangle = mxGetScalar(prhs[3]);
  
  plhs[0]     = mxCreateDoubleMatrix(1, 1, mxREAL);
  Volume      = mxGetPr(plhs[0]);
  
  TetrahedronSphereVolume(Vertex, Triangle, numVertex, numTriangle, Volume);
  
}

