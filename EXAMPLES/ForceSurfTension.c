/*
    ForceSurfTension.c
 
 
    Created by yuanxun bao on 7/15/14.

*/

#include <mex.h>
#include <math.h>
#include <stdio.h>

double e1[3] = {1, 0, 0};
double e2[3] = {0, 1, 0};
double e3[3] = {0, 0, 1};

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

void ForceSurfTension(double *V, double *T, int nV, int nT, double *F){
  
  double X1[3], X2[3], X3[3];
  int t[3];
  double v1[3], v2[3], v3[3], w[3];
  double u1[3], u2[3], u3[3];
  double areaCross;
  int l,j;

  
  for(l = 0; l < nT; l++){
    
    /* indices of the 3 vertices of the l-th trianle */
    t[0] = (int) T[l] - 1;
    t[1] = (int) T[l + nT] - 1;
    t[2] = (int) T[l + 2*nT] - 1;
    
    /* printf("%d, %d, %d\n",t[0],t[1],t[2]); */
    
    /* fill the 3 vertices from V */
    for (j = 0; j<3; j++){
      X1[j] = V[ j*nV + t[0]];
      X2[j] = V[ j*nV + t[1]];
      X3[j] = V[ j*nV + t[2]];
    }
    
    /* printf("%f, %f, %f\n",X1[0],X1[1],X1[2]); */
    
    /* get 3 edges of a triangle */
    for (j = 0; j<3; j++){
      v1[j] = X2[j] - X3[j];
      v2[j] = X3[j] - X1[j];
      v3[j] = X1[j] - X2[j];
    }
    
    /*  |(X1-X2) x (X2-X3)| */
    crossProduct(v3,v1,w);
    areaCross = sqrt( w[0]*w[0] + w[1]*w[1] + w[2] * w[2] );
    
    /* contribution for F from X1 */
    crossProduct(e1,v1,u1);
    crossProduct(e2,v1,u2);
    crossProduct(e3,v1,u3);
    
    F[ t[0]        ] += 0.5 / areaCross * dotProduct(w,u1);
    F[ t[0] + nV   ] += 0.5 / areaCross * dotProduct(w,u2);
    F[ t[0] + 2*nV ] += 0.5 / areaCross * dotProduct(w,u3);
    
    /* contribution for F from X2 */
    crossProduct(e1,v2,u1);
    crossProduct(e2,v2,u2);
    crossProduct(e3,v2,u3);
    
    F[ t[1]        ] += 0.5 / areaCross * dotProduct(w,u1);
    F[ t[1] + nV   ] += 0.5 / areaCross * dotProduct(w,u2);
    F[ t[1] + 2*nV ] += 0.5 / areaCross * dotProduct(w,u3);
    
    /* contribution for F from X3 */
    crossProduct(e1,v3,u1);
    crossProduct(e2,v3,u2);
    crossProduct(e3,v3,u3);
    
    F[ t[2]        ] += 0.5 / areaCross * dotProduct(w,u1);
    F[ t[2] + nV   ] += 0.5 / areaCross * dotProduct(w,u2);
    F[ t[2] + 2*nV ] += 0.5 / areaCross * dotProduct(w,u3);

    
  }
  
  /* minus sign for -dE/dX */
  for(j = 0; j < 3*nV; j++){
    F[j] = -F[j];
  }

}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  
  int numVertex, numTriangle;
  double *Vertex, *Triangle;
  double *Force;
  
  Vertex      = mxGetPr(prhs[0]);
  Triangle    = mxGetPr(prhs[1]);
  numVertex   = mxGetScalar(prhs[2]);
  numTriangle = mxGetScalar(prhs[3]);
  
  plhs[0]     = mxCreateDoubleMatrix((mwSize)numVertex, 3, mxREAL);
  Force       = mxGetPr(plhs[0]);
  
  /*
  for(int j = 0; j < numVertex; j++){
    printf("%f, %f %f \n", Vertex[j], Vertex[j+numVertex], Vertex[j+2*numVertex]);
  }
  */

  ForceSurfTension(Vertex, Triangle,numVertex, numTriangle, Force);


}