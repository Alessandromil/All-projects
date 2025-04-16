#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Code is an algorithm that, givem a time serie Xn, codes it into
 * a new serie Yn of same lenght N.
 *
 * This is a MEX-file for MATLAB.
 */



char* code(double *Sn,  int n, int s, double p)
{
	int j;
	double A,b;
  char *coded = calloc(n +1, sizeof(char));
  memset(coded, '\0', n);
	coded[0]='1';

	if (s==2){
		for (j=1;j<n;j++){
		b=(double)(Sn[j-1]*p);
		A =(double) b < 0 ? -b : b; /*absolute value*/
		if (Sn[j]<=Sn[j-1]+A){  /* code the Sn time series in 0 and 1 */
			coded[j]='0';}
		if (Sn[j]>Sn[j-1]+A){
			coded[j]='1';}
		}
	}
	else{
		for (j=1;j<n;j++){
			b=(double)(Sn[j-1]*p);
			A =(double) b < 0 ? -b : b; /*absolute value*/
			if (Sn[j]<=Sn[j-1]+A && Sn[j]>=Sn[j-1]-A){  /* code the Sn time series in 0 1 and 2*/
				coded[j]='2';}
			if (Sn[j]>Sn[j-1]+A){
				coded[j]='1';}
			if (Sn[j]<Sn[j-1]-A){
				coded[j]='0';}
		}
	}
  coded[j] = '\0';
  return coded;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double  *S,*v,perc;
  int N,symbol;


  /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgTxt("Three input required: Data Vector, number of symbol (2 or 3) and the percentange. ");
  } 
    //else if(nlhs!=1) {
    //mexErrMsgTxt("One output arguments required.");
  //}

  /* The first input must be a noncomplex vector double.*/

  if( mxIsComplex(prhs[0]) || !(mxGetM(prhs[0])>1 && mxGetN(prhs[0])==1) )
  {
    mexErrMsgTxt("Input must be a noncomplex  double vector!.");
  }

  /* The second impust must be a scalar value: 2 or 3*/
  if( mxIsComplex(prhs[1]) || !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1))
    {
      mexErrMsgTxt("Input must be a noncomplex  scalar (2 or 3)!.");
  }

  /* The second impust must be a scalar value: 2 or 3*/
    if( mxIsComplex(prhs[1]) || !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1))
      {
        mexErrMsgTxt("Input must be a noncomplex  scalar (2 or 3)!.");
  }
  /* The third impust must be a noncomplex  double number*/
      if( mxIsComplex(prhs[1]) || !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1) )
        {
          mexErrMsgTxt("Input must be a noncomplex  scalar >0 and <1.");
  }

  S=mxGetPr(prhs[0]); /*DATA VECTOR*/
  N=(int)mxGetM(prhs[0]); /*lenght of Data Vector*/
  symbol=(int)*mxGetPr(prhs[1]); /*number of symbols for coding*/
  perc=(double)*mxGetPr(prhs[2]);/* variation percentage accepted for the eveness*/

  if (!(symbol==2 || symbol==3)){
  	  mexErrMsgTxt("The second imput must be a noncomplex  scalar (2 or 3)!.");
  }
  if (!(perc>=0.0 || perc<1.0)){
  	  mexErrMsgTxt("The third imput must be a noncomplex  scalar (2 or 3)!.");
  }
  char *coded = code(S,N,symbol,perc);
  plhs[0]= mxCreateString(coded);


  return;
}
