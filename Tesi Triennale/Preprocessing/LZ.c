#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* LZ is an algorithm that, givem an encoded time series Xn (either binary or ternary) 
 * calculate the Lempel - Ziv coefficient c.
 *
 * This is a MEX-file for MATLAB.
 */

double LZ(char *Y, int k)
{
    size_t n = strlen(Y);

    // s array
    char *s = malloc(n * n);
    if (s == NULL)
    {
        free(s);
        return 1;
    }
    memset(s, '\0', sizeof(n*n));
    s[0] = Y[0];

    // q array
    char *q = malloc(n * n);
    if (q == NULL)
    {
        free(q);
        return 1;
    }
    memset(q, '\0', sizeof(n * n));
    double c = 1.0;
    int indexQ = 0;
    for (int i = 1; i < n; i++)
    {

        //  q=[q Sn(1,i)];
        q[indexQ] = Y[i];
        q[indexQ + 1] = '\0';
        indexQ++;
        
        char *temp = calloc(strlen(s) + 1, sizeof *s);
        if (temp == NULL)
        {
            free(temp);
            return 1;
        }
        memset(temp, '\0', strlen(s)+1);
        memcpy(temp, s, strlen(s));

        //         sq=[s q];
        char *seq = calloc(strlen(q) + strlen(s) + 1, sizeof *s);
        if (seq == NULL)
        {
            free(seq);
            return 1;
        }

        //         sq1=sq(1,1:end-1);
        memset(seq, '\0', strlen(q) + strlen(s) + 1);
        memcpy(seq, s, strlen(s));
        seq = strncat(seq, q, strlen(seq) + strlen(q));
        seq[strlen(seq) - 1] = '\0';

        //         if isempty(strfind(sq1,q))
        //             s=[s q];
        //             c=c+1;
        //             q=[];
        memset(s, '\0', strlen(s));
        memcpy(s, temp, strlen(temp));
        free(temp);
        if (!strstr(seq, q))
        {
            s = strncat(s, q, strlen(q) + strlen(s));
            memset(q, '\0', strlen(q));
            indexQ = 0;
            c++;
        }
    }
    if(strlen(s) != n) c++;
    free(q);
    free(s);
    //c++;
    double b = n / (log10(n) / log10(k));
    c = c / b;
    return c;
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    char *Y = calloc(100, sizeof(char));
    int k;
    double c;
    mxArray *mx;
    mexCallMATLAB(1,&mx,1,(mxArray **)prhs,"char");
    /* Check for proper number of arguments. */
    if (nrhs != 2)
    {
        mexErrMsgTxt("Two input required: Data Vector (string) and number of symbol (2 or 3).");
    }
    else if (nlhs > 1)
    {
        mexErrMsgTxt("At most one output arguments required.");
    }

    /* The first input must be a noncomplex vector string.*/
    if (!mxIsClass(prhs[0], "string"))
    {
        mexErrMsgTxt("Input must be a string vector!.");
    }

    /* The second impust must be a scalar value: 2 or 3*/
    if (mxIsComplex(prhs[1]) || !(mxGetM(prhs[1]) == 1 && mxGetN(prhs[1]) == 1))
    {
        mexErrMsgTxt("Input must be a noncomplex  scalar (2 or 3)!.");
    }

    Y = mxArrayToString(mx);  /*DATA VECTOR*/       
    mxDestroyArray(mx);
    k = mxGetScalar(prhs[1]);  /*number of coding symbols*/

    if (!(k == 2 || k == 3))
    {
        mexErrMsgTxt("The second imput must be a noncomplex  scalar (2 or 3)!.");
    }

    /* call the subroutine */
    c = LZ(Y, k);
    plhs[0] = mxCreateDoubleScalar(c);
    mxFree(Y);

    return;
}
