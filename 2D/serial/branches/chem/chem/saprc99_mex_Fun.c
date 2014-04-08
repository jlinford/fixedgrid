/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Matlab Gateway for the Derivative Function Fun 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "mex.h"
#define min( x, y ) (x) < (y) ? (x) : (y)
#define max( x, y ) (x) > (y) ? (x) : (y)

void mexFunction( int nlhs, mxArray *plhs[], 
                     int nrhs, const mxArray *prhs[] )
{
 int mrows, mcols;
 double *V, *F, *RCT, *Vdot;

/* Check for the right number and size of input arguments */
 if ( nrhs != 3 ) {
   mexErrMsgTxt("saprc99_Fun requires 3 input vectors: V(74), F(5), RCT(211)");
 }
 mrows =  mxGetM(prhs[0]); mcols = mxGetN(prhs[0]);
 if ( ( mrows != 74 )||( mcols != 1 ) ) {
   mexPrintf("First saprc99_Fun input argument is of size V(%d,%d).",  
               mrows, mcols);
   mexErrMsgTxt("First saprc99_Fun input argument should be a column vector V(74,1)");
 }
 mrows =  mxGetM(prhs[1]); mcols = mxGetN(prhs[1]);
 if ( ( mrows != 5 )||( mcols != 1 ) ) {
   mexPrintf("Second saprc99_Fun input argument is of size F(%d,%d).",  
               mrows, mcols);
   mexErrMsgTxt("Second saprc99_Fun input argument should be a column vector F(5,1)");
 }
 mrows =  mxGetM(prhs[2]); mcols = mxGetN(prhs[2]);
 if ( (  mrows != 211 )||( mcols != 1 ) ) {
   mexPrintf("Third saprc99_Fun input argument is of size RCT(%d,%d).",  
               mrows, mcols);
   mexErrMsgTxt("Third saprc99_Fun input argument should be a column vector RCT(211,1)");
 }
 
/* Check for the right number of output arguments */
 if ( nlhs != 1 ) {
   mexErrMsgTxt("saprc99_Fun requires 1 output column vector: Vdot(74)");
 }

 V   = mxGetPr(prhs[0]);
 F   = mxGetPr(prhs[1]);
 RCT = mxGetPr(prhs[2]);

 plhs[0] = mxCreateDoubleMatrix(74,1,mxREAL);
 Vdot = mxGetPr(plhs[0]);

 Fun( V, F, RCT, Vdot );

}
