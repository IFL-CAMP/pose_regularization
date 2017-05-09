/* MATLAB wrapper for pose denoising
 * 
 * by Maximilian Baust on
 *
 * 05-09-2017
 * 
 */

#include "mex.h"
#include "../include/manifold.h"

using namespace manifold;

/* encapsulation of algorithm */
void denoiseSE3prod( double* inputArray,
                     double* outputArray,
                     int gridDimX,
                     int gridDimY,
                     int gridDimZ,
                     int p,
                     int q,
                     int r,
					 double inner_factor,
                     double alpha,
                     double beta,
					 double gamma,
                     int steps,
					 int inner_steps)
{
    proximalPointAlgorithm( inputArray, outputArray, 
                            gridDimX, gridDimY, gridDimZ,
                            p, q, r, inner_factor, alpha, beta,gamma, steps, inner_steps );
    
} /* end encapsulation of algorithm */

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

	int p = 1;					/* exponent for data term */
	int q = 1;                  /* exponent for 1st order regularizer */
	int r = 1;                  /* exponent for 2nd order regularizer */
	double inner_factor;        /* scale for the tau - second order */
	double alpha;				/* weighting for 1st order regularization */
	double beta;				/* weighting for 2nd order regularization */
	double gamma;				/* weighting for mixed 2nd order regularization */
	int gridDimX = 1;           /* number of columns of data grid */
	int gridDimY = 1;			/* number of fows of data grid */
	int gridDimZ = 1;           /* number of slices of data grid */
	int steps = 1;              /* number of iteration steps */
	int inner_steps = 1;        /* number of inner iteration steps */
	
    double* inMatrix;           /* MxN input matrix */
    double* outMatrix;          /* MxN output matrix */

    /* check number of ouput arguments */
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:nlhs","One output required!");
    }
    
    /* get element length */
    int elementLength = (int)mxGetM(prhs[0]);
    
    /* check number of input arguments */
    switch ( nrhs )
    {
        /* 3D case */
        case 13:
            
            /* create a pointer to the real data in the input matrix  */
            inMatrix = mxGetPr(prhs[0]);
            
            /* check element length */
            if ( wrongElementLength( elementLength ) ) {
                mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:elementLength","Element length is not correct!");
            }
            
			/* get all parameters */
			gridDimX = (int)mxGetScalar(prhs[1]);
			gridDimY = (int)mxGetScalar(prhs[2]);
			gridDimZ = (int)mxGetScalar(prhs[3]);
			p = (int)mxGetScalar(prhs[4]);
			q = (int)mxGetScalar(prhs[5]);
			r = (int)mxGetScalar(prhs[6]);
			inner_factor = (double)mxGetScalar(prhs[7]);
			alpha = (double)mxGetScalar(prhs[8]);
			beta = (double)mxGetScalar(prhs[9]);
			gamma = (double)mxGetScalar(prhs[10]);
			steps = (int)mxGetScalar(prhs[11]);
			inner_steps = (int)mxGetScalar(prhs[12]);
            
            /* create the output matrix */
            plhs[0] = mxCreateDoubleMatrix((mwSize)elementLength,(mwSize)(gridDimX*gridDimY*gridDimZ), mxREAL);

            /* get a pointer to the real data in the output matrix */
            outMatrix = mxGetPr(plhs[0]);
            
            break;
        
        /* 2D case */
        case 12:
            
            /* create a pointer to the real data in the input matrix  */
            inMatrix = mxGetPr(prhs[0]);
            
            /* check element length */
            if ( wrongElementLength( elementLength ) ) {
                mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:elementLength","Element length is not correct!");
            }
            
			/* get all parameters */
			gridDimX = (int)mxGetScalar(prhs[1]);
			gridDimY = (int)mxGetScalar(prhs[2]);
			p = (int)mxGetScalar(prhs[3]);
			q = (int)mxGetScalar(prhs[4]);
			r = (int)mxGetScalar(prhs[5]);
			inner_factor = (double)mxGetScalar(prhs[6]);
			alpha = (double)mxGetScalar(prhs[7]);
			beta = (double)mxGetScalar(prhs[8]);
			gamma = (double)mxGetScalar(prhs[9]);
			steps = (int)mxGetScalar(prhs[10]);
			inner_steps = (int)mxGetScalar(prhs[11]);
            
            /* create the output matrix */
            plhs[0] = mxCreateDoubleMatrix((mwSize)elementLength,(mwSize)(gridDimX*gridDimY), mxREAL);

            /* get a pointer to the real data in the output matrix */
            outMatrix = mxGetPr(plhs[0]);
            
            break;
            
        /* 1D case */
        case 11:
            
            /* create a pointer to the real data in the input matrix  */
            inMatrix = mxGetPr(prhs[0]);
            
            /* check element length */
            if ( wrongElementLength( elementLength ) ) {
                mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:elementLength","Element length is not correct!");
            }
            
			/* get all parameters */
			gridDimX = (int)mxGetScalar(prhs[1]);
			p = (int)mxGetScalar(prhs[2]);
			q = (int)mxGetScalar(prhs[3]);
			r = (int)mxGetScalar(prhs[4]);
			inner_factor = (double)mxGetScalar(prhs[5]);
			alpha = (double)mxGetScalar(prhs[6]);
			beta = (double)mxGetScalar(prhs[7]);
			gamma = (double)mxGetScalar(prhs[8]);
			steps = (int)mxGetScalar(prhs[9]);
			inner_steps = (int)mxGetScalar(prhs[10]);
            
            /* create the output matrix */
            plhs[0] = mxCreateDoubleMatrix((mwSize)elementLength,(mwSize)(gridDimX), mxREAL);

            /* get a pointer to the real data in the output matrix */
            outMatrix = mxGetPr(plhs[0]);
            
            break;
       
        /* else */
        default:
            mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:nrhs","At least nine input parameters are required!");

    } /* end check number of input arguments */

    /* call encapsulation of proximal point algorithm */
	denoiseSE3prod( inMatrix, outMatrix, gridDimX, gridDimY, gridDimZ, p, q, r, inner_factor, alpha, beta, gamma, steps, inner_steps );
    
} /* end gateway function */
