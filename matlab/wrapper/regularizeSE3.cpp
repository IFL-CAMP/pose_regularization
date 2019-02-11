/* MATLAB wrapper for pose denoising
 * 
 * by Maximilian Baust on
 *
 * 05-09-2017
 * 
 */

#include "pose_regularization.h"
#include "mex.h"
#include <cstring>

using namespace pose_regularization;

/* encapsulation of algorithm */
void regularizeSE3(double *inputArray,
                   double *outputArray,
                   size_t gridDimX,
                   Regularization p,
                   Regularization q,
                   Regularization r,
                   double inner_factor,
                   double alpha,
                   double beta,
                   int steps,
                   int inner_steps)
{
    std::vector<SE3_ELEMENT> input(gridDimX), output(gridDimX);
    for (int i = 0; i < gridDimX; i++) {
        std::memcpy(input[i].data(), inputArray + (i * 16), 16 * sizeof(double));
    }

    proximalPointAlgorithm(input, output,
                           p, q, r, inner_factor, alpha, beta, steps, inner_steps);

    for (int i = 0; i < gridDimX; i++) {
        std::memcpy(outputArray + (i * 16), output[i].data(), 16 * sizeof(double));
    }

} /* end encapsulation of algorithm */

/* gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    Regularization p = L1;        /* exponent for data term */
    Regularization q = L1;      /* exponent for 1st order regularizer */
    Regularization r = L1;      /* exponent for 2nd order regularizer */
    double inner_factor;        /* scale for the tau - second order */
    double alpha;                /* weighting for 1st order regularization */
    double beta;                /* weighting for 2nd order regularization */
    size_t gridDimX = 1;        /* number of columns of data grid */
    int steps = 1;              /* number of iteration steps */
    int inner_steps = 1;        /* number of inner iteration steps */

    double *inMatrix;           /* MxN input matrix */
    double *outMatrix;          /* MxN output matrix */

    /* check number of ouput arguments */
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:nlhs", "One output required!");
    }

    /* get element length */
    int elementLength = (int) mxGetM(prhs[0]);

    /* check number of input arguments */
    switch (nrhs) {

        /* 1D case */
        case 10:

            /* create a pointer to the real data in the input matrix  */
            inMatrix = mxGetPr(prhs[0]);

            /* check element length */
            if (elementLength != SE3_ELEMENT_LENGTH) {
                mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:elementLength", "Element length is not correct!");
            }

            /* get all parameters */
            gridDimX = (size_t) mxGetScalar(prhs[1]);
            p = (Regularization) mxGetScalar(prhs[2]);
            q = (Regularization) mxGetScalar(prhs[3]);
            r = (Regularization) mxGetScalar(prhs[4]);
            inner_factor = (double) mxGetScalar(prhs[5]);
            alpha = (double) mxGetScalar(prhs[6]);
            beta = (double) mxGetScalar(prhs[7]);
            steps = (int) mxGetScalar(prhs[8]);
            inner_steps = (int) mxGetScalar(prhs[9]);

            /* create the output matrix */
            plhs[0] = mxCreateDoubleMatrix((mwSize) elementLength, (mwSize)(gridDimX), mxREAL);

            /* get a pointer to the real data in the output matrix */
            outMatrix = mxGetPr(plhs[0]);

            break;

            /* else */
        default:
            mexErrMsgIdAndTxt("manifoldtoolbox:denoiseSE3prod:nrhs", "Nine input parameters are required!");

    } /* end check number of input arguments */

    /* call encapsulation of proximal point algorithm */
    regularizeSE3(inMatrix, outMatrix, gridDimX, p, q, r, inner_factor, alpha, beta, steps, inner_steps);

} /* end gateway function */
