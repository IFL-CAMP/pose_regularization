/* lightweight implementation for processing manifold valued data
 *
 * written by Maximilian Baust
 *
 * e-mail: maximilian.baust@tum.de
 *
 * 08-16-2014
 *
 * This file has been simplified for pose denoising
 *
 * by Maximilian Baust on
 *
 * 05-09-2017
 *
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <manifold.h>

constexpr double SE3_EPSILON = 0.00000001;

namespace manifold {

/***************** A*B matrix multiplication for SO(3) ******************/
    void AmultB( double * A, double * B, double * C )
    {
        C[0] = A[0]*B[0] + A[4]*B[1] + A[8]*B[2];
        C[1] = A[1]*B[0] + A[5]*B[1] + A[9]*B[2];
        C[2] = A[2]*B[0] + A[6]*B[1] + A[10]*B[2];
        C[3] = 0.0;
 
        C[4] = A[0]*B[4] + A[4]*B[5] + A[8]*B[6];
        C[5] = A[1]*B[4] + A[5]*B[5] + A[9]*B[6];
        C[6] = A[2]*B[4] + A[6]*B[5] + A[10]*B[6];
        C[7] = 0.0;
        
        C[8] = A[0]*B[8] + A[4]*B[9] + A[8]*B[10];
        C[9] = A[1]*B[8] + A[5]*B[9] + A[9]*B[10];
        C[10] = A[2]*B[8] + A[6]*B[9] + A[10]*B[10];
        C[11] = 0.0;
        
        C[12] = 0.0;
        C[13] = 0.0;
        C[14] = 0.0;
        C[15] = 1.0;
    };
    /************** end A*B matrix multiplication *****************/
    
    /**************** A^-1*B matrix multiplication for SO(3) ****************/
    void invAmultB( double * A, double * B, double * C )
    {
        C[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
        C[1] = A[4]*B[0] + A[5]*B[1] + A[6]*B[2];
        C[2] = A[8]*B[0] + A[9]*B[1] + A[10]*B[2];
        C[3] = 0.0;
        
        C[4] = A[0]*B[4] + A[1]*B[5] + A[2]*B[6];
        C[5] = A[4]*B[4] + A[5]*B[5] + A[6]*B[6];
        C[6] = A[8]*B[4] + A[9]*B[5] + A[10]*B[6];
        C[7] =  0.0;
        
        C[8] = A[0]*B[8] + A[1]*B[9] + A[2]*B[10];
        C[9] = A[4]*B[8] + A[5]*B[9] + A[6]*B[10];
        C[10] = A[8]*B[8] + A[9]*B[9] + A[10]*B[10];
        C[11] = 0.0;
        
        C[12] = 0.0;
        C[13] = 0.0;
        C[14] = 0.0;
        C[15] = 1.0;
    };
    /************** end A^-1*B matrix multiplication **************/


	/********************* matrix exponential *********************/
	void matrixExponential(double * Ai, double * Ao)
	{

		/* get omegas */
		double omega_x = Ai[6];
		double omega_y = Ai[8];
		double omega_z = Ai[1];

		/* calculate norm of omega */
		double no = std::sqrt(omega_x*omega_x + omega_y*omega_y + omega_z*omega_z);

		if (std::abs(no) > SE3_EPSILON)
		{

			/* calculate factors */
			double f1 = std::sin(no)/no;
			double f2 = (1.0 - std::cos(no)) / (no*no);


			/* calculate matrix exponential of destination element (tangential element, called v in our conversions) */
			Ao[0] = 1.0 - f2*(omega_y*omega_y + omega_z*omega_z);
			Ao[1] = f1*omega_z + f2*omega_x*omega_y;
			Ao[2] = f2*omega_x*omega_z - f1*omega_y;
			Ao[3] = 0.0;

			Ao[4] = f2*omega_x*omega_y - f1*omega_z;
			Ao[5] = 1.0 - f2*(omega_x*omega_x + omega_z*omega_z);
			Ao[6] = f1*omega_x + f2*omega_y*omega_z;
			Ao[7] = 0.0;

			Ao[8] = f1*omega_y + f2*omega_x*omega_z;
			Ao[9] = f2*omega_y*omega_z - f1*omega_x;
			Ao[10] = 1.0 - f2*(omega_x*omega_x + omega_y*omega_y);
			Ao[11] = 0.0;

			Ao[12] = 0.0;
			Ao[13] = 0.0;
			Ao[14] = 0.0;
			Ao[15] = 1.0;
		}
		else
		{
			Ao[0] = 1.0;
			Ao[1] = 0.0;
			Ao[2] = 0.0;
			Ao[3] = 0.0;

			Ao[4] = 0.0;
			Ao[5] = 1.0;
			Ao[6] = 0.0;
			Ao[7] = 0.0;

			Ao[8] = 0.0;
			Ao[9] = 0.0;
			Ao[10] = 1.0;
			Ao[11] = 0.0;

			Ao[12] = 0.0;
			Ao[13] = 0.0;
			Ao[14] = 0.0;
			Ao[15] = 1.0;

		} /* end if ( abs(no) > SE3_EPSILON ) */

	};
	/******************* end matrix exponential *******************/
    
    /********************* Schur decomposition *********************/
	double schurDecomposition(double * Ai, double * U)
	{

        /* get a, b and c of matrix
         * A = [0 -a b; a 0 -c ; -b c 0]
         */
        double a = Ai[1];
        double b = Ai[8];
        double c = Ai[6];
        
        /* compute beta and gamma */
        double beta = std::sqrt( a*a + b*b + c*c );
        double gamma = std::sqrt( a*a + b*b );
        
        /* compose U */
        if ( beta > SE3_EPSILON && gamma > SE3_EPSILON )
        {
            U[0] = c/beta;
            U[1] = b/beta;
            U[2] = a/beta;

            U[3] = 0.0;
            U[4] = a/gamma;
            U[5] = -b/gamma;

            U[6] = -gamma/beta;
            U[7] = (b*c)/(beta*gamma);
            U[8] = (a*c)/(beta*gamma);
        }
        else
        {
            
            U[0] = 1.0;
            U[1] = 0.0;
            U[2] = 0.0;
        
            U[3] = 0.0;
            U[4] = 1.0;
            U[5] = 0.0;
                
            U[6] = 0.0;
            U[7] = 0.0;
            U[8] = 1.0;
        } /* end of beta > eps */
        
		return beta;
		
	};
	/******************* end Schur decompositon *******************/
    
    /************* logarithmMap::operator() *************/
    double LogarithmMap::operator() ( double * baseElPtr,
                                      double * destElPtr,
                                      double * resultElPtr ) 
    {
        
        /* allocate memory */
        double C[_elementLength];
        
        /* compute inverted base element * dest element (only for SO(3)-component) */
        invAmultB( baseElPtr, destElPtr, C );
        
        /* calculate phi */
//         double phi =  acos( 0.5*(C[0] + C[5] + C[10] - 1.0) );
        double phi = (acos((std::complex<double>)(0.5*(C[0] + C[5] + C[10] - 1.0)))).real();
		
        /* calculate result for SO(3) part */
        if ( (std::abs(phi) > SE3_EPSILON) && (std::abs(phi) < M_PI-SE3_EPSILON) )
        {
            /* get omegas */
            double omega_x = 0.5*phi*(C[6] - C[9])/std::sin(phi);
            double omega_y = 0.5*phi*(C[8] - C[2])/std::sin(phi);
            double omega_z = 0.5*phi*(C[1] - C[4])/std::sin(phi);
            
            /* first column of result */
            resultElPtr[0] =  0.0;
            resultElPtr[1] =  omega_z;
            resultElPtr[2] = -omega_y;
            resultElPtr[3] =  0.0;
            
            /* second column of result */
            resultElPtr[4] = -omega_z;
            resultElPtr[5] =  0.0;
            resultElPtr[6] =  omega_x;
            resultElPtr[7] =  0.0;
            
            /* third column of result */
            resultElPtr[8] =  omega_y;
            resultElPtr[9] = -omega_x;
            resultElPtr[10] = 0.0;
            resultElPtr[11] = 0.0;
            
            /* fourth column of result is calculated later */

			
        }
        else
        {
            for ( int i = 0; i < _elementLength - 4; ++i )
            {
                resultElPtr[i] = 0.0;
            }
        } // end if ( abs(phi) > 0 )
        
               
		resultElPtr[12] = destElPtr[12] - baseElPtr[12];
		resultElPtr[13] = destElPtr[13] - baseElPtr[13];
		resultElPtr[14] = destElPtr[14] - baseElPtr[14];
		resultElPtr[15] = 0.0;
        
        /* compute length of logarithm without considering last element
         * (note that translational component does not have factor 1/2)
         */
        double dist = 0.0;
        for ( int i=0; i < _elementLength-4; ++i ) {
            dist += 0.5*resultElPtr[i]*resultElPtr[i];
        }
        for ( int i=_elementLength-4; i < _elementLength-1; ++i ) {
            dist += resultElPtr[i]*resultElPtr[i];
        }
        
        return std::sqrt( dist );
    };
    /*********** end logarithmMap::operator() ***********/
    
    /************* exponentialMap::operator() *************/
    void ExponentialMap::operator() ( double * baseElPtr,
                                      double * destElPtr,
                                      double * resultElPtr ) 
    {
       
		/* compute matrix exponential of destination element */
		double expTangElPtr[_elementLength];
		matrixExponential(destElPtr, expTangElPtr);

        /* multiply base element with matrix exponential of tangential element */
        AmultB( baseElPtr, expTangElPtr, resultElPtr );
        
        /* change fourth column of result */
		resultElPtr[12] = baseElPtr[12] + destElPtr[12];
		resultElPtr[13] = baseElPtr[13] + destElPtr[13];
		resultElPtr[14] = baseElPtr[14] + destElPtr[14];
		resultElPtr[15] = 1.0;

        
    };
    /*********** end exponentialMap::operator() ***********/
    
    /****************** innerProduct::operator() ******************/
    double InnerProduct::operator() ( double * elementPtr1,
                                      double * elementPtr2,
                                      double * baseElPtr )
    {
        
		double prod = 0.0;
		for (int i = 0; i < _elementLength - 4; ++i) {
			prod += 0.5*elementPtr1[i] * elementPtr2[i];
		}
        for (int i = _elementLength - 4; i < _elementLength - 1; ++i) {
			prod += elementPtr1[i] * elementPtr2[i];
		}
		return prod;
    };
    /**************** end innerProduct::operator() ****************/
    
    /*************** adjustCoefficients::operator() ***************/
    void AdjustCoefficients::operator() ( double * inputCoefficients,
                                          double * outputCoefficients,
                                          double * eigenvalues,
                                          double t )
    {  

		/* generic initialization */
		outputCoefficients[0] = 0.5*inputCoefficients[0];
		outputCoefficients[1] = 0.5*inputCoefficients[1];
		outputCoefficients[2] = 0.5*inputCoefficients[2];
		outputCoefficients[3] = 0.5*inputCoefficients[3];
		outputCoefficients[4] = 0.5*inputCoefficients[4];
		outputCoefficients[5] = 0.5*inputCoefficients[5];

		/* modify if necessary */
		if (t > SE3_EPSILON && eigenvalues[1] > SE3_EPSILON)
		{
			outputCoefficients[1] = inputCoefficients[1] * std::sin(std::sqrt(eigenvalues[1])*t / 2.0) / (std::sin(std::sqrt(eigenvalues[1])*t) + 0.000000000001);
		}
		if (t > SE3_EPSILON && eigenvalues[2] > SE3_EPSILON)
		{
			outputCoefficients[2] = inputCoefficients[2] * std::sin(std::sqrt(eigenvalues[2])*t / 2.0) / (std::sin(std::sqrt(eigenvalues[2])*t) + 0.000000000001);
		}
        
       
    };
    /************* end adjustCoefficients::operator() *************/
    
    /****************** computeFrame::operator() ******************/
    void ComputeFrame::operator() ( double * tagentSpaceElement,
                                    double * basisElements,
                                    double * eigenvalues,
									double * p) // p not used
    {
		/* Based on the tangential element this functor computes
         * the moving frame, i.e., the basis as well as the 
         * eigenvalues of the corresponding Jacobi field and stroes
         * them in the arrays:
         * 
         * basisElements (of size _elementLength*_dimension) and
         * 
         * eigenvalues (of size _dimension).
         *
         */
        
        /* perform Schur decompositon */
        double beta = 0.0;
        double Um[9];
        beta = schurDecomposition(tagentSpaceElement, Um);

		/* compute eigenvalues */
		eigenvalues[0] = 0.0;
		eigenvalues[1] = 0.25*beta*beta;
		eigenvalues[2] = 0.25*beta*beta;
		eigenvalues[3] = 0.0;
		eigenvalues[4] = 0.0;
		eigenvalues[5] = 0.0;
        
		/* first basis element */
		int index = 0;

		basisElements[index] = 0.0;
		basisElements[index + 1] = Um[3]*Um[7] - Um[6]*Um[4];
		basisElements[index + 2] = Um[3]*Um[8] - Um[6]*Um[5];
		basisElements[index + 3] = 0.0;

		basisElements[index + 4] = Um[6]*Um[4] - Um[3]*Um[7];
		basisElements[index + 5] = 0.0;
		basisElements[index + 6] = Um[4]*Um[8] - Um[7]*Um[5];
		basisElements[index + 7] = 0.0;

		basisElements[index + 8] = Um[6]*Um[5] - Um[3]*Um[8];
		basisElements[index + 9] = Um[7]*Um[5] - Um[4]*Um[8];
		basisElements[index + 10] = 0.0;
		basisElements[index + 11] = 0.0;

		basisElements[index + 12] = 0.0;
		basisElements[index + 13] = 0.0;
		basisElements[index + 14] = 0.0;
		basisElements[index + 15] = 0.0;


		/* second basis element */
		index = _elementLength;
		basisElements[index] = 0.0;
		basisElements[index + 1] = Um[0]*Um[4] - Um[3]*Um[1];
		basisElements[index + 2] = Um[0]*Um[5] - Um[3]*Um[2];
		basisElements[index + 3] = 0.0;

		basisElements[index + 4] = Um[3]*Um[1] - Um[0]*Um[4];
		basisElements[index + 5] = 0.0;
		basisElements[index + 6] = Um[1]*Um[5] - Um[4]*Um[2];
		basisElements[index + 7] = 0.0;

		basisElements[index + 8] = Um[3]*Um[2] - Um[0]*Um[5];
		basisElements[index + 9] = Um[4]*Um[2] - Um[1]*Um[5];
		basisElements[index + 10] = 0.0;
		basisElements[index + 11] = 0.0;

		basisElements[index + 12] = 0.0;
		basisElements[index + 13] = 0.0;
		basisElements[index + 14] = 0.0;
		basisElements[index + 15] = 0.0;


		/* third basis element */
		index = 2 * _elementLength;
		basisElements[index] = 0.0;
		basisElements[index + 1] = Um[0]*Um[7] - Um[6]*Um[1];
		basisElements[index + 2] = Um[0]*Um[8] - Um[6]*Um[2];
		basisElements[index + 3] = 0.0;

		basisElements[index + 4] = Um[6]*Um[1] - Um[0]*Um[7];
		basisElements[index + 5] = 0.0;
		basisElements[index + 6] = Um[1]*Um[8] - Um[7]*Um[2];
		basisElements[index + 7] = 0.0;

		basisElements[index + 8] = Um[6]*Um[2] - Um[0]*Um[8];
		basisElements[index + 9] = Um[7]*Um[2] - Um[1]*Um[8];
		basisElements[index + 10] = 0.0;
		basisElements[index + 11] = 0.0;

		basisElements[index + 12] = 0.0;
		basisElements[index + 13] = 0.0;
		basisElements[index + 14] = 0.0;
		basisElements[index + 15] = 0.0;

		for (int ind = 3 * _elementLength; ind < 6 * _elementLength; ++ind)
        {
			basisElements[ind] = 0.0;
        }

        /* eigenvectors corresponding to translation */
        
        basisElements[60] = 1.0;
        basisElements[61] = 0.0;
        basisElements[62] = 0.0;
        
        basisElements[76] = 0.0;
        basisElements[77] = 1.0;
        basisElements[78] = 0.0;
        
        basisElements[76] = 0.0;
        basisElements[77] = 0.0;
        basisElements[78] = 1.0;
        
    };
    /**************** end computeFrame::operator() ****************/
    
    /*************** parallelTransport::operator() ****************/
    void ParallelTransport::operator() ( double * inputFrameVectors,
                                         double * outputFrameVectors,
                                         double * baseElPtr,
                                         double * destElPtr )
    {	
		/* The parallel transport takes the inputFrameVectors 
         * (and array of size _elementLength*_dimension ), 
         * moves them into the direction log( baseElPtr, destElPtr ), 
         * and stores them in inputFrameVectors 
         * (of size _elementLength*_dimension ). 
         * The implementation goes along with the paper 
         * 
         * Q. Reentmesters 
         * A gradient method for geodesic data fitting on some
         * symmetric Riemannian manifolds
         * IEEE Conference on Decision and Control
         * and European Control Conference (2011)
         */
        
        /* at first we implement the SO(3) parallel transport */
                
		/* allocate memory */
		double v[_elementLength];
		double exp_aux[_elementLength];
		double inv_dest_base[_elementLength];
		double left_factor[_elementLength];
		double aux[_elementLength];
        
        /* initialize log map */
		LogarithmMap logMap(_elementLength);
        
        /* compute */
        logMap(baseElPtr, destElPtr, v);
		for (int i = 0; i < _elementLength; ++i)
		{
			v[i] *= 0.5;
		}

		/* keep only SO(3) part */
		v[12] = 0.0;
        v[13] = 0.0;
        v[14] = 0.0;
        v[15] = 0.0;
        
        /* compute matrix exponential */
		matrixExponential(v, exp_aux);

		/* compute inv(dest)*base */
		invAmultB(destElPtr, baseElPtr, inv_dest_base);
        
        /* compute left factor of product */
        AmultB(inv_dest_base, exp_aux, left_factor);

		/* compute parallel transport for all basis elements */
		for (int i = 0; i < _dimension-3; ++i)
		{
			int index = i*_elementLength;
            AmultB(left_factor, &inputFrameVectors[index], aux);
			AmultB(aux, exp_aux, &outputFrameVectors[index]);
            
            outputFrameVectors[index + _elementLength - 4] = 0.0; 
            outputFrameVectors[index + _elementLength - 3] = 0.0; 
            outputFrameVectors[index + _elementLength - 2] = 0.0; 
            outputFrameVectors[index + _elementLength - 1] = 0.0; 
            
		}

      
        /* compute parallel transport for translational elements (flat geometry!) */
		for (int i = 3 * _elementLength; i < 6 * _elementLength; ++i)
        {
			outputFrameVectors[i] = inputFrameVectors[i];
        }
 
    };
    /************** end parallelTransport::operator() *************/
    
    /*************** conjugateMap::operator() ****************/
    void ConjugateMap::operator() ( double * baseElPtr,
                                    double * destElPtr,
                                    double * resultElPtr )
    {
        /* not implemented */
    };
    /************** end conjugateMap::operator() *************/
    
} /* end namespace */