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
#include <limits>
#include <manifold.h>

constexpr double EPSILON = 0.000001;

#ifdef SE3PROD
    constexpr int ELEMENT_LENGTH = 16;
    constexpr int DIMENSION = 6;
#endif

namespace manifold {
    
    /********** function for asserting the element length **********/
    bool wrongElementLength( int value )
    {
        bool flag = true;
        if ( value == ELEMENT_LENGTH ) { flag = false; }
        return flag;
    }
    
    /************* geodesic path length for data term *************/
    double dataPathLength( double lambda,
                           double dist,
                           int p )
    {
        /* initialize tau */
        double tau = 0.0;
        
        /* change tau in case the geodesic distance is positive */
        if ( dist > 0.0 )
        {
            //tau determines how far we shoot !
            switch (p)
            {
                /* L^2-norm */
                case 2:
                    tau = lambda/(1.0+lambda);
                    break;
                    
                    /* L^1-norm => soft thresholding */
                case 1:
                    if(lambda<dist)
                    {
                        tau = lambda/dist;
                    }
                    else
                    {
                        tau = 1.0;
                    }
                    break;
                    
                    /* Huber-norm */
                case 0:
                    
                    /* old implementation of Laurent, should be checked in a future version */
                    double tau_Huber = 1.0;
                    double omega = 1.0;
                    
                    double l_t2_huber = lambda*tau_Huber*tau_Huber;
                    
                    if(dist*sqrt(2.0)*tau_Huber < omega*(1.0+2.0*l_t2_huber))
                    {
                        tau = (2.0*l_t2_huber)/(1.0+2.0*l_t2_huber);
                    }
                    else
                    {
                        if(sqrt(2.0)*lambda*omega*tau_Huber  < dist)
                        {
                            tau = sqrt(2.0)*lambda*omega*tau_Huber/dist;
                        }
                        else
                        {
                            tau = 1.0;
                        }
                    } /* end if(dist*sqrt(2.0)...*/
                    
                    break;
                    
            } /* end switch (p) */
        } /* end if ( dist > 0 ) */
        
        return tau;
    }
    
    /************* geodesic path length for regularizer *************/
    double regPathLength( double alpha,
            double lambda,
            double dist,
            int q )
    {
        /* initialize tau */
        double tau = 0.0;
        
        /* change tau in case the geodesic distance is positive */
        if ( dist > 0.0 )
        {
            switch (q)
            {
                /* L^2-norm */
                case 2:
                    tau = lambda*alpha/(2.0*alpha*lambda+1.0);
                    break;
                    
                    /* L^1-norm => soft thresholding */
                case 1:
                    tau = lambda*alpha/dist;
                    if(tau > 0.5)
                    {
                        tau = 0.5;
                    }
                    break;
                    
                    /* Huber-norm */
                case 0:
                    
                    /* old implementation of Laurent, should be checked in a future version */
                    double tau_Huber = 1.0;
                    double omega = 1.0;
                    
                    if ( dist < omega*(1.0+4.0*lambda*alpha*tau_Huber*tau_Huber)/(sqrt(2.0)*tau_Huber) )
                    {
                        tau = (2.0*alpha*lambda*tau_Huber*tau_Huber)/(1.0+4.0*alpha*lambda*tau_Huber*tau_Huber);
                    }
                    else
                    {
                        tau = sqrt(2.0)*alpha*lambda*omega*tau_Huber/dist;
                        if(tau > 0.5)
                        {
                            tau = 0.5;
                        }
                    } /* end if ( dist < omega*... */
                    break;
                    
            } /* end switch (q) */
        } /* end if ( dist > 0.0 ) */
        
        return tau;
    }
    
    /*************** proximalMapFirstOrder::operator() ***************/
    void ProximalMapFirstOrder::operator() ( double * startElPtr,
                                             double * destElPtr,
                                             double lambda,
                                             double * startElResPtr )
    {
        /* initialize memory for intermediate results */
        double * cache = new double[_elementLength]; 
        
        /* compute logarithm and distance */
        double dist = _log( startElPtr, destElPtr, cache );
        
        /* calculate time for data term geodesic */
        double t = dataPathLength( lambda, dist, _regularization );
        
        /* multiply computed logarithm with t */
        for ( int i=0; i < _elementLength; ++i ) {
            cache[i] *= t;
        }
        
        /* compute proximal mapping for data term */
        _exp( startElPtr, cache, startElResPtr );
        
        /* free memory */
        delete[] cache;
        
    };
    
    void ProximalMapFirstOrder::operator() ( double * startElPtr,
                                             double * destElPtr,
                                             double alpha,
                                             double lambda,
                                             double * startElResPtr )
    {
        /* initialize memory for intermediate results */
        double * cache = new double[_elementLength];
        
        /* compute logarithm and distance */
        double dist = _log( startElPtr, destElPtr, cache );
        
        /* calculate time for data term geodesic */
        double t = regPathLength( alpha, lambda, dist, _regularization );
        
        /* multiply computed logarithm with t */
        for ( int i=0; i < _elementLength; ++i ) {
            cache[i] *= t;
        }
        
        /* compute proximal mapping for data term */
        _exp( startElPtr, cache, startElResPtr );
        
        
        /* free memory */
        delete[] cache;
        
    };
    
    void ProximalMapFirstOrder::operator() ( double * startElPtr,
                                             double * destElPtr,
                                             double alpha,
                                             double lambda,
                                             double * startElResPtr,
                                             double * destElResPtr )
    {
        /* initialize memory for intermediate results */
        double * cache1 = new double[_elementLength];
        double * cache2 = new double[_elementLength];
        
        /* compute logarithm and distance */
        double dist = _log( startElPtr, destElPtr, cache1 );
        _log( destElPtr, startElPtr, cache2 );
        
        /* calculate times for data term geodesics */
        double t = regPathLength( alpha, lambda, dist, _regularization );
        
        /* multiply computed logarithm with t */
        for ( int i=0; i < _elementLength; ++i ) {
            cache1[i] *= t;
            cache2[i] *= t;
        }
        
        /* compute proximal mapping for data term */
        _exp( startElPtr, cache1, startElResPtr );
        _exp( destElPtr, cache2, destElResPtr );
        
        /* free memory */
        delete[] cache1;
        delete[] cache2;
    };
    /************* end proximalMapFirstOrder::operator() *************/
        
    /*************** proximalMapSecondOrder::operator() ***************/
    void ProximalMapSecondOrder::operator() ( double * x1,
                                              double * x2,
                                              double * x3,
                                              double beta,
                                              double lambda,
											  double inner_factor,
											  int inner_steps,
                                              double * y1,
                                              double * y2,
                                              double * y3 )
    {
        /* allocate memory */
        double * aux1 = new double[_elementLength];
        double * aux2 = new double[_elementLength];
        double * aux3 = new double[_elementLength];
        double * ym = new double[_elementLength];
        double * cache = new double[_elementLength];
        double * log_ym_y2 = new double[_elementLength];
        double * log_y2_ym = new double[_elementLength];
        double * log_y1_y3 = new double[_elementLength];
        double * log_y3_y1 = new double[_elementLength];
        double * log_y1_x1 = new double[_elementLength];
        double * log_y2_x2 = new double[_elementLength];
        double * log_y3_x3 = new double[_elementLength];
        double * grad_1 = new double[_elementLength];
        double * grad_2 = new double[_elementLength];
        double * grad_3 = new double[_elementLength];
        double * basis_x_coeffs = new double[_elementLength];
        double * basisElements1 = new double[_dimension*_elementLength];
        double * basisElements2 = new double[_dimension*_elementLength];
        double * eigenvalues = new double[_dimension];
        double * coeffs1 = new double[_dimension];
        double * coeffs2 = new double[_dimension];
        double length2ndDiff = 0.0;
        double factor = 0.0;
        double tau = 0.0;
        double t = 0.0;
        
        /* initialization */
        for ( int i = 0; i < _elementLength; ++i )
        {
            y1[i] = x1[i];
            y2[i] = x2[i];
            y3[i] = x3[i];
        }
        
        /* adjust beta according to damping parameters */
        beta = beta*lambda;

        /* main iteration (gradient descent) */
        for ( int i = 0; i < inner_steps; ++i )
        {
               
            /* compute mean of y1 and y3, i.e., [y1, y3]_(d/2) */
            _log( y1, y3, cache );
            for ( int j = 0; j < _elementLength; ++j )
            {
                cache[j] *= 0.5;
            }
            _exp( y1, cache, ym );
            
            /* compute direction at xm and length of 2nd difference */
            length2ndDiff = _log( ym, y2, log_ym_y2 );
            
            /* compute tau */
            //tau = (inner_factor/beta)/std::pow((double)(i)+1.0,0.95 + 0.5*std::pow(i+1,-0.18));
            tau = inner_factor/std::pow((double)(i)+1.0,0.95 + 0.5*std::pow(i+1,-0.18));
            
            /******** gradient for y1 ********/
            
            /* compute direction of geodesic */
            _log( y3, y1, log_y3_y1 );
            
            /* compute basis of eigenvectors at y3 */
            _frame( log_y3_y1, basisElements1, eigenvalues, NULL );
            
            /* move basis to xm via parallel transport */
            _par( basisElements1, basisElements2, y3, ym );
            
            /* compute coefficients at ym */
            for ( int k = 0; k < _dimension; ++k )
            {
                int index = k*_elementLength;
                coeffs1[k] = _prod( &basisElements2[index], log_ym_y2, ym );
            }
            
            /* move basis towards y1 */
            _par( basisElements2, basisElements1, ym, y1 );
            
            /* adjust coefficients at y1 */
            t = _prod( log_y3_y1, log_y3_y1, y1 );
            _adjust( coeffs1, coeffs2, eigenvalues, t );
            
            /* compute basis x coefficients */
            for ( int j = 0; j < _elementLength; ++j )
            {
                basis_x_coeffs[j] = 0.0;
            }
            
            for ( int k = 0; k < _dimension; ++k )
            {
                int index = k*_elementLength;
                for ( int j = 0; j < _elementLength; ++j )
                {
                    basis_x_coeffs[j] += coeffs2[k]*basisElements1[index+j];
                }
            }
            
            /* check for regularization */
            if ( _regularization == 2 )
            {
                factor = beta;
            }
            if ( _regularization == 1 )
            {
                //factor = std::sqrt(_prod( basis_x_coeffs, basis_x_coeffs, y1 ));
                
				if (length2ndDiff > EPSILON)
                {
                    factor = beta/length2ndDiff;
                }
                else
                {
                    factor = 0.0;
                }
            }
            
            /* adjust length of basis_x_coeffs and compose gradient at y1 */
            _log( y1, x1, log_y1_x1 );
            for ( int j = 0; j < _elementLength; ++j )
            {
                grad_1[j] = tau*( basis_x_coeffs[j]*factor - log_y1_x1[j] );
            }
            
            /******** gradient for y2 ********/
            
            /* compute log( y2, ym ) */
            _log( y2, ym, log_y2_ym );
                        
            /* check for regularization */
            if ( _regularization == 2 )
            {
                factor = beta;
            }
            if ( _regularization == 1 )
            {
                //factor = std::sqrt(_prod( log_y2_ym, log_y2_ym, y2 ));
				if (length2ndDiff > EPSILON)
                {
                    factor = beta/length2ndDiff;
                }
                else
                {
                    factor = 0.0;
                }
            }
            
            /* adjust length of basis_x_coeffs and compose gradient at y2 */
            _log( y2, x2, log_y2_x2 );
            for ( int j = 0; j < _elementLength; ++j )
            {
                grad_2[j] = tau*( log_y2_ym[j]*factor - log_y2_x2[j] );
            }
            
            /******** gradient for y3 ********/
            
            /* compute direction of geodesic */
            _log( y1, y3, log_y1_y3 );
            
            /* compute basis of eigenvectors at y1 */
            _frame( log_y1_y3, basisElements1, eigenvalues, NULL );
            
            /* move basis to xm via parallel transport */
            _par( basisElements1, basisElements2, y1, ym );
            
            /* compute coefficients at ym */
            for ( int k = 0; k < _dimension; ++k )
            {
                int index = k*_elementLength;
                coeffs1[k] = _prod( &basisElements2[index], log_ym_y2, ym );
            }
            
            /* move basis towards y3 */
            _par( basisElements2, basisElements1, ym, y3 );
            
            /* adjust coefficients at y3 */
            t = _prod( log_y1_y3, log_y1_y3, y3 );
            _adjust( coeffs1, coeffs2, eigenvalues, t );
            
            /* compute basis x coefficients */
            for ( int j = 0; j < _elementLength; ++j )
            {
                basis_x_coeffs[j] = 0.0;
            }
            
            for ( int k = 0; k < _dimension; ++k )
            {
                int index = k*_elementLength;
                for ( int j = 0; j < _elementLength; ++j )
                {
                    basis_x_coeffs[j] += coeffs2[k]*basisElements1[index+j];
                }
            }
            
            /* check for regularization */
            if ( _regularization == 2 )
            {
                factor = beta;
            }
            if ( _regularization == 1 )
            {
                //factor = std::sqrt(_prod( basis_x_coeffs, basis_x_coeffs, y3 ));
                
				if (length2ndDiff > EPSILON)
                {
                    factor = beta/length2ndDiff;
                }
                else
                {
                    factor = 0.0;
                }
            }
            
            /* adjust length of basis_x_coeffs and compose gradient at y1 */
            _log( y3, x3, log_y3_x3 );
            for ( int j = 0; j < _elementLength; ++j )
            {
                grad_3[j] = tau*( basis_x_coeffs[j]*factor - log_y3_x3[j] );
            }
                    
            /* update step */
            for ( int j = 0; j < _elementLength; ++j )
            {
                aux1[j] = y1[j];
                aux2[j] = y2[j];
                aux3[j] = y3[j];
            }
            _exp( aux1, grad_1, y1 );
            _exp( aux2, grad_2, y2 );
            _exp( aux3, grad_3, y3 );

        } /* end main iteration */
        
        /* free memory */
        delete[] aux1;
        delete[] aux2;
        delete[] aux3;
        delete[] cache;
        delete[] ym;
        delete[] log_ym_y2;
        delete[] log_y2_ym;
        delete[] log_y3_y1;
        delete[] log_y1_y3;
        delete[] log_y1_x1;
        delete[] log_y2_x2;
        delete[] log_y3_x3;
        delete[] grad_1;
        delete[] grad_2;
        delete[] grad_3;
        delete[] basis_x_coeffs;
        delete[] basisElements1;
        delete[] basisElements2;
        delete[] eigenvalues;
        delete[] coeffs1;
        delete[] coeffs2;
       
    };
    /************* end proximalMapSecondOrder::operator() *************/

    /****************** proximal point algorithm 3D *******************/
	void proximalPointAlgorithm( double* f,
                                 double* x,
                                 int m,
                                 int n,
                                 int s,
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
        /* initialize functors for proximal mappings */
        ProximalMapFirstOrder proxData( p, ELEMENT_LENGTH );
        ProximalMapFirstOrder proxReg1stOrder( q, ELEMENT_LENGTH );
        ProximalMapSecondOrder proxReg2ndOrder( r, DIMENSION, ELEMENT_LENGTH );
		
        /* compute total length of array */
        int arrayLength = m*n*s*ELEMENT_LENGTH;
        
        /* initialize x */
        for ( int i = 0; i < arrayLength; ++i ) { x[i] = f[i]; }
        
        /* initialize caches */
        double * cache1 = new double[ELEMENT_LENGTH];
        double * cache2 = new double[ELEMENT_LENGTH];
        double * cache3 = new double[ELEMENT_LENGTH];
		double * cache4 = new double[ELEMENT_LENGTH];

        /* main iteration */
        for ( int ll = 0; ll < steps; ++ll ) {
            
            /* compute lambda */
            double lambda = 1.0/pow((double)(ll+1),0.95 + 0.5*pow((double)(ll+1),-0.18));
            
            /* proximal mapping of data term */
            for ( int k = 0; k < s; ++k ) {
                for ( int j = 0; j < n; ++j ) {
                    for ( int i = 0; i < m; ++i ) {
                        int index = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                        proxData( &x[index], &f[index], lambda, cache1 );
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                            x[index+l] = cache1[l];
                        }                        
                    }
                }
            } /* end proximal mapping of data term */
            
            /* 1st order proximal mappings */
            if ( alpha > 0)
            {
                /* 1st order proximal mapping w.r.t. x-direction */
				
                for ( int k = 0; k < s; ++k ) {
                    for ( int j = 0; j < n-1; ++j ) {
                        for ( int i = 0; i < m; ++i ) {
                            int index1 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index2 = (k*m*n + (j+1)*m + i)*ELEMENT_LENGTH;
                            proxReg1stOrder( &x[index1], &x[index2], alpha, lambda, cache1, cache2 );
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                                x[index1+l] = cache1[l];
                                x[index2+l] = cache2[l];
                            }
                        }
                    }
                } /* end 1st order proximal mapping w.r.t. x-direction */

                /* 1st order proximal mapping w.r.t. y-direction */
				
                for ( int k = 0; k < s; ++k ) {
                    for ( int j = 0; j < n; ++j ) {
                        for ( int i = 0; i < m-1; ++i ) {
                            int index1 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index2 = (k*m*n + j*m + (i+1))*ELEMENT_LENGTH;
                            proxReg1stOrder( &x[index1], &x[index2], alpha, lambda, cache1, cache2 );
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                                x[index1+l] = cache1[l];
                                x[index2+l] = cache2[l];
                            }
                        }
                    }
                } /* end 1st order proximal mapping w.r.t. y-direction */
				
                /* 1st order proximal mapping w.r.t. z-direction */
                for ( int k = 0; k < s-1; ++k ) {
                    for ( int j = 0; j < n; ++j ) {
                        for ( int i = 0; i < m; ++i ) {
                            int index1 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index2 = ((k+1)*m*n + j*m + i)*ELEMENT_LENGTH;
                            proxReg1stOrder( &x[index1], &x[index2], alpha, lambda, cache1, cache2 );
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                                x[index1+l] = cache1[l];
                                x[index2+l] = cache2[l];
                            }
                        }
                    }
                } /* end 1st order proximal mapping w.r.t. z-direction */
            }
            /* end 1st order proximal mappings */
            
            /* 2nd order proximal mappings */
            if ( beta > 0)
            {
                
                /* 2nd order proximal mapping w.r.t. y-direction */
                for ( int k = 0; k < s; ++k ) {
                    for ( int j = 0; j < n; ++j ) {
                        for ( int i = 1; i < m-1; ++i ) {
                            
                            int index1 = (k*m*n + j*m + (i-1))*ELEMENT_LENGTH;
                            int index2 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index3 = (k*m*n + j*m + (i+1))*ELEMENT_LENGTH;

                            proxReg2ndOrder( &x[index1], &x[index2], &x[index3], beta, lambda, inner_factor, inner_steps, cache1, cache2, cache3 );
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                                x[index1+l] = cache1[l];
                                x[index2+l] = cache2[l];
                                x[index3+l] = cache3[l];
                            }                      
                        }
                    }
                } /* end 2nd order proximal mapping w.r.t. y-direction */

                /* 2nd order proximal mapping w.r.t. x-direction */
                for (int k = 0; k < s; ++k) {
                    for (int j = 1; j < n-1; ++j) {
                        for (int i = 1; i < m-1; ++i) {
                            int index1 = (k*m*n + (j - 1)*m + i)*ELEMENT_LENGTH;
                            int index2 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index3 = (k*m*n + (j + 1)*m + i)*ELEMENT_LENGTH;
                            proxReg2ndOrder(&x[index1], &x[index2], &x[index3], beta, lambda, inner_factor, inner_steps, cache1, cache2, cache3);
                            for (int l = 0; l < ELEMENT_LENGTH; ++l) {
                                x[index1 + l] = cache1[l];
                                x[index2 + l] = cache2[l];
                                x[index3 + l] = cache3[l];
                            }
                        }
                    }
                } /* end 2nd order proximal mapping w.r.t. x-direction */
                
                /* 2nd order proximal mapping w.r.t. z-direction */
                for ( int k = 1; k < s-1; ++k ) {
                    for ( int j = 0; j < n; ++j ) {
                        for ( int i = 0; i < m; ++i ) {
                            int index1 = ((k-1)*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index2 = (k*m*n + j*m + i)*ELEMENT_LENGTH;
                            int index3 = ((k+1)*m*n + j*m + i)*ELEMENT_LENGTH;
							proxReg2ndOrder(&x[index1], &x[index2], &x[index3], beta, lambda, inner_factor, inner_steps, cache1, cache2, cache3);
                            for ( int l = 0; l < ELEMENT_LENGTH; ++l ) {
                                x[index1+l] = cache1[l];
                                x[index2+l] = cache2[l];
                                x[index3+l] = cache3[l];
                            }
                        }
                    }
                } /* end 2nd order proximal mapping w.r.t. z-direction */
            }
            /* end 2nd order proximal mappings */

            /* REMARK: Usually this algorithm would require the implementation of mixed
             *         proximal mappings (diagonal directions). However, these have been
             *         removed here, because pose denoising is currently only done in 1D!
             */            
            
        } /* end main iteration */
        
        /* free memory */
        delete[] cache1;
        delete[] cache2;
        delete[] cache3;
		delete[] cache4;
		
    };
    /************** end proximal point algorithm 3D *************/
    
 
    /************** iterative Karcher mean **************/
    /* This part is not used either,
     * but left for future use by Marco!
     */
    void iterativeKarcherMean( double* f,
                               double* x,
                               double* result,
                               int numberOfElements,
                               int steps )
    {
        
        /* initialize functors for geodesic distances */
        LogarithmMap logMap( ELEMENT_LENGTH );
        ExponentialMap expMap( ELEMENT_LENGTH );
        
        /* initialize cache memory */
        double * cache1 = new double[ELEMENT_LENGTH];
        double * cache2 = new double[ELEMENT_LENGTH];
        
        /* initialize constant for efficient normalization */
        double c = 1.0/( (double)numberOfElements );
		
		/* store result */
        for ( int k = 0; k < ELEMENT_LENGTH; ++k ) { result[k] = x[k]; }
        
        /* main iteration */
        for ( int i = 0; i < steps; ++i )
        {
            /* initialize cache for update */
            for ( int k = 0; k < ELEMENT_LENGTH; ++k ) { cache1[k] = 0.0; }
            
            /* sum up all tangential elements */
            for ( int j = 0; j < numberOfElements; ++j )
            {
                /* compute logarithm */
                int index = j*ELEMENT_LENGTH;
                logMap( &result[0], &f[index], &cache2[0] );
                
                /* add computed element */
                for ( int k = 0; k < ELEMENT_LENGTH; ++k) { cache1[k] += c*cache2[k]; }
                
            } /* end of summ of all tangential contributions */
            
            /* compute exponential map */
            expMap( &result[0], &cache1[0], &cache2[0] );
            
            /* update result */
            for ( int k = 0; k < ELEMENT_LENGTH; ++k ) { result[k] = cache2[k]; }
            
        } /* end main iteration */

        /* free memory */
        delete[] cache1;
        delete[] cache2;
    }
    /************ end iterative Karcher mean ************/
    
} /* end namespace */



