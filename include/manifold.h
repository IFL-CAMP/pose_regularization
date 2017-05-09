/* lightweight template class for processing manifold valued data
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

namespace manifold {
    
    /* function for asserting the element length */
    bool wrongElementLength( int value );
    
    /* define regularization types */
    enum regularization { HUBER=0, L1, L2 };
    
    /* functor for logarithm map */
    class LogarithmMap
    {
    public:
        LogarithmMap ( int elementLength ) :
            _elementLength( elementLength ) {}
            
        double operator() ( double * baseElPtr,
                            double * destElPtr,
                            double * resultElPtr );
    private:
        int _elementLength;
    };
    
    /* functor for exponential map */
    class ExponentialMap
    {
    public:
        ExponentialMap ( int elementLength ) :
            _elementLength( elementLength ) {}
            
        void operator() ( double * baseElPtr,
                          double * dirElPtr,
                          double * resultElPtr );
    private:
        int _elementLength;
    };

	/* functor for inner product */
    class InnerProduct
    {
    public:
        InnerProduct ( int elementLength ) :
            _elementLength( elementLength ) {}
            
        double operator() ( double * elementPtr1,
                            double * elementPtr2,
                            double * baseElPtr );

    private:
        int _elementLength;
    };
    
    /* functor for adjusting coefficients of moving frame */
    class AdjustCoefficients
    {
    public:
        AdjustCoefficients ( int dimension ) :
            _dimension( dimension ) {}
            
        void operator() ( double * inputCoefficients,
                          double * oputputCoefficients,
                          double * eigenvalues,
                          double t );

    private:
        int _dimension;
    };
    
    /* functor for computing the moving frame */
    class ComputeFrame
    {
    public:
        ComputeFrame ( int dimension, int elementLength ) :
            _dimension( dimension ),
            _elementLength( elementLength ) {}
            
        void operator() ( double * tagentSpaceElement,
                          double * basisElements,
                          double * eigenvalues, 
						  double * p);

    private:
        int _dimension;
        int _elementLength;     
    };
    
    /* functor for computing the parallel transport */
    class ParallelTransport
    {
    public:
        ParallelTransport ( int dimension, int elementLength ) :
            _dimension( dimension ),
            _elementLength( elementLength ) {}
            
        void operator() ( double * inputFrameVectors,
                          double * outputFrameVectors,
                          double * baseElPtr,
                          double * destElPtr );

    private:
        int _dimension;
        int _elementLength;
    };
    
    class ConjugateMap
    {
    public:
        ConjugateMap ( int elementLength ) :
        _elementLength( elementLength ) {}
        
        void operator() ( double * baseElPtr,
                           double * destElPtr,
                           double * resultElPtr );
    private:
        int _elementLength;
    };

    
    /* functor for computing first order proximal mappings */
    class ProximalMapFirstOrder
    {
        
    public:
        
        ProximalMapFirstOrder ( int regularization, int elementLength ) :
            _regularization( regularization ),
            _elementLength( elementLength ),
            _log( elementLength ),
            _exp( elementLength ){}
            
        void operator() ( double * startElPtr, 
                          double * destElPtr, 
                          double lambda, 
                          double * resElStartPtr );
        
        void operator() ( double * startElPtr, 
                          double * destElPtr,
                          double alpha,
                          double lambda, 
                          double * resElStartPtr );
        
        void operator() ( double * startElPtr, 
                          double * destElPtr,
                          double alpha, 
                          double lambda,
                          double * resElStartPtr,
                          double * resElDestPtr );
    private:
        
        int _elementLength;
        int _regularization;
        LogarithmMap _log;
        ExponentialMap _exp;
        
    };
    /* end functor for computing first order proximal mappings */

    /* functor for computing second order proximal mappings */
    class ProximalMapSecondOrder
    {
        
    public:
            
        ProximalMapSecondOrder( int regularization, int dimension, int elementLength ) :
            _regularization( regularization ),
            _dimension( dimension ),
            _elementLength( elementLength ),
            _log( elementLength ),
            _exp( elementLength ),
            _prod( elementLength ),
            _adjust( dimension ),
            _frame( dimension, elementLength ),
            _par( dimension, elementLength ){}
        
            void operator() ( double * x1,
                              double * x2,
                              double * x3,
                              double beta,
                              double lambda,
							  double inner_factor,
							  int inner_steps,
                              double * y1,
                              double * y2,
                              double * y3);
    private:
        
        int _elementLength;
        int _dimension;
        int _regularization;
        LogarithmMap _log;
        ExponentialMap _exp;
        InnerProduct _prod;
        AdjustCoefficients _adjust;
        ComputeFrame _frame;
        ParallelTransport _par;
        
    };
    /* end functor for computing second order proximal mappings */

	/* proximal point algorithm for 3D */
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
								 int inner_steps);
    
    
    /* iterative Karcher mean */
    /* This part is not used,
     * but left for future use by Marco!
     */
    void iterativeKarcherMean( double* f,
                               double* x,
                               double* result,
                               int numberOfElements,
                               int steps );
    
} /* end namespace manifold */

