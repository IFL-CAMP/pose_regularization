/* lightweight template class for processing pose_regularization valued data
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

#ifndef POSE_REGULARIZATION_H
#define POSE_REGULARIZATION_H

#include <array>
#include <vector>

namespace pose_regularization {
    constexpr size_t SE3_DIMENSIONS = 6;
    constexpr size_t SE3_ELEMENT_LENGTH = 16;
    using SE3_ELEMENT = std::array<double, SE3_ELEMENT_LENGTH>;
    // TODO: add row-major conversion operator from vec<vec>?
    // TODO: add span<double> when available...

    /// Loss functions available for regularization
    enum Regularization
    {
        HUBER, L1, L2
    };

    /*!
     * Regularize a 6-DoF pose stream, expressed as a vector of transformation matrices
     *
     * The input sequence of 4x4 Euclidean transformation matrices is regularized
     * according to the provided parameters, and the result is returned.
     *
     * The action of each parameter is discussed below; for best results,
     * they should be chosen carefully (for example, with a grid search
     * study). However, a good starting point for the choice of the parameters
     * in case of strong jitter noise is the following:
     * 
     * proximalPointAlgorithm(poses_noisy, poses_regularized, L2, HUBER, L2, 50,
     * 1, 0.25, 100, 50)
     *
     * If faster execution is desired, 2nd order regularization can be
     * disabled (it is worthwhile to study its effect anyway):
     * proximalPointAlgorithm(poses_noisy, poses_regularized, L2, HUBER, HUBER, 50,
     * 0, 0, 100, 0)
     *
     * poses_noisy is expected to be a stack of 4x4 matrices, representing
     * poses in the Euclidean space with respect to a single coordinate
     * system. If N is the number of poses, this Numpy array should have shape
     * Nx4x4.
     *
     * p, q, and r define the regularization method for the data term, and
     * for the 1st order and 2nd order output regularization terms
     * respectively. 0 defines the HUBER norm, 1 the L1 regularization
     * (linear) and 2 the L2 regularization method (quadratic).
     *
     * alpha and beta are dampening coefficients for the 1st and 2nd order
     * regularization, respectively. Either of the two regularizations can be
     * disabled by setting the relative coefficient to zero.
     *
     * inner_factor is a further parameter for the action of the 2nd order
     * regularizer.
     *
     * steps is the number of global regularization steps performed.
     *
     * inner_steps is the number of steps performed within the 2nd order
     * regularization process.
     *
     * @param f input 6-DoF pose stream, composed of row-major flattened 4x4 transformation matrices
     * @param x output pose stream, in same format as input
     * @param p regularization method for the similarity term
     * @param q regularization method for the 1st order forward difference
     * @param r regularization method for the 2nd order difference
     * @param inner_factor further parameter for the 2nd order regularization
     * @param alpha dampening coefficient for the 1st order regularization
     * @param beta dampening coefficient for the 2nd order regularization
     * @param steps number of iterations 
     * @param inner_steps number of iterations of the 2nd order regularization
     */
    void proximalPointAlgorithm(const std::vector<SE3_ELEMENT> &f,
                                std::vector<SE3_ELEMENT> &x,
                                Regularization p,
                                Regularization q,
                                Regularization r,
                                double inner_factor,
                                double alpha,
                                double beta,
                                int steps,
                                int inner_steps);

} /* end namespace pose_regularization */

#endif
