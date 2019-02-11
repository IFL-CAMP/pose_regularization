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

    /* define regularization types */
    enum Regularization
    {
        HUBER, L1, L2
    };

    /* proximal point algorithm for 3D */
    // f, x: column-major flattened out 4x4 matrices
    /*!
     * Regularize a 6-DoF pose stream.
     *
     * //TODO
     *
     * @param f input 6-DoF pose stream, composed of row-major flattened 4x4 transformation matrices
     * @param x output pose stream, in same format as input
     * @param p regularization method for the similarity term
     * @param q regularization method for the 1st order forward difference
     * @param r regularization method for the 2nd order difference
     * @param inner_factor // TODO
     * @param alpha
     * @param beta
     * @param steps
     * @param inner_steps
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