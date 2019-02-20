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

#include <pose_regularization.h>
#include <iostream>
#include <cmath>

constexpr double EPSILON = 0.000001;
constexpr double SE3_EPSILON = 0.00000001;

namespace pose_regularization {

    /*!
     * Hand-optimized product of 4x4 rotation matrices.
     *
     * Given matrices A and B represented as flat arrays of double of length 16,
     * returns the matrix product of A and B.
     *
     * Only the rotational component of the matrix will be considered (the 3x3 top-left submatrix).
     *
     * @param A First input matrix
     * @param B Second input matrix
     * @param C Output matrix: product of A and B
     */
    void SO3Product_AB(const SE3_ELEMENT &A, const SE3_ELEMENT &B, SE3_ELEMENT &C)
    {
        C[0] = A[0] * B[0] + A[4] * B[1] + A[8] * B[2];
        C[1] = A[1] * B[0] + A[5] * B[1] + A[9] * B[2];
        C[2] = A[2] * B[0] + A[6] * B[1] + A[10] * B[2];
        C[3] = 0.0;

        C[4] = A[0] * B[4] + A[4] * B[5] + A[8] * B[6];
        C[5] = A[1] * B[4] + A[5] * B[5] + A[9] * B[6];
        C[6] = A[2] * B[4] + A[6] * B[5] + A[10] * B[6];
        C[7] = 0.0;

        C[8] = A[0] * B[8] + A[4] * B[9] + A[8] * B[10];
        C[9] = A[1] * B[8] + A[5] * B[9] + A[9] * B[10];
        C[10] = A[2] * B[8] + A[6] * B[9] + A[10] * B[10];
        C[11] = 0.0;

        C[12] = 0.0;
        C[13] = 0.0;
        C[14] = 0.0;
        C[15] = 1.0;
    };

    /*!
     * Hand-optimized product of 4x4 rotation matrices.
     *
     * Given matrices A and B represented as flat arrays of double of length 16,
     * returns the matrix product of the inverse of A and B.
     *
     * Only the rotational component of the matrix will be considered (the 3x3 top-left submatrix).
     *
     * @param A First input matrix
     * @param B Second input matrix
     * @param C Output matrix: product of inv(A) and B
     */
    void SO3InvProduct_invAB(const SE3_ELEMENT &A, const SE3_ELEMENT &B, SE3_ELEMENT &C)
    {
        C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
        C[1] = A[4] * B[0] + A[5] * B[1] + A[6] * B[2];
        C[2] = A[8] * B[0] + A[9] * B[1] + A[10] * B[2];
        C[3] = 0.0;

        C[4] = A[0] * B[4] + A[1] * B[5] + A[2] * B[6];
        C[5] = A[4] * B[4] + A[5] * B[5] + A[6] * B[6];
        C[6] = A[8] * B[4] + A[9] * B[5] + A[10] * B[6];
        C[7] = 0.0;

        C[8] = A[0] * B[8] + A[1] * B[9] + A[2] * B[10];
        C[9] = A[4] * B[8] + A[5] * B[9] + A[6] * B[10];
        C[10] = A[8] * B[8] + A[9] * B[9] + A[10] * B[10];
        C[11] = 0.0;

        C[12] = 0.0;
        C[13] = 0.0;
        C[14] = 0.0;
        C[15] = 1.0;
    };

    /*!
     * Hand-optimized exponential of 4x4 rotation matrices.
     *
     * Given a matrix A represented as a flat array of double of length 16,
     * returns the matrix exponential of A.
     *
     * Only the rotational component of the matrix will be considered (the 3x3 top-left submatrix).
     *
     * @param Ai Input matrix
     * @param Ao Output matrix: exp(A)
     */
    void SO3Exponential(const SE3_ELEMENT &Ai, SE3_ELEMENT &Ao)
    {
        double omega_x = Ai[6];
        double omega_y = Ai[8];
        double omega_z = Ai[1];

        /* calculate norm of omega */
        double no = std::sqrt(std::pow(omega_x, 2) + std::pow(omega_y, 2) + std::pow(omega_z, 2));

        if (std::abs(no) > SE3_EPSILON) {
            /* calculate factors */
            double f1 = std::sin(no) / no;
            double f2 = (1.0 - std::cos(no)) / std::pow(no, 2);

            /* calculate matrix exponential of destination element (tangential element, called v in our conversions) */
            Ao[0] = 1.0 - f2 * (std::pow(omega_y, 2) + omega_z * omega_z);
            Ao[1] = f1 * omega_z + f2 * omega_x * omega_y;
            Ao[2] = f2 * omega_x * omega_z - f1 * omega_y;
            Ao[3] = 0.0;

            Ao[4] = f2 * omega_x * omega_y - f1 * omega_z;
            Ao[5] = 1.0 - f2 * (omega_x * omega_x + std::pow(omega_z, 2));
            Ao[6] = f1 * omega_x + f2 * omega_y * omega_z;
            Ao[7] = 0.0;

            Ao[8] = f1 * omega_y + f2 * omega_x * omega_z;
            Ao[9] = f2 * omega_y * omega_z - f1 * omega_x;
            Ao[10] = 1.0 - f2 * (std::pow(omega_x, 2) + std::pow(omega_y, 2));
            Ao[11] = 0.0;

            Ao[12] = 0.0;
            Ao[13] = 0.0;
            Ao[14] = 0.0;
            Ao[15] = 1.0;
        } else {
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

    /*!
     * Schur decomposition of a 4x4 rotation matrix
     *
     * Only the rotational component of the matrix will be considered (the 3x3 top-left submatrix).
     *
     * @param Ai Input matrix
     * @param U // TODO
     * @return
     */
    double SO3SchurDecomposition(const SE3_ELEMENT &Ai, std::array<double, 9> &U)
    {
        /* get a, b and c of matrix
         * A = [0 -a b; a 0 -c ; -b c 0]
         */
        double a = Ai[1];
        double b = Ai[8];
        double c = Ai[6];

        double beta = std::sqrt(std::pow(a, 2) + std::pow(b, 2) + std::pow(c, 2));
        double gamma = std::sqrt(std::pow(a, 2) + std::pow(b, 2));

        /* compose U */
        if (beta > SE3_EPSILON && gamma > SE3_EPSILON) {
            U[0] = c / beta;
            U[1] = b / beta;
            U[2] = a / beta;

            U[3] = 0.0;
            U[4] = a / gamma;
            U[5] = -b / gamma;

            U[6] = -gamma / beta;
            U[7] = (b * c) / (beta * gamma);
            U[8] = (a * c) / (beta * gamma);
        } else {
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

    /*!
     * Compute the logarithmic mapping for the vector (x1,x2) at x1
     *
     * @param x1 Starting point and tangential element
     * @param x2 Tail of the direction vector at x1
     * @param y1 Result of the projection of (x1,x2) on the tangential space
     * @return Length of the logarithm
     */
    double SE3LogarithmMap(const SE3_ELEMENT &x1, const SE3_ELEMENT &x2, SE3_ELEMENT &y1)
    {
        SE3_ELEMENT C;

        /* compute inverted base element * dest element (only for SO(3)-component) */
        SO3InvProduct_invAB(x1, x2, C);

        /* calculate phi */
         double phi =  acos( 0.5*(C[0] + C[5] + C[10] - 1.0) );
//        double phi = (acos((std::complex<double>) (0.5 * (C[0] + C[5] + C[10] - 1.0)))).real();

        /* calculate result for SO(3) part */
        if ((std::abs(phi) > SE3_EPSILON) && (std::abs(phi) < M_PI - SE3_EPSILON)) {
            /* get omegas */
            double omega_x = 0.5 * phi * (C[6] - C[9]) / std::sin(phi);
            double omega_y = 0.5 * phi * (C[8] - C[2]) / std::sin(phi);
            double omega_z = 0.5 * phi * (C[1] - C[4]) / std::sin(phi);

            /* first column of result */
            y1[0] = 0.0;
            y1[1] = omega_z;
            y1[2] = -omega_y;
            y1[3] = 0.0;

            /* second column of result */
            y1[4] = -omega_z;
            y1[5] = 0.0;
            y1[6] = omega_x;
            y1[7] = 0.0;

            /* third column of result */
            y1[8] = omega_y;
            y1[9] = -omega_x;
            y1[10] = 0.0;
            y1[11] = 0.0;

            /* fourth column of result is calculated later */
        } else {
            for (auto &resEl: y1) {
                resEl = 0.0;
            }
        } // end if ( abs(phi) > 0 )


        y1[12] = x2[12] - x1[12];
        y1[13] = x2[13] - x1[13];
        y1[14] = x2[14] - x1[14];
        y1[15] = 0.0;

        /* compute length of logarithm without considering last element
         * (note that translational component does not have factor 1/2)
         */
        double dist = 0.0;
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH - 4; ++i) {
            dist += 0.5 * y1[i] * y1[i];
        }
        for (size_t i = SE3_ELEMENT_LENGTH - 4; i < SE3_ELEMENT_LENGTH - 1; ++i) {
            dist += y1[i] * y1[i];
        }

        return std::sqrt(dist);
    };

    /*!
     * Compute the exponential mapping for the vector (x1,x2) at x1
     *
     * @param x1 Starting point and tangential element
     * @param x2 Tail of the direction vector at x1 in tangent space
     * @param y1 Result of the projection of (x1,x2) on the SE3 manifold
     */
    void SE3ExponentialMap(const SE3_ELEMENT &x1, const SE3_ELEMENT &x2, SE3_ELEMENT &y1)
    {

        /* compute matrix exponential of destination element */
        SE3_ELEMENT expTangEl;
        SO3Exponential(x2, expTangEl);

        /* multiply base element with matrix exponential of tangential element */
        SO3Product_AB(x1, expTangEl, y1);

        /* change fourth column of result */
        y1[12] = x1[12] + x2[12];
        y1[13] = x1[13] + x2[13];
        y1[14] = x1[14] + x2[14];
        y1[15] = 1.0;
    };

    /*!
     * Compute the inner product of two 4x4 transformation matrices
     *
     * @param A First input matrix
     * @param B Second input matrix
     * @return Inner product of the two inputs (scalar)
     */
    double SE3InnerProduct(const SE3_ELEMENT &A, const SE3_ELEMENT &B)
    {

        double prod = 0.0;
        // rotational components
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH - 4; ++i) {
            prod += 0.5 * A[i] * B[i];
        }
        // translational components
        for (size_t i = SE3_ELEMENT_LENGTH - 4; i < SE3_ELEMENT_LENGTH - 1; ++i) {
            prod += A[i] * B[i];
        }
        return prod;
    };

    /*!
     * // TODO
     * @param inputCoefficients
     * @param outputCoefficients
     * @param eigenvalues
     * @param t
     */
    void SE3AdjustCoefficientsMovingFrame(const std::array<double, SE3_DIMENSIONS> &inputCoefficients,
                                          std::array<double, SE3_DIMENSIONS> &outputCoefficients,
                                          const std::array<double, SE3_DIMENSIONS> &eigenvalues,
                                          const double &t)
    {

        /* generic initialization */
        for (size_t i = 0; i < SE3_DIMENSIONS; i++) {
            outputCoefficients[i] = 0.5 * inputCoefficients[i];
        }

        /* modify if necessary */
        if (t > SE3_EPSILON && eigenvalues[1] > SE3_EPSILON) {
            outputCoefficients[1] = inputCoefficients[1] * std::sin(std::sqrt(eigenvalues[1]) * t / 2.0) /
                                    (std::sin(std::sqrt(eigenvalues[1]) * t) + 0.000000000001);
        }
        if (t > SE3_EPSILON && eigenvalues[2] > SE3_EPSILON) {
            outputCoefficients[2] = inputCoefficients[2] * std::sin(std::sqrt(eigenvalues[2]) * t / 2.0) /
                                    (std::sin(std::sqrt(eigenvalues[2]) * t) + 0.000000000001);
        }
    };

    /*!
     * Compute the basis and eigenvalues of Jacobi field at tangential element
     * 
     * @param tangentSpaceElement Input tangential element
     * @param basisElements Output Jacobi field basis vectors
     * @param eigenvalues Output corresponding eigenvalues
     */
    void SE3ComputeFrame(const SE3_ELEMENT &tangentSpaceElement,
                         std::array<SE3_ELEMENT, SE3_DIMENSIONS> &basisElements,
                         std::array<double, SE3_DIMENSIONS> &eigenvalues)
    {
        /* perform Schur decompositon */
        std::array<double, 9> Um{};
        double beta = SO3SchurDecomposition(tangentSpaceElement, Um);

        /* compute eigenvalues */
        eigenvalues[0] = 0.0;
        eigenvalues[1] = 0.25 * beta * beta;
        eigenvalues[2] = 0.25 * beta * beta;
        eigenvalues[3] = 0.0;
        eigenvalues[4] = 0.0;
        eigenvalues[5] = 0.0;

        /* compute basis elements */
        basisElements[0][0] = 0.0;
        basisElements[0][1] = Um[3] * Um[7] - Um[6] * Um[4];
        basisElements[0][2] = Um[3] * Um[8] - Um[6] * Um[5];
        basisElements[0][3] = 0.0;

        basisElements[0][4] = Um[6] * Um[4] - Um[3] * Um[7];
        basisElements[0][5] = 0.0;
        basisElements[0][6] = Um[4] * Um[8] - Um[7] * Um[5];
        basisElements[0][7] = 0.0;

        basisElements[0][8] = Um[6] * Um[5] - Um[3] * Um[8];
        basisElements[0][9] = Um[7] * Um[5] - Um[4] * Um[8];
        basisElements[0][10] = 0.0;
        basisElements[0][11] = 0.0;

        basisElements[0][12] = 0.0;
        basisElements[0][13] = 0.0;
        basisElements[0][14] = 0.0;
        basisElements[0][15] = 0.0;


        basisElements[1][0] = 0.0;
        basisElements[1][1] = Um[0] * Um[4] - Um[3] * Um[1];
        basisElements[1][2] = Um[0] * Um[5] - Um[3] * Um[2];
        basisElements[1][3] = 0.0;

        basisElements[1][4] = Um[3] * Um[1] - Um[0] * Um[4];
        basisElements[1][5] = 0.0;
        basisElements[1][6] = Um[1] * Um[5] - Um[4] * Um[2];
        basisElements[1][7] = 0.0;

        basisElements[1][8] = Um[3] * Um[2] - Um[0] * Um[5];
        basisElements[1][9] = Um[4] * Um[2] - Um[1] * Um[5];
        basisElements[1][10] = 0.0;
        basisElements[1][11] = 0.0;

        basisElements[1][12] = 0.0;
        basisElements[1][13] = 0.0;
        basisElements[1][14] = 0.0;
        basisElements[1][15] = 0.0;


        basisElements[2][0] = 0.0;
        basisElements[2][1] = Um[0] * Um[7] - Um[6] * Um[1];
        basisElements[2][2] = Um[0] * Um[8] - Um[6] * Um[2];
        basisElements[2][3] = 0.0;

        basisElements[2][4] = Um[6] * Um[1] - Um[0] * Um[7];
        basisElements[2][5] = 0.0;
        basisElements[2][6] = Um[1] * Um[8] - Um[7] * Um[2];
        basisElements[2][7] = 0.0;

        basisElements[2][8] = Um[6] * Um[2] - Um[0] * Um[8];
        basisElements[2][9] = Um[7] * Um[2] - Um[1] * Um[8];
        basisElements[2][10] = 0.0;
        basisElements[2][11] = 0.0;

        basisElements[2][12] = 0.0;
        basisElements[2][13] = 0.0;
        basisElements[2][14] = 0.0;
        basisElements[2][15] = 0.0;

        /* eigenvectors corresponding to translation */
        for (size_t basisElementInd = 3; basisElementInd < SE3_DIMENSIONS; basisElementInd++) {
            for (size_t ind = 0; ind < SE3_ELEMENT_LENGTH; ++ind) {
                basisElements[basisElementInd][ind] = 0.0;
            }
        }

        basisElements[3][12] = 1.0; // x
        basisElements[4][13] = 1.0; // y
        basisElements[5][14] = 1.0; // z
    };

    /*!
     * Compute the Parallel Transport along the direction specified by two points
     *
     * The parallel transport takes the inputFrameVectors
     * (and array of size SE3_ELEMENT_LENGTH*SE3_DIMENSIONS ),
     * moves them into the direction log( baseEl, destEl ),
     * and stores them in inputFrameVectors
     * (of size SE3_ELEMENT_LENGTH*SE3_DIMENSIONS ).
     * The implementation goes along with the paper
     *
     * Q. Reentmesters
     * A gradient method for geodesic data fitting on some
     * symmetric Riemannian manifolds
     * IEEE Conference on Decision and Control
     * and European Control Conference (2011)
     * 
     * @param inputFrameVectors Input set of vectors constituting an orthogonal basis
     * @param outputFrameVectors Resulting transported vectors
     * @param baseEl First point of the direction vector
     * @param destEl Second point of the direction vector
     */
    void SE3ParallelTransport(const std::array<SE3_ELEMENT, SE3_DIMENSIONS> &inputFrameVectors,
                              std::array<SE3_ELEMENT, SE3_DIMENSIONS> &outputFrameVectors,
                              const SE3_ELEMENT &baseEl,
                              const SE3_ELEMENT &destEl)
    {
        /* at first we implement the SO(3) parallel transport */

        /* compute logarithm map*/
        SE3_ELEMENT v;
        SE3LogarithmMap(baseEl, destEl, v);
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH; ++i) {
            v[i] *= 0.5;
        }

        /* keep only SO(3) part */
        for (size_t i = 12; i < SE3_ELEMENT_LENGTH; i++) {
            v[i] = 0.0;
        }

        /* compute matrix exponential */
        SE3_ELEMENT exp_aux;
        SO3Exponential(v, exp_aux);

        /* compute inv(dest)*base */
        SE3_ELEMENT inv_dest_base;
        SO3InvProduct_invAB(destEl, baseEl, inv_dest_base);

        /* compute left factor of product */
        SE3_ELEMENT left_factor;
        SO3Product_AB(inv_dest_base, exp_aux, left_factor);

        /* compute parallel transport for all basis elements */
        SE3_ELEMENT aux;
        for (size_t i = 0; i < SE3_DIMENSIONS - 3; ++i) {
            SO3Product_AB(left_factor, inputFrameVectors[i], aux);
            SO3Product_AB(aux, exp_aux, outputFrameVectors[i]);

            for (size_t j = 12; j < SE3_ELEMENT_LENGTH; j++) {
                outputFrameVectors[i][j] = 0.0;
            }
        }

        /* compute parallel transport for translational elements (flat geometry!) */
        for (size_t i = 3; i < 6; ++i) {
            outputFrameVectors[i] = inputFrameVectors[i];
        }

    };

    // TODO: find way to avoid branches in following functions

    /*!
     * Adjust the geodesic distance for the data term
     *
     * The computed distance passed per input is adjusted according
     * to the chosen regularization method and the lambda dampening parameter.
     *
     * @param lambda Dampening (0 <= lambda <= 1)
     * @param dist Input geodesic distance
     * @param p Regularization method
     * @return Updated geodesic distance
     */
    double dataPathLength(double lambda,
                          double dist,
                          Regularization p)
    {
        /* initialize tau */
        double tau = 0.0;

        /* change tau in case the geodesic distance is positive */
        if (dist > 0.0) {
            //tau determines how far we shoot !
            switch (p) {
                /* L^2-norm */
                case L2:
                    tau = lambda / (1.0 + lambda);
                    break;

                    /* L^1-norm => soft thresholding */
                case L1:
                    if (lambda < dist) {
                        tau = lambda / dist;
                    } else {
                        tau = 1.0;
                    }
                    break;

                    /* Huber-norm */
                case HUBER:

                    /* old implementation of Laurent, should be checked in a future version */
                    double tau_Huber = 1.0;
                    double omega = 1.0;

                    double l_t2_huber = lambda * tau_Huber * tau_Huber;

                    if (dist * sqrt(2.0) * tau_Huber < omega * (1.0 + 2.0 * l_t2_huber)) {
                        tau = (2.0 * l_t2_huber) / (1.0 + 2.0 * l_t2_huber);
                    } else {
                        if (sqrt(2.0) * lambda * omega * tau_Huber < dist) {
                            tau = sqrt(2.0) * lambda * omega * tau_Huber / dist;
                        } else {
                            tau = 1.0;
                        }
                    } /* end if(dist*sqrt(2.0)...*/

                    break;

            } /* end switch (p) */
        } /* end if ( dist > 0 ) */

        return tau;
    }

    /*!
     * Adjust the geodesic distance for the regularization term
     *
     * The computed distance passed per input is adjusted according
     * to the chosen regularization method and the alpha and lambda
     * dampening parameters.
     *
     * @param alpha Dampening (0 <= alpha <= 1)
     * @param lambda Dampening (0 <= lambda <= 1)
     * @param dist Input geodesic distance
     * @param p Regularization method
     * @return Updated geodesic distance
     */
    double regPathLength(double alpha,
                         double lambda,
                         double dist,
                         Regularization q)
    {
        /* initialize tau */
        double tau = 0.0;

        /* change tau in case the geodesic distance is positive */
        if (dist > 0.0) {
            switch (q) {
                /* L^2-norm */
                case L2:
                    tau = lambda * alpha / (2.0 * alpha * lambda + 1.0);
                    break;

                    /* L^1-norm => soft thresholding */
                case L1:
                    tau = lambda * alpha / dist;
                    if (tau > 0.5) {
                        tau = 0.5;
                    }
                    break;

                    /* Huber-norm */
                case HUBER:

                    /* old implementation of Laurent, should be checked in a future version */
                    double tau_Huber = 1.0;
                    double omega = 1.0;

                    if (dist < omega * (1.0 + 4.0 * lambda * alpha * tau_Huber * tau_Huber) / (sqrt(2.0) * tau_Huber)) {
                        tau = (2.0 * alpha * lambda * tau_Huber * tau_Huber) /
                              (1.0 + 4.0 * alpha * lambda * tau_Huber * tau_Huber);
                    } else {
                        tau = sqrt(2.0) * alpha * lambda * omega * tau_Huber / dist;
                        if (tau > 0.5) {
                            tau = 0.5;
                        }
                    } /* end if ( dist < omega*... */
                    break;

            } /* end switch (q) */
        } /* end if ( dist > 0.0 ) */

        return tau;
    }

    /*!
     * Compute the 1st order proximal map
     *
     * @param r Regularization method
     * @param x1 First element
     * @param x2 Second element
     * @param lambda Dampening parameter
     * @param y1 Resulting updated first element
     */
    void SE3ProximalMapFirstOrder(Regularization r, const SE3_ELEMENT &x1, const SE3_ELEMENT &x2,
                                  const double &lambda, SE3_ELEMENT &y1)
    {
        /* initialize memory for intermediate results */
        SE3_ELEMENT cache;

        /* compute logarithm and distance */
        double dist = SE3LogarithmMap(x1, x2, cache);

        /* calculate time for data term geodesic */
        double t = dataPathLength(lambda, dist, r);

        /* multiply computed logarithm with t */
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH; ++i) {
            cache[i] *= t;
        }

        /* compute proximal mapping for data term */
        SE3ExponentialMap(x1, cache, y1);
    };

//    void SE3ProximalMapFirstOrder(Regularization r, const SE3_ELEMENT &startEl, const SE3_ELEMENT &destEl,
//                                  const double &alpha, const double &lambda, SE3_ELEMENT &resElStart)
//    {
//        /* initialize memory for intermediate results */
//        SE3_ELEMENT cache;
//
//        /* compute logarithm and distance */
//        double dist = SE3LogarithmMap(startEl, destEl, cache);
//
//        /* calculate time for data term geodesic */
//        double t = regPathLength(alpha, lambda, dist, r);
//
//        /* multiply computed logarithm with t */
//        for (size_t i = 0; i < SE3_ELEMENT_LENGTH; ++i) {
//            cache[i] *= t;
//        }
//
//        /* compute proximal mapping for data term */
//        SE3ExponentialMap(startEl, cache, resElStart);
//    };

    /*!
     * Compute the 1st order proximal map
     *
     * @param r Regularization method
     * @param x1 First element
     * @param x2 Second element
     * @param alpha Dampening parameter
     * @param lambda Dampening parameter
     * @param y1 Resulting updated first element
     * @param y2 Resulting updated second element
     */
    void SE3ProximalMapFirstOrder(Regularization r, const SE3_ELEMENT &x1, const SE3_ELEMENT &x2,
                                  const double &alpha, const double &lambda,
                                  SE3_ELEMENT &y1, SE3_ELEMENT &y2)
    {
        /* initialize memory for intermediate results */
        SE3_ELEMENT cache1;
        SE3_ELEMENT cache2;

        /* compute logarithm and distance */
        double dist = SE3LogarithmMap(x1, x2, cache1);
        SE3LogarithmMap(x2, x1, cache2);

        /* calculate times for data term geodesics */
        double t = regPathLength(alpha, lambda, dist, r);

        /* multiply computed logarithm with t */
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH; ++i) {
            cache1[i] *= t;
            cache2[i] *= t;
        }

        /* compute proximal mapping for data term */
        SE3ExponentialMap(x1, cache1, y1);
        SE3ExponentialMap(x2, cache2, y2);
    };

    /*!
     * Compute the 2nd order proximal map
     *
     * @param r Regularization method
     * @param x1 First element
     * @param x2 Second element
     * @param x3 Third element
     * @param beta Dampening parameter
     * @param lambda Dampening parameter
     * @param y1 Resulting updated first element
     * @param y2 Resulting updated second element
     * @param y3 Resulting updated third element
     */
    void SE3ProximalMapSecondOrder(Regularization r,
                                   const SE3_ELEMENT &x1,
                                   const SE3_ELEMENT &x2,
                                   const SE3_ELEMENT &x3,
                                   double beta,
                                   double lambda,
                                   double inner_factor,
                                   int inner_steps,
                                   SE3_ELEMENT &y1,
                                   SE3_ELEMENT &y2,
                                   SE3_ELEMENT &y3)
    {
        /* initialization */
        for (size_t i = 0; i < SE3_ELEMENT_LENGTH; ++i) {
            y1[i] = x1[i];
            y2[i] = x2[i];
            y3[i] = x3[i];
        }

        /* adjust beta according to damping parameters */
        beta = beta * lambda;

        /* main iteration (gradient descent) */
        SE3_ELEMENT cache;
        for (int i = 0; i < inner_steps; ++i) {

            /* compute mean of y1 and y3, i.e., [y1, y3]_(d/2) */
            SE3_ELEMENT ym;
            SE3LogarithmMap(y1, y3, cache);
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                cache[j] *= 0.5;
            }
            SE3ExponentialMap(y1, cache, ym);

            /* compute direction at xm and length of 2nd difference */
            SE3_ELEMENT log_ym_y2;
            double length2ndDiff = SE3LogarithmMap(ym, y2, log_ym_y2);

            /* compute tau */
            //tau = (inner_factor/beta)/std::pow((double)(i)+1.0,0.95 + 0.5*std::pow(i+1,-0.18));
            double tau = inner_factor / std::pow((double) (i) + 1.0, 0.95 + 0.5 * std::pow(i + 1, -0.18));

            /******** gradient for y1 ********/
            SE3_ELEMENT log_y1_y3;
            SE3_ELEMENT log_y3_y1;
            SE3_ELEMENT log_y1_x1;
            SE3_ELEMENT log_y2_x2;
            SE3_ELEMENT log_y3_x3;

            /* compute direction of geodesic */
            SE3LogarithmMap(y3, y1, log_y3_y1);

            /* compute basis of eigenvectors at y3 */
            std::array<SE3_ELEMENT, SE3_DIMENSIONS> basisElements1{};
            std::array<SE3_ELEMENT, SE3_DIMENSIONS> basisElements2{};
            std::array<double, SE3_DIMENSIONS> eigenvalues{};
            SE3ComputeFrame(log_y3_y1, basisElements1, eigenvalues);

            /* move basis to xm via parallel transport */

            SE3ParallelTransport(basisElements1, basisElements2, y3, ym);

            /* compute coefficients at ym */
            std::array<double, SE3_DIMENSIONS> coeffs1{};
            std::array<double, SE3_DIMENSIONS> coeffs2{};

            for (size_t k = 0; k < SE3_DIMENSIONS; ++k) {
                coeffs1[k] = SE3InnerProduct(basisElements2[k], log_ym_y2);
            }

            /* move basis towards y1 */
            SE3ParallelTransport(basisElements2, basisElements1, ym, y1);

            /* adjust coefficients at y1 */
            double t = SE3InnerProduct(log_y3_y1, log_y3_y1);
            SE3AdjustCoefficientsMovingFrame(coeffs1, coeffs2, eigenvalues, t);

            /* compute basis x coefficients */
            SE3_ELEMENT basis_x_coeffs;
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                basis_x_coeffs[j] = 0.0;
            }

            for (size_t k = 0; k < SE3_DIMENSIONS; ++k) {
                for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                    basis_x_coeffs[j] += coeffs2[k] * basisElements1[k][j];
                }
            }

            /* check for regularization */
            double factor = 0.0;
            if (r == L2) {
                factor = beta;
            }
            if (r == L1) {
                //factor = std::sqrt(SE3InnerProduct( basis_x_coeffs, basis_x_coeffs, y1 ));

                if (length2ndDiff > EPSILON) {
                    factor = beta / length2ndDiff;
                } else {
                    factor = 0.0;
                }
            }

            /* adjust length of basis_x_coeffs and compose gradient at y1 */
            SE3_ELEMENT grad_1;
            SE3LogarithmMap(y1, x1, log_y1_x1);
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                grad_1[j] = tau * (basis_x_coeffs[j] * factor - log_y1_x1[j]);
            }

            /******** gradient for y2 ********/

            /* compute log( y2, ym ) */
            SE3_ELEMENT log_y2_ym;
            SE3LogarithmMap(y2, ym, log_y2_ym);

            /* check for regularization */
            if (r == L2) {
                factor = beta;
            }
            if (r == L1) {
                //factor = std::sqrt(SE3InnerProduct( log_y2_ym, log_y2_ym, y2 ));
                if (length2ndDiff > EPSILON) {
                    factor = beta / length2ndDiff;
                } else {
                    factor = 0.0;
                }
            }

            /* adjust length of basis_x_coeffs and compose gradient at y2 */
            SE3_ELEMENT grad_2;
            SE3LogarithmMap(y2, x2, log_y2_x2);
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                grad_2[j] = tau * (log_y2_ym[j] * factor - log_y2_x2[j]);
            }

            /******** gradient for y3 ********/

            /* compute direction of geodesic */
            SE3LogarithmMap(y1, y3, log_y1_y3);

            /* compute basis of eigenvectors at y1 */
            SE3ComputeFrame(log_y1_y3, basisElements1, eigenvalues);

            /* move basis to xm via parallel transport */
            SE3ParallelTransport(basisElements1, basisElements2, y1, ym);

            /* compute coefficients at ym */
            for (size_t k = 0; k < SE3_DIMENSIONS; ++k) {
                coeffs1[k] = SE3InnerProduct(basisElements2[k], log_ym_y2);
            }

            /* move basis towards y3 */
            SE3ParallelTransport(basisElements2, basisElements1, ym, y3);

            /* adjust coefficients at y3 */
            t = SE3InnerProduct(log_y1_y3, log_y1_y3);
            SE3AdjustCoefficientsMovingFrame(coeffs1, coeffs2, eigenvalues, t);

            /* compute basis x coefficients */
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                basis_x_coeffs[j] = 0.0;
            }

            for (size_t k = 0; k < SE3_DIMENSIONS; ++k) {
                for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                    basis_x_coeffs[j] += coeffs2[k] * basisElements1[k][j];
                }
            }

            /* check for regularization */
            if (r == L2) {
                factor = beta;
            }
            if (r == L1) {
                //factor = std::sqrt(SE3InnerProduct( basis_x_coeffs, basis_x_coeffs, y3 ));

                if (length2ndDiff > EPSILON) {
                    factor = beta / length2ndDiff;
                } else {
                    factor = 0.0;
                }
            }

            /* adjust length of basis_x_coeffs and compose gradient at y1 */
            SE3_ELEMENT grad_3;
            SE3LogarithmMap(y3, x3, log_y3_x3);
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                grad_3[j] = tau * (basis_x_coeffs[j] * factor - log_y3_x3[j]);
            }

            /* update step */
            SE3_ELEMENT aux1;
            SE3_ELEMENT aux2;
            SE3_ELEMENT aux3;
            for (size_t j = 0; j < SE3_ELEMENT_LENGTH; ++j) {
                aux1[j] = y1[j];
                aux2[j] = y2[j];
                aux3[j] = y3[j];
            }
            SE3ExponentialMap(aux1, grad_1, y1);
            SE3ExponentialMap(aux2, grad_2, y2);
            SE3ExponentialMap(aux3, grad_3, y3);

        } /* end main iteration */
    };

    // documented in header
    void proximalPointAlgorithm(const std::vector<SE3_ELEMENT> &f,
                                std::vector<SE3_ELEMENT> &x,
                                Regularization p,
                                Regularization q,
                                Regularization r,
                                double inner_factor,
                                double alpha,
                                double beta,
                                int steps,
                                int inner_steps)
    {
        /* initialize x */
        x = f;

        /* main iteration */
        for (int ll = 0; ll < steps; ++ll) {

            /* compute lambda */
            double lambda = 1.0 / pow((double) (ll + 1), 0.95 + 0.5 * pow((double) (ll + 1), -0.18));

            /* proximal mapping of data term */
            for (size_t i = 0; i < f.size(); i++) {
                SE3_ELEMENT cache;
                SE3ProximalMapFirstOrder(p, x[i], f[i], lambda, cache);
                x[i] = cache;
            }

            /* 1st order proximal mappings */
            if (alpha > 0) {
                /* 1st order proximal mapping w.r.t. x-direction */
                for (size_t i = 0; i < f.size() - 1; i++) {
                    SE3_ELEMENT cache1, cache2;
                    SE3ProximalMapFirstOrder(q, x[i], x[i + 1], alpha, lambda, cache1, cache2);
                    x[i] = cache1;
                    x[i + 1] = cache2;
                }
            }

            /* 2nd order proximal mappings */
            if (beta > 0) {

                /* 2nd order proximal mapping w.r.t. x-direction */
                for (size_t i = 1; i < f.size() - 1; i++) {
                    SE3_ELEMENT cache1, cache2, cache3;
                    SE3ProximalMapSecondOrder(r, x[i - 1], x[i], x[i + 1], beta, lambda, inner_factor, inner_steps,
                                              cache1, cache2, cache3);
                    x[i - 1] = cache1;
                    x[i] = cache2;
                    x[i + 1] = cache3;
                }
            }

            /* REMARK: Usually this algorithm would require the implementation of proximal mappings in more directions,
             *         and mixed proximal mappings (diagonal directions). However, these have been
             *         removed here, because pose denoising is currently only done in 1D (e.g. time)!
             */

        } /* end main iteration */
    };


} /* end namespace */