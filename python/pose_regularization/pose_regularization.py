
from enum import IntEnum
import numpy as np
from . import pypose_regularization


class Regularization(IntEnum):
    HUBER = 0
    L1 = 1
    L2 = 2

def regularize_matrices_4x4(pose_matrices,
    p, q, r, inner_factor,
    alpha, beta, steps, inner_steps):
    # type: (np.ndarray, Regularization, Regularization, Regularization, float, float, float, int, int) -> np.ndarray
    """Regularize a series of 6 DoF poses
    
    The input sequence of 4x4 Euclidean transformation matrices is regularized
    according to the provided parameters, and the result is returned.
    
    The action of each parameter is discussed below; for best results,
    they should be chosen carefully (for example, with a grid search
    study). However, a good starting point for the choice of the parameters
    in case of strong jitter noise is the following:
    poses_regularized = regularize_matrices_4x4(poses_noisy, L2, HUBER, L2, 50,
    1, 0.25, 100, 50)

    If faster execution is desired, 2nd order regularization can be
    disabled (it is worthwhile to study its effect anyway):
    poses_regularized = regularize_matrices_4x4(poses_noisy, L2, HUBER, HUBER, 50,
    0, 0, 100, 0)

    poses_noisy is expected to be a stack of 4x4 matrices, representing
    poses in the Euclidean space with respect to a single coordinate
    system. If N is the number of poses, this Numpy array should have shape
    Nx4x4.

    p, q, and r define the regularization method for the data term, and
    for the 1st order and 2nd order output regularization terms
    respectively. 0 defines the HUBER norm, 1 the L1 regularization
    (linear) and 2 the L2 regularization method (quadratic).

    alpha and beta are dampening coefficients for the 1st and 2nd order
    regularization, respectively. Either of the two regularizations can be
    disabled by setting the relative coefficient to zero.

    inner_factor is a further parameter for the action of the 2nd order
    regularizer.

    steps is the number of global regularization steps performed.

    inner_steps is the number of steps performed within the 2nd order
    regularization process.

    :param pose_matrices: input pose signal (Nx4x4 numpy array)
    :param p: regularizer for data fidelity term (see the Regularization class)
    :param q: regularizer for 1st order regularization term (see the Regularization class)
    :param r: regularizer for 2nd order regularization term (see the Regularization class)
    :param inner_factor: further parameter for the action of the 2nd order
    regularizer
    :param alpha: damping coefficient for 1st order regularization
    :param beta: damping coefficient for 2nd order regularization
    :param steps: number of global regularization steps performed
    :param inner_steps: number of steps performed within the 2nd order
    regularization process
    :return: output pose signal (Nx4x4 numpy array)
    """
    return pypose_regularization.proximal_point_algorithm(pose_matrices,
                                                          int(p), int(q), int(r), inner_factor,
                                                          alpha, beta, steps, inner_steps)