
from . import pypose_regularization

def regularize_matrices_4x4(pose_matrices,
    p, q, r, inner_factor,
    alpha, beta, steps, inner_steps):
    return pypose_regularization.proximal_point_algorithm(pose_matrices,
                                                          p, q, r, inner_factor,
                                                          alpha, beta, steps, inner_steps)