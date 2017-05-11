
cimport numpy as np

np.import_array()

cdef extern from '../include/manifold.h' namespace 'manifold':
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
                                 int inner_steps
    );

def proximal_point_algorithm(
        np.ndarray in_array not None,
        np.ndarray out_array not None,
        grid_shape,
        p, q, r, inner_factor,
        alpha, beta, gamma, steps, inner_steps
    ):

    proximalPointAlgorithm(
        <double*> np.PyArray_DATA(in_array),
        <double*> np.PyArray_DATA(out_array),
        grid_shape[0],
        grid_shape[1],
        grid_shape[2],
        p,
        q,
        r,
        inner_factor,
        alpha,
        beta,
        gamma,
        steps,
        inner_steps)
