
cimport numpy as np

np.import_array()

cdef extern from 'manifold.h' namespace 'manifold' nogil:
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

    cdef double* c_in_array = <double*> np.PyArray_DATA(in_array)
    cdef double* c_out_array = <double*> np.PyArray_DATA(out_array)
    cdef int c_grid_shape_x = grid_shape[0]
    cdef int c_grid_shape_y = grid_shape[1]
    cdef int c_grid_shape_z = grid_shape[2]
    cdef int c_p = p
    cdef int c_q = q
    cdef int c_r = r
    cdef double c_inner_factor = inner_factor
    cdef double c_alpha = alpha
    cdef double c_beta = beta
    cdef double c_gamma = gamma
    cdef int c_steps = steps
    cdef int c_inner_steps = inner_steps

    with nogil:
        proximalPointAlgorithm(
            c_in_array,
            c_out_array,
            c_grid_shape_x,
            c_grid_shape_y,
            c_grid_shape_z,
            c_p,
            c_q,
            c_r,
            c_inner_factor,
            c_alpha,
            c_beta,
            c_gamma,
            c_steps,
            c_inner_steps)
