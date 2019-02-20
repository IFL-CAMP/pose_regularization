from libcpp.vector cimport vector
cimport numpy as np
import numpy as np

np.import_array()

cdef extern from "<array>" namespace "std" nogil:
  cdef cppclass se3_element "std::array<double, 16>":
    se3_element() except+
    double& operator[](size_t)

cdef extern from "pose_regularization.h" namespace "pose_regularization":
    cdef enum Regularization:
            HUBER = 0
            L1 = 1
            L2 = 2
    void proximalPointAlgorithm( const vector[se3_element] &f,
                                 vector[se3_element] &x,
                                 Regularization p,
                                 Regularization q,
                                 Regularization r,
                                 double inner_factor,
                                 double alpha,
                                 double beta,
                                 int steps,
                                 int inner_steps
    );
cdef extern from "../src/pose_regularization.cpp":
    pass

def proximal_point_algorithm(
        np.ndarray in_array not None,
        p, q, r, inner_factor,
        alpha, beta, steps, inner_steps
    ):

    # TODO: find a better way
    cdef vector[se3_element] in_vec = vector[se3_element](in_array.shape[0])
    # numpy is row-major
    for i in range(in_array.shape[0]):
        for j in range(4):
            for k in range(4):
                in_vec[i][j*4+k] = in_array[i,k,j]
    cdef vector[se3_element] out_vec

    proximalPointAlgorithm(
        in_vec,
        out_vec,
        p,
        q,
        r,
        inner_factor,
        alpha,
        beta,
        steps,
        inner_steps)

    # TODO: reformat out_array
    out_array = np.zeros_like(in_array)
    for i in range(out_array.shape[0]):
        for j in range(4):
            for k in range(4):
                out_array[i,k,j] = out_vec[i][j*4+k]
    return out_array