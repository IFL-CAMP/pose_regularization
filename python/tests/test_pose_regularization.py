# /usr/bin/env python2

import os
import numpy as np
from scipy.io import loadmat
import pose_regularization

np.set_printoptions(precision=5)


def replicate_regularization(filename):

    data = loadmat(filename)

    p_n = np.rollaxis(data['Poses_Noisy'], 2, 0)
    p_r = np.rollaxis(data['Poses_Regularized'], 2, 0)

    b = pose_regularization.regularize_matrices_4x4(
        p_n, data['p'], data['q'], data['r'],
        data['inner_factor'], data['alpha'],
        data['beta'], data['steps'], data['inner_steps'])

    ATOL = 1e-07 # relaxed for portability

    assert(np.isclose(b, p_r, atol=ATOL).all())


TEST_DATA_DIR = '../test/data'


def test_enums():
    assert(int(pose_regularization.Regularization.HUBER) == 0)
    assert(int(pose_regularization.Regularization.L1) == 1)
    assert(int(pose_regularization.Regularization.L2) == 2)


def test_optical():
    replicate_regularization(os.path.join(TEST_DATA_DIR,'case5_preset_optical.mat'))


def test_em1():
    replicate_regularization(os.path.join(TEST_DATA_DIR,'case5_preset_em1.mat'))


def test_em2():
    replicate_regularization(os.path.join(TEST_DATA_DIR,'case5_preset_em2.mat'))

if __name__ == '__main__':
    test_enums()
    test_optical()
    test_em1()
    test_em2()
    print('passed!')
