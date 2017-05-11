import numpy as np
from pose_denoising import proximal_point_algorithm


def test_trivial():
    identities = np.array([np.eye(4)] * 10)
    column_identities = np.reshape(identities, (10, 16))
    smoothed_identities = np.zeros(column_identities.shape)
    proximal_point_algorithm(column_identities, smoothed_identities,
                             (identities.shape[2], 1, 1),
                             1, 0, 0, 5.0,
                             2, 0.0, 0.0, 100, 0)
    assert np.array_equal(smoothed_identities, column_identities)
