from __future__ import division

import h5py
import numpy as np
from pose_denoising import proximal_point_algorithm
import utils

data = h5py.File('../../data/testData.mat')

sensor_in_world = data['sensorInWorld']

number_of_poses = len(sensor_in_world)

column_poses = np.reshape(sensor_in_world, (number_of_poses, 16))
smoothed_column_poses = np.zeros(column_poses.shape)

proximal_point_algorithm(column_poses, smoothed_column_poses,
                         (number_of_poses, 1, 1),
                         1, 0, 0, 5.0,
                         2, 0.0, 0.0, 100, 0)

smoothed_poses = np.reshape(smoothed_column_poses, (number_of_poses, 4, 4))

utils.plot_transformation_matrix_array(smoothed_poses)