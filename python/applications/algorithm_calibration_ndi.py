from __future__ import division

import h5py
import numpy as np
from matplotlib import pyplot as plt
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


fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

smoothed_poses = np.reshape(smoothed_column_poses, (number_of_poses, 4, 4))

utils.plot_transformation_matrix_array(sensor_in_world, ax1)
utils.plot_transformation_matrix_array(smoothed_poses, ax2)


def on_move(event):
    if event.inaxes == ax1:
        ax2.view_init(elev=ax1.elev, azim=ax1.azim)
    elif event.inaxes == ax2:
        ax1.view_init(elev=ax2.elev, azim=ax2.azim)
    else:
        return
    fig.canvas.draw_idle()

c1 = fig.canvas.mpl_connect('motion_notify_event', on_move)

plt.show()
