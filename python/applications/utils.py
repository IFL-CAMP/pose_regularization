from __future__ import division

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_transformation_matrix_array(ma, ax=None, fig=None):

    translations = ma[:, 3, :3]
    x, y, z = translations[:, 0], translations[:, 1], translations[:, 2]

    # TODO: vectorize this stuff!
    x_vec = [1, 0, 0, 1]
    y_vec = [0, 1, 0, 1]
    z_vec = [0, 0, 1, 1]

    x_axis_comp = [np.dot(m, x_vec) for m in ma]
    y_axis_comp = [np.dot(m, y_vec) for m in ma]
    z_axis_comp = [np.dot(m, z_vec) for m in ma]

    if ax is None:
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    if matplotlib.__version__ < '2':
        x_axis_comp_norm = np.array([(xc/xc[-1]) for xc in x_axis_comp])
        y_axis_comp_norm = np.array([(xc/xc[-1]) for xc in y_axis_comp])
        z_axis_comp_norm = np.array([(xc/xc[-1]) for xc in z_axis_comp])
        ax.quiver(x, y, z, x_axis_comp_norm[:, 0], x_axis_comp_norm[:, 1], x_axis_comp_norm[:, 2], color='r')
        ax.quiver(x, y, z, y_axis_comp_norm[:, 0], y_axis_comp_norm[:, 1], y_axis_comp_norm[:, 2], color='b')
        ax.quiver(x, y, z, z_axis_comp_norm[:, 0], z_axis_comp_norm[:, 1], z_axis_comp_norm[:, 2], color='g')
    else:
        ax.quiver(x, y, z, x_axis_comp[:, 0], x_axis_comp[:, 1], x_axis_comp[:, 2], normalize=True, color='r')
        ax.quiver(x, y, z, y_axis_comp[:, 0], y_axis_comp[:, 1], y_axis_comp[:, 2], normalize=True, color='b')
        ax.quiver(x, y, z, z_axis_comp[:, 0], z_axis_comp[:, 1], z_axis_comp[:, 2], normalize=True, color='g')
