# Pose Regularization

This directory contains a Cython wrapper for the pose regularization library.

### Building the Python extension

Numpy and Cython must be installed in order to build this extension. It is recommended 
to use a virtual environment, as usual.

Once in the `python` directory, you can build and install the extension into your environment:
`pip install .`

### Testing the Python extension

After building and installing the extension as in the previous section, you can go into the 
`python/pose_regularization/test` directory, and execute the tests:

`python test_pose_regularization.py`