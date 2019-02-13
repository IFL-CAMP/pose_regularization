# Pose Regularization

This directory contains a Cython wrapper for the pose regularization library.

### Building the Python extension

Numpy and Cython must be installed in order to build this extension. It is recommended 
to use a virtual environment, as usual.
```$bash
virtualenv -p <PATH_TO_YOUR_INTERPRETER> <DESIRED_VIRTUALENV_LOCATION>
source <DESIRED_VIRTUALENV_LOCATION>/bin/activate
```

Once in the `python` directory, you can build and install the extension into your environment:
```
cd <REPO_ROOT>/python
pip install numpy cython
python setup.py bdist_wheel
pip install dist/pose_regularization-*.whl
```

The resulting wheel file in `python/dist` can be saved for later use or redistributed 
(to compatible systems).

### Using the package

The algorithm is implemented as a function `regularize_matrices_4x4`, which takes as input a 
stack of 4x4 Euclidean transformation matrices and a set of parameters. The returned output is 
the regularization result, in the same format as the input.

The following example shows the application of the algorithm to dummy data. For the action of 
each parameter, please refer to the docstrings.

```$python
import numpy as np
from pose_regularization import regularize_matrices_4x4, Regularization as R

poses_noisy = np.array([np.eye(4)]*5)  # replace with your input signal

poses_regularized = regularize_matrices_4x4(poses_noisy, 
    R.L1, R.L2, R.HUBER, 0.25,
    50, 1, 100, 50) 

``` 

### Testing the Python extension

After building and installing the extension as in the previous section, you can go into the 
`python/pose_regularization/test` directory, and execute the tests:
```$bash
cd <REPO_ROOT>/python/pose_regularization/test
python test_pose_regularization.py
```