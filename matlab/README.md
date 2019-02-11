# Pose Regularization Toolbox

This directory contains a MATLAB MEX wrapper for the pose regularization library. 
MATLAB interface scripts are included for convenience.


### Building the MEX file

A MATLAB compile script is included. It is sufficient to go into the `matlab/wrapper` 
directory and call the MATLAB `compile` function.

No particular dependencies are needed. If compilation fails, please check that 
you are able to compile example MEX files.


### Using the package

The library allows to regularize an input stream of poses in Euclidean space, 
thus removing jitter. Employing an explicit fidelity term allows to maintain 
features of the original signal, such as sharp variations in the direction of 
movement.

The input signal, for example the output of a tracking system, should be expressed 
as a series of 4x4 transformation matrices (where the top-left 3x3 submatrix is a 
rotation matrix, and the 3x1 vector next to it is the translation vector). 

After adding the `matlab` directory to your path (recursively), you can 
use the `regularizeMatrices4x4` function as follows:

`poses_regularized = regularizeMatrices4x4(poses_noisy, p, q, r, alpha,
 beta, inner_factor, steps, inner_steps);`
 
Please refer to the global documentation of the library for the effect of each 
parameter. 