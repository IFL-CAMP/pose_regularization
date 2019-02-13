# Total Variation Regularization of Pose Streams

This repository contains the implementation of the method proposed in 
```
@article{esposito2019total,
  title={Total Variation Regularization of Pose Signals
         with an Application to 3D Freehand Ultrasound},
  author={Esposito, Marco and Hennersperger, Christoph and Göbl, Rüdiger and Demaret, Laurent and Storath, Martin and Navab, Nassir
           and Baust, Maximilian and Weinmann, Andreas},
  journal={IEEE Transactions on Medical Imaging},
  volume={TBD},
  number={TBD},
  pages={TBD},
  year={2019},
  publisher={ IEEE }
}
```
(PREPRINT available [here](https://ieeexplore.ieee.org/document/8638838)) for the regularization of 6-degrees-of-freedom pose streams, such as the output of a tracking system. The method has 
the ability to compensate for significant jitter, while maintaining the features of the original signal (e.g. sudden 
and sharp variations in the movement direction). 

This is achieved by minimizing the following functional:

![equation](https://latex.codecogs.com/svg.latex?E%28%5Cmathbf%7Bx%7D%29%20%3D%20D%28%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7Bp%7D%29%20&plus;%20%5Calpha%20R%28%5Cmathbf%7Bx%7D%29)

where the ![equation](https://latex.codecogs.com/svg.latex?D%28%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7Bp%7D%29) term penalizes 
the deviation of the output signal ![equation](https://latex.codecogs.com/svg.latex?%5Cmathbf%7Bx%7D) from the original 
signal ![equation](https://latex.codecogs.com/svg.latex?%5Cmathbf%7Bp%7D), and 
![equation](https://latex.codecogs.com/svg.latex?R%28%5Cmathbf%7Bx%7D%29) penalizes large variations within consecutive 
poses of the output stream. The joint action of the two functionals, balanced by the coefficient 
![equation](https://latex.codecogs.com/svg.latex?%5Calpha%20%3E%200), allows to find the desired trade-off between the 
fidelity to the original signal and the regularity of the output.

This implementation includes first- and second-order regularization for the 
![equation](https://latex.codecogs.com/svg.latex?R%28%5Cmathbf%7Bx%7D%29) term, and a first-order regularizer for 
![equation](https://latex.codecogs.com/svg.latex?D%28%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7Bp%7D%29). Each regularizer can 
be combined with an L1, L2 or HUBER penalization.

For further details on the algorithm, the reader is referred to the published paper. <!-- TODO: add link --> 
Practical instructions to use the algorithm, including a brief description of the action of each parameter, is provided 
in the documentation of each interface (Doxygen for C++, MATLAB help, Python docstrings). <!-- TODO: add shared section with more detail --> 

## License

This software is released under the LGPLv3 license. Hence, it can be used and included into proprietary software; 
however, any modification to the code of the library itself must be released publicly.

## Usage

The software is implemented in C++. A MATLAB and a Python interface are also provided.

### C++

The library can be built using CMake. There are no external dependencies. C++11 is required.

### MATLAB

A `compile.m` script is provided in the `matlab/wrapper` folder. This script can be executed to compile a `mex` 
extension.

After running the compile script, the `matlab` directory must be added recursively to your path to use the contained 
functions. <!-- TODO: document functions again? -->

Alternatively, the `mex` file can be copied together with the scripts into your path. 

More details are available in the respective [page](matlab/README.md).

### Python

Coming soon as a PyPI package.

Meanwhile, the Cython extension can be built from source and installed as detailed in the respective 
[instructions](python/README.md).

<!-- TODO: comment tests? -->