function [ regularizedMatrices ] = regularizeMatrices4x4(matrices, p, q, r, alpha, beta, inner_factor, steps, inner_steps)
%regularizeMatrices4x4 Regularize a series of 6 DoF poses.
%   poses_regularized = regularizeMatrices4x4(poses_noisy, p, q, r, alpha,
%   beta, inner_factor, steps, inner_steps) regularizes the input pose
%   stream according to the defined parameters.
%
%   The action of each parameter is discussed below; for best results, 
%   they should be chosen carefully (for example, with a grid search 
%   study). However, a good starting point for the choice of the parameters 
%   in case of strong jitter noise is the following:
%   poses_regularized = regularizeMatrices4x4(poses_noisy, 2, 0, 2, 50,
%   1, 0.25, 100, 50);
%
%   If faster execution is desired, 2nd order regularization can be
%   disabled (it is worthwhile to study its effect anyway):
%   poses_regularized = regularizeMatrices4x4(poses_noisy, 2, 0, 0, 50,
%   0, 0, 100, 0);
%   
%   
%   poses_noisy is expected to be a stack of 4x4 matrices, representing
%   poses in the Euclidean space with respect to a single coordinate
%   system. If N is the number of poses, this array should have shape
%   4x4xN.
%
%   p, q, and r define the regularization method for the data term, and 
%   for the 1st order and 2nd order output regularization terms
%   respectively. 0 defines the HUBER norm, 1 the L1 regularization
%   (linear) and 2 the L2 regularization method (quadratic).
%
%   alpha and beta are dampening coefficients for the 1st and 2nd order 
%   regularization, respectively. Either of the two regularizations can be
%   disabled by setting the relative coefficient to zero.
%   
%   inner_factor is a further parameter for the action of the 2nd order
%   regularizer.
%
%   steps is the number of global regularization steps performed.
%
%   inner_steps is the number of steps performed within the 2nd order
%   regularization process.
%
%   See also regularizeMatricesColumn.

    assert(size(matrices, 1) == 4, 'expected stack of 4x4 matrices (shape 4x4xN)');
    assert(size(matrices, 2) == 4, 'expected stack of 4x4 matrices (shape 4x4xN)');

    matricesColumn = reshape(matrices, 16, []);
    
    regularizedMatricesColumn = regularizeMatricesColumn(matricesColumn, p, q, r, alpha, beta, inner_factor, steps, inner_steps);

    regularizedMatrices = reshape(regularizedMatricesColumn, 4, 4, []);
end

