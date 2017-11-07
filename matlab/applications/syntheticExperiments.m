close all;
clear;
clc;

addpath('..')

% generate test data
[Pgt, Pnoise] = genTestData(100);

% get number of poses
N = size(Pgt,3);

% only 1st order huber regularization with 100 iterations
Pden1 = denoiseSE3prod(reshape(Pnoise,16,N),N,1,1,...
                        1,0,0,0,...
                        10,0,0, 1000,0);

Pden1 = reshape(Pden1,4,4,N);


% l1 penalty for data term, 1st-order TV and 2nd order TV
Pden2 = denoiseSE3prod(reshape(Pnoise,16,N),N,1,1,...
                       1,1,1,0.25,...
                       2,1,0,500,100);

Pden2 = reshape(Pden2,4,4,N);

% visulization
subplot(2,2,1)
plotPoseMatrix(Pnoise,2.0)
title('noisy poses')

subplot(2,2,2)
plotPoseMatrix(Pgt,2.0)
title('ground truth')

subplot(2,2,3)
plotPoseMatrix(Pden1,2.0)
title('only 1st order Huber')

subplot(2,2,4)
plotPoseMatrix(Pden2,2.0)
title('1st and 2nd order TV')
