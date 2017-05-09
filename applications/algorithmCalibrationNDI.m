close all;
clear;
clc;

addpath('..')

load('../data/testData.mat')

% get number of poses
N = size(sensorInWorld,3);

% convert pose tracks for algorithm
NDI = zeros(16,N);
for i = 1:N
    P = sensorInWorld(:,:,i);
    NDI(:,i) = P(:);
end

% run algorithm without 2nd order regularization
NDInew = denoiseSE3prod(NDI,size(sensorInWorld,3),1,1,...
                        1,0,0,5.0,...
                        2,0.0,0.0,100,0);


% re-format poses
trajectoryNDInew  = sensorInWorld;
for i = 1:N
    P = NDInew(:,i);
    trajectoryNDInew(:,:,i) = reshape(P,4,4);
end


subplot(1,2,1)
plotPoseMatrix(sensorInWorld,2.0)
title('NDI Poses')

subplot(1,2,2)
plotPoseMatrix(trajectoryNDInew,2.0)
title('NDI Poses Denoised')