close all;
clear;
clc;

addpath('..')

load('../../data/algorithmCalibrationNDI')

% get number of poses
N = size(trajectoryNDI,3);

% convert pose tracks for algorithm
NDI = zeros(16,N);
for i = 1:N
    P = trajectoryNDI(:,:,i);
    NDI(:,i) = P(:);
end

% run algorithm
NDInew = denoiseSE3prod(NDI,size(trajectoryNDI,3),1,1,...
                        2,2,0,1.0,...
                        12,0.0,0.0, 100,0);


% re-format poses
trajectoryNDInew  = trajectoryNDI;
for i = 1:N
    P = NDInew(:,i);
    trajectoryNDInew(:,:,i) = reshape(P,4,4);
end

% convert poses to cell an dcompute errors
robotCell{1} = [];
ndiCell{1} = [];
ndiCellnew{1} = [];
RMSbefore = 0;
RMSafter = 0;
for i = 1:N
    robotCell{i} = trajectoryRob(:,:,i);
    ndiCell{i} = trajectoryNDI(:,:,i);
    ndiCellnew{i} = trajectoryNDInew(:,:,i);
    
    E = trajectoryRob(:,:,i) - trajectoryNDI(:,:,i);
    RMSbefore = RMSbefore + sum(E(:).^2);
    E = trajectoryRob(:,:,i) - trajectoryNDInew(:,:,i);
    RMSafter = RMSafter + sum(E(:).^2);
    
end

RMSbefore = sqrt( RMSbefore/N )
RMSafter = sqrt( RMSafter/N )


offset = 50;

subplot(1,3,1)
% plotPose(robotCell(offset+1:offset+100),2.0)
plotPoseMatrix(trajectoryRob(:,:,offset+1:offset+400),2.0)
title('KUKA Poses (ground truth)')

subplot(1,3,2)
% plotPose(ndiCell(offset+1:offset+100),2.0)
plotPoseMatrix(trajectoryNDI(:,:,offset+1:offset+400),2.0)
title('NDI Poses')

subplot(1,3,3)
% plotPose(ndiCellnew(offset+1:offset+100),2.0)
plotPoseMatrix(trajectoryNDInew(:,:,offset+1:offset+400),2.0)
title('NDI Poses Denoised')