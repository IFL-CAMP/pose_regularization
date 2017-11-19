close all;
clear;
clc;

addpath('..')

load('paramsFirstOrder.mat')
load('paramsFirstAndSecondOrder.mat')

file_prefix = '/home/marco/.ros/pose_denoising/logs/test_square_out_of_plane_2017-11-02-15-00-32';

ground_truth = 'robot';
comparison_target = 'em_integrated';

USE_ICP = true;
if USE_ICP
    file_suffix = '_probeviz_registered.mat';
else
    file_suffix = '_probeviz.mat';
end

ground_truth_data = load(strcat(file_prefix, '_', ground_truth, '_probeviz.mat'));
ground_truth_data = ground_truth_data.tracking;

target_data = load(strcat(file_prefix, '_', comparison_target, file_suffix));
target_data = target_data.tracking;

% compute reference error
err_ref = comp_dist(ground_truth_data, target_data);

RESULTS3 = {};
for bestResult2Index = 1:5
    results = zeros(5, 7);
        
    % find optimal parameters
    e = RESULTS2{bestResult2Index}(:,2);
    [sortede, sortedeindices] = sort(e);
    eIndicesToTry = sortedeindices(1:5);
    for i = 1:length(eIndicesToTry)
        index = eIndicesToTry(i);

        p = RESULTS2{bestResult2Index}(index,3);
        q = RESULTS2{bestResult2Index}(index,4);
        r = RESULTS2{bestResult2Index}(index,5);
        alpha = RESULTS2{bestResult2Index}(index,6);
        beta = RESULTS2{bestResult2Index}(index,7);

        denoised_target_data = denoiseMatrices(target_data, p, q, r, alpha, beta, 0, 0.25, 100, 50);

        err = comp_dist(ground_truth_data, denoised_target_data);
        results(i, :) = [err_ref, err, p, q, r, alpha, beta];
    end
    RESULTS3{bestResult2Index} = results;
end

save('paramsTracking.mat','RESULTS3');
