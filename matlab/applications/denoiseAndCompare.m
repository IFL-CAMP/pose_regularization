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

denoised_target_data = denoiseMatrices(target_data);

DO_COMPARATIVE_PLOTS = true;
DO_PLOTS = false;

if DO_PLOTS
    subplot(2,2,1)
    plotPoseMatrix(ground_truth_data,2.0)
    title(ground_truth, 'Interpreter', 'none')

    subplot(2,2,2)
    plotPoseMatrix(target_data,2.0)
    title(comparison_target, 'Interpreter', 'none')

    subplot(2,2,3)
    plotPoseMatrix(denoised_target_data,2.0)
    title(strcat(comparison_target, ' denoised'), 'Interpreter', 'none')
end

if DO_COMPARATIVE_PLOTS
    subplot(1,2,1)
    plotPoseMatrix(ground_truth_data(:,:,1:5:end),2.0)
    plotPoseMatrix(target_data(:,:,1:5:end),2.0)
    title('original')

    subplot(1,2,2)
    plotPoseMatrix(ground_truth_data(:,:,1:5:end),2.0)
    plotPoseMatrix(denoised_target_data(:,:,1:5:end),2.0)
    title('denoised')
end

fprintf('Distance to ground truth: %6.4f\n', comp_dist(ground_truth_data, target_data))
fprintf('Distance to ground truth after denoising: %6.4f\n', comp_dist(ground_truth_data, denoised_target_data))