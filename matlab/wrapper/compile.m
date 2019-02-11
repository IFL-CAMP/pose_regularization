close all;
clear;
clc;

% compile SE3prod manifold as release
mex -I../../include regularizeSE3.cpp ../../src/pose_regularization.cpp

% % compile SE3prod manifold as debug
% mex -g -I../../include regularizeSE3.cpp ../../src/pose_regularization.cpp


