close all;
clear;
clc;

% compile SE3prod manifold as release
mex -DSE3PROD -I./src/generic ./src/SE3prod/denoiseSE3prod.cpp ./src/generic/manifoldGeneric.cpp ./src/SE3prod/manifoldSE3prod.cpp

% % compile SE3prod manifold as debug
% mex -g -DSE3PROD -I./src/generic ./src/SE3prod/denoiseSE3prod.cpp ./src/generic/manifoldGeneric.cpp ./src/SE3prod/manifoldSE3prod.cpp


