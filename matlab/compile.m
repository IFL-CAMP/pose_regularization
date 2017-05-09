close all;
clear;
clc;

% compile SE3prod manifold as release
mex -DSE3PROD -I../include denoiseSE3prod.cpp ../src/manifoldGeneric.cpp ../src/manifoldSE3prod.cpp

% % compile SE3prod manifold as debug
% mex -g -DSE3PROD -I../include denoiseSE3prod.cpp ../src/manifoldGeneric.cpp ../src/manifoldSE3prod.cpp


