close all;
clear;
%clc;

dbstop('error')

% define tolerance for relative error
TOL = 1e-4;

% define number of steps
steps = 1;
innerSteps = 1;
innerFactor = 1;

% define manifold to be tested
manifold = 'SE3prod';

%% prepare evaluation

% add main path
addpath('..')

if strcmp(manifold,'SE3prod')
    
    % define handles
    logHandle = @SE3prod.log;
    expHandle = @SE3prod.exp;
    prodHandle = @SE3prod.prod;
    frameHandle = @SE3prod.frame;
    parHandle = @SE3prod.par;
    adjustHandle = @SE3prod.adjust;
    
    % cpp handle
    cppHandle = @denoiseSE3prod;
    
    % load data
    load('../data/testDataSE3.mat')
    
    % restrict data to 11x11x11
    f = f(:,1:11^3);
end

%% compare the implementations

flag = true;

%% initialize parameters for 1st order case

alpha = 0;
q = -1;
gamma = 0; %param for mixed second order

% 1D tests
for beta = [0.5 1 2]
    for p = [0 1 2]
        for r = [1]
            
            % run MATLAB implementation
            x_m = proximalPointAlgorithm( f(:,1:11),11, 1, 1,...
                                          p, q, r, alpha, beta, gamma, steps,...
                                          logHandle, expHandle, prodHandle,...
                                          frameHandle, parHandle, adjustHandle );
            
            % run C++ implementation
            x_c = cppHandle( f(:,1:11) ,11, p, q, r,innerFactor, alpha, beta,gamma, steps, innerSteps);
            
            
            % break if error is too large
            error = norm(x_m(:) - x_c(:),2)/11;
            if  error > TOL
                flag = false;
                disp(['1D ' manifold ' unit test failed with error ' num2str(error)...
                      ', at p = ' num2str(p)...
                      ', q = ' num2str(q)...
                      ',f r = ' num2str(r)...
                      ', alpha = ' num2str(alpha)...
                      ', and beta = ' num2str(beta)]);
            else
                disp(['1D ' manifold ' unit test passed with error ' num2str(error)...
                      ', at p = ' num2str(p)...
                      ', q = ' num2str(q)...
                      ', r = ' num2str(r)...
                      ', alpha = ' num2str(alpha)...
                      ', and beta = ' num2str(beta)]);
            end
        end
    end
end

 
% final check
% clc
% if flag
%     disp(['Unit tests for ' manifold ' passed!']);
% else
%     disp(['Unit tests for ' manifold ' failed!']);
% end

