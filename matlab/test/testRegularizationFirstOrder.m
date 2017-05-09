close all;
clear;
clc;

% dbstop('error')

% define tolerance for relative error
TOL = 1e-2;

% define number of steps
steps = 10;
innerSteps = 1;
innerFactor  = 1;
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
    
    % cpp handle
    cppHandle = @denoiseSE3prod;
    
    % load data
    load('../data/testDataSE3.mat')
    
    % restrict data to 11x11x11
    f = f(:,1:11^3);
end
%% compare the implementations

flag = true;

%% initialize parameters for 2nd order case
r = -1;
beta = 0;
gamma = 0; %mixed second order param
frameHandle = [];
parHandle = [];
adjustHandle = [];

%% 1D tests
for alpha = [0.5 1 2]
    for p = [0 1 2]  
        for q = [0 1 2 ]
                        
            % run MATLAB implementation
            x_m = proximalPointAlgorithm( f(:,1:11), 11, 1, 1,...
                                          p, q, r, alpha, beta, gamma, steps,...
                                          logHandle, expHandle, prodHandle,...
                                          frameHandle, parHandle, adjustHandle );
            
            % run C++ implementation
            tic
            x_c = cppHandle( f(:,1:11), 11, p, q, r, innerFactor, alpha, beta,gamma, steps,innerSteps );
            toc
            
            % break if error is too large
            error = norm(x_m(:) - x_c(:),2)/11;
            if  error > TOL
                flag = false;
                disp(['1D ' manifold ' unit test failed with error ' num2str(error)...
                      ', at p = ' num2str(p)...
                      ', q = ' num2str(q)...
                      ', and alpha = ' num2str(alpha)]);
            else
                disp(['1D ' manifold ' unit test passed with error ' num2str(error)...
                      ', at p = ' num2str(p)...
                      ', q = ' num2str(q)...
                      ', and alpha = ' num2str(alpha)]);
            end

        end
    end
end



%% final check
% clc
if flag
    disp(['Unit tests for ' manifold ' passed!']);
else
    disp(['Unit tests for ' manifold ' failed!']);
end
%     
    