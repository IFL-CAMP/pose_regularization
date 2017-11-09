close all;
clear;
clc;

addpath('..')

iteration_steps = 1000;
count_kappa = 1;
RESULTS = {};
DATA_GT = {};
DATA_NOISE = {};
for testCase = 1:5
    
    % generate test data
    [Pgt, Pnoise] = genTestData(testCase);
    DATA_GT{count_kappa} = Pgt;
    DATA_NOISE{count_kappa} = Pnoise;
    
    % get number of poses
    N = size(Pgt,3);
    
    % reshape
    Pnoise = reshape(Pnoise,16,N);
    Pgt = reshape(Pgt,16,N);
    
    % compute reference error
    err_ref = comp_dist(Pgt,Pnoise);
    
    % loop over p, q and alpha
    count = 1;
    results = zeros(3*3*8,5);
    for p = [0,1,2]
        for q = [0,1,2]
            for alpha = [0.01, 0.1, 0.5, 1, 2, 5, 10, 20]
                
                Pden = denoiseSE3prod(Pnoise,N,1,1,...
                    p,q,0,0,... % r = 0, inner_factor = 0
                    alpha,0,0, iteration_steps,0); % beta = gamma = 0, inner_steps = 0
                
                err = comp_dist(Pgt,Pden);
                results(count,:) = [err_ref, err, p, q, alpha];
                count = count + 1;
                
            end
        end
    end
    
    RESULTS{count_kappa} = results;
    count_kappa = count_kappa + 1;
      
end

save('paramsFirstOrder.mat','RESULTS','DATA_GT','DATA_NOISE');
