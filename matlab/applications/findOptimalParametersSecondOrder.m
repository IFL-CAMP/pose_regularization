% close all;
% clear;
% clc;

addpath('..')

load paramsFirstOrder.mat

iteration_steps = 1000;
inner_steps = 50;
inner_factor = 0.25;
count_kappa = 1;
RESULTS2 = {};
for kappa = [10,50,100,500,1000]
    
    % load test data
    Pgt = DATA_GT{count_kappa};
    Pnoise = DATA_NOISE{count_kappa};
    
    % get number of poses
    N = size(Pgt,3);
    
    % reshape
    Pnoise = reshape(Pnoise,16,N);
    Pgt = reshape(Pgt,16,N);
    
    % find optimal parameters
    e = RESULTS{count_kappa}(:,1);
    [~, index] = min(e);
    p = RESULTS{count_kappa}(index,2);
    q = RESULTS{count_kappa}(index,3);
    alpha = RESULTS{count_kappa}(index,4);
    
    % loop over r and beta
    count = 1;
    results = zeros(12,6);
    for r = [1 2]
        for beta = [0, 0.1, 0.5, 1, 2, 5]
            Pden = denoiseSE3prod(Pnoise,N,1,1,...
                p,q,r,inner_factor,... % r = 0
                alpha,beta,0, iteration_steps,inner_steps); % gamma = 0
             
            err = comp_dist(Pgt,Pden);
            results(count,:) = [err, p, q, r, alpha, beta];
            count = count + 1;
            
        end
    end

    RESULTS2{count_kappa} = results;
    count_kappa = count_kappa + 1;
      
end

save('paramsFirstAndSecondOrder.mat','RESULTS2','DATA_GT','DATA_NOISE');

