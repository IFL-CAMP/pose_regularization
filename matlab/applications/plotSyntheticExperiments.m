close; %all
clear; %all
clc;

% load parameter study for 1st order
load paramsFirstOrder.mat

% define plot routine
pltRoutine = @(data) plotPoseMatrix(data,2.0);

% write png for GT data
plotInMaxScreenResolution(pltRoutine,DATA_GT{1})
f = getframe(gcf);
imwrite(f.cdata,'data_gt.png')
close

% generate images for noisy data
for i = 1:5
    plotInMaxScreenResolution(pltRoutine,DATA_NOISE{i})
    f = getframe(gcf);
    imwrite(f.cdata,['data_noisy_' num2str(i) '.png'])
    close
end

% generate images for 1st order denoised data
iteration_steps = 1000;
inner_steps = 100;
inner_factor = 0.25;
r = 0;
gamma = 0;
beta = 0;

for i = 1:5
    
    % load test data
    Pgt = DATA_GT{i};
    Pnoise = DATA_NOISE{i};
    
    % get number of poses
    N = size(Pgt,3);
    
    % reshape
    Pnoise = reshape(Pnoise,16,N);
    Pgt = reshape(Pgt,16,N);
    
    % compute reference error
    err_ref = comp_dist(Pgt,Pnoise);
    
    % find optimal parameters
    e = RESULTS{i}(:,2);
    [~, index] = min(e);
    p = RESULTS{i}(index,3);
    q = RESULTS{i}(index,4);
    alpha = RESULTS{i}(index,5);
    
    % denoise data
    Pden = denoiseSE3prod(Pnoise,N,1,1,...
                p,q,r,inner_factor,... % r = 0
                alpha,beta,gamma, iteration_steps,inner_steps); % gamma = 0
             
    plotInMaxScreenResolution(pltRoutine,reshape(Pden,4,4,N))
    f = getframe(gcf);
    imwrite(f.cdata,['data_den1st_' num2str(i) '.png'])
    close
end

% load parameter study for 2nd order
load paramsFirstAndSecondOrder.mat

% generate images for 2nd order denoised data
iteration_steps = 1000;
inner_steps = 100;
inner_factor = 0.25;
r = 0;
gamma = 0;

for i = 1:5
    
    % load test data
    Pgt = DATA_GT{i};
    Pnoise = DATA_NOISE{i};
    
    % get number of poses
    N = size(Pgt,3);
    
    % reshape
    Pnoise = reshape(Pnoise,16,N);
    Pgt = reshape(Pgt,16,N);
    
    % compute reference error
    err_ref = comp_dist(Pgt,Pnoise);
    
    % find optimal parameters
    e = RESULTS2{i}(:,2);
    [~, index] = min(e);
    p = RESULTS2{i}(index,3);
    q = RESULTS2{i}(index,4);
    alpha = RESULTS2{i}(index,5);
    beta = RESULTS2{i}(index,6);
    
    % denoise data
    Pden = denoiseSE3prod(Pnoise,N,1,1,...
                p,q,r,inner_factor,... % r = 0
                alpha,beta,gamma, iteration_steps,inner_steps); % gamma = 0
             
    plotInMaxScreenResolution(pltRoutine,reshape(Pden,4,4,N))
    f = getframe(gcf);
    imwrite(f.cdata,['data_den2nd_' num2str(i) '.png'])
    close
end


