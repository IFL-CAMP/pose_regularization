function [ denoisedMatrices ] = denoiseMatrices( matrices, p, q, r, alpha, beta, gamma, steps )
%DENOISEMATRICES Summary of this function goes here
%   Detailed explanation goes here

    % get number of poses
    N = size(matrices,3);

    % convert pose tracks for algorithm
    matricesColumn = zeros(16,N);
    for i = 1:N
        P = matrices(:,:,i);
        matricesColumn(:,i) = P(:);
    end

    % run algorithm without 2nd order regularization
    denoisedMatricesColumn = denoiseSE3prod(matricesColumn,...
                            N,1,1,...
                            p,q,r,5.0,...
                            alpha,beta,gamma,steps,50);


    % re-format poses
    for i = 1:N
        P = denoisedMatricesColumn(:,i);
        denoisedMatrices(:,:,i) = reshape(P,4,4);
    end

end

