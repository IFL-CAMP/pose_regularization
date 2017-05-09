function [f,g] = gammaNegLikelihood(x,k,theta)
% [f,g] = gammaNegLikelihood(x,k,theta)
%
% Evaluate the negative log likelihood of x under a gamma distribution 
% with shape parameter k and scale parameter theta.
%
% E[X] = k*theta
% E[(X-E[X])^2] = k*theta^2
%
% f = -log( (x^(k-1)/(theta.^k*gamma(k))) * exp(-x/theta) )
% g = df/dx
    if nargin == 2 && length(k) == 2
        theta = k(2);
        k = k(1);
    end

    f = zeros(size(x));
    posX = x > 0;
    f(posX) = -(k-1)*log(x(posX)) + k*log(theta) + gammaln(k) + x(posX)/theta;
    f(~posX) = inf;
    if nargout > 1
        g = zeros(size(x));
        g(posX) = -(k-1) ./ x(posX) + 1./theta;
%         g(~posX) = 0; % do something more clever here maybe?
    end
end
