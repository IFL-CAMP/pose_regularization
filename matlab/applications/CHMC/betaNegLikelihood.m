function [f g] = betaNegLikelihood(x,a,b,normalized)
%[f g] = betaNegLikelihood(x,a,b)
%
% Evaluate the negative log likelihood of x under a beta distribution 
% with parameters a and b.
%
% E[X] = a/(a+b)
% E[(X-E[X])^2] = a*b/((a+b)^2*(a+b+1))
%
% f = -log( x^(a-1) * (1-x)^(b-1) / beta(a,b) )
% g = df/dx
    if nargin == 2 && length(a) == 2
        b = a(2);
        a = a(1);
    end

    f = zeros(size(x));
    posX = x > 0 & x < 1;
    f(posX) = -(a-1)*log( x(posX) ) - (b-1)*log( (1-x(posX)) );
    if normalized
        f(posX) = f(posX) + betaln(a,b);
    end
    f(~posX) = inf;
    if nargout > 1
        g = zeros(size(x));
        g(posX) = -(a-1)./x(posX) + (b-1)./(1-x(posX));
%         g(~posX) = 0; % do something more clever here maybe?
    end
end
