function x = gammaRand(k,theta,sz)
% x = gammaRand(params,sz)
% x = gammaRand(k,theta,sz)
%
% If sz is not specificied, it is assumed that sz = 1.
    if length(k) == 2 && nargin < 3
        if nargin < 2
            sz = 1;
        else
            sz = theta;
        end
        theta = k(2);
        k = k(1);
    elseif nargin < 3
        sz = 1;
    end
    if length(sz) == 1
        sz = [sz 1];
    end
    x = theta*mrandg(k,sz);
end
