function x = betaRand(a,b,sz)
    if length(a) == 2 && nargin < 3
        if nargin < 2
            sz = 1;
        else
            sz = b;
        end
        b = a(2);
        a = a(1);
    elseif nargin < 3
        sz = 1;
    end
    if length(sz) == 1
        sz = [sz 1];
    end
    theta = 0.5*(1/a + 1/b);
    x = gammaRand(a,theta,sz);
    x2 = gammaRand(b,theta,sz);
    x = x ./ (x + x2);
end
