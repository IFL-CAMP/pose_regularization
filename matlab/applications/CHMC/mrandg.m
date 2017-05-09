function x = mrandg(beta,sz)
% x = mrandg(beta,sz)
%
% Generate Gamma random variates with shape parameter beta and unit scale.
% This function has been checked against randg in the Statistics Toolbox.

    if isscalar(beta)
        if (beta > 1.0)
            x = grand_cheng_gb(beta,sz);
        elseif (beta == 1.0)
            x = -log(rand(sz));
        else
            x = grand_johnk(beta,sz);
        end
    else
        if ndims(beta) == 2 && any(size(beta) == 1)
            outSize = [length(beta),sz];
        else
            outSize = [size(beta),sz];
        end
        nSamples = prod(sz);
        x = zeros([numel(beta),nSamples]);

        k = find(beta > 1);
        for i = 1:length(k)
            kk = k(i);
            x(kk,:) = grand_cheng_gb(beta(kk),[1,nSamples]);
        end

        k = find(beta < 1);
        for i = 1:length(k)
            kk = k(i);
            x(kk,:) = grand_johnk(beta(kk),[1,nSamples]);
        end

        k = find(beta == 1);
        x(k,:) = -log(rand([length(k),nSamples]));

        x = reshape(x,outSize);
    end
end

function xr = grand_cheng_gb(alpha,sz)
    theta = 4.5;
    onePlusLogTheta = 1 + log(theta);

    a = 1.0/sqrt(2*alpha - 1);
    b = alpha - log(4);
    c = alpha + (1.0/a);

    xr = zeros(sz);
    Nxr = numel(xr);
    Nsampled = 0;
    while Nsampled < Nxr
        Nsample = 3*(Nxr - Nsampled);
        % Step 1
        u1 = rand([1,Nsample]);
        u2 = rand([1,Nsample]);

        % Step 2
        v = a*log(u1./(1.0-u1));
        x = alpha*exp(v);

        % Step 3'
        z = (u1.^2).*u2;
        r = b + (c*v) - x;
        accepted = find(((r + onePlusLogTheta - (theta*z)) >= 0) | (r >= log(z)));

        Naccepted = min([length(accepted),Nxr - Nsampled]);
        accepted = accepted(1:Naccepted);

        xr(Nsampled+1:Nsampled+Naccepted) = x(accepted);
        Nsampled = Nsampled + Naccepted;
    end
end

function xr = grand_johnk(a,sz)
    xr = zeros(sz);

    Nxr = numel(xr);
    Nsampled = 0;
    while Nsampled < Nxr
        Nsample = 3*(Nxr - Nsampled);
        u = rand([1,Nsample]);
        v = rand([1,Nsample]);
        x = u.^(1.0/a);
        y = v.^(1.0/(1.0-a));
        accepted = find((x+y) <= 1.0);

        Naccepted = min([length(accepted),Nxr - Nsampled]);
        accepted = accepted(1:Naccepted);
        e = -log(rand([1,Naccepted]));
        xr(Nsampled+1:Nsampled+Naccepted) = (x(1,accepted).*e)./(x(1,accepted)+y(1,accepted));
        Nsampled = Nsampled + Naccepted;
    end
end
