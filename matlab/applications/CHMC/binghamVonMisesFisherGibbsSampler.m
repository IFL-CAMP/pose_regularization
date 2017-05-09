function [x0,xs,stats] = binghamVonMisesFisherGibbsSampler(A,c,nSamples,x0,rejSample)
% [x0 xs] = binghamVonMisesFisherGibbsSampler(A,c,nSamples,x0)
%
% Do Gibbs sampling for a vector Bingham-von Mises-Fisher distribution.
% Generates nSamples starting from x0.  If x0 is not given, a default is
% used.  If nSamples isn't provided it defaults to 1.
%
% The vector Bingham-von Mises-Fisher distribution has pdf:
%    p(x) \propto exp(c.'*x + 0.5*x.'(A + A.')x) for ||x|| = 1
%
% This function implements the Gibbs sampler described by Peter Hoff,
% "Simulation of the matrix Bingham-von Mises-Fisher distribution..."
% (2007).
    nDims = size(A,1);
    if nargin < 4 || isempty(x0)
        x0 = randn([nDims,1]);
        x0 = x0/norm(x0);
    end
    if nargin < 3
        nSamples = 1;
    end
    if nargin < 4
        rejSample = false;
    end

    if nargout > 1
        xs = zeros(length(x0),nSamples-1);
    end
    startT = tic;
    
    A = 0.5*(A + A.');
    [E,L] = eig(A);
    Et = E.';
    L = diag(L);
%     L = L - median(L);
    d = Et*c;

    for sNum = 1:nSamples
        y = Et*x0;
        s = sign(y);
        y2 = y.^2;
%         prevy = y;
        for i = randperm(nDims)
%             prevyy = y;
            q = y2./(1-y2(i));
            subi = [1:(i-1) (i+1):nDims];
            aa = L(i) - (q(subi).'*L(subi));
            bb = sum(s(subi).*sqrt(q(subi)).*d(subi));
            if rejSample
                theta = sampleTheta((nDims-3)/2, aa, bb, d(i));
            else
                theta = sampleThetaGibbs(y2(i), (nDims-3)/2, aa, bb, d(i));
            end
            y2(i) = theta;
            y2(subi) = (1-theta).*q(subi);
            lp = sqrt(theta)*[ -d(i), d(i) ];
            lp = lp - logsum(lp);
            if log(rand) < lp(1)
                s(i) = -1;
            else
                s(i) = 1;
            end
%             assert(all(isreal(y)));
%             assert(abs(sum(y.^2)-1) < 1e-6);
%             y = y./norm(y);
            y2 = y2./sum(y2); % prevent drift
        end
%         assert(all(y2 >= 0));
        x0 = Et.'*(s.*sqrt(y2));
%         assert(all(isreal(x0)));
        if nargout > 1 && sNum < nSamples
            xs(:,sNum) = x0;
        end
    end
    
    if nargout > 2
        stats.N = nSamples;
        stats.t = toc(startT);
    end
end

function theta = sampleThetaGibbs(theta,k,a,b,c)
    d1=1/2;
    if k > min([a,b,c])
        d2 = k;
    else
        d2=1+min([k, max( [k-a,-1/2] )] );
    end
    
    betaBlockSize = 20;
    th = betaRand(d1,d2,[betaBlockSize,1]);
%     logptheta = -0.5*log(th) + k*log(th) + th.*a + (1-th).^(1/2) .* b + log(exp(-c/2) + exp(b/2));
%     logptheta = -0.5*log(th) + k*log(th) + th.*a + (1-th).^(1/2) .* b;
    logpth = computeLogPTheta(th,k,a,b,c);
    logptheta = computeLogPTheta(theta,k,a,b,c);
    logbth = -betaNegLikelihood(th,d1,d2,false);
    logbtheta = -betaNegLikelihood(theta,d1,d2,false);
    for i = 1:betaBlockSize
        acc = log(rand) < (logpth(i) - logptheta) + (logbtheta - logbth(i));
        if acc
            theta = th(i);
            logptheta = logpth(i);
            logbtheta = logbth(i);
        end
    end
end

function logptheta = computeLogPTheta(th,k,a,b,c)
    logptheta = -0.5*log(th) + k*log(1-th) + th.*a + (1-th).^(1/2) .* b + log(exp(-c*th.^0.5) + exp(c*th.^0.5));
end

function theta = sampleTheta(k,a,b,c)
    d1=1/2;
    if k > min([a,b,c])
        d2 = k;
    else
        d2 = 1+min([k, max( [k-a,-1/2] )] );
    end
    lmx = 0;
    if a > (k - d2 + 1) 
        lmx = (k-d2+1)*(log(k-d2+1)-log(a))+ a-(k-d2+1);
    end
  
    bb = max([0,b]);
    thmx = c^2/(c^2+bb^2);
    if abs(c)==0 && bb==0
        thmx = 1;
    end
  
    lfgmx = lmx + (b*sqrt(1-thmx)+abs(c)*sqrt(thmx)) + log(1+exp(-2*abs(c)));

    betaBlockSize = 20;
    while true
        th = betaRand(d1,d2,[betaBlockSize,1]);
%         ptheta = th.^(-1/2) .* (1-th).^k .* exp(th.*a + (1-th).^(1/2).*b) .* (exp(-c/2) + exp(b/2));
%         logptheta = -0.5*log(th) + k*log(th) + th.*a + (1-th).^(1/2) .* b + log(exp(-c/2) + exp(b/2));
%         logbtheta = -betaNegLikelihood(th,d1,d2);
        lfgth = (k-d2+1)*log(1-th) + a*th + (b*sqrt(1-th)+abs(c)*sqrt(th)) + log(1+exp(-2*sqrt(th)*abs(c)));
        tst = ( log(rand([betaBlockSize,1])) < lfgth - lfgmx );
        firstI = find(tst,1,'first');
        if ~isempty(firstI)
            theta = th(firstI);
            break;
        end
    end
    
    assert(theta >= 0);
    assert(theta <= 1);
    assert(isreal(theta));
end
