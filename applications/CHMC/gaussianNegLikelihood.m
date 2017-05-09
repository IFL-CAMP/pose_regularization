function [f, dfdx, dfdmu] = gaussianNegLikelihood(x,mu,covar,normalized)
% [f, dfdx, dfdmu] = gaussianNegLikelihood(x,mu,covar,normalized)
% 
% Gaussian negative log likelihood with covariance matrix covar
    if nargin < 4
        normalized = true;
    end
    
    [nDim nPts] = size(x);
    if nPts == 0
        f = [];
        dfdx = [];
        dfdmu = [];
        return;
    end

    if nDim == 1
%         covar = reshape(covar,1,[]);
        f = 0.5*(x - mu).^2./covar;
        if normalized
            f = f + 0.5*log(2*pi*covar);
        end
        if nargout > 1
            dfdx = (x - mu)./covar;
        end
        if nargout > 2
            dfdmu = -dfdx;
        end

    else
%         if nargout > 2
%             [f dfdx dfdmu] = nll_gaussian(x,mu,covar,normalized,0);
%         elseif nargout > 1
%             [f dfdx] = nll_gaussian(x,mu,covar,normalized,0);
%         else 
%             f = nll_gaussian(x,mu,covar,normalized,0);
%         end
        if nDim ~= size(x,1)
            error('stats:gaussianNegLikelihood:BadXSize','size(x,1) ~= numel(mu)');
        end
        
        [c1,c2,c3] = size(covar);
        nSigma = c1*c2;
        if size(x,2) == size(mu,2), err = x-mu; else err = bsxfun(@minus,x,mu); end
        nErr = size(err,2);
        if c3 ~= 1 && c3 ~= nErr
            error('stats:gaussianNegLikelihood:BadCovarSize','size(covar,3) ~= 1 or size(x,2)');
        end
        if nSigma == 1 % scalar times identity covariance
            covar = reshape(covar,[1,c3]);
            if normalized
                logdetSigma = nDim*log(covar);
            end
            covarInvErr = bsxfun(@rdivide,err,covar);
        elseif nDim == nSigma % diagonal covariance
            covar = reshape(covar,[nDim,c3]);
            if normalized
                logdetSigma = sum(log(covar));
            end
            covarInvErr = bsxfun(@rdivide,err,covar);
        elseif nDim^2 == nSigma
            if c3 == 1
                cholCovar = chol(covar,'lower');
                if normalized
                    logdetSigma = 2*sum(log(diag(cholCovar)));
                end
                covarInvErr = cholCovar.'\(cholCovar\err);
            else
                covarInvErr = zeros(size(err));
                logdetSigma = zeros([1,c3]);
                for i = 1:c3
                    cholCovar = chol(covar(:,:,i),'lower');
                    if normalized
                        logdetSigma(1,i) = 2*sum(log(diag(cholCovar)));
                    end
                    covarInvErr(:,i) = cholCovar.'\(cholCovar\err(:,i));
                end
            end
        else
            error('stats:gaussianNegLikelihood:BadParamSize','mu and covar have mismatched sizes');
        end
        
        f = 0.5*sum(covarInvErr.*err,1);
        if normalized
            f = f + 0.5*(nDim*log(2*pi) + logdetSigma);
        end
        if nargout > 1
%             dfdx = bsxfun(@minus,covarInvX,covarInvMu);
            dfdx = covarInvErr;
        end
        if nargout > 2
            dfdmu = -dfdx;
        end
    end
end
