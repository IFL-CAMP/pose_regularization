function x = gaussianRand(mu,covar,sz)
% x = gaussianRand(mu,covar,sz)
    if nargin < 3
        sz = 1;
    end
    [nDim nMu] = size(mu);

    if length(sz) == 1 && nDim == 1
        sz = [1 sz];
    elseif length(sz) ~= 1 && nDim > 1
        error('stats:gaussianRand:BadSize','Bad sz parameter specified.  Must be a scalar for ND Gaussians.');
    end

    if nDim == 1
        x = mu + sqrt(covar)*randn(sz);
    else
        % TODO: If we ever feel like optimizing this we could start by
        % exploiting these cases for covar instead of defaulting to the
        % full covariance matrix case.
        nSigma = numel(covar);
        if nSigma == 1
            if any(covar < 0)
                error('stats:gaussianRand:BadScalarCovar','Covar must be symmetric, positive semi-definite.');
            end
            cholCovar = sqrt(covar);
            nnzDim = nDim;
        elseif nDim == nSigma
            if any(covar < 0)
                error('stats:gaussianRand:BadDiagCovar','Covar must be symmetric, positive semi-definite.');
            end
            cholCovar = diag(sqrt(covar));
            nnzDim = nDim;
        elseif nDim^2 == nSigma && ndims(covar) == 2 && all(size(covar) == [nDim nDim])
%             cholCovar = chol(covar).';
            cholCovar = covarFactor(covar).';
            nnzDim = size(cholCovar,2);
            if size(cholCovar,1) ~= nDim
                error('stats:gaussianRand:BadFullCovar','Covar must be symmetric, positive semi-definite.');
            end
        else
            error('stats:gaussianRand:BadParamSize','mu and covar have mismatched sizes');
        end

        if nMu == sz
            x = mu + cholCovar*randn([nnzDim,sz]);
        else
            x = bsxfun(@plus,mu,cholCovar*randn([nnzDim,sz]));
        end
    end
end
