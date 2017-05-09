function [f g] = binghamVonMisesFisherNegLikelihood(x,A,c)
    nDims = length(c);
    c = reshape(c,nDims,1);
    if size(x,1) ~= nDims
        if size(x,2) == nDims && size(x,1) == 1
            x = x.';
        else
            error('stats:binghamVonMisesFisherNegLikelihood:BadParameterSize', ...
                'The sizes of x and c do not match.');
        end
    end
    nData = size(x,2);
%     f = -(x.'*c + diag(x.'*A*x));
    f = -(x.'*c + sum(x.*(A*x),1).');
    if nargout > 1
%         g = -(repmat(c,[1,nData]) + 2*A*x);
        g = bsxfun(@minus,-c,2*A*x);
    end
end
