function ls = logsum(xx,dim)
% ls = logsum(x,dim)
%
% returns the log of sum of logs, summing over dimension dim
% computes ls = log(sum(exp(x),dim))
% but in a way that tries to avoid underflow/overflow
%
% basic idea: shift before exp and reshift back
% log(sum(exp(x))) = alpha + log(sum(exp(x-alpha)));
%

if numel(xx) <= 1 || (nargin > 1 && size(xx,dim)<=1)
    ls=xx;
    return;
end

xdims=size(xx);
if(nargin<2) 
    nonsingletons=find(xdims>1);
    dim=nonsingletons(1);
end

infxx = isinf(xx);
if any(infxx)
    xx(infxx) = -inf;
end

alpha = max(xx,[],dim)-log(realmax)/2+2*log(xdims(dim));
% repdims=ones(size(xdims)); repdims(dim)=xdims(dim);
% xx = xx-repmat(alpha,repdims);
xx = bsxfun(@minus,xx,alpha);
xx(infxx) = 0;
xx(~infxx) = exp(xx(~infxx));
ls = alpha+log(sum(xx,dim));


