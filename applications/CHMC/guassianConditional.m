function [mupa,sigmapa] = guassianConditional(mu,sigma,aInds,bInds,c)
% compute the conditional distribution p(x(aInds)|x(bInds) = c) where
% p(x) = N(x|mu,sigma)

    sigmaa = sigma(aInds,aInds);
    sigmab = sigma(bInds,bInds);
    sigmac = sigma(aInds,bInds);
    mupa = mu(aInds) + sigmac * (sigmab\(c - mu(bInds)));
    sigmapa = sigmaa - sigmac * (sigmab \ sigmac.');
end