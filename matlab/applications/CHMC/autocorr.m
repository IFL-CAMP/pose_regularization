function autocorr = autocorr(x,dim,doMeanSub)
  if nargin < 2
      dim = 1;
  end
  if nargin < 3
      doMeanSub = true;
  end
  if doMeanSub
      x = bsxfun(@minus,x,mean(x,dim));
  end
  n = size(x,dim);
  m = 2^(ceil(log2(n)));
  f = fft(x,m,dim);
  a = real(ifft(f.*conj(f),m,dim));
  permDims = 1:ndims(x);
  permDims = [dim permDims(1:(dim-1)) permDims((dim+1):end)];
  a = permute(a,permDims);
  autocorr = bsxfun(@rdivide,a(1:n,:),a(1,:));
  sza = size(a);
  sza(1) = n;
  autocorr = ipermute(reshape(autocorr,sza),permDims);
end