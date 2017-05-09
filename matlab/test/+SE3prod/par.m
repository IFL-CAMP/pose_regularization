function vectors = par(vectors, base, dest)

base(13:15) = 0;
dest(13:15) = 0;
v = reshape(SE3prod.log( base, dest ),4,4);
base = reshape(base,4,4);
dest = reshape(dest,4,4);
% base(1:3,4) = 0;
% dest(1:3,4) = 0;
% v(1:3,4) = 0;

expv2 = expSO3xR3matrix(v/2);
% expv2(1:3,4) = base(1:3,4) + v(1:3,4)/2;
for i = 1:size(vectors,2)-3;
    w = reshape(vectors(:,i),4,4);
    gamma = dest'*base*expv2*w*expv2;
    vectors(:,i) = gamma(:);
end

end

function result = expSO3xR3matrix(A)

% use the same constants
EPSILON = 0.000001;

% initialize solution
result = zeros(4,4);
result(4,4) = 1;

% get omega and omega hat
omega = [A(3,2) A(1,3) A(2,1)];
omegaHat = A(1:3,1:3);

% compute norm of omega
normOmega = norm(omega,2);

if ( normOmega > EPSILON )
    
    % compute factors
    factor1 = real(sin(normOmega)/normOmega);
    factor2 = real((1-cos(normOmega))/power(normOmega,2));
    
    % compute A
    result(1:3,1:3) = eye(3) + factor1*omegaHat + factor2*omegaHat*omegaHat;
    
else
    
    result(1:3,1:3) = eye(3);
    
end

% result(1:3,4) = A(1:3,4);

end
