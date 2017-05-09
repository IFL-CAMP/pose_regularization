function result = exp( base, dest )
    base = reshape(base,4,4);
    dest = reshape(dest,4,4);
    result = expSO3xR3matrix( dest );
    result = base*result;
    result(1:3,4) = base(1:3,4) + dest(1:3,4);
    result(4,4) = 1;
    result = result(:);
end

function result = expSO3xR3matrix(A)

% use the same constants
EPSILON = 0.000001;

% initialize solution
result = zeros(4,4);
% result(4,4) = 1;

% get omega and omega hat
omega = [A(3,2) A(1,3) A(2,1)];
omegaHat = A(1:3,1:3);

% compute norm of omega
normOmega = norm(omega,2);

if ( normOmega > EPSILON )
    
    % compute factors
    factor1 = sin(normOmega)/normOmega;
    factor2 = (1-cos(normOmega))/power(normOmega,2);

    
    
    % compute A
    result(1:3,1:3) = eye(3) + factor1*omegaHat + factor2*omegaHat*omegaHat;
    
else
    
    result(1:3,1:3) = eye(3);
    
end

% result(1:3,4) = A(1:3,4); %% check this out  - gets multiplied with base

end
