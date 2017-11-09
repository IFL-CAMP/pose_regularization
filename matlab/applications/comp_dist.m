function err = comp_dist(Pgt,Pden)

    err = 0;
    if size(Pgt, 1) == 16
        for i = 1:size(Pgt,2)
            err = err + comp_elem_dist(Pgt(:,i),Pden(:,i));
        end
    else
        for i = 1:size(Pgt,3)
            err = err + comp_elem_dist(Pgt(:,:,i),Pden(:,:,i));
        end
    end

end

function val = comp_elem_dist(base, dest)
    if size(base, 1) == 16
        base = reshape(base,4,4);
        dest = reshape(dest,4,4);
    end
    base(1:3,1:3) = base(1:3,1:3)';
    result = logSO3xR3matrix( base*dest );
    result(1:3,4) = dest(1:3,4)-base(1:3,4);
    
    %val = 0*(trace(result(1:3,1:3)'*result(1:3,1:3))) + result(1:3,4)' * result(1:3,4);
     val = 0.5*(trace(result(1:3,1:3)'*result(1:3,1:3))) + result(1:3,4)' * result(1:3,4);

end


function L = logSO3xR3matrix(A)

% use the same constants
PI = 3.141592653589793;
EPSILON = 0.000001;

R = A(1:3,1:3);

% phi = acos(0.5*(R(1,1)+R(2,2)+R(3,3) - 1));
phi = real(acos(0.5*(R(1,1)+R(2,2)+R(3,3) - 1))); 

L = zeros(4);
L(1:3,4) = A(1:3,4);

if ( abs(phi)>EPSILON && abs(phi)<PI-EPSILON )
    
    omegaHat = 0.5*phi*(R-R')/(sin(phi));
    
    L(1:3,1:3) = omegaHat;

end

end