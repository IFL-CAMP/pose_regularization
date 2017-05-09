function [basis, eigvals] = frame(v, base)

% call Schur decomposition from Eigen as Matlab one is different 
v = reshape(v,4,4); 
v(:,4) = 0; %% CHECK this
[U, T] = schur(v);
e1 = [1; 0; 0; 0];
e2 = [0; 1; 0; 0];
e3 = [0; 0; 1; 0];
beta = max(T(:));

b1 = U*(e3*e2' - e2*e3')*U';
b2 = U*(e2*e1' - e1*e2')*U';
b3 = U*(e3*e1' - e1*e3')*U';
b4 = zeros(4,4);
b4(1,4) =1;
b5 = zeros(4,4);
b5(2,4) =1;
b6 = zeros(4,4);
b6(3,4) =1;
basis = [b1(:) b2(:) b3(:) b4(:) b5(:) b6(:)];
% eigvals = [0; beta*beta/4; beta*beta/4 ];
eigvals = [0; beta*beta/4; beta*beta/4 ;  0 ;  0;  0];
% basis = ones(size(basis));
end