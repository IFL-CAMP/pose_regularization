function value = prod( v, w, base )

% base not relevant due to definition of log as tangential element at Id!

% to be changed
a = 0.5;
b = 0.5;
% v = reshape(v,4,4);
% w = reshape(w,4,4);
% 
% omega_v = [v(3,2); v(1,3); v(2,1)];
% omega_w = [w(3,2); w(1,3); w(2,1)];
% t_v = v(1:3,4);
% t_w = w(1:3,4);
% 
% value = a*omega_v'*omega_w + b*t_v'*t_w;
%     value = 0;
%     for i=1:15
%     value = value + v(i)*w(i);   
%     end
    v = reshape(v,4,4);
    w = reshape(w,4,4);

    value = a*trace(v(1:3,1:3)'*w(1:3,1:3)) + b*v(1:3,4)' *w(1:3,4);

end 