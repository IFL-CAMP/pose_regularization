function [f,A] = linearConstraint(x,A,b)
    f = A*x-b*ones([1,size(x,2)]);
end