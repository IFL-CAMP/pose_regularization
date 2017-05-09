function [c,dc] = quatConst(q)
%     q = reshape(q,[],1);
    c = 0.5*(q.'*q - 1);
    if nargout > 1
        dc = q.';
    end
end